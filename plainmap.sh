#!/usr/bin/env bash
# PlainMap v0.1 (development)
# Transparent, failure-aware mapping pipeline for ancient and modern DNA
#
# This iteration:
#   - Mixed SE + PE manifests supported (header-based pairing):
#       * Illumina/CASAVA: " 1:" / " 2:"
#       * ENA-style: "/1" / "/2"
#   - FASTP strategy:
#       * For each ORIGINAL unit (SE file or PE file-pair):
#           - If unit read count > max_reads_per_chunk: chunk BEFORE fastp (fastp safety valve for huge datasets)
#           - Else: run fastp once on the unit
#       * Append trimmed outputs into pooled trimmed FASTQs (gzip members concatenated)
#       * Map ONCE from pooled trimmed outputs (avoids many BWA/sort invocations)
#   - Prefer pigz where possible (parallel compression/decompression)
#   - fastp HTML is skipped (JSON only)
#   - Stats use unambiguous “fragments” terminology and fragment-aware mapping counts
#   - avg_depth computed correctly from samtools coverage; pct_covered reported

set -euo pipefail

###############################################################################
# REQUIRE BASH
###############################################################################
if [[ -z "${BASH_VERSION:-}" ]]; then
  echo "ERROR: PlainMap must be run with bash (not sh)."
  exit 1
fi

###############################################################################
# DEFAULTS
###############################################################################
VERSION="0.1"

THREADS=1
MINLEN=30
MISMATCH=0.01
LIBTYPE="modern"            # modern|ancient
MAX_READS_PER_CHUNK=0       # 0 disables pre-fastp chunking
MAPQ=20

ADAPTER_R1=""
ADAPTER_R2=""
TRIM_ONLY=0

DRYRUN=0
VALIDATE_ONLY=0
RESUME=1
KEEP_INTERMEDIATE=0
TMPDIR_USER=""
RESET=0

FASTP=${FASTP:-fastp}
BWA=${BWA:-bwa}
SAMTOOLS=${SAMTOOLS:-samtools}
MAPDAMAGE=${MAPDAMAGE:-mapDamage}
PYTHON=${PYTHON:-python3}

###############################################################################
# HELP
###############################################################################
usage() {
cat <<EOF
PlainMap v$VERSION (development)

Required:
  -manifest FILE
  -prefix STRING
  -ref FILE
  -outdir DIR

Optional:
  -library-type modern|ancient    (default: modern)
  -t, --threads INT               (default: 1)
  -minlength INT                  (default: 30)
  -mismatch FLOAT                 (ancient only; default: 0.01)
  -max-reads-per-chunk INT        (default: 0; disabled)
      Pre-fastp safety valve: if a single input unit has >INT reads, it is split into chunks
      (<=INT reads each) and fastp runs per chunk. Small units skip chunking automatically.
  --adapter-r1 SEQ                (explicit adapter; otherwise detection is used)
  --adapter-r2 SEQ
  --trim-only                     Trim only (fastp); exit before mapping
  --tmpdir DIR                    Custom temp dir (e.g. node-local scratch)
  --keep-intermediate             Keep <outdir>/<prefix>/work
  --resume | --no-resume
  --dry-run                       Print plan only (no gzip/pigz -t, no work)
  --validate                      Check tools + gzip/pigz -t all manifest FASTQs, then exit
  --reset                         Clear ALL checkpoints for this sample

Tool overrides:
  --fastp CMD
  --bwa CMD
  --samtools CMD
  --mapdamage CMD
  --python CMD
EOF
exit 0
}

###############################################################################
# PARSE ARGS
###############################################################################
MANIFEST=""
SAMPLE=""
REF=""
OUT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -manifest) MANIFEST="$2"; shift 2 ;;
    -prefix) SAMPLE="$2"; shift 2 ;;
    -ref) REF="$2"; shift 2 ;;
    -outdir) OUT="$2"; shift 2 ;;
    -library-type) LIBTYPE="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    -minlength) MINLEN="$2"; shift 2 ;;
    -mismatch) MISMATCH="$2"; shift 2 ;;
    -max-reads-per-chunk) MAX_READS_PER_CHUNK="$2"; shift 2 ;;
    --adapter-r1) ADAPTER_R1="$2"; shift 2 ;;
    --adapter-r2) ADAPTER_R2="$2"; shift 2 ;;
    --trim-only) TRIM_ONLY=1; shift ;;
    --tmpdir) TMPDIR_USER="$2"; shift 2 ;;
    --keep-intermediate) KEEP_INTERMEDIATE=1; shift ;;
    --resume) RESUME=1; shift ;;
    --no-resume) RESUME=0; shift ;;
    --dry-run) DRYRUN=1; shift ;;
    --validate) VALIDATE_ONLY=1; shift ;;
    --reset) RESET=1; shift ;;
    --fastp) FASTP="$2"; shift 2 ;;
    --bwa) BWA="$2"; shift 2 ;;
    --samtools) SAMTOOLS="$2"; shift 2 ;;
    --mapdamage) MAPDAMAGE="$2"; shift 2 ;;
    --python) PYTHON="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown option $1"; usage ;;
  esac
done

[[ -z "$MANIFEST" || -z "$SAMPLE" || -z "$REF" || -z "$OUT" ]] && usage
[[ "$LIBTYPE" == "modern" || "$LIBTYPE" == "ancient" ]] || { echo "ERROR: -library-type must be modern|ancient"; exit 1; }

###############################################################################
# PATH RESOLUTION (ABSOLUTE, SYMLINK-SAFE)
###############################################################################
resolve_path() {
  "$PYTHON" - "$1" <<'PY'
import os,sys
p=sys.argv[1]
p=os.path.expanduser(p)
print(os.path.realpath(os.path.abspath(p)))
PY
}

mkdir -p "$OUT"
MANIFEST="$(resolve_path "$MANIFEST")"
REF="$(resolve_path "$REF")"
OUT="$(resolve_path "$OUT")"
MANIFEST_DIR="$(dirname "$MANIFEST")"

###############################################################################
# LOGGING (ROTATE LOG PER RUN)
###############################################################################
LOG="$OUT/${SAMPLE}_plainmap.log"
if [[ -f "$LOG" ]]; then
  ts=$(date '+%Y%m%d-%H%M%S')
  mv "$LOG" "$OUT/${SAMPLE}_plainmap.${ts}.log"
fi
exec > >(tee "$LOG") 2>&1

log() { echo "[$(date '+%F %T')] $*"; }
die() { log "ERROR: $*"; exit 1; }

PIPE_T0=$(date +%s)
SORT_THREADS=$(( THREADS > 2 ? THREADS / 2 : 1 ))

###############################################################################
# PICK gzip vs pigz (for compression/decompression/tests)
###############################################################################
HAVE_PIGZ=0
if command -v pigz >/dev/null 2>&1; then
  HAVE_PIGZ=1
fi

gz_test() {
  local f="$1"
  if [[ $HAVE_PIGZ -eq 1 ]]; then
    pigz -p "$THREADS" -t "$f" >/dev/null
  else
    gzip -t "$f" >/dev/null
  fi
}

# Used for split --filter. Must write to "$FILE.fastq.gz".
gz_filter_cmd() {
  if [[ $HAVE_PIGZ -eq 1 ]]; then
    echo "pigz -p $THREADS > \$FILE.fastq.gz"
  else
    echo "gzip > \$FILE.fastq.gz"
  fi
}

# Decompress command for process substitution (mapping).
decomp_cmd() {
  local f="$1"
  if [[ $HAVE_PIGZ -eq 1 ]]; then
    echo "pigz -dc -p $THREADS \"$f\""
  else
    echo "zcat \"$f\""
  fi
}

###############################################################################
# WORK DIRS
###############################################################################
WORK="$OUT/${SAMPLE}/work"
RAW="$WORK/raw"              # symlinks to original inputs (per unit)
CHUNKS="$WORK/chunks"        # raw chunks (only for big units)
POOL="$WORK/pool"            # pooled trimmed outputs (.fastq.gz)
FINAL="$WORK/final"
CKPT="$WORK/.checkpoints"
REPORTS="$OUT/${SAMPLE}/fastp_reports"   # JSON only

STATS="$OUT/${SAMPLE}.plainmap.stats.tsv"
COV_TSV="$OUT/${SAMPLE}.plainmap.coverage.tsv"

if [[ -n "$TMPDIR_USER" ]]; then
  TMP="$TMPDIR_USER/plainmap_${SAMPLE}"
else
  TMP="$WORK/tmp"
fi

mkdir -p "$RAW" "$CHUNKS" "$POOL" "$FINAL" "$CKPT" "$TMP" "$REPORTS"

###############################################################################
# FILE CHECKS (SYMLINK-SAFE)
###############################################################################
file_nonempty() {
  local f="$1"
  [[ -e "$f" ]] || return 1
  local t="$f"
  if command -v readlink >/dev/null 2>&1; then
    for _ in 1 2 3 4 5 6 7 8 9 10; do
      [[ -L "$t" ]] || break
      t=$(readlink -f "$t" 2>/dev/null || readlink "$t" 2>/dev/null || echo "$t")
    done
  fi
  [[ -s "$t" ]]
}

###############################################################################
# CHECKPOINTS
###############################################################################
ckpt_file() { echo "$CKPT/${SAMPLE}.$1.done"; }
ckpt_mark() { [[ $RESUME -eq 1 ]] && date -Is > "$(ckpt_file "$1")"; }
ckpt_clear() { rm -f "$(ckpt_file "$1")" 2>/dev/null || true; }
files_ok() { local f; for f in "$@"; do file_nonempty "$f" || return 1; done; return 0; }
ckpt_ok() { local step="$1"; shift; [[ $RESUME -eq 1 ]] && [[ -f "$(ckpt_file "$step")" ]] && files_ok "$@"; }

if [[ $RESET -eq 1 ]]; then
  log "Reset requested: clearing all checkpoints for sample $SAMPLE"
  rm -f "$CKPT/${SAMPLE}."*.done 2>/dev/null || true
fi

###############################################################################
# PREFLIGHT
###############################################################################
require_cmd() { command -v "$1" >/dev/null 2>&1 || die "Tool not found in PATH: $1"; }
check_split_filter_support() { split --help 2>/dev/null | grep -q -- '--filter' || die "Need GNU split with --filter support"; }

quick_preflight() {
  [[ -r "$MANIFEST" ]] || die "Manifest not readable: $MANIFEST"
  [[ -r "$REF" ]] || die "Reference not readable: $REF"
  require_cmd zcat
  require_cmd split
  require_cmd "$FASTP"; require_cmd "$BWA"; require_cmd "$SAMTOOLS"; require_cmd "$PYTHON"
  check_split_filter_support
  [[ "$LIBTYPE" != "ancient" ]] || require_cmd "$MAPDAMAGE"
  if [[ $HAVE_PIGZ -eq 1 ]]; then
    log "Compression: pigz (threads=$THREADS)"
  else
    require_cmd gzip
    log "Compression: gzip"
  fi
}

if [[ $VALIDATE_ONLY -eq 1 ]]; then
  log "Validation mode: checking tools + gzip integrity for all manifest FASTQs"
  quick_preflight
  n=0
  while read -r f; do
    [[ -z "$f" || "$f" =~ ^# ]] && continue
    n=$((n+1))
    [[ "$f" = /* ]] || f="$MANIFEST_DIR/$f"
    f="$(resolve_path "$f")"
    [[ -r "$f" ]] || die "FASTQ not readable: $f"
    gz_test "$f" || die "Corrupt/truncated gzip: $f"
  done < "$MANIFEST"
  [[ $n -gt 0 ]] || die "Manifest empty: $MANIFEST"
  log "Validation OK"
  exit 0
fi

if [[ $DRYRUN -eq 1 ]]; then
  log "Dry-run mode: printing plan only"
  quick_preflight
  log "Resolved paths:"
  log "  manifest: $MANIFEST"
  log "  ref:      $REF"
  log "  outdir:   $OUT"
  log "Work directory: $WORK"
  log "Plan:"
  log "  1) Classify manifest into PE keys and SE keys (Illumina/ENA header detection)"
  log "  2) For each original unit: optional chunking BEFORE fastp if >max_reads_per_chunk"
  log "  3) fastp per unit/chunk (JSON only), append trimmed outputs into pooled .fastq.gz"
  log "  4) Map once from pooled trimmed inputs"
  log "  5) Merge -> dedup -> RG -> mapDamage (ancient) -> coverage/depth -> stats -> cleanup"
  exit 0
fi

quick_preflight

###############################################################################
# BWA INDEX (SAFE)
###############################################################################
ensure_bwa_index() {
  local ref="$1"
  local lock="${ref}.bwa.lock"
  local idx=( "${ref}.bwt" "${ref}.pac" "${ref}.ann" "${ref}.amb" "${ref}.sa" )
  local missing=0
  for f in "${idx[@]}"; do [[ -f "$f" ]] || missing=1; done
  [[ "$missing" -eq 0 ]] && { log "Reference index: present"; return; }

  if mkdir "$lock" 2>/dev/null; then
    trap 'rm -rf "$lock"' EXIT
    log "Indexing reference"
    "$BWA" index "$ref"
    rm -rf "$lock"; trap - EXIT
  else
    log "Waiting for BWA index to finish"
    while :; do
      sleep 10
      missing=0
      for f in "${idx[@]}"; do [[ -f "$f" ]] || missing=1; done
      [[ "$missing" -eq 0 ]] && break
    done
  fi
}
ensure_bwa_index "$REF"

###############################################################################
# FASTP WRAPPER + ARGS (JSON ONLY)
###############################################################################
run_fastp() {
  local stderr_file="$1"; shift
  log "$FASTP $*"
  : > "$stderr_file"
  set +e
  "$FASTP" "$@" 2> >(tee -a "$stderr_file" >&2)
  local status=$?
  set -e
  [[ $status -eq 0 ]] || die "fastp exited non-zero ($status)"
  if grep -qiE "igzip|invalid gzip|premature|truncated|error" "$stderr_file"; then
    die "fastp reported gzip/igzip error (see $stderr_file)"
  fi
}

# Base args; PE adapter detection appended automatically unless explicit adapters are supplied.
FASTP_ARGS=(-l "$MINLEN" -g -w "$THREADS")
[[ -n "$ADAPTER_R1" ]] && FASTP_ARGS+=(--adapter_sequence "$ADAPTER_R1")
[[ -n "$ADAPTER_R2" ]] && FASTP_ARGS+=(--adapter_sequence_r2 "$ADAPTER_R2")

###############################################################################
# STATS HELPERS
###############################################################################
count_reads_fastq_gz() { zcat "$1" | awk 'END{print NR/4}'; }

json_get_int() {
  "$PYTHON" - "$1" "$2" <<'PY'
import json,sys
path=sys.argv[1]; keys=sys.argv[2].split(".")
with open(path) as f: d=json.load(f)
cur=d
try:
    for k in keys: cur=cur[k]
except Exception:
    cur=0
try: print(int(cur))
except Exception: print(0)
PY
}

pct() {
  "$PYTHON" - "$1" "$2" <<'PY'
import sys
n=float(sys.argv[1]); d=float(sys.argv[2])
print("0" if d==0 else f"{(100.0*n/d):.6f}")
PY
}

ratio() {
  "$PYTHON" - "$1" "$2" <<'PY'
import sys
n=float(sys.argv[1]); d=float(sys.argv[2])
print("0" if d==0 else f"{(n/d):.6f}")
PY
}

hash8() {
  "$PYTHON" - "$1" <<'PY'
import sys,hashlib
s=sys.argv[1].encode()
print(hashlib.md5(s).hexdigest()[:8])
PY
}

###############################################################################
# HEADER PARSING / UNIT CLASSIFICATION
###############################################################################
# Key is derived from read header token (first token) with /1 or /2 stripped.
declare -A KEY_R1_FILES
declare -A KEY_R2_FILES

is_r1_header() { local hdr="$1"; [[ "$hdr" =~ [[:space:]]1: ]] || [[ "$hdr" =~ /1([[:space:]]|$) ]]; }
is_r2_header() { local hdr="$1"; [[ "$hdr" =~ [[:space:]]2: ]] || [[ "$hdr" =~ /2([[:space:]]|$) ]]; }
read_token() { local hdr="$1"; printf "%s" "${hdr%%[[:space:]]*}"; }

read_key_from_hdr() {
  local hdr="$1"
  local tok
  tok="$(read_token "$hdr")"
  tok="${tok%/1}"
  tok="${tok%/2}"
  printf "%s" "$tok"
}

append_file() {
  local -n A="$1"
  local key="$2"
  local file="$3"
  if [[ -n "${A[$key]:-}" ]]; then
    A["$key"]+=$'\n'"$file"
  else
    A["$key"]="$file"
  fi
}

val_to_array() {
  local v="$1"
  mapfile -t _ARR <<<"$v"
}

log "STEP: reading manifest and classifying FASTQs (supports mixed SE/PE)"
n_fastq=0
while read -r f; do
  [[ -z "$f" || "$f" =~ ^# ]] && continue
  [[ "$f" = /* ]] || f="$MANIFEST_DIR/$f"
  f="$(resolve_path "$f")"
  [[ -r "$f" ]] || die "FASTQ not readable: $f"
  n_fastq=$((n_fastq+1))

  hdr=$(zcat "$f" | head -n1 || true)
  [[ -n "$hdr" ]] || die "Could not read FASTQ header (empty?): $f"

  key="$(read_key_from_hdr "$hdr")"
  [[ -n "$key" ]] || die "Could not derive read key from header: $hdr (file: $f)"

  if is_r1_header "$hdr"; then
    append_file KEY_R1_FILES "$key" "$f"
  elif is_r2_header "$hdr"; then
    append_file KEY_R2_FILES "$key" "$f"
  else
    die "Could not determine read direction (R1/R2) from FASTQ header: $hdr (file: $f)"
  fi
done < "$MANIFEST"

[[ $n_fastq -gt 0 ]] || die "Manifest empty: $MANIFEST"

# Disallow R2-only keys
for key in "${!KEY_R2_FILES[@]}"; do
  [[ -n "${KEY_R1_FILES[$key]:-}" ]] || die "Found R2 with no matching R1 key: $key"
done

KEYS=( "${!KEY_R1_FILES[@]}" )
[[ ${#KEYS[@]} -gt 0 ]] || die "No R1 reads found in manifest"
IFS=$'\n' KEYS=( $(printf "%s\n" "${KEYS[@]}" | sort) ); unset IFS

N_SE_KEYS=0
N_PE_KEYS=0
for key in "${KEYS[@]}"; do
  if [[ -n "${KEY_R2_FILES[$key]:-}" ]]; then
    N_PE_KEYS=$((N_PE_KEYS+1))
  else
    N_SE_KEYS=$((N_SE_KEYS+1))
  fi
done
log "Detected keys: total=${#KEYS[@]} (PE keys=$N_PE_KEYS, SE keys=$N_SE_KEYS)"

SEQ_MODE="SE"
[[ $N_PE_KEYS -gt 0 ]] && SEQ_MODE="PE"
[[ $N_PE_KEYS -gt 0 && $N_SE_KEYS -gt 0 ]] && SEQ_MODE="MIX"

###############################################################################
# RAW CHUNKING (PER UNIT, ADAPTIVE)
###############################################################################
split_fastq() {
  local fq="$1" prefix="$2" reads="$3"
  local lines=$((reads * 4))
  local filter
  filter="$(gz_filter_cmd)"
  zcat "$fq" | split -l "$lines" -d -a 4 --filter="$filter" - "$prefix"
}

###############################################################################
# POOLED TRIMMED OUTPUTS (gzip concatenation)
###############################################################################
POOL_SE_GZ="$POOL/${SAMPLE}.pool.SE.fastq.gz"                 # SE inputs (true SE)
POOL_PE_R1_GZ="$POOL/${SAMPLE}.pool.PE.R1.fastq.gz"           # modern PE trimmed R1
POOL_PE_R2_GZ="$POOL/${SAMPLE}.pool.PE.R2.fastq.gz"           # modern PE trimmed R2
POOL_SE_RESCUE_GZ="$POOL/${SAMPLE}.pool.SE.rescue.fastq.gz"   # modern PE merged+unpaired
POOL_ANCIENT_MERGED_GZ="$POOL/${SAMPLE}.pool.ancient.merged.fastq.gz"  # ancient PE merged reads

# Initialize pools as empty gzip streams (we can just truncate; appending concatenated gzip members is valid)
: > "$POOL_SE_GZ"
: > "$POOL_PE_R1_GZ"
: > "$POOL_PE_R2_GZ"
: > "$POOL_SE_RESCUE_GZ"
: > "$POOL_ANCIENT_MERGED_GZ"

###############################################################################
# GLOBAL ACCUMULATORS (FRAGMENT-SPACE)
###############################################################################
RAW_R1_READS=0
RAW_R2_READS=0
RAW_FRAGMENTS=0

TRIMMED_FRAGMENTS=0
MERGED_READS=0
UNPAIRED_READS=0

###############################################################################
# PROCESS: per-unit optional chunking BEFORE fastp; fastp; append to pools
###############################################################################
log "STEP: fastp per unit (adaptive chunking before fastp), append trimmed outputs to pooled FASTQs"

# Helper: decide chunking for a given unit by R1 read count
unit_make_raw_chunks_se() {
  local in_r1="$1" unit_tag="$2"
  local reads
  reads=$(count_reads_fastq_gz "$in_r1")
  if [[ "$MAX_READS_PER_CHUNK" -gt 0 && "$reads" -gt "$MAX_READS_PER_CHUNK" ]]; then
    log "  [$unit_tag] Pre-fastp chunking enabled (reads=$reads > $MAX_READS_PER_CHUNK)"
    rm -f "$CHUNKS/${unit_tag}.R1_"*.fastq.gz 2>/dev/null || true
    split_fastq "$in_r1" "$CHUNKS/${unit_tag}.R1_" "$MAX_READS_PER_CHUNK"
  else
    log "  [$unit_tag] No pre-fastp chunking (reads=$reads)"
    rm -f "$CHUNKS/${unit_tag}.R1_0000.fastq.gz" 2>/dev/null || true
    ln -sf "$in_r1" "$CHUNKS/${unit_tag}.R1_0000.fastq.gz"
  fi
}

unit_make_raw_chunks_pe() {
  local in_r1="$1" in_r2="$2" unit_tag="$3"
  local reads1 reads2
  reads1=$(count_reads_fastq_gz "$in_r1")
  reads2=$(count_reads_fastq_gz "$in_r2")
  [[ "$reads1" -eq "$reads2" ]] || die "  [$unit_tag] R1/R2 read count mismatch (R1=$reads1 R2=$reads2). Refusing to chunk unsynchronized PE."

  if [[ "$MAX_READS_PER_CHUNK" -gt 0 && "$reads1" -gt "$MAX_READS_PER_CHUNK" ]]; then
    log "  [$unit_tag] Pre-fastp chunking enabled (reads=$reads1 > $MAX_READS_PER_CHUNK)"
    rm -f "$CHUNKS/${unit_tag}.R1_"*.fastq.gz "$CHUNKS/${unit_tag}.R2_"*.fastq.gz 2>/dev/null || true
    split_fastq "$in_r1" "$CHUNKS/${unit_tag}.R1_" "$MAX_READS_PER_CHUNK"
    split_fastq "$in_r2" "$CHUNKS/${unit_tag}.R2_" "$MAX_READS_PER_CHUNK"
  else
    log "  [$unit_tag] No pre-fastp chunking (reads=$reads1)"
    rm -f "$CHUNKS/${unit_tag}.R1_0000.fastq.gz" "$CHUNKS/${unit_tag}.R2_0000.fastq.gz" 2>/dev/null || true
    ln -sf "$in_r1" "$CHUNKS/${unit_tag}.R1_0000.fastq.gz"
    ln -sf "$in_r2" "$CHUNKS/${unit_tag}.R2_0000.fastq.gz"
  fi
}

# Determine whether to use PE adapter detection
FASTP_PE_EXTRA=()
if [[ -z "$ADAPTER_R1" && -z "$ADAPTER_R2" ]]; then
  FASTP_PE_EXTRA+=(--detect_adapter_for_pe)
fi

for key in "${KEYS[@]}"; do
  kid="$(hash8 "$key")"
  has_r2=0
  [[ -n "${KEY_R2_FILES[$key]:-}" ]] && has_r2=1

  val_to_array "${KEY_R1_FILES[$key]}"; R1_FILES=( "${_ARR[@]}" )
  if [[ $has_r2 -eq 1 ]]; then
    val_to_array "${KEY_R2_FILES[$key]}"; R2_FILES=( "${_ARR[@]}" )
  else
    R2_FILES=()
  fi

  log "Key K${kid} mode=$([[ $has_r2 -eq 1 ]] && echo PE || echo SE) nR1=${#R1_FILES[@]} nR2=${#R2_FILES[@]}"

  # Raw counts
  for f in "${R1_FILES[@]}"; do
    c=$(count_reads_fastq_gz "$f")
    RAW_R1_READS=$("$PYTHON" - <<PY
print(int($RAW_R1_READS) + int($c))
PY
)
    RAW_FRAGMENTS=$("$PYTHON" - <<PY
print(int($RAW_FRAGMENTS) + int($c))
PY
)
  done
  if [[ $has_r2 -eq 1 ]]; then
    for f in "${R2_FILES[@]}"; do
      c=$(count_reads_fastq_gz "$f")
      RAW_R2_READS=$("$PYTHON" - <<PY
print(int($RAW_R2_READS) + int($c))
PY
)
    done
  fi

  if [[ $has_r2 -eq 1 ]]; then
    [[ ${#R1_FILES[@]} -eq ${#R2_FILES[@]} ]] || die "Key mismatch: key=$key has ${#R1_FILES[@]} R1 files but ${#R2_FILES[@]} R2 files. Provide matching lane splits."

    for i in "${!R1_FILES[@]}"; do
      r1="${R1_FILES[$i]}"
      r2="${R2_FILES[$i]}"
      uid="$(hash8 "${key}|${r1}|${r2}")"
      unit_tag="${SAMPLE}.K${kid}.U${uid}"

      # Stable links
      in_r1="$RAW/${unit_tag}.R1.fastq.gz"
      in_r2="$RAW/${unit_tag}.R2.fastq.gz"
      ln -sf "$r1" "$in_r1"
      ln -sf "$r2" "$in_r2"
      file_nonempty "$in_r1" || die "Missing/empty R1: $in_r1"
      file_nonempty "$in_r2" || die "Missing/empty R2: $in_r2"

      # Make raw chunks (adaptive)
      unit_make_raw_chunks_pe "$in_r1" "$in_r2" "$unit_tag"

      shopt -s nullglob
      R1_CHUNKS=( "$CHUNKS/${unit_tag}.R1_"*.fastq.gz )
      shopt -u nullglob
      [[ ${#R1_CHUNKS[@]} -gt 0 ]] || die "No raw chunks found for $unit_tag"

      for r1c in "${R1_CHUNKS[@]}"; do
        base=$(basename "$r1c")
        cid="${base#${unit_tag}.R1_}"
        cid="${cid%.fastq.gz}"
        r2c="$CHUNKS/${unit_tag}.R2_${cid}.fastq.gz"
        [[ -f "$r2c" ]] || die "Missing raw R2 chunk: $r2c"

        fastp_err="$TMP/fastp_${unit_tag}_${cid}.stderr.log"
        fp_json="$REPORTS/${unit_tag}.${cid}.fastp.json"

        if [[ "$LIBTYPE" == "ancient" ]]; then
          # aDNA PE: merge and map merged reads only
          merged_out="$TMP/${unit_tag}.${cid}.merged.fastq.gz"
          run_fastp "$fastp_err" "${FASTP_ARGS[@]}" "${FASTP_PE_EXTRA[@]}" \
            -i "$r1c" -I "$r2c" -m --merged_out "$merged_out" \
            -j "$fp_json"

          file_nonempty "$merged_out" || die "fastp produced empty merged output: $merged_out"
          cat "$merged_out" >> "$POOL_ANCIENT_MERGED_GZ"

          r1_after=$(json_get_int "$fp_json" "read1_after_filtering.total_reads")
          merged_reads=$(json_get_int "$fp_json" "merged_and_filtered.total_reads")

          TRIMMED_FRAGMENTS=$("$PYTHON" - <<PY
print(int($TRIMMED_FRAGMENTS) + int($r1_after))
PY
)
          MERGED_READS=$("$PYTHON" - <<PY
print(int($MERGED_READS) + int($merged_reads))
PY
)

          rm -f "$merged_out"

        else
          # modern PE: keep paired + rescue (merged+unpaired)
          trim1="$TMP/${unit_tag}.${cid}.R1.trim.fastq.gz"
          trim2="$TMP/${unit_tag}.${cid}.R2.trim.fastq.gz"
          merged="$TMP/${unit_tag}.${cid}.merged.fastq.gz"
          up1="$TMP/${unit_tag}.${cid}.unpaired1.fastq.gz"
          up2="$TMP/${unit_tag}.${cid}.unpaired2.fastq.gz"

          run_fastp "$fastp_err" "${FASTP_ARGS[@]}" "${FASTP_PE_EXTRA[@]}" \
            -i "$r1c" -I "$r2c" \
            --out1 "$trim1" --out2 "$trim2" \
            -m --merged_out "$merged" \
            --unpaired1 "$up1" --unpaired2 "$up2" \
            -j "$fp_json"

          file_nonempty "$trim1" || die "fastp produced empty trim1: $trim1"
          file_nonempty "$trim2" || die "fastp produced empty trim2: $trim2"

          cat "$trim1" >> "$POOL_PE_R1_GZ"
          cat "$trim2" >> "$POOL_PE_R2_GZ"

          # Rescue pool (only if present)
          [[ -s "$merged" ]] && cat "$merged" >> "$POOL_SE_RESCUE_GZ"
          [[ -s "$up1" ]] && cat "$up1" >> "$POOL_SE_RESCUE_GZ"
          [[ -s "$up2" ]] && cat "$up2" >> "$POOL_SE_RESCUE_GZ"

          r1_after=$(json_get_int "$fp_json" "read1_after_filtering.total_reads")
          r2_after=$(json_get_int "$fp_json" "read2_after_filtering.total_reads")
          merged_reads=$(json_get_int "$fp_json" "merged_and_filtered.total_reads")

          tp=$("$PYTHON" - <<PY
print(min(int($r1_after), int($r2_after)))
PY
)
          unp=$("$PYTHON" - <<PY
tp=min(int($r1_after), int($r2_after))
print(max(0, (int($r1_after)+int($r2_after)) - 2*tp))
PY
)

          TRIMMED_FRAGMENTS=$("$PYTHON" - <<PY
print(int($TRIMMED_FRAGMENTS) + int($tp))
PY
)
          UNPAIRED_READS=$("$PYTHON" - <<PY
print(int($UNPAIRED_READS) + int($unp))
PY
)
          MERGED_READS=$("$PYTHON" - <<PY
print(int($MERGED_READS) + int($merged_reads))
PY
)

          rm -f "$trim1" "$trim2" "$merged" "$up1" "$up2"
        fi
      done

      # Cleanup raw chunks for this unit (save disk)
      rm -f "$CHUNKS/${unit_tag}.R1_"*.fastq.gz "$CHUNKS/${unit_tag}.R2_"*.fastq.gz 2>/dev/null || true
    done

  else
    # SE: each R1 file is its own unit
    for r1 in "${R1_FILES[@]}"; do
      uid="$(hash8 "${key}|${r1}")"
      unit_tag="${SAMPLE}.K${kid}.U${uid}"

      in_r1="$RAW/${unit_tag}.R1.fastq.gz"
      ln -sf "$r1" "$in_r1"
      file_nonempty "$in_r1" || die "Missing/empty SE R1: $in_r1"

      unit_make_raw_chunks_se "$in_r1" "$unit_tag"

      shopt -s nullglob
      R1_CHUNKS=( "$CHUNKS/${unit_tag}.R1_"*.fastq.gz )
      shopt -u nullglob
      [[ ${#R1_CHUNKS[@]} -gt 0 ]] || die "No raw chunks found for $unit_tag"

      for r1c in "${R1_CHUNKS[@]}"; do
        base=$(basename "$r1c")
        cid="${base#${unit_tag}.R1_}"
        cid="${cid%.fastq.gz}"

        fastp_err="$TMP/fastp_${unit_tag}_${cid}.stderr.log"
        fp_json="$REPORTS/${unit_tag}.${cid}.fastp.json"
        trim_out="$TMP/${unit_tag}.${cid}.SE.trim.fastq.gz"

        run_fastp "$fastp_err" "${FASTP_ARGS[@]}" \
          -i "$r1c" --out1 "$trim_out" \
          -j "$fp_json"

        file_nonempty "$trim_out" || die "fastp produced empty SE trim: $trim_out"
        cat "$trim_out" >> "$POOL_SE_GZ"

        r1_after=$(json_get_int "$fp_json" "read1_after_filtering.total_reads")
        TRIMMED_FRAGMENTS=$("$PYTHON" - <<PY
print(int($TRIMMED_FRAGMENTS) + int($r1_after))
PY
)

        rm -f "$trim_out"
      done

      rm -f "$CHUNKS/${unit_tag}.R1_"*.fastq.gz 2>/dev/null || true
    done
  fi
done

###############################################################################
# TRIM-ONLY EXIT
###############################################################################
if [[ $TRIM_ONLY -eq 1 ]]; then
  log "Trim-only mode enabled; writing stats and exiting"
  PCT_REMAIN=$(pct "$TRIMMED_FRAGMENTS" "$RAW_FRAGMENTS")
  RATIO_TRIM=$(ratio "$TRIMMED_FRAGMENTS" "$RAW_FRAGMENTS")

  {
    printf "sample\tlibrary_type\tseq_mode\traw_R1_reads\traw_R2_reads\traw_fragments\ttrimmed_fragments\tmerged_reads\tunpaired_reads\tpct_fragments_remaining\ttrimmed_over_raw\n"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$SAMPLE" "$LIBTYPE" "$SEQ_MODE" \
      "$RAW_R1_READS" "$RAW_R2_READS" "$RAW_FRAGMENTS" \
      "$TRIMMED_FRAGMENTS" "$MERGED_READS" "$UNPAIRED_READS" \
      "$PCT_REMAIN" "$RATIO_TRIM"
  } > "$STATS"

  PIPE_T1=$(date +%s)
  log "STATS written: $STATS"
  log "Total wall time: $((PIPE_T1 - PIPE_T0)) seconds"
  exit 0
fi

###############################################################################
# MAPPING ONCE FROM POOLED TRIMMED OUTPUTS
###############################################################################
MERGED_BAM="$FINAL/${SAMPLE}.merged.bam"
DEDUP_BAM="$FINAL/${SAMPLE}.dedup.bam"
FINAL_BAM="$FINAL/${SAMPLE}.final.bam"
OUT_BAM="$OUT/${SAMPLE}.final.RG.bam"

TMP_BAMS=()

log "STEP: mapping pooled trimmed inputs (map once)"

map_se_modern() {
  local fq_gz="$1" outbam="$2"
  local cmd
  cmd=$(decomp_cmd "$fq_gz")
  log "  Mapping modern SE: $fq_gz"
  "$BWA" mem -t "$THREADS" "$REF" <(eval "$cmd") > "$TMP/${SAMPLE}.pooled.SE.sam"
  "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.pooled.SE.sam" > "$TMP/${SAMPLE}.pooled.SE.bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$outbam" "$TMP/${SAMPLE}.pooled.SE.bam"
  rm -f "$TMP/${SAMPLE}.pooled.SE.sam" "$TMP/${SAMPLE}.pooled.SE.bam"
}

map_pe_modern() {
  local r1_gz="$1" r2_gz="$2" outbam="$3"
  local c1 c2
  c1=$(decomp_cmd "$r1_gz")
  c2=$(decomp_cmd "$r2_gz")
  log "  Mapping modern PE: $r1_gz + $r2_gz"
  "$BWA" mem -t "$THREADS" "$REF" <(eval "$c1") <(eval "$c2") > "$TMP/${SAMPLE}.pooled.PE.sam"
  "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.pooled.PE.sam" > "$TMP/${SAMPLE}.pooled.PE.bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$outbam" "$TMP/${SAMPLE}.pooled.PE.bam"
  rm -f "$TMP/${SAMPLE}.pooled.PE.sam" "$TMP/${SAMPLE}.pooled.PE.bam"
}

map_se_ancient_aln() {
  local fq_gz="$1" outbam="$2" tag="$3"
  local cmd
  cmd=$(decomp_cmd "$fq_gz")
  log "  Mapping ancient SE ($tag): $fq_gz"
  "$BWA" aln -l 999 -n "$MISMATCH" -t "$THREADS" "$REF" <(eval "$cmd") > "$TMP/${SAMPLE}.${tag}.sai"
  "$BWA" samse "$REF" "$TMP/${SAMPLE}.${tag}.sai" <(eval "$cmd") > "$TMP/${SAMPLE}.${tag}.sam"
  "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.${tag}.sam" > "$TMP/${SAMPLE}.${tag}.bam.tmp"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$outbam" "$TMP/${SAMPLE}.${tag}.bam.tmp"
  rm -f "$TMP/${SAMPLE}.${tag}.sai" "$TMP/${SAMPLE}.${tag}.sam" "$TMP/${SAMPLE}.${tag}.bam.tmp"
}

# Only map if pool non-empty
if [[ "$LIBTYPE" == "ancient" ]]; then
  # aDNA:
  #  - SE pool maps as SE
  #  - ancient merged pool maps as SE (merged reads)
  if [[ -s "$POOL_SE_GZ" ]]; then
    out="$FINAL/${SAMPLE}.map.SE.bam"
    map_se_ancient_aln "$POOL_SE_GZ" "$out" "SE"
    TMP_BAMS+=( "$out" )
  fi
  if [[ -s "$POOL_ANCIENT_MERGED_GZ" ]]; then
    out="$FINAL/${SAMPLE}.map.MERGED.bam"
    map_se_ancient_aln "$POOL_ANCIENT_MERGED_GZ" "$out" "MERGED"
    TMP_BAMS+=( "$out" )
  fi
else
  # modern:
  #  - PE pool maps as PE
  #  - SE pool + rescue pool map as SE and are merged with PE
  if [[ -s "$POOL_PE_R1_GZ" && -s "$POOL_PE_R2_GZ" ]]; then
    out="$FINAL/${SAMPLE}.map.PE.bam"
    map_pe_modern "$POOL_PE_R1_GZ" "$POOL_PE_R2_GZ" "$out"
    TMP_BAMS+=( "$out" )
  fi

  # Combine SE sources into one mapping (to avoid multiple BWAs)
  POOL_ALL_SE_GZ="$POOL/${SAMPLE}.pool.SE.all.fastq.gz"
  : > "$POOL_ALL_SE_GZ"
  [[ -s "$POOL_SE_GZ" ]] && cat "$POOL_SE_GZ" >> "$POOL_ALL_SE_GZ"
  [[ -s "$POOL_SE_RESCUE_GZ" ]] && cat "$POOL_SE_RESCUE_GZ" >> "$POOL_ALL_SE_GZ"

  if [[ -s "$POOL_ALL_SE_GZ" ]]; then
    out="$FINAL/${SAMPLE}.map.SE.bam"
    map_se_modern "$POOL_ALL_SE_GZ" "$out"
    TMP_BAMS+=( "$out" )
  fi
fi

[[ ${#TMP_BAMS[@]} -gt 0 ]] || die "No pooled mapping inputs were non-empty; nothing to map."

log "STEP: merge mapped BAMs"
"$SAMTOOLS" merge -@ "$SORT_THREADS" "$MERGED_BAM" "${TMP_BAMS[@]}"
file_nonempty "$MERGED_BAM" || die "Merge failed: $MERGED_BAM"

###############################################################################
# DEDUP + RG
###############################################################################
log "STEP: dedup"
if [[ "$LIBTYPE" == "modern" && $N_PE_KEYS -gt 0 ]]; then
  "$SAMTOOLS" sort -n -@ "$SORT_THREADS" -o "$FINAL/${SAMPLE}.namesort.bam" "$MERGED_BAM"
  "$SAMTOOLS" fixmate -m -@ "$SORT_THREADS" "$FINAL/${SAMPLE}.namesort.bam" "$FINAL/${SAMPLE}.fixmate.bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$FINAL/${SAMPLE}.coordsort.bam" "$FINAL/${SAMPLE}.fixmate.bam"
  "$SAMTOOLS" markdup -r -@ "$SORT_THREADS" "$FINAL/${SAMPLE}.coordsort.bam" "$DEDUP_BAM"
  rm -f "$FINAL/${SAMPLE}.namesort.bam" "$FINAL/${SAMPLE}.fixmate.bam" "$FINAL/${SAMPLE}.coordsort.bam"
else
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$FINAL/${SAMPLE}.coordsort.bam" "$MERGED_BAM"
  "$SAMTOOLS" markdup -r -s -@ "$SORT_THREADS" "$FINAL/${SAMPLE}.coordsort.bam" "$DEDUP_BAM"
  rm -f "$FINAL/${SAMPLE}.coordsort.bam"
fi
file_nonempty "$DEDUP_BAM" || die "Dedup failed: $DEDUP_BAM"

log "STEP: final sort+index"
"$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$FINAL_BAM" "$DEDUP_BAM"
"$SAMTOOLS" index "$FINAL_BAM"
file_nonempty "$FINAL_BAM" || die "Final BAM missing/empty"

log "STEP: read groups"
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina\tLB:${SAMPLE}"
"$SAMTOOLS" addreplacerg -r "$RG" -o "$OUT_BAM" "$FINAL_BAM"
"$SAMTOOLS" index "$OUT_BAM"
file_nonempty "$OUT_BAM" || die "RG BAM missing/empty"

if [[ "$LIBTYPE" == "ancient" ]]; then
  log "STEP: mapDamage"
  "$MAPDAMAGE" -i "$OUT_BAM" --merge-reference-sequences --no-stats -r "$REF" -d "$OUT/${SAMPLE}_mapdamage"
fi

###############################################################################
# COVERAGE/DEPTH + SUMMARY STATS (FRAGMENT-SPACE)
###############################################################################
log "STEP: coverage/depth statistics"
"$SAMTOOLS" coverage -q "$MAPQ" "$OUT_BAM" > "$COV_TSV"

AVG_DEPTH=$(
  awk 'BEGIN{sum=0;len=0}
       NR==1{next}
       {L=($3-$2+1); sum += $7*L; len += L}
       END{if(len>0) printf "%.6f", sum/len; else print 0}' \
    "$COV_TSV"
)

PCT_COVERED=$(
  awk 'BEGIN{cov=0;len=0}
       NR==1{next}
       {L=($3-$2+1); cov += $5; len += L}
       END{if(len>0) printf "%.6f", 100.0*cov/len; else print 0}' \
    "$COV_TSV"
)

log "STEP: summary statistics"
EXCL=2308  # 4(unmapped)+256(secondary)+2048(supplementary)

MAPPED_READS_ALL=$("$SAMTOOLS" view -c -F "$EXCL" "$MERGED_BAM")
MAPPED_READS_UNIQUE=$("$SAMTOOLS" view -c -F "$EXCL" "$DEDUP_BAM")

PE_FRAG_ALL=$("$SAMTOOLS" view -c -f 1 -f 64 -F "$EXCL" "$MERGED_BAM")
SE_FRAG_ALL=$("$SAMTOOLS" view -c -F "$EXCL" -F 1 "$MERGED_BAM")
MAPPED_FRAGMENTS_ALL=$("$PYTHON" - <<PY
print(int($PE_FRAG_ALL) + int($SE_FRAG_ALL))
PY
)

PE_FRAG_UNIQUE=$("$SAMTOOLS" view -c -f 1 -f 64 -F "$EXCL" "$DEDUP_BAM")
SE_FRAG_UNIQUE=$("$SAMTOOLS" view -c -F "$EXCL" -F 1 "$DEDUP_BAM")
MAPPED_FRAGMENTS_UNIQUE=$("$PYTHON" - <<PY
print(int($PE_FRAG_UNIQUE) + int($SE_FRAG_UNIQUE))
PY
)

AVG_READLEN=$(
  "$SAMTOOLS" view -F "$EXCL" "$DEDUP_BAM" \
  | awk '{l=length($10); if(l>0){sum+=l;n++}} END{ if(n>0) printf "%.2f", sum/n; else printf "0.00"}'
)

MAPPED_BP=$(awk -v n="$MAPPED_READS_ALL" -v l="$AVG_READLEN" 'BEGIN{printf "%.0f", n*l}')

PCT_REMAIN=$(pct "$TRIMMED_FRAGMENTS" "$RAW_FRAGMENTS")
RATIO_TRIM=$(ratio "$TRIMMED_FRAGMENTS" "$RAW_FRAGMENTS")

ENDOG=$(ratio "$MAPPED_FRAGMENTS_UNIQUE" "$RAW_FRAGMENTS")

DUPRATE=$(
  echo "1 - $(ratio "$MAPPED_FRAGMENTS_UNIQUE" "$MAPPED_FRAGMENTS_ALL")" \
  | bc -l \
  | awk '{printf "%.6f", $0}'
)

{
  printf "sample\tlibrary_type\tseq_mode\traw_R1_reads\traw_R2_reads\traw_fragments\ttrimmed_fragments\tmerged_reads\tunpaired_reads\tpct_fragments_remaining\ttrimmed_over_raw\tmapped_reads_all\tmapped_reads_unique\tmapped_fragments_all\tmapped_fragments_unique\tendog\tduprate\tavg_readlen\tmapped_bp\tavg_depth\tpct_covered\n"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$SAMPLE" "$LIBTYPE" "$SEQ_MODE" \
    "$RAW_R1_READS" "$RAW_R2_READS" "$RAW_FRAGMENTS" \
    "$TRIMMED_FRAGMENTS" "$MERGED_READS" "$UNPAIRED_READS" \
    "$PCT_REMAIN" "$RATIO_TRIM" \
    "$MAPPED_READS_ALL" "$MAPPED_READS_UNIQUE" \
    "$MAPPED_FRAGMENTS_ALL" "$MAPPED_FRAGMENTS_UNIQUE" \
    "$ENDOG" "$DUPRATE" \
    "$AVG_READLEN" "$MAPPED_BP" \
    "$AVG_DEPTH" "$PCT_COVERED"
} > "$STATS"

###############################################################################
# CLEANUP
###############################################################################
if [[ $KEEP_INTERMEDIATE -eq 0 ]]; then
  log "Cleaning intermediates (work dir only)"
  rm -rf "$WORK"
  [[ -n "$TMPDIR_USER" ]] && rm -rf "$TMPDIR_USER/plainmap_${SAMPLE}" || true
fi

PIPE_T1=$(date +%s)
log "STATUS: SUCCESS"
log "Total wall time: $((PIPE_T1 - PIPE_T0)) seconds"
log "Final BAM: $OUT_BAM"
log "Coverage table: $COV_TSV"
log "Stats table: $STATS"
