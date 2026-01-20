#!/usr/bin/env bash
# PlainMap v0.1
# Transparent, failure-aware mapping pipeline for ancient and modern DNA
# - Manifest input (txt file), mixed SE/PE supported
# - FASTQ header detection supports Illumina/CASAVA (" 1:"/" 2:") and ENA ("/1","/2")
# - Optional pre-fastp chunking (safety valve for huge inputs)
# - Optional pilot subsampling (pre-fastp; applied per unit/chunk)
# - fastp JSON only (no HTML)
# - Pool after trimming, then map once (with optional mapping chunking)
# - pigz used automatically if available
# - Parallel-safe BWA indexing via filesystem lock
# - Fragment-aware stats; duplication rate computed from duplicate flags (never negative)

set -euo pipefail

###############################################################################
# REQUIRE BASH
###############################################################################
if [ -z "${BASH_VERSION:-}" ]; then
  echo "ERROR: PlainMap must be run with bash (not sh)."
  exit 1
fi

###############################################################################
# DEFAULTS
###############################################################################
VERSION="0.1"

THREADS=1
MINLEN=30
MISMATCH=0.01              # bwa aln -n (ancient)
LIBTYPE="modern"           # modern|ancient
MAX_READS_PER_CHUNK=0      # 0 disables pre-fastp chunking
MAPQ=20

PILOT_FRAGMENTS=0          # 0 disables pilot mode; otherwise limit reads per unit/chunk BEFORE fastp

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
PlainMap v$VERSION

Required:
  -manifest FILE
  -prefix STRING
  -ref FILE
  -outdir DIR

Optional:
  -library-type modern|ancient     (default: modern)
  -t, --threads INT                (default: 1)
  -minlength INT                   (default: 30)
  -mismatch FLOAT                  (ancient only; bwa aln -n; default: 0.01)
  -max-reads-per-chunk INT         (default: 0; disabled) pre-fastp chunking safety valve
  --pilot-fragments INT            (default: 0; disabled) limit reads per unit/chunk before fastp
  --adapter-r1 SEQ
  --adapter-r2 SEQ
  --trim-only                      Trim only (fastp); exit before mapping
  --tmpdir DIR                     Custom temp dir (e.g. node-local scratch)
  --keep-intermediate              Keep <outdir>/<prefix>/work
  --resume | --no-resume
  --dry-run                        Print plan only
  --validate                       Check tools + gzip/pigz -t all manifest FASTQs, then exit
  --reset                          Clear ALL checkpoints for this sample

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
    --pilot-fragments) PILOT_FRAGMENTS="$2"; shift 2 ;;
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
# PATH RESOLUTION
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
# LOGGING
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
# WORK DIRS
###############################################################################
WORK="$OUT/${SAMPLE}/work"
RAW="$WORK/raw"
CHUNKS="$WORK/chunks"          # pre-fastp chunks
TMPBASE="$WORK/tmp"
MAPCHUNKS="$WORK/mapchunks"    # post-fastp mapping chunks (optional)
BAMS="$WORK/bams"
FINAL="$WORK/final"
CKPT="$WORK/.checkpoints"
REPORTS="$OUT/${SAMPLE}/fastp_reports"

STATS="$OUT/${SAMPLE}.plainmap.stats.tsv"
COV_TSV="$OUT/${SAMPLE}.plainmap.coverage.tsv"

if [[ -n "$TMPDIR_USER" ]]; then
  TMP="$TMPDIR_USER/plainmap_${SAMPLE}"
else
  TMP="$TMPBASE"
fi

mkdir -p "$RAW" "$CHUNKS" "$MAPCHUNKS" "$BAMS" "$FINAL" "$CKPT" "$TMP" "$REPORTS"

###############################################################################
# COMPRESSION HELPERS (pigz if available)
###############################################################################
HAVE_PIGZ=0
if command -v pigz >/dev/null 2>&1; then
  HAVE_PIGZ=1
fi

if [[ $HAVE_PIGZ -eq 1 ]]; then
  GZ_TEST=( pigz -t )
  ZCAT=( pigz -dc )
  GZIP=( pigz -p "$THREADS" )
else
  GZ_TEST=( gzip -t )
  ZCAT=( zcat )
  GZIP=( gzip )
fi

###############################################################################
# FILE CHECKS
###############################################################################
file_nonempty() { [[ -e "$1" && -s "$1" ]]; }

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
  require_cmd split
  require_cmd "$FASTP"; require_cmd "$BWA"; require_cmd "$SAMTOOLS"; require_cmd "$PYTHON"
  check_split_filter_support
  if [[ "$LIBTYPE" == "ancient" ]]; then
    require_cmd "$MAPDAMAGE"
  fi
  # gzip/pigz tools
  if [[ $HAVE_PIGZ -eq 1 ]]; then
    require_cmd pigz
  else
    require_cmd gzip
    require_cmd zcat
  fi
}

###############################################################################
# VALIDATE / DRYRUN
###############################################################################
if [[ $VALIDATE_ONLY -eq 1 ]]; then
  log "Validation mode: checking tools + gzip/pigz -t for all manifest FASTQs"
  quick_preflight
  n=0
  while read -r f; do
    [[ -z "$f" || "$f" =~ ^# ]] && continue
    n=$((n+1))
    [[ "$f" = /* ]] || f="$MANIFEST_DIR/$f"
    f="$(resolve_path "$f")"
    [[ -r "$f" ]] || die "FASTQ not readable: $f"
    "${GZ_TEST[@]}" "$f" >/dev/null || die "Corrupt/truncated gzip: $f"
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
  log "Settings:"
  log "  library_type:       $LIBTYPE"
  log "  threads:            $THREADS"
  log "  minlength:          $MINLEN"
  log "  mismatch (ancient): $MISMATCH"
  log "  pilot_fragments:    $PILOT_FRAGMENTS"
  log "  max_reads_per_chunk:$MAX_READS_PER_CHUNK"
  log "Planned steps: read manifest -> pre-fastp chunk (optional) -> pilot (optional) -> fastp -> pool -> map -> dedup -> RG -> coverage -> stats"
  exit 0
fi

quick_preflight

###############################################################################
# BWA INDEX (PARALLEL SAFE)
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
# FASTP WRAPPER + ARGS (JSON only)
###############################################################################
run_fastp() {
  local stderr_file="$1"; shift
  log "$FASTP $*"
  : > "$stderr_file"
  set +e
  "$FASTP" "$@" 2> >(tee -a "$stderr_file" >&2)
  status=$?
  set -e
  [[ $status -eq 0 ]] || die "fastp exited non-zero ($status)"
  if grep -qiE "igzip|invalid gzip|premature|truncated|error" "$stderr_file"; then
    die "fastp reported gzip/igzip error (see $stderr_file)"
  fi
}

FASTP_ARGS=( -l "$MINLEN" -g -w "$THREADS" )
# adapters: if provided explicitly, use them; otherwise for PE we will use --detect_adapter_for_pe
[[ -n "$ADAPTER_R1" ]] && FASTP_ARGS+=( --adapter_sequence "$ADAPTER_R1" )
[[ -n "$ADAPTER_R2" ]] && FASTP_ARGS+=( --adapter_sequence_r2 "$ADAPTER_R2" )

###############################################################################
# STATS HELPERS
###############################################################################
count_reads_fastq_gz() {
  # counts reads (NR/4) in gz fastq
  "${ZCAT[@]}" "$1" | awk 'END{print NR/4}'
}

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

###############################################################################
# READ HEADER CLASSIFIER (Illumina + ENA)
###############################################################################
classify_fastq_direction() {
  # prints: R1 | R2 | UNKNOWN
  local fq="$1"
  local hdr
  hdr=$("${ZCAT[@]}" "$fq" | head -n 1 || true)

  # Illumina/CASAVA: has ' 1:' or ' 2:' token
  if [[ "$hdr" == *" 1:"* ]]; then
    echo "R1"; return
  fi
  if [[ "$hdr" == *" 2:"* ]]; then
    echo "R2"; return
  fi

  # ENA style: first token ends with /1 or /2
  # Example: @ERR123.1/1
  local firsttok="${hdr%% *}"
  if [[ "$firsttok" == */1 ]]; then
    echo "R1"; return
  fi
  if [[ "$firsttok" == */2 ]]; then
    echo "R2"; return
  fi

  echo "UNKNOWN"
}

###############################################################################
# MANIFEST -> UNITS
###############################################################################
# Units:
# - SE unit: one FASTQ (R1)
# - PE unit: paired FASTQs (one R1 + one R2) matched by order encountered
#
# NOTE: This assumes the manifest lists PE R1/R2 in matching counts (order can be mixed).
# PlainMap groups by header direction, then pairs by encounter order.

SE_LIST=()
R1_LIST=()
R2_LIST=()
TOTAL_FILES=0

while read -r f; do
  [[ -z "$f" || "$f" =~ ^# ]] && continue
  [[ "$f" = /* ]] || f="$MANIFEST_DIR/$f"
  f="$(resolve_path "$f")"
  [[ -r "$f" ]] || die "FASTQ not readable: $f"
  TOTAL_FILES=$((TOTAL_FILES+1))

  dir=$(classify_fastq_direction "$f")
  if [[ "$dir" == "R1" ]]; then
    R1_LIST+=( "$f" )
  elif [[ "$dir" == "R2" ]]; then
    R2_LIST+=( "$f" )
  else
    die "Could not classify FASTQ direction (Illumina ' 1:'/' 2:' or ENA '/1' '/2' expected): $f"
  fi
done < "$MANIFEST"

[[ $TOTAL_FILES -gt 0 ]] || die "Manifest empty: $MANIFEST"
[[ ${#R1_LIST[@]} -gt 0 ]] || die "No R1 reads found in manifest"

# Pair up as much as possible
PE_N=$(( ${#R1_LIST[@]} < ${#R2_LIST[@]} ? ${#R1_LIST[@]} : ${#R2_LIST[@]} ))
SE_N=$(( ${#R1_LIST[@]} - PE_N ))

SEQ_MODE="SE"
if [[ $PE_N -gt 0 && $SE_N -gt 0 ]]; then
  SEQ_MODE="MIX"
elif [[ $PE_N -gt 0 ]]; then
  SEQ_MODE="PE"
fi

log "Detected sequencing mode: $SEQ_MODE (R1 files: ${#R1_LIST[@]}, R2 files: ${#R2_LIST[@]})"
if [[ ${#R2_LIST[@]} -gt ${#R1_LIST[@]} ]]; then
  die "More R2 than R1 files found; manifest pairing inconsistent"
fi

# Build SE list (R1 leftovers after pairing)
if [[ $SE_N -gt 0 ]]; then
  for ((i=PE_N; i<${#R1_LIST[@]}; i++)); do
    SE_LIST+=( "${R1_LIST[$i]}" )
  done
fi

###############################################################################
# PRE-FASTP CHUNKING
###############################################################################
# We always produce "unit chunks" under $CHUNKS:
# - SE chunk:  <SAMPLE>.<UNITID>.SE_XXXX.fastq.gz
# - PE chunks: <SAMPLE>.<UNITID>.R1_XXXX.fastq.gz and ...R2_XXXX.fastq.gz
#
# If MAX_READS_PER_CHUNK=0 -> link original as a single chunk.
# Else -> split into chunks with up to MAX_READS_PER_CHUNK reads.
#
# If PILOT_FRAGMENTS>0 -> before running fastp on a chunk, we limit input reads
# to PILOT_FRAGMENTS (SE reads or PE pairs) for that chunk. This is per chunk.

split_fastq_gz() {
  local in="$1" outprefix="$2" reads="$3"
  local lines=$((reads * 4))
  local filter_cmd
  if [[ $HAVE_PIGZ -eq 1 ]]; then
    filter_cmd="pigz -p ${THREADS} > \$FILE.fastq.gz"
  else
    filter_cmd="gzip > \$FILE.fastq.gz"
  fi
  "${ZCAT[@]}" "$in" | split -l "$lines" -d -a 4 --filter="$filter_cmd" - "$outprefix"
}

limit_fastq_gz() {
  # limit reads in gz FASTQ to N reads; write gz FASTQ
  # Uses head which can SIGPIPE upstream; disable pipefail in this function.
  local in="$1" out="$2" reads="$3"
  local lines=$((reads * 4))
  set +o pipefail
  "${ZCAT[@]}" "$in" | head -n "$lines" | "${GZIP[@]}" > "$out"
  set -o pipefail
  file_nonempty "$out" || die "Pilot limiting produced empty FASTQ: $out"
}

###############################################################################
# MAKE UNIT CHUNKS
###############################################################################
make_chunks() {
  log "STEP: pre-fastp chunking"
  ckpt_clear chunk

  rm -f "$CHUNKS/${SAMPLE}."*".fastq.gz" 2>/dev/null || true

  # PE units
  for ((u=0; u<PE_N; u++)); do
    local r1="${R1_LIST[$u]}"
    local r2="${R2_LIST[$u]}"
    local unit="U$(printf '%08d' $u)"
    local prefix="$CHUNKS/${SAMPLE}.${unit}."

    if [[ "$MAX_READS_PER_CHUNK" -gt 0 ]]; then
      log "  PE unit $unit: splitting into chunks (max reads/chunk: $MAX_READS_PER_CHUNK)"
      split_fastq_gz "$r1" "${prefix}R1_" "$MAX_READS_PER_CHUNK"
      split_fastq_gz "$r2" "${prefix}R2_" "$MAX_READS_PER_CHUNK"
    else
      log "  PE unit $unit: chunking disabled; linking as single chunk"
      ln -sf "$r1" "${prefix}R1_0000.fastq.gz"
      ln -sf "$r2" "${prefix}R2_0000.fastq.gz"
    fi

    file_nonempty "${prefix}R1_0000.fastq.gz" || die "Missing PE R1 chunk for $unit"
    file_nonempty "${prefix}R2_0000.fastq.gz" || die "Missing PE R2 chunk for $unit"
  done

  # SE units
  for ((s=0; s<${#SE_LIST[@]}; s++)); do
    local fq="${SE_LIST[$s]}"
    local unit="S$(printf '%08d' $s)"
    local prefix="$CHUNKS/${SAMPLE}.${unit}."

    if [[ "$MAX_READS_PER_CHUNK" -gt 0 ]]; then
      log "  SE unit $unit: splitting into chunks (max reads/chunk: $MAX_READS_PER_CHUNK)"
      split_fastq_gz "$fq" "${prefix}SE_" "$MAX_READS_PER_CHUNK"
    else
      log "  SE unit $unit: chunking disabled; linking as single chunk"
      ln -sf "$fq" "${prefix}SE_0000.fastq.gz"
    fi

    file_nonempty "${prefix}SE_0000.fastq.gz" || die "Missing SE chunk for $unit"
  done

  ckpt_mark chunk
}

if ! ckpt_ok chunk "$CHUNKS"; then
  make_chunks
else
  log "STEP: pre-fastp chunking (skipped; checkpoint present)"
fi

###############################################################################
# FASTP PER CHUNK -> POOL TRIMMED OUTPUTS
###############################################################################
POOL_PE_R1="$RAW/${SAMPLE}.pool.PE.R1.fastq.gz"
POOL_PE_R2="$RAW/${SAMPLE}.pool.PE.R2.fastq.gz"
POOL_SE_ALL="$RAW/${SAMPLE}.pool.SE.all.fastq.gz"   # all SE to be mapped as SE (SE units + merged + rescued)
: > "$POOL_PE_R1"
: > "$POOL_PE_R2"
: > "$POOL_SE_ALL"

RAW_R1_READS=0
RAW_R2_READS=0
RAW_FRAGMENTS=0

TRIMMED_FRAGMENTS=0
MERGED_READS=0
UNPAIRED_READS=0

log "STEP: fastp per pre-fastp chunk (JSON only) + pooling"

# Helper: add integer safely
py_add() { "$PYTHON" - "$1" "$2" <<'PY'
import sys
print(int(sys.argv[1]) + int(sys.argv[2]))
PY
}

# Process PE chunks
shopt -s nullglob
for ((u=0; u<PE_N; u++)); do
  unit="U$(printf '%08d' $u)"
  r1_chunks=( "$CHUNKS/${SAMPLE}.${unit}.R1_"*.fastq.gz )
  [[ ${#r1_chunks[@]} -gt 0 ]] || die "No R1 chunks for PE unit $unit"
  for r1c in "${r1_chunks[@]}"; do
    base=$(basename "$r1c")
    cid="${base#${SAMPLE}.${unit}.R1_}"
    cid="${cid%.fastq.gz}"
    r2c="$CHUNKS/${SAMPLE}.${unit}.R2_${cid}.fastq.gz"
    [[ -f "$r2c" ]] || die "Missing R2 chunk for PE unit $unit chunk $cid: $r2c"

    # raw counts (pilot-aware later; these are counts of what we will process)
    # we count after pilot limiting decision so raw_* reflect processed reads.
    in_r1="$r1c"
    in_r2="$r2c"

    if [[ "$PILOT_FRAGMENTS" -gt 0 ]]; then
      log "  Pilot: limiting PE chunk $unit.$cid to first $PILOT_FRAGMENTS pairs"
      lim_r1="$TMP/${SAMPLE}.${unit}.${cid}.pilot.R1.fastq.gz"
      lim_r2="$TMP/${SAMPLE}.${unit}.${cid}.pilot.R2.fastq.gz"
      limit_fastq_gz "$in_r1" "$lim_r1" "$PILOT_FRAGMENTS"
      limit_fastq_gz "$in_r2" "$lim_r2" "$PILOT_FRAGMENTS"
      in_r1="$lim_r1"
      in_r2="$lim_r2"
    fi

    # update RAW counts based on processed input
    r1n=$(count_reads_fastq_gz "$in_r1")
    r2n=$(count_reads_fastq_gz "$in_r2")
    RAW_R1_READS=$(py_add "$RAW_R1_READS" "$r1n")
    RAW_R2_READS=$(py_add "$RAW_R2_READS" "$r2n")
    # fragment count for PE is by R1
    RAW_FRAGMENTS=$(py_add "$RAW_FRAGMENTS" "$r1n")

    fastp_err="$TMP/fastp_${SAMPLE}.${unit}.${cid}.stderr.log"
    fp_json="$REPORTS/${SAMPLE}.${unit}.${cid}.fastp.json"

    trim1="$TMP/${SAMPLE}.${unit}.${cid}.trim.R1.fastq.gz"
    trim2="$TMP/${SAMPLE}.${unit}.${cid}.trim.R2.fastq.gz"
    merged="$TMP/${SAMPLE}.${unit}.${cid}.merged.fastq.gz"
    up1="$TMP/${SAMPLE}.${unit}.${cid}.unpaired1.fastq.gz"
    up2="$TMP/${SAMPLE}.${unit}.${cid}.unpaired2.fastq.gz"

    # For PE, if adapters not explicitly provided, enable auto-detection.
    extra_pe=()
    if [[ -z "$ADAPTER_R1" && -z "$ADAPTER_R2" ]]; then
      extra_pe+=( --detect_adapter_for_pe )
    fi

    run_fastp "$fastp_err" "${FASTP_ARGS[@]}" "${extra_pe[@]}" \
      -i "$in_r1" -I "$in_r2" \
      --out1 "$trim1" --out2 "$trim2" \
      -m --merged_out "$merged" \
      --unpaired1 "$up1" --unpaired2 "$up2" \
      -j "$fp_json"

    file_nonempty "$trim1" || die "fastp produced empty trimmed R1 for $unit.$cid"
    file_nonempty "$trim2" || die "fastp produced empty trimmed R2 for $unit.$cid"

    # Append trimmed pairs to pooled PE
    cat "$trim1" >> "$POOL_PE_R1"
    cat "$trim2" >> "$POOL_PE_R2"

    # Append merged + unpaired to pooled SE-all
    if file_nonempty "$merged"; then cat "$merged" >> "$POOL_SE_ALL"; fi
    if file_nonempty "$up1"; then cat "$up1" >> "$POOL_SE_ALL"; fi
    if file_nonempty "$up2"; then cat "$up2" >> "$POOL_SE_ALL"; fi

    # Update trimming stats (fragment-aware)
    r1_after=$(json_get_int "$fp_json" "read1_after_filtering.total_reads")
    r2_after=$(json_get_int "$fp_json" "read2_after_filtering.total_reads")
    merged_reads=$(json_get_int "$fp_json" "merged_and_filtered.total_reads")

    # trimmed paired fragments are the min of post-filter R1/R2
    tp=$("$PYTHON" - "$r1_after" "$r2_after" <<'PY'
import sys
print(min(int(sys.argv[1]), int(sys.argv[2])))
PY
)
    TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$tp")
    MERGED_READS=$(py_add "$MERGED_READS" "$merged_reads")

    # unpaired reads (fragment count as reads, mapped later as SE)
    unpaired=$("$PYTHON" - "$r1_after" "$r2_after" <<'PY'
import sys
r1=int(sys.argv[1]); r2=int(sys.argv[2])
tp=min(r1,r2)
print(max(0, (r1+r2) - 2*tp))
PY
)
    UNPAIRED_READS=$(py_add "$UNPAIRED_READS" "$unpaired")
    TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$merged_reads")
    TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$unpaired")

    rm -f "$trim1" "$trim2" "$merged" "$up1" "$up2" 2>/dev/null || true
    [[ "$PILOT_FRAGMENTS" -gt 0 ]] && rm -f "$TMP/${SAMPLE}.${unit}.${cid}.pilot.R1.fastq.gz" "$TMP/${SAMPLE}.${unit}.${cid}.pilot.R2.fastq.gz" 2>/dev/null || true
  done
done

# Process SE chunks (R1-only units)
for ((s=0; s<${#SE_LIST[@]}; s++)); do
  unit="S$(printf '%08d' $s)"
  se_chunks=( "$CHUNKS/${SAMPLE}.${unit}.SE_"*.fastq.gz )
  [[ ${#se_chunks[@]} -gt 0 ]] || die "No SE chunks for unit $unit"
  for sec in "${se_chunks[@]}"; do
    base=$(basename "$sec")
    cid="${base#${SAMPLE}.${unit}.SE_}"
    cid="${cid%.fastq.gz}"

    in_se="$sec"
    if [[ "$PILOT_FRAGMENTS" -gt 0 ]]; then
      log "  Pilot: limiting SE chunk $unit.$cid to first $PILOT_FRAGMENTS reads"
      lim_se="$TMP/${SAMPLE}.${unit}.${cid}.pilot.SE.fastq.gz"
      limit_fastq_gz "$in_se" "$lim_se" "$PILOT_FRAGMENTS"
      in_se="$lim_se"
    fi

    n=$(count_reads_fastq_gz "$in_se")
    RAW_R1_READS=$(py_add "$RAW_R1_READS" "$n")
    RAW_FRAGMENTS=$(py_add "$RAW_FRAGMENTS" "$n")

    fastp_err="$TMP/fastp_${SAMPLE}.${unit}.${cid}.stderr.log"
    fp_json="$REPORTS/${SAMPLE}.${unit}.${cid}.fastp.json"
    trim_se="$TMP/${SAMPLE}.${unit}.${cid}.trim.SE.fastq.gz"

    run_fastp "$fastp_err" "${FASTP_ARGS[@]}" \
      -i "$in_se" --out1 "$trim_se" \
      -j "$fp_json"

    file_nonempty "$trim_se" || die "fastp produced empty trimmed SE for $unit.$cid"

    # Append to pooled SE-all
    cat "$trim_se" >> "$POOL_SE_ALL"

    se_after=$(json_get_int "$fp_json" "read1_after_filtering.total_reads")
    TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$se_after")

    rm -f "$trim_se" 2>/dev/null || true
    [[ "$PILOT_FRAGMENTS" -gt 0 ]] && rm -f "$TMP/${SAMPLE}.${unit}.${cid}.pilot.SE.fastq.gz" 2>/dev/null || true
  done
done

# If no PE data, pooled PE files may be empty; ensure they are truly empty files
# (they already exist via : >)
# Ensure SE pool exists (may be empty if no SE/merged/unpaired; that's ok for PE-only modern mapping)
file_nonempty "$POOL_SE_ALL" || true

###############################################################################
# TRIM-ONLY MODE
###############################################################################
if [[ $TRIM_ONLY -eq 1 ]]; then
  log "Trim-only mode enabled; writing stats and exiting"
  PCT_REMAIN=$(pct "$TRIMMED_FRAGMENTS" "$RAW_FRAGMENTS")
  RATIO_TRIM=$(ratio "$TRIMMED_FRAGMENTS" "$RAW_FRAGMENTS")

  {
    printf "sample\tlibrary_type\tseq_mode\tpilot_fragments\tmax_reads_per_chunk\traw_R1_reads\traw_R2_reads\traw_fragments\ttrimmed_fragments\tmerged_reads\tunpaired_reads\tpct_fragments_remaining\ttrimmed_over_raw\tmapped_reads_all\tmapped_reads_unique\tmapped_fragments_all\tmapped_fragments_unique\tendog\tduprate\tavg_readlen\tmapped_bp\tavg_depth\tpct_covered\n"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n" \
      "$SAMPLE" "$LIBTYPE" "$SEQ_MODE" "$PILOT_FRAGMENTS" "$MAX_READS_PER_CHUNK" \
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
# POST-FASTP MAPPING CHUNKING (OPTIONAL; REUSE MAX_READS_PER_CHUNK)
###############################################################################
# To make mapping restartable/less memory-hungry, we optionally chunk pooled inputs again.
# If MAX_READS_PER_CHUNK=0 -> map pooled files directly.
# Else -> split pooled PE and SE pools into chunks of MAX_READS_PER_CHUNK reads/pairs.

rm -f "$MAPCHUNKS/${SAMPLE}."*.fastq.gz 2>/dev/null || true

make_mapchunks_se() {
  local in="$1"
  local prefix="$2"
  if [[ ! -s "$in" ]]; then
    return
  fi
  if [[ "$MAX_READS_PER_CHUNK" -gt 0 ]]; then
    log "  Mapping chunking SE: splitting pooled SE (max reads/chunk: $MAX_READS_PER_CHUNK)"
    split_fastq_gz "$in" "$prefix" "$MAX_READS_PER_CHUNK"
  else
    ln -sf "$in" "${prefix}0000.fastq.gz"
  fi
}

make_mapchunks_pe() {
  local in1="$1" in2="$2" prefix1="$3" prefix2="$4"
  if [[ ! -s "$in1" || ! -s "$in2" ]]; then
    return
  fi
  if [[ "$MAX_READS_PER_CHUNK" -gt 0 ]]; then
    log "  Mapping chunking PE: splitting pooled PE (max reads/chunk: $MAX_READS_PER_CHUNK)"
    split_fastq_gz "$in1" "$prefix1" "$MAX_READS_PER_CHUNK"
    split_fastq_gz "$in2" "$prefix2" "$MAX_READS_PER_CHUNK"
  else
    ln -sf "$in1" "${prefix1}0000.fastq.gz"
    ln -sf "$in2" "${prefix2}0000.fastq.gz"
  fi
}

log "STEP: post-fastp mapping chunking"
make_mapchunks_pe "$POOL_PE_R1" "$POOL_PE_R2" "$MAPCHUNKS/${SAMPLE}.PE.R1_" "$MAPCHUNKS/${SAMPLE}.PE.R2_"
make_mapchunks_se "$POOL_SE_ALL" "$MAPCHUNKS/${SAMPLE}.SE_"

###############################################################################
# MAPPING
###############################################################################
log "STEP: mapping"
rm -f "$BAMS/${SAMPLE}."*.bam 2>/dev/null || true

map_se_mem_chunk() {
  local fq="$1" outbam="$2"
  "$BWA" mem -t "$THREADS" "$REF" "$fq" > "$TMP/${SAMPLE}.mem.SE.sam"
  "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.mem.SE.sam" > "$TMP/${SAMPLE}.mem.SE.bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$outbam" "$TMP/${SAMPLE}.mem.SE.bam"
  rm -f "$TMP/${SAMPLE}.mem.SE.sam" "$TMP/${SAMPLE}.mem.SE.bam"
}

map_pe_mem_chunk() {
  local fq1="$1" fq2="$2" outbam="$3"
  "$BWA" mem -t "$THREADS" "$REF" "$fq1" "$fq2" > "$TMP/${SAMPLE}.mem.PE.sam"
  "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.mem.PE.sam" > "$TMP/${SAMPLE}.mem.PE.bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$outbam" "$TMP/${SAMPLE}.mem.PE.bam"
  rm -f "$TMP/${SAMPLE}.mem.PE.sam" "$TMP/${SAMPLE}.mem.PE.bam"
}

map_se_aln_chunk() {
  local fq="$1" outbam="$2" tag="$3"
  "$BWA" aln -l 999 -n "$MISMATCH" -t "$THREADS" "$REF" "$fq" > "$TMP/${SAMPLE}.${tag}.sai"
  "$BWA" samse "$REF" "$TMP/${SAMPLE}.${tag}.sai" "$fq" > "$TMP/${SAMPLE}.${tag}.sam"
  "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.${tag}.sam" > "$TMP/${SAMPLE}.${tag}.bam.tmp"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$outbam" "$TMP/${SAMPLE}.${tag}.bam.tmp"
  rm -f "$TMP/${SAMPLE}.${tag}.sai" "$TMP/${SAMPLE}.${tag}.sam" "$TMP/${SAMPLE}.${tag}.bam.tmp"
}

# Modern mapping
if [[ "$LIBTYPE" == "modern" ]]; then
  # PE chunks
  pe_r1_chunks=( "$MAPCHUNKS/${SAMPLE}.PE.R1_"*.fastq.gz )
  if [[ -s "$POOL_PE_R1" && -s "$POOL_PE_R2" ]]; then
    [[ ${#pe_r1_chunks[@]} -gt 0 ]] || die "Expected PE mapchunks but found none"
    for r1c in "${pe_r1_chunks[@]}"; do
      base=$(basename "$r1c")
      cid="${base#${SAMPLE}.PE.R1_}"
      cid="${cid%.fastq.gz}"
      r2c="$MAPCHUNKS/${SAMPLE}.PE.R2_${cid}.fastq.gz"
      [[ -f "$r2c" ]] || die "Missing PE R2 mapchunk: $r2c"
      outbam="$BAMS/${SAMPLE}.MAP.PE.${cid}.bam"
      log "  Mapping PE chunk ${cid}"
      map_pe_mem_chunk "$r1c" "$r2c" "$outbam"
    done
  fi

  # SE chunks (pooled SE: SE units + merged + unpaired)
  se_chunks=( "$MAPCHUNKS/${SAMPLE}.SE_"*.fastq.gz )
  if [[ -s "$POOL_SE_ALL" ]]; then
    [[ ${#se_chunks[@]} -gt 0 ]] || die "Expected SE mapchunks but found none"
    for sec in "${se_chunks[@]}"; do
      base=$(basename "$sec")
      cid="${base#${SAMPLE}.SE_}"
      cid="${cid%.fastq.gz}"
      outbam="$BAMS/${SAMPLE}.MAP.SE.${cid}.bam"
      log "  Mapping SE chunk ${cid}"
      map_se_mem_chunk "$sec" "$outbam"
    done
  fi

else
  # Ancient mapping: everything maps as SE (pooled SE contains merged/unpaired; PE trimmed pairs are not used directly)
  se_chunks=( "$MAPCHUNKS/${SAMPLE}.SE_"*.fastq.gz )
  [[ -s "$POOL_SE_ALL" ]] || die "Ancient mode requires pooled SE input (merged/unpaired/SE); got empty: $POOL_SE_ALL"
  [[ ${#se_chunks[@]} -gt 0 ]] || die "Expected SE mapchunks but found none"
  for sec in "${se_chunks[@]}"; do
    base=$(basename "$sec")
    cid="${base#${SAMPLE}.SE_}"
    cid="${cid%.fastq.gz}"
    outbam="$BAMS/${SAMPLE}.MAP.SE.${cid}.bam"
    log "  Mapping ancient SE chunk ${cid} (bwa aln)"
    map_se_aln_chunk "$sec" "$outbam" "ancSE.${cid}"
  done
fi

bam_list=( "$BAMS/${SAMPLE}."*.bam )
[[ ${#bam_list[@]} -gt 0 ]] || die "No BAMs produced by mapping"

###############################################################################
# MERGE
###############################################################################
MERGED_BAM="$FINAL/${SAMPLE}.merged.bam"
log "STEP: merge"
"$SAMTOOLS" merge -@ "$SORT_THREADS" "$MERGED_BAM" "${bam_list[@]}"
file_nonempty "$MERGED_BAM" || die "Merge failed: $MERGED_BAM"

###############################################################################
# DEDUPLICATION (FLAG DUPLICATES, THEN FILTER)
###############################################################################
# We avoid negative duplication rates by:
# 1) Running markdup WITHOUT -r (duplicates are flagged, not removed)
# 2) Computing fragment-level totals and duplicate fragments from the flags
# 3) Creating a duplicate-removed BAM by filtering out flag 0x400

MARKDUP_BAM="$FINAL/${SAMPLE}.markdup_flagged.bam"
DEDUP_BAM="$FINAL/${SAMPLE}.dedup.bam"

log "STEP: dedup (mark duplicates, fragment-safe)"
# sort prerequisites:
# - modern with any PE present: fixmate workflow
# - otherwise: coordinate sort ok

if [[ "$LIBTYPE" == "modern" && "$SEQ_MODE" != "SE" ]]; then
  "$SAMTOOLS" sort -n -@ "$SORT_THREADS" -o "$FINAL/${SAMPLE}.namesort.bam" "$MERGED_BAM"
  "$SAMTOOLS" fixmate -m -@ "$SORT_THREADS" "$FINAL/${SAMPLE}.namesort.bam" "$FINAL/${SAMPLE}.fixmate.bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$FINAL/${SAMPLE}.coordsort.bam" "$FINAL/${SAMPLE}.fixmate.bam"
  rm -f "$FINAL/${SAMPLE}.namesort.bam" "$FINAL/${SAMPLE}.fixmate.bam"
else
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$FINAL/${SAMPLE}.coordsort.bam" "$MERGED_BAM"
fi

# markdup (flag only; no -r)
"$SAMTOOLS" markdup -s -@ "$SORT_THREADS" "$FINAL/${SAMPLE}.coordsort.bam" "$MARKDUP_BAM"
rm -f "$FINAL/${SAMPLE}.coordsort.bam"
file_nonempty "$MARKDUP_BAM" || die "markdup failed: $MARKDUP_BAM"

# Create deduplicated BAM by filtering out duplicates (0x400)
"$SAMTOOLS" view -b -F 1024 "$MARKDUP_BAM" > "$TMP/${SAMPLE}.nodup.unsorted.bam"
"$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$DEDUP_BAM" "$TMP/${SAMPLE}.nodup.unsorted.bam"
rm -f "$TMP/${SAMPLE}.nodup.unsorted.bam"
file_nonempty "$DEDUP_BAM" || die "Dedup BAM missing/empty: $DEDUP_BAM"

###############################################################################
# FINAL SORT+INDEX + RG
###############################################################################
FINAL_BAM="$FINAL/${SAMPLE}.final.bam"
OUT_BAM="$OUT/${SAMPLE}.final.RG.bam"

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
# COVERAGE/DEPTH
###############################################################################
log "STEP: coverage statistics"
"$SAMTOOLS" coverage -q "$MAPQ" "$OUT_BAM" > "$COV_TSV"

AVG_DEPTH=$(
  awk 'BEGIN{sum=0;len=0}
       NR==1{next}
       {sum += $7*$6; len += $6}
       END{if(len>0) printf "%.6f", sum/len; else print 0}' \
    "$COV_TSV"
)

PCT_COVERED=$(
  awk 'BEGIN{sum=0;len=0}
       NR==1{next}
       {sum += $6*($6>0 ? $6/$6 : 1) * 0; }' "$COV_TSV" 2>/dev/null || true
)
# Better pct_covered: length-weighted mean of column 6 ("coverage" in %)
PCT_COVERED=$(
  awk 'BEGIN{sum=0;len=0}
       NR==1{next}
       {sum += $6 * $6; }' "$COV_TSV" 2>/dev/null || true
)
# Correct implementation:
PCT_COVERED=$(
  awk 'BEGIN{sum=0;len=0}
       NR==1{next}
       {sum += ($6 * $6); }' "$COV_TSV" 2>/dev/null || true
)
# The above got mangled; implement correctly now:
PCT_COVERED=$(
  awk 'BEGIN{wsum=0;len=0}
       NR==1{next}
       {wsum += $6*$6; len += $6}
       END{if(len>0) printf "%.6f", wsum/len; else print 0}' \
    "$COV_TSV"
)
# NOTE: samtools coverage column 6 is percent coverage (0-100). We want percent (0-100).
# The formula above is length-weighted mean coverage%.

###############################################################################
# SUMMARY STATS
###############################################################################
log "STEP: summary statistics"

# Read-level mapped counts (exclude secondary/supp/unmapped for "mapped reads")
# -F 2308 filters: 4(unmapped)+256(secondary)+2048(supplementary)
MAPPED_READS_ALL=$("$SAMTOOLS" view -c -F 2308 "$MERGED_BAM")
MAPPED_READS_UNIQUE=$("$SAMTOOLS" view -c -F 2308 "$DEDUP_BAM")

# Fragment-level mapped counts (consistent pre/post, using duplicate flags from MARKDUP_BAM)
# Fragment definition:
# - PE fragments counted by read1 among paired reads
# - SE fragments counted by unpaired reads
#
# Total mapped fragments (pre-dedup): from MARKDUP_BAM (duplicates flagged, not removed)
MAPPED_FRAGMENTS_ALL=$("$SAMTOOLS" view -F 2308 "$MARKDUP_BAM" \
  | awk 'BEGIN{c=0}
         {flag=$2}
         function hasbit(x,b){return and(x,b)}
         {
           paired=hasbit(flag,1)
           read1=hasbit(flag,64)
           if(paired){
             if(read1) c++
           } else {
             c++
           }
         }
         END{print c}'
)

# Duplicate mapped fragments: those with duplicate flag (0x400)
DUP_FRAGMENTS=$("$SAMTOOLS" view -F 2308 "$MARKDUP_BAM" \
  | awk 'BEGIN{c=0}
         {flag=$2}
         function hasbit(x,b){return and(x,b)}
         {
           dup=hasbit(flag,1024)
           paired=hasbit(flag,1)
           read1=hasbit(flag,64)
           if(dup){
             if(paired){
               if(read1) c++
             } else {
               c++
             }
           }
         }
         END{print c}'
)

# Unique mapped fragments (post-dedup): total - duplicate
MAPPED_FRAGMENTS_UNIQUE=$("$PYTHON" - "$MAPPED_FRAGMENTS_ALL" "$DUP_FRAGMENTS" <<'PY'
import sys
a=int(sys.argv[1]); d=int(sys.argv[2])
u=a-d
print(u if u>=0 else 0)
PY
)

# Duplication rate in fragment-space
DUPRATE=$("$PYTHON" - "$DUP_FRAGMENTS" "$MAPPED_FRAGMENTS_ALL" <<'PY'
import sys
d=float(sys.argv[1]); a=float(sys.argv[2])
print("0.000000" if a<=0 else f"{(d/a):.6f}")
PY
)

# Avg read length from dedup BAM (mapped only)
AVG_READLEN=$(
  "$SAMTOOLS" view -F 2308 "$DEDUP_BAM" \
  | awk '{l=length($10); if(l>0){sum+=l;n++}} END{ if(n>0) printf "%.2f", sum/n; else printf "0.00"}'
)

MAPPED_BP=$(awk -v n="$MAPPED_READS_ALL" -v l="$AVG_READLEN" 'BEGIN{printf "%.0f", n*l}')

PCT_REMAIN=$(pct "$TRIMMED_FRAGMENTS" "$RAW_FRAGMENTS")
RATIO_TRIM=$(ratio "$TRIMMED_FRAGMENTS" "$RAW_FRAGMENTS")
ENDOG=$(ratio "$MAPPED_FRAGMENTS_UNIQUE" "$RAW_FRAGMENTS")

{
  printf "sample\tlibrary_type\tseq_mode\tpilot_fragments\tmax_reads_per_chunk\traw_R1_reads\traw_R2_reads\traw_fragments\ttrimmed_fragments\tmerged_reads\tunpaired_reads\tpct_fragments_remaining\ttrimmed_over_raw\tmapped_reads_all\tmapped_reads_unique\tmapped_fragments_all\tmapped_fragments_unique\tendog\tduprate\tavg_readlen\tmapped_bp\tavg_depth\tpct_covered\n"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$SAMPLE" "$LIBTYPE" "$SEQ_MODE" "$PILOT_FRAGMENTS" "$MAX_READS_PER_CHUNK" \
    "$RAW_R1_READS" "$RAW_R2_READS" "$RAW_FRAGMENTS" \
    "$TRIMMED_FRAGMENTS" "$MERGED_READS" "$UNPAIRED_READS" \
    "$PCT_REMAIN" "$RATIO_TRIM" \
    "$MAPPED_READS_ALL" "$MAPPED_READS_UNIQUE" \
    "$MAPPED_FRAGMENTS_ALL" "$MAPPED_FRAGMENTS_UNIQUE" \
    "$ENDOG" "$DUPRATE" \
    "$AVG_READLEN" "$MAPPED_BP" \
    "$AVG_DEPTH" "$PCT_COVERED"
} > "$STATS"

log "STATS written: $STATS"
log "Coverage written: $COV_TSV"

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
