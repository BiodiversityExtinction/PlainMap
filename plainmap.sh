#!/usr/bin/env bash
# PlainMap v0.1
# Transparent, failure-aware mapping pipeline for ancient and modern DNA
#
# Key properties (current dev version):
# - Manifest can contain mixed SE + PE, mixed platforms, mixed naming conventions
# - SE/PE classification is derived from FASTQ headers:
#     * Illumina/CASAVA: header contains ' 1:' / ' 2:'
#     * ENA/older Illumina: first token ends with /1 or /2
#     * SRA-export FASTQ (often lacks mate markers): R1/R2 headers are identical;
#       PlainMap pairs files by identical normalized first-read name (key) and
#       assigns mates deterministically by manifest order (first=R1, second=R2),
#       with strict sanity checks.
# - fastp runs per unit (or per unit-chunk), then pools:
#     * pooled unmerged PE (R1+R2) kept separate
#     * pooled SE-like reads (merged + unpaired + SE) kept separate
# - mapping is restartable: optional mapping chunking after pooling, with strict mate-synchrony checks
# - pigz is used automatically if available (decompress/test/compress)
# - fragment-aware stats + duplication rate derived from duplicate flags (samtools markdup)
#
# Library types:
#   modern     : bwa mem; maps unmerged PE as PE and SE-like as SE
#   ancient    : bwa aln/samse; maps SE-like only (merged + unpaired + SE)
#   historical : bwa aln; maps unmerged PE with aln+sampe AND SE-like with aln+samse
#
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
MISMATCH=0.01              # bwa aln -n (ancient/historical)
LIBTYPE="modern"           # modern|ancient|historical
MAX_READS_PER_CHUNK=0      # 0 disables pre-fastp chunking safety valve
PILOT_FRAGMENTS=0          # 0 disables pilot; otherwise GLOBAL cap on fragments BEFORE fastp
MAPQ=20                    # mapping quality threshold (samtools view -q)

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
  --manifest FILE
  --prefix STRING
  --ref FILE
  --outdir DIR

Optional:
  --library-type modern|ancient|historical   (default: modern)
  --threads INT                              (default: 1)
  --minlength INT                            (default: 30)
  --mismatch FLOAT                           (ancient/historical; bwa aln -n; default: 0.01)
  --mapq INT                                 Mapping quality threshold (default: 20)
  --max-reads-per-chunk INT                  (default: 0; disabled) pre-fastp chunking safety valve
  --pilot-fragments INT                      (default: 0; disabled) GLOBAL cap on fragments before fastp
  --adapter-r1 SEQ
  --adapter-r2 SEQ
  --trim-only                                Trim only (fastp); exit before mapping
  --tmpdir DIR                               Custom temp dir (e.g. node-local scratch)
  --keep-intermediate                        Keep <outdir>/<prefix>/work
  --resume | --no-resume
  --dry-run                                  Print plan only
  --validate                                 Check tools + gzip/pigz -t all manifest FASTQs, then exit
  --reset                                    Clear ALL checkpoints for this sample

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
# PARSE ARGS (double-dash only)
###############################################################################
MANIFEST=""
SAMPLE=""
REF=""
OUT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --prefix) SAMPLE="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --outdir) OUT="$2"; shift 2 ;;
    --library-type) LIBTYPE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --minlength) MINLEN="$2"; shift 2 ;;
    --mismatch) MISMATCH="$2"; shift 2 ;;
    --mapq) MAPQ="$2"; shift 2 ;;
    --max-reads-per-chunk) MAX_READS_PER_CHUNK="$2"; shift 2 ;;
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
[[ "$LIBTYPE" == "modern" || "$LIBTYPE" == "ancient" || "$LIBTYPE" == "historical" ]] || {
  echo "ERROR: --library-type must be modern|ancient|historical"
  exit 1
}
[[ "$MAPQ" =~ ^[0-9]+$ ]] || { echo "ERROR: --mapq must be an integer"; exit 1; }
[[ "$THREADS" =~ ^[0-9]+$ ]] || { echo "ERROR: --threads must be an integer"; exit 1; }
[[ "$MAX_READS_PER_CHUNK" =~ ^[0-9]+$ ]] || { echo "ERROR: --max-reads-per-chunk must be an integer"; exit 1; }
[[ "$PILOT_FRAGMENTS" =~ ^[0-9]+$ ]] || { echo "ERROR: --pilot-fragments must be an integer"; exit 1; }

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
MAPCHUNKS="$WORK/mapchunks"    # post-fastp mapping chunks
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
  if [[ $HAVE_PIGZ -eq 1 ]]; then
    require_cmd pigz
  else
    require_cmd gzip
    require_cmd zcat
  fi
  [[ "$LIBTYPE" != "ancient" ]] || require_cmd "$MAPDAMAGE"
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
  log "  library_type:        $LIBTYPE"
  log "  threads:             $THREADS"
  log "  minlength:           $MINLEN"
  log "  mismatch (aln):      $MISMATCH"
  log "  mapq:                $MAPQ"
  log "  pilot_fragments:     $PILOT_FRAGMENTS (GLOBAL cap)"
  log "  max_reads_per_chunk: $MAX_READS_PER_CHUNK"
  log "Plan: manifest -> build SE/PE units -> pre-fastp chunk (optional) -> pilot (optional) -> fastp -> pool (SE-like vs unmerged PE) -> map (chunked optional) -> dedup -> RG -> coverage -> stats"
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
[[ -n "$ADAPTER_R1" ]] && FASTP_ARGS+=( --adapter_sequence "$ADAPTER_R1" )
[[ -n "$ADAPTER_R2" ]] && FASTP_ARGS+=( --adapter_sequence_r2 "$ADAPTER_R2" )

###############################################################################
# JSON HELPERS (robust across fastp schemas)
###############################################################################
json_get_int_any() {
  # usage: json_get_int_any <json> <key1> <key2> ...
  "$PYTHON" - "$@" <<'PY'
import json,sys
path=sys.argv[1]
keys_list=sys.argv[2:]
with open(path) as f:
    d=json.load(f)

def get_path(obj, path):
    cur=obj
    for k in path.split("."):
        if isinstance(cur, dict) and k in cur:
            cur=cur[k]
        else:
            return None
    return cur

for k in keys_list:
    v=get_path(d,k)
    if v is None:
        continue
    try:
        iv=int(v)
    except Exception:
        continue
    if iv != 0:
        print(iv)
        sys.exit(0)

for k in keys_list:
    v=get_path(d,k)
    if v is None:
        continue
    try:
        print(int(v))
        sys.exit(0)
    except Exception:
        pass

print(0)
PY
}

###############################################################################
# STATS HELPERS
###############################################################################
count_reads_fastq_gz() { "${ZCAT[@]}" "$1" | awk 'END{print NR/4}'; }

py_add() { "$PYTHON" - "$1" "$2" <<'PY'
import sys
print(int(sys.argv[1]) + int(sys.argv[2]))
PY
}

py_sub_nonneg() {
  "$PYTHON" - "$1" "$2" <<'PY'
import sys
a=int(sys.argv[1]); b=int(sys.argv[2])
x=a-b
print(x if x>=0 else 0)
PY
}

###############################################################################
# FASTQ HEADER PARSING / DIRECTION / PAIR KEYS (HEADER-ONLY, SRA-SAFE)
###############################################################################
first_header_line() {
  local fq="$1"
  "${ZCAT[@]}" "$fq" | head -n 1 || true
}

normalize_qname() {
  # normalize a FASTQ header line into a "mate key":
  # - take first token
  # - strip leading @
  # - strip trailing /1 or /2 if present
  local hdr="$1"
  local first="${hdr%% *}"
  first="${first#@}"
  first="${first%/1}"
  first="${first%/2}"
  echo "$first"
}

classify_fastq_direction() {
  # prints: R1 | R2 | UNKNOWN
  local fq="$1"
  local hdr firsttok rest secondtok
  hdr="$(first_header_line "$fq")"
  [[ -n "$hdr" ]] || { echo "UNKNOWN"; return; }

  # Illumina/CASAVA token " 1:" / " 2:"
  if [[ "$hdr" == *" 1:"* ]]; then echo "R1"; return; fi
  if [[ "$hdr" == *" 2:"* ]]; then echo "R2"; return; fi

  # Tokenize
  firsttok="${hdr%% *}"
  rest="${hdr#"$firsttok"}"; rest="${rest# }"
  secondtok="${rest%% *}"

  # Mate marker sometimes appears on token2 (common in many public FASTQs)
  if [[ "$secondtok" == */1 ]]; then echo "R1"; return; fi
  if [[ "$secondtok" == */2 ]]; then echo "R2"; return; fi

  # ENA/older Illumina style: token1 ends with /1 or /2
  if [[ "$firsttok" == */1 ]]; then echo "R1"; return; fi
  if [[ "$firsttok" == */2 ]]; then echo "R2"; return; fi

  # SRA-style (occasionally): token2 is literally "1" or "2"
  if [[ "$secondtok" == "1" ]]; then echo "R1"; return; fi
  if [[ "$secondtok" == "2" ]]; then echo "R2"; return; fi

  echo "UNKNOWN"
}

check_headers_identical_sra() {
  # Compare normalized qname roots (strips /1,/2 and keeps stable part)
  local f1="$1" f2="$2"
  local h1 h2 k1 k2
  h1="$(first_header_line "$f1")"
  h2="$(first_header_line "$f2")"
  k1="$(normalize_qname "$h1")"
  k2="$(normalize_qname "$h2")"
  [[ "$k1" == "$k2" ]] || die "SRA/unknown-PE sanity check failed: normalized first read names differ:
$f1 -> $k1
$f2 -> $k2"
}

###############################################################################
# BUILD UNITS SAFELY (NO ZIP-BY-INDEX ACROSS DIFFERENT KEYS; SRA UNKNOWN handled)
###############################################################################
declare -A R1_BY_KEY
declare -A R2_BY_KEY
declare -A UNK_BY_KEY        # newline-separated "idx:path" entries
declare -A DIR_BY_PATH
declare -A HDR_BY_PATH

ALL_FILES=()
TOTAL_FILES=0
IDX=0

while read -r f; do
  [[ -z "$f" || "$f" =~ ^# ]] && continue
  [[ "$f" = /* ]] || f="$MANIFEST_DIR/$f"
  f="$(resolve_path "$f")"
  [[ -r "$f" ]] || die "FASTQ not readable: $f"

  TOTAL_FILES=$((TOTAL_FILES+1))
  ALL_FILES+=( "$f" )
  IDX=$((IDX+1))

  hdr="$(first_header_line "$f")"
  [[ -n "$hdr" ]] || die "Could not read FASTQ header for: $f"
  HDR_BY_PATH["$f"]="$hdr"

  dir="$(classify_fastq_direction "$f")"
  DIR_BY_PATH["$f"]="$dir"

  key="$(normalize_qname "$hdr")"
  [[ -n "$key" ]] || die "Could not derive header key for: $f"

  if [[ "$dir" == "R1" ]]; then
    [[ -z "${R1_BY_KEY[$key]:-}" ]] || die "Multiple R1 files share the same first-read key '$key' (unexpected). Offending: $f and ${R1_BY_KEY[$key]}"
    R1_BY_KEY["$key"]="$f"
  elif [[ "$dir" == "R2" ]]; then
    [[ -z "${R2_BY_KEY[$key]:-}" ]] || die "Multiple R2 files share the same first-read key '$key' (unexpected). Offending: $f and ${R2_BY_KEY[$key]}"
    R2_BY_KEY["$key"]="$f"
  else
    # UNKNOWN: store with manifest order index to deterministically assign mates later
    if [[ -z "${UNK_BY_KEY[$key]:-}" ]]; then
      UNK_BY_KEY["$key"]="${IDX}:${f}"
    else
      UNK_BY_KEY["$key"]+=$'\n'"${IDX}:${f}"
    fi
  fi
done < "$MANIFEST"

[[ $TOTAL_FILES -gt 0 ]] || die "Manifest empty: $MANIFEST"

PE_KEYS=()
SE_FILES=()

# 1) If we have explicit R2 without R1 for a key, error (safer than guessing)
for k in "${!R2_BY_KEY[@]}"; do
  [[ -n "${R1_BY_KEY[$k]:-}" ]] || die "Found R2 without matching R1 for key '$k' (check manifest): ${R2_BY_KEY[$k]}"
done

# 2) Explicitly marked R1 keys: either PE (if matching R2) or SE
for k in "${!R1_BY_KEY[@]}"; do
  if [[ -n "${R2_BY_KEY[$k]:-}" ]]; then
    PE_KEYS+=( "$k" )
  else
    SE_FILES+=( "${R1_BY_KEY[$k]}" )
  fi
done

# 3) UNKNOWN keys: if exactly 1 file -> SE; if exactly 2 -> SRA-style PE; else error
for k in "${!UNK_BY_KEY[@]}"; do
  mapfile -t entries <<< "${UNK_BY_KEY[$k]}"
  if [[ "${#entries[@]}" -eq 1 ]]; then
    # single unknown file => treat as SE
    f="${entries[0]#*:}"
    SE_FILES+=( "$f" )
  elif [[ "${#entries[@]}" -eq 2 ]]; then
    # two unknown files => treat as SRA-style PE; assign by manifest order
    # sort by idx (small, so simple compare)
    e1="${entries[0]}"; e2="${entries[1]}"
    i1="${e1%%:*}"; f1="${e1#*:}"
    i2="${e2%%:*}"; f2="${e2#*:}"
    if [[ "$i2" -lt "$i1" ]]; then
      tmpi="$i1"; tmpf="$f1"
      i1="$i2"; f1="$f2"
      i2="$tmpi"; f2="$tmpf"
    fi

    # sanity: headers should match (SRA-style), and read counts should match
    check_headers_identical_sra "$f1" "$f2"
    n1=$(count_reads_fastq_gz "$f1")
    n2=$(count_reads_fastq_gz "$f2")
    [[ "$n1" -eq "$n2" ]] || die "SRA-style PE sanity check failed: read counts differ for key '$k' (R1?=$n1 R2?=$n2):
$f1
$f2"

    # Assign R1/R2 deterministically (manifest order); store as explicit R1/R2 for downstream
    [[ -z "${R1_BY_KEY[$k]:-}" && -z "${R2_BY_KEY[$k]:-}" ]] || die "Internal error: UNKNOWN key '$k' already classified"
    R1_BY_KEY["$k"]="$f1"
    R2_BY_KEY["$k"]="$f2"
    PE_KEYS+=( "$k" )
  else
    die "Ambiguous SRA-style grouping for key '$k': found ${#entries[@]} files with no mate markers. PlainMap refuses to guess."
  fi
done

# Determine seq_mode summary
SEQ_MODE="SE"
if [[ ${#PE_KEYS[@]} -gt 0 && ${#SE_FILES[@]} -gt 0 ]]; then
  SEQ_MODE="MIX"
elif [[ ${#PE_KEYS[@]} -gt 0 ]]; then
  SEQ_MODE="PE"
fi

log "Detected sequencing mode: $SEQ_MODE (PE keys: ${#PE_KEYS[@]}, SE files: ${#SE_FILES[@]})"

###############################################################################
# PRE-FASTP CHUNKING HELPERS
###############################################################################
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
  local in="$1" out="$2" reads="$3"
  local lines=$((reads * 4))
  set +o pipefail
  "${ZCAT[@]}" "$in" | head -n "$lines" | "${GZIP[@]}" > "$out"
  set -o pipefail
  file_nonempty "$out" || die "Pilot limiting produced empty FASTQ: $out"
}

pilot_take() {
  "$PYTHON" - "$1" "$2" <<'PY'
import sys
c=int(sys.argv[1]); p=int(sys.argv[2])
if p<=0: print(0)
else: print(min(c,p))
PY
}

###############################################################################
# CREATE PRE-FASTP CHUNKS PER UNIT (PE stays paired)
###############################################################################
make_chunks() {
  log "STEP: pre-fastp chunking"
  ckpt_clear chunk
  rm -f "$CHUNKS/${SAMPLE}."*.fastq.gz 2>/dev/null || true

  # PE units: one R1 and one R2 per key
  local idx=0
  for k in "${PE_KEYS[@]}"; do
    local r1="${R1_BY_KEY[$k]}"
    local r2="${R2_BY_KEY[$k]}"
    [[ -n "$r1" && -n "$r2" ]] || die "Internal error: PE key '$k' missing mate"

    local unit="PE$(printf '%04d' $idx)"
    local prefix="$CHUNKS/${SAMPLE}.${unit}."

    if [[ "$MAX_READS_PER_CHUNK" -gt 0 ]]; then
      log "  PE unit $unit (key=$k): splitting into chunks (max reads/chunk: $MAX_READS_PER_CHUNK)"
      split_fastq_gz "$r1" "${prefix}R1_" "$MAX_READS_PER_CHUNK"
      split_fastq_gz "$r2" "${prefix}R2_" "$MAX_READS_PER_CHUNK"
    else
      log "  PE unit $unit (key=$k): chunking disabled; linking as single chunk"
      ln -sf "$r1" "${prefix}R1_0000.fastq.gz"
      ln -sf "$r2" "${prefix}R2_0000.fastq.gz"
    fi

    file_nonempty "${prefix}R1_0000.fastq.gz" || die "Missing PE R1 chunk for $unit"
    file_nonempty "${prefix}R2_0000.fastq.gz" || die "Missing PE R2 chunk for $unit"
    idx=$((idx+1))
  done

  # SE units: each file independent
  local sidx=0
  for fq in "${SE_FILES[@]}"; do
    local unit="SE$(printf '%04d' $sidx)"
    local prefix="$CHUNKS/${SAMPLE}.${unit}."
    if [[ "$MAX_READS_PER_CHUNK" -gt 0 ]]; then
      log "  SE unit $unit: splitting into chunks (max reads/chunk: $MAX_READS_PER_CHUNK)"
      split_fastq_gz "$fq" "${prefix}SE_" "$MAX_READS_PER_CHUNK"
    else
      log "  SE unit $unit: chunking disabled; linking as single chunk"
      ln -sf "$fq" "${prefix}SE_0000.fastq.gz"
    fi
    file_nonempty "${prefix}SE_0000.fastq.gz" || die "Missing SE chunk for $unit"
    sidx=$((sidx+1))
  done

  ckpt_mark chunk
}

if ! ckpt_ok chunk "$CHUNKS"; then
  make_chunks
else
  log "STEP: pre-fastp chunking (skipped; checkpoint present)"
fi

###############################################################################
# FASTP PER PRE-FASTP CHUNK -> POOL TRIMMED OUTPUTS
# - unmerged PE kept separate
# - SE-like (merged + unpaired + SE) pooled together
###############################################################################
POOL_PE_R1="$RAW/${SAMPLE}.pool.PE.R1.fastq.gz"
POOL_PE_R2="$RAW/${SAMPLE}.pool.PE.R2.fastq.gz"
POOL_SE_ALL="$RAW/${SAMPLE}.pool.SE.all.fastq.gz"

: > "$POOL_PE_R1"
: > "$POOL_PE_R2"
: > "$POOL_SE_ALL"

RAW_R1_READS=0
RAW_R2_READS=0
RAW_FRAGMENTS=0

TRIMMED_FRAGMENTS=0
MERGED_READS=0
UNPAIRED_READS=0

PILOT_LEFT="$PILOT_FRAGMENTS"

log "STEP: fastp per unit-chunk (JSON only) + pooling (unmerged PE separate from SE-like)"
shopt -s nullglob

# Process PE unit-chunks
pe_r1_chunks=( "$CHUNKS/${SAMPLE}.PE"*".R1_"*.fastq.gz )
for r1c in "${pe_r1_chunks[@]}"; do
  [[ "$PILOT_FRAGMENTS" -gt 0 && "$PILOT_LEFT" -le 0 ]] && break

  base=$(basename "$r1c")
  r2c="${r1c/.R1_/.R2_}"
  [[ -f "$r2c" ]] || die "Missing PE mate chunk for: $r1c (expected $r2c)"

  in_r1="$r1c"
  in_r2="$r2c"

  pairs_in_chunk=$(count_reads_fastq_gz "$in_r1")
  take_pairs="$pairs_in_chunk"

  if [[ "$PILOT_FRAGMENTS" -gt 0 ]]; then
    take_pairs=$(pilot_take "$pairs_in_chunk" "$PILOT_LEFT")
    [[ "$take_pairs" -le 0 ]] && break
    if [[ "$take_pairs" -lt "$pairs_in_chunk" ]]; then
      log "  Pilot: limiting PE chunk $(basename "$r1c") to first $take_pairs pairs (global remaining before: $PILOT_LEFT)"
    fi
    lim_r1="$TMP/${SAMPLE}.$(basename "${r1c%.fastq.gz}").pilot.R1.fastq.gz"
    lim_r2="$TMP/${SAMPLE}.$(basename "${r2c%.fastq.gz}").pilot.R2.fastq.gz"
    limit_fastq_gz "$in_r1" "$lim_r1" "$take_pairs"
    limit_fastq_gz "$in_r2" "$lim_r2" "$take_pairs"
    in_r1="$lim_r1"
    in_r2="$lim_r2"
    PILOT_LEFT="$(py_sub_nonneg "$PILOT_LEFT" "$take_pairs")"
  fi

  RAW_R1_READS=$(py_add "$RAW_R1_READS" "$take_pairs")
  RAW_R2_READS=$(py_add "$RAW_R2_READS" "$take_pairs")
  RAW_FRAGMENTS=$(py_add "$RAW_FRAGMENTS" "$take_pairs")

  tag="$(basename "${r1c%.fastq.gz}")"
  tag="${tag//\//_}"

  fastp_err="$TMP/fastp_${SAMPLE}.${tag}.stderr.log"
  fp_json="$REPORTS/${SAMPLE}.${tag}.fastp.json"

  trim1="$TMP/${SAMPLE}.${tag}.trim.R1.fastq.gz"
  trim2="$TMP/${SAMPLE}.${tag}.trim.R2.fastq.gz"
  merged="$TMP/${SAMPLE}.${tag}.merged.fastq.gz"
  up1="$TMP/${SAMPLE}.${tag}.unpaired1.fastq.gz"
  up2="$TMP/${SAMPLE}.${tag}.unpaired2.fastq.gz"

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

  file_nonempty "$trim1" || die "fastp produced empty trimmed R1 for $tag"
  file_nonempty "$trim2" || die "fastp produced empty trimmed R2 for $tag"

  cat "$trim1" >> "$POOL_PE_R1"
  cat "$trim2" >> "$POOL_PE_R2"

  file_nonempty "$merged" && cat "$merged" >> "$POOL_SE_ALL"
  file_nonempty "$up1"   && cat "$up1"   >> "$POOL_SE_ALL"
  file_nonempty "$up2"   && cat "$up2"   >> "$POOL_SE_ALL"

  out1_reads=$(json_get_int_any "$fp_json" "read1_after_filtering.total_reads")
  [[ "$out1_reads" -le 0 ]] && out1_reads=$(count_reads_fastq_gz "$trim1")
  TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$out1_reads")

  merged_reads=$(json_get_int_any "$fp_json" "merged_and_filtered.total_reads")
  if [[ "$merged_reads" -le 0 && -s "$merged" ]]; then
    merged_reads=$(count_reads_fastq_gz "$merged")
  fi
  MERGED_READS=$(py_add "$MERGED_READS" "$merged_reads")
  TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$merged_reads")

  up1n=0; up2n=0
  [[ -s "$up1" ]] && up1n=$(count_reads_fastq_gz "$up1")
  [[ -s "$up2" ]] && up2n=$(count_reads_fastq_gz "$up2")
  UNPAIRED_READS=$(py_add "$UNPAIRED_READS" "$up1n")
  UNPAIRED_READS=$(py_add "$UNPAIRED_READS" "$up2n")
  TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$up1n")
  TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$up2n")

  rm -f "$trim1" "$trim2" "$merged" "$up1" "$up2" 2>/dev/null || true
  [[ "$PILOT_FRAGMENTS" -gt 0 ]] && rm -f "$lim_r1" "$lim_r2" 2>/dev/null || true
done

# Process SE unit-chunks
se_chunks=( "$CHUNKS/${SAMPLE}.SE"*".SE_"*.fastq.gz )
for sec in "${se_chunks[@]}"; do
  [[ "$PILOT_FRAGMENTS" -gt 0 && "$PILOT_LEFT" -le 0 ]] && break

  in_se="$sec"
  reads_in_chunk=$(count_reads_fastq_gz "$in_se")
  take_reads="$reads_in_chunk"

  if [[ "$PILOT_FRAGMENTS" -gt 0 ]]; then
    take_reads=$(pilot_take "$reads_in_chunk" "$PILOT_LEFT")
    [[ "$take_reads" -le 0 ]] && break
    if [[ "$take_reads" -lt "$reads_in_chunk" ]]; then
      log "  Pilot: limiting SE chunk $(basename "$sec") to first $take_reads reads (global remaining before: $PILOT_LEFT)"
    fi
    lim_se="$TMP/${SAMPLE}.$(basename "${sec%.fastq.gz}").pilot.SE.fastq.gz"
    limit_fastq_gz "$in_se" "$lim_se" "$take_reads"
    in_se="$lim_se"
    PILOT_LEFT="$(py_sub_nonneg "$PILOT_LEFT" "$take_reads")"
  fi

  RAW_R1_READS=$(py_add "$RAW_R1_READS" "$take_reads")
  RAW_FRAGMENTS=$(py_add "$RAW_FRAGMENTS" "$take_reads")

  tag="$(basename "${sec%.fastq.gz}")"
  tag="${tag//\//_}"

  fastp_err="$TMP/fastp_${SAMPLE}.${tag}.stderr.log"
  fp_json="$REPORTS/${SAMPLE}.${tag}.fastp.json"
  trim_se="$TMP/${SAMPLE}.${tag}.trim.SE.fastq.gz"

  run_fastp "$fastp_err" "${FASTP_ARGS[@]}" \
    -i "$in_se" --out1 "$trim_se" \
    -j "$fp_json"

  file_nonempty "$trim_se" || die "fastp produced empty trimmed SE for $tag"

  cat "$trim_se" >> "$POOL_SE_ALL"

  se_after=$(json_get_int_any "$fp_json" "read1_after_filtering.total_reads" "summary.after_filtering.total_reads")
  [[ "$se_after" -le 0 ]] && se_after=$(count_reads_fastq_gz "$trim_se")
  TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$se_after")

  rm -f "$trim_se" 2>/dev/null || true
  [[ "$PILOT_FRAGMENTS" -gt 0 ]] && rm -f "$lim_se" 2>/dev/null || true
done

###############################################################################
# TRIM-ONLY MODE
###############################################################################
if [[ $TRIM_ONLY -eq 1 ]]; then
  log "Trim-only mode enabled; writing stats and exiting"

  {
    printf "sample\tlibrary_type\tseq_mode\tpilot_fragments\tmax_reads_per_chunk\tmapq\traw_R1_reads\traw_R2_reads\traw_fragments\ttrimmed_fragments\tmerged_reads\tunpaired_reads\tmapped_reads_all\tmapped_reads_unique\tmapped_fragments_all\tmapped_fragments_unique\tendog\tduprate\tavg_readlen\tmapped_bp\tavg_depth\tpct_covered\n"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n" \
      "$SAMPLE" "$LIBTYPE" "$SEQ_MODE" "$PILOT_FRAGMENTS" "$MAX_READS_PER_CHUNK" "$MAPQ" \
      "$RAW_R1_READS" "$RAW_R2_READS" "$RAW_FRAGMENTS" \
      "$TRIMMED_FRAGMENTS" "$MERGED_READS" "$UNPAIRED_READS"
  } > "$STATS"

  PIPE_T1=$(date +%s)
  log "STATS written: $STATS"
  log "Total wall time: $((PIPE_T1 - PIPE_T0)) seconds"
  exit 0
fi

###############################################################################
# POST-FASTP MAPPING CHUNKING (OPTIONAL; REUSE MAX_READS_PER_CHUNK)
# CRITICAL: PE pools must remain synchronized for sampe; we enforce this.
###############################################################################
rm -f "$MAPCHUNKS/${SAMPLE}."*.fastq.gz 2>/dev/null || true

make_mapchunks_se() {
  local in="$1" prefix="$2"
  [[ -s "$in" ]] || return 0
  if [[ "$MAX_READS_PER_CHUNK" -gt 0 ]]; then
    log "  Mapping chunking SE: splitting pooled SE-like (max reads/chunk: $MAX_READS_PER_CHUNK)"
    split_fastq_gz "$in" "$prefix" "$MAX_READS_PER_CHUNK"
  else
    ln -sf "$in" "${prefix}0000.fastq.gz"
  fi
}

make_mapchunks_pe() {
  local in1="$1" in2="$2" prefix1="$3" prefix2="$4"
  [[ -s "$in1" && -s "$in2" ]] || return 0

  local n1 n2
  n1=$(count_reads_fastq_gz "$in1")
  n2=$(count_reads_fastq_gz "$in2")
  [[ "$n1" -eq "$n2" ]] || die "Pooled PE R1/R2 read counts differ (R1=$n1 R2=$n2). This would break sampe."

  if [[ "$MAX_READS_PER_CHUNK" -gt 0 ]]; then
    log "  Mapping chunking PE: splitting pooled PE (max reads/chunk: $MAX_READS_PER_CHUNK)"
    split_fastq_gz "$in1" "$prefix1" "$MAX_READS_PER_CHUNK"
    split_fastq_gz "$in2" "$prefix2" "$MAX_READS_PER_CHUNK"
  else
    ln -sf "$in1" "${prefix1}0000.fastq.gz"
    ln -sf "$in2" "${prefix2}0000.fastq.gz"
  fi

  shopt -s nullglob
  local r1_chunks=( "${prefix1}"*.fastq.gz )
  local r2_chunks=( "${prefix2}"*.fastq.gz )
  [[ "${#r1_chunks[@]}" -eq "${#r2_chunks[@]}" ]] || die "PE mapping chunk count mismatch after splitting"

  for ((i=0; i<${#r1_chunks[@]}; i++)); do
    local f1="${r1_chunks[$i]}"
    local f2="${r2_chunks[$i]}"
    local h1 h2 n1x n2x
    h1="$(first_header_line "$f1")"; h2="$(first_header_line "$f2")"
    n1x="$(normalize_qname "$h1")"; n2x="$(normalize_qname "$h2")"
    [[ "$n1x" == "$n2x" ]] || die "PE mapping chunks out-of-sync (first read name mismatch): $(basename "$f1")='$n1x' vs $(basename "$f2")='$n2x'"
  done
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
  local fq="$1" outbam="$2" tag="$3"
  "$BWA" mem -t "$THREADS" "$REF" "$fq" > "$TMP/${SAMPLE}.${tag}.sam"
  "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.${tag}.sam" > "$TMP/${SAMPLE}.${tag}.bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$outbam" "$TMP/${SAMPLE}.${tag}.bam"
  rm -f "$TMP/${SAMPLE}.${tag}.sam" "$TMP/${SAMPLE}.${tag}.bam"
}

map_pe_mem_chunk() {
  local fq1="$1" fq2="$2" outbam="$3" tag="$4"
  "$BWA" mem -t "$THREADS" "$REF" "$fq1" "$fq2" > "$TMP/${SAMPLE}.${tag}.sam"
  "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.${tag}.sam" > "$TMP/${SAMPLE}.${tag}.bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$outbam" "$TMP/${SAMPLE}.${tag}.bam"
  rm -f "$TMP/${SAMPLE}.${tag}.sam" "$TMP/${SAMPLE}.${tag}.bam"
}

map_se_aln_chunk() {
  local fq="$1" outbam="$2" tag="$3"
  "$BWA" aln -l 999 -n "$MISMATCH" -t "$THREADS" "$REF" "$fq" > "$TMP/${SAMPLE}.${tag}.sai"
  "$BWA" samse "$REF" "$TMP/${SAMPLE}.${tag}.sai" "$fq" > "$TMP/${SAMPLE}.${tag}.sam"
  "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.${tag}.sam" > "$TMP/${SAMPLE}.${tag}.bam.tmp"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$outbam" "$TMP/${SAMPLE}.${tag}.bam.tmp"
  rm -f "$TMP/${SAMPLE}.${tag}.sai" "$TMP/${SAMPLE}.${tag}.sam" "$TMP/${SAMPLE}.${tag}.bam.tmp"
}

map_pe_aln_chunk() {
  local fq1="$1" fq2="$2" outbam="$3" tag="$4"

  local h1 h2 n1x n2x
  h1="$(first_header_line "$fq1")"; h2="$(first_header_line "$fq2")"
  n1x="$(normalize_qname "$h1")"; n2x="$(normalize_qname "$h2")"
  [[ "$n1x" == "$n2x" ]] || die "Refusing to run sampe: PE chunks first read names differ: $(basename "$fq1")='$n1x' vs $(basename "$fq2")='$n2x'"

  "$BWA" aln -l 999 -n "$MISMATCH" -t "$THREADS" "$REF" "$fq1" > "$TMP/${SAMPLE}.${tag}.1.sai"
  "$BWA" aln -l 999 -n "$MISMATCH" -t "$THREADS" "$REF" "$fq2" > "$TMP/${SAMPLE}.${tag}.2.sai"
  "$BWA" sampe "$REF" "$TMP/${SAMPLE}.${tag}.1.sai" "$TMP/${SAMPLE}.${tag}.2.sai" "$fq1" "$fq2" > "$TMP/${SAMPLE}.${tag}.sam"
  "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.${tag}.sam" > "$TMP/${SAMPLE}.${tag}.bam.tmp"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$outbam" "$TMP/${SAMPLE}.${tag}.bam.tmp"
  rm -f "$TMP/${SAMPLE}.${tag}.1.sai" "$TMP/${SAMPLE}.${tag}.2.sai" "$TMP/${SAMPLE}.${tag}.sam" "$TMP/${SAMPLE}.${tag}.bam.tmp"
}

# Map unmerged PE (if any)
if [[ -s "$POOL_PE_R1" && -s "$POOL_PE_R2" ]]; then
  pe_r1_chunks=( "$MAPCHUNKS/${SAMPLE}.PE.R1_"*.fastq.gz )
  [[ ${#pe_r1_chunks[@]} -gt 0 ]] || die "Expected PE mapchunks but found none"
  for r1c in "${pe_r1_chunks[@]}"; do
    base=$(basename "$r1c")
    cid="${base#${SAMPLE}.PE.R1_}"
    cid="${cid%.fastq.gz}"
    r2c="$MAPCHUNKS/${SAMPLE}.PE.R2_${cid}.fastq.gz"
    [[ -f "$r2c" ]] || die "Missing PE R2 mapchunk: $r2c"
    outbam="$BAMS/${SAMPLE}.MAP.PE.${cid}.bam"

    if [[ "$LIBTYPE" == "modern" ]]; then
      log "  Mapping PE chunk ${cid} (mem)"
      map_pe_mem_chunk "$r1c" "$r2c" "$outbam" "mem.PE.${cid}"
    elif [[ "$LIBTYPE" == "historical" ]]; then
      log "  Mapping PE chunk ${cid} (aln/sampe)"
      map_pe_aln_chunk "$r1c" "$r2c" "$outbam" "aln.PE.${cid}"
    else
      :
    fi
  done
fi

# Map SE-like pool (merged + unpaired + SE)
if [[ -s "$POOL_SE_ALL" ]]; then
  se_chunks=( "$MAPCHUNKS/${SAMPLE}.SE_"*.fastq.gz )
  [[ ${#se_chunks[@]} -gt 0 ]] || die "Expected SE mapchunks but found none"
  for sec in "${se_chunks[@]}"; do
    base=$(basename "$sec")
    cid="${base#${SAMPLE}.SE_}"
    cid="${cid%.fastq.gz}"
    outbam="$BAMS/${SAMPLE}.MAP.SE.${cid}.bam"

    if [[ "$LIBTYPE" == "modern" ]]; then
      log "  Mapping SE chunk ${cid} (mem)"
      map_se_mem_chunk "$sec" "$outbam" "mem.SE.${cid}"
    else
      log "  Mapping SE chunk ${cid} (aln/samse)"
      map_se_aln_chunk "$sec" "$outbam" "aln.SE.${cid}"
    fi
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
MARKDUP_BAM="$FINAL/${SAMPLE}.markdup_flagged.bam"
DEDUP_BAM="$FINAL/${SAMPLE}.dedup.bam"

log "STEP: dedup (mark duplicates, then filter; duprate from duplicate flags)"

if [[ ( "$LIBTYPE" == "modern" || "$LIBTYPE" == "historical" ) && "$SEQ_MODE" != "SE" ]]; then
  "$SAMTOOLS" sort -n -@ "$SORT_THREADS" -o "$FINAL/${SAMPLE}.namesort.bam" "$MERGED_BAM"
  "$SAMTOOLS" fixmate -m -@ "$SORT_THREADS" "$FINAL/${SAMPLE}.namesort.bam" "$FINAL/${SAMPLE}.fixmate.bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$FINAL/${SAMPLE}.coordsort.bam" "$FINAL/${SAMPLE}.fixmate.bam"
  rm -f "$FINAL/${SAMPLE}.namesort.bam" "$FINAL/${SAMPLE}.fixmate.bam"
else
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$FINAL/${SAMPLE}.coordsort.bam" "$MERGED_BAM"
fi

"$SAMTOOLS" markdup -s -@ "$SORT_THREADS" "$FINAL/${SAMPLE}.coordsort.bam" "$MARKDUP_BAM"
rm -f "$FINAL/${SAMPLE}.coordsort.bam"
file_nonempty "$MARKDUP_BAM" || die "markdup failed: $MARKDUP_BAM"

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
  awk 'BEGIN{wsum=0;len=0}
       NR==1{next}
       {ref_len=($3-$2+1); wsum += $7*ref_len; len += ref_len}
       END{if(len>0) printf "%.6f", wsum/len; else print 0}' \
    "$COV_TSV"
)

PCT_COVERED=$(
  awk 'BEGIN{wsum=0;len=0}
       NR==1{next}
       {ref_len=($3-$2+1); wsum += $6*ref_len; len += ref_len}
       END{if(len>0) printf "%.6f", wsum/len; else print 0}' \
    "$COV_TSV"
)

###############################################################################
# SUMMARY STATS
###############################################################################
log "STEP: summary statistics"

MAPPED_READS_ALL=$("$SAMTOOLS" view -c -F 2308 "$MERGED_BAM")
MAPPED_READS_UNIQUE=$("$SAMTOOLS" view -c -F 2308 "$DEDUP_BAM")

MAPPED_FRAGMENTS_ALL=$("$SAMTOOLS" view -F 2308 "$MARKDUP_BAM" \
  | awk 'BEGIN{c=0}
         function hasbit(x,b){return and(x,b)}
         {flag=$2;
          paired=hasbit(flag,1);
          read1=hasbit(flag,64);
          if(paired){ if(read1) c++ } else { c++ }
         }
         END{print c}'
)

DUP_FRAGMENTS=$("$SAMTOOLS" view -F 2308 "$MARKDUP_BAM" \
  | awk 'BEGIN{c=0}
         function hasbit(x,b){return and(x,b)}
         {flag=$2;
          dup=hasbit(flag,1024);
          paired=hasbit(flag,1);
          read1=hasbit(flag,64);
          if(dup){
            if(paired){ if(read1) c++ } else { c++ }
          }
         }
         END{print c}'
)

MAPPED_FRAGMENTS_UNIQUE="$(py_sub_nonneg "$MAPPED_FRAGMENTS_ALL" "$DUP_FRAGMENTS")"

DUPRATE=$("$PYTHON" - "$DUP_FRAGMENTS" "$MAPPED_FRAGMENTS_ALL" <<'PY'
import sys
d=float(sys.argv[1]); a=float(sys.argv[2])
print("0.000000" if a<=0 else f"{(d/a):.6f}")
PY
)

AVG_READLEN=$(
  "$SAMTOOLS" view -F 2308 "$DEDUP_BAM" \
  | awk '{l=length($10); if(l>0){sum+=l;n++}} END{ if(n>0) printf "%.2f", sum/n; else printf "0.00"}'
)

MAPPED_BP=$(awk -v n="$MAPPED_READS_ALL" -v l="$AVG_READLEN" 'BEGIN{printf "%.0f", n*l}')

ENDOG=$("$PYTHON" - "$MAPPED_FRAGMENTS_UNIQUE" "$RAW_FRAGMENTS" <<'PY'
import sys
n=float(sys.argv[1]); d=float(sys.argv[2])
print("0" if d==0 else f"{(n/d):.6f}")
PY
)

{
  printf "sample\tlibrary_type\tseq_mode\tpilot_fragments\tmax_reads_per_chunk\tmapq\traw_R1_reads\traw_R2_reads\traw_fragments\ttrimmed_fragments\tmerged_reads\tunpaired_reads\tmapped_reads_all\tmapped_reads_unique\tmapped_fragments_all\tmapped_fragments_unique\tendog\tduprate\tavg_readlen\tmapped_bp\tavg_depth\tpct_covered\n"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$SAMPLE" "$LIBTYPE" "$SEQ_MODE" "$PILOT_FRAGMENTS" "$MAX_READS_PER_CHUNK" "$MAPQ" \
    "$RAW_R1_READS" "$RAW_R2_READS" "$RAW_FRAGMENTS" \
    "$TRIMMED_FRAGMENTS" "$MERGED_READS" "$UNPAIRED_READS" \
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
