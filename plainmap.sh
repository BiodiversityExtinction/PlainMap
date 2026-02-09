#!/usr/bin/env bash
# PlainMap v0.1
# Transparent, failure-aware mapping pipeline for ancient and modern DNA
#
# Resume/robustness + space-saving (implemented):
# - Run signature guard: refuses to "resume" with mismatched settings/manifest unless --no-resume or --reset
# - Pre-fastp chunking checkpoint
# - fastp is chunk-resumable via per-chunk "parts" outputs (no re-trimming on restart)
# - pooling is deterministic + atomic (tmp then rename)
# - mapping is chunk-resumable: existing BAMs are skipped only if samtools quickcheck passes
# - mapchunks are preserved under resume and deleted "as you go" after each chunk maps successfully (Option B)
# - pre-fastp chunks are deleted once pooling is done (space)
# - per-chunk BAMs are deleted once merge succeeds (unless --keep-intermediate)
# - samtools merge uses -f and is atomic (survives partial merge on crash)
# - atomic outputs for BAMs/TSVs (tmp then rename)
# - BWA index lock has stale detection + timeout + writable-dir check
# - temp cleanup trap removes *.tmp.* in TMP only (safe)
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
MAX_READS_PER_CHUNK=0      # 0 disables chunking safety valve
PILOT_FRAGMENTS=0          # 0 disables pilot; otherwise PER-UNIT cap on fragments BEFORE fastp
MAPQ=20                    # mapping quality threshold (samtools view -q)

ADAPTER_R1=""
ADAPTER_R2=""
TRIM_ONLY=0

RUN_MAPDAMAGE=0

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
  --max-reads-per-chunk INT                  (default: 0; disabled)
  --pilot-fragments INT                      (default: 0; disabled) PER-UNIT cap
  --adapter-r1 SEQ
  --adapter-r2 SEQ
  --trim-only                                Trim only (fastp); exit before mapping
  --run-mapdamage                            Run mapDamage (any library type; default OFF)
  --no-mapdamage                             Do not run mapDamage (default)
  --tmpdir DIR                               Custom temp dir (e.g. node-local scratch)
  --keep-intermediate                        Keep <outdir>/<prefix>/work (disables some space-saving deletions)
  --resume | --no-resume                     (default: resume)
  --dry-run                                  Print plan only
  --validate                                 Check tools + gzip/pigz -t all manifest FASTQs, then exit
  --reset                                    Clear checkpoints AND remove intermediate outputs for this sample

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
    --run-mapdamage) RUN_MAPDAMAGE=1; shift ;;
    --no-mapdamage) RUN_MAPDAMAGE=0; shift ;;
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
PARTS="$RAW/parts"             # resumable per-chunk fastp outputs

STATS="$OUT/${SAMPLE}.plainmap.stats.tsv"
COV_TSV="$OUT/${SAMPLE}.plainmap.coverage.tsv"

if [[ -n "$TMPDIR_USER" ]]; then
  TMP="$TMPDIR_USER/plainmap_${SAMPLE}"
else
  TMP="$TMPBASE"
fi

mkdir -p "$RAW" "$CHUNKS" "$MAPCHUNKS" "$BAMS" "$FINAL" "$CKPT" "$TMP" "$REPORTS" "$PARTS"

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
# TEMP CLEANUP TRAP (SAFE: TMP ONLY)
###############################################################################
cleanup_tmp() {
  rm -f "$TMP"/*.tmp.* 2>/dev/null || true
  rm -f "$TMP"/*.sam "$TMP"/*.bam "$TMP"/*.bam.tmp "$TMP"/*.sai 2>/dev/null || true
}
trap cleanup_tmp EXIT

###############################################################################
# FILE CHECKS
###############################################################################
file_nonempty() { [[ -e "$1" && -s "$1" ]]; }
gz_ok() { "${GZ_TEST[@]}" "$1" >/dev/null 2>&1; }

###############################################################################
# CHECKPOINTS
###############################################################################
ckpt_file() { echo "$CKPT/${SAMPLE}.$1.done"; }
ckpt_mark() { [[ $RESUME -eq 1 ]] && date -Is > "$(ckpt_file "$1")"; }
ckpt_clear() { rm -f "$(ckpt_file "$1")" 2>/dev/null || true; }
files_ok() { local f; for f in "$@"; do file_nonempty "$f" || return 1; done; return 0; }
ckpt_ok() { local step="$1"; shift; [[ $RESUME -eq 1 ]] && [[ -f "$(ckpt_file "$step")" ]] && files_ok "$@"; }

###############################################################################
# RUN SIGNATURE (GUARD AGAINST RESUMING WITH DIFFERENT SETTINGS)
###############################################################################
SIG_FILE="$CKPT/${SAMPLE}.run_signature.txt"
compute_signature() {
  "$PYTHON" - "$MANIFEST" "$REF" <<'PY'
import hashlib,sys,os
manifest=sys.argv[1]
ref=sys.argv[2]
h=hashlib.sha256()

def upd(s):
    h.update(s.encode('utf-8',errors='ignore'))
    h.update(b'\n')

# Include key settings from environment passed via bash heredoc substitution:
# We'll read them from environment variables (set in bash before calling).
keys = [
  "PLAINMAP_VERSION","LIBTYPE","THREADS","MINLEN","MISMATCH","MAPQ",
  "MAX_READS_PER_CHUNK","PILOT_FRAGMENTS","ADAPTER_R1","ADAPTER_R2","RUN_MAPDAMAGE",
  "FASTP","BWA","SAMTOOLS","MAPDAMAGE","PYTHON"
]
for k in keys:
    upd(f"{k}={os.environ.get(k,'')}")
upd(f"MANIFEST_REALPATH={os.path.realpath(os.path.abspath(manifest))}")
upd(f"REF_REALPATH={os.path.realpath(os.path.abspath(ref))}")

# Hash manifest contents (paths order matters!)
with open(manifest,'rb') as f:
    h.update(f.read())

print(h.hexdigest())
PY
}

export PLAINMAP_VERSION="$VERSION"
export LIBTYPE THREADS MINLEN MISMATCH MAPQ MAX_READS_PER_CHUNK PILOT_FRAGMENTS ADAPTER_R1 ADAPTER_R2 RUN_MAPDAMAGE FASTP BWA SAMTOOLS MAPDAMAGE PYTHON

CURRENT_SIG="$(compute_signature)"

if [[ $RESET -eq 1 ]]; then
  log "Reset requested: clearing checkpoints and removing intermediate outputs for sample $SAMPLE"
  rm -f "$CKPT/${SAMPLE}."*.done 2>/dev/null || true
  rm -f "$SIG_FILE" 2>/dev/null || true
  rm -rf "$WORK" "$OUT/${SAMPLE}/fastp_reports" 2>/dev/null || true
  mkdir -p "$RAW" "$CHUNKS" "$MAPCHUNKS" "$BAMS" "$FINAL" "$CKPT" "$TMP" "$REPORTS" "$PARTS"
fi

if [[ $RESUME -eq 1 ]]; then
  if [[ -f "$SIG_FILE" ]]; then
    PREV_SIG="$(cat "$SIG_FILE" || true)"
    if [[ -n "$PREV_SIG" && "$PREV_SIG" != "$CURRENT_SIG" ]]; then
      die "Resume refused: run signature mismatch (settings/manifest/ref changed). Use --no-resume or --reset."
    fi
  else
    echo "$CURRENT_SIG" > "$SIG_FILE"
  fi
else
  # no-resume: overwrite signature so future resumes use new settings
  echo "$CURRENT_SIG" > "$SIG_FILE"
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
  if [[ "$RUN_MAPDAMAGE" -eq 1 ]]; then
    require_cmd "$MAPDAMAGE"
  fi
  # quickcheck exists in samtools (new-ish but widely available)
  "$SAMTOOLS" quickcheck -h >/dev/null 2>&1 || die "samtools quickcheck not available; need a newer samtools"
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
  log "  pilot_fragments:     $PILOT_FRAGMENTS (PER-UNIT cap)"
  log "  max_reads_per_chunk: $MAX_READS_PER_CHUNK"
  log "  run_mapdamage:       $RUN_MAPDAMAGE"
  log "  resume:              $RESUME"
  log "Plan: manifest -> build units -> pre-fastp chunk (ckpt) -> fastp parts (resumable) -> pool (ckpt) -> mapchunk -> map (chunk-resume + delete-as-you-go) -> merge -f (atomic) -> dedup -> RG -> mapDamage (opt) -> coverage -> stats"
  exit 0
fi

quick_preflight

###############################################################################
# BWA INDEX (PARALLEL SAFE + STALE LOCK + TIMEOUT + WRITABLE CHECK)
###############################################################################
ensure_bwa_index() {
  local ref="$1"
  local lock="${ref}.bwa.lock"
  local idx=( "${ref}.bwt" "${ref}.pac" "${ref}.ann" "${ref}.amb" "${ref}.sa" )
  local missing=0
  for f in "${idx[@]}"; do [[ -f "$f" ]] || missing=1; done
  [[ "$missing" -eq 0 ]] && { log "Reference index: present"; return; }

  local refdir
  refdir="$(dirname "$ref")"
  [[ -w "$refdir" ]] || die "Reference directory not writable: $refdir. Copy reference to writable location or pre-index."

  if mkdir "$lock" 2>/dev/null; then
    trap 'rm -rf "$lock"' EXIT
    log "Indexing reference"
    "$BWA" index "$ref"
    rm -rf "$lock"; trap - EXIT
  else
    log "Waiting for BWA index to finish (lock exists: $lock)"

    local lock_mtime=0
    lock_mtime=$(stat -c %Y "$lock" 2>/dev/null || echo 0)
    local now
    now=$(date +%s)
    local age=$(( now - lock_mtime ))
    if [[ "$lock_mtime" -gt 0 && "$age" -gt 7200 ]]; then
      die "BWA index lock appears stale (>2h): $lock. If no other job is indexing, remove it: rm -rf '$lock'"
    fi

    # Wait up to 2 hours (720 * 10s)
    for _i in {1..720}; do
      sleep 10
      missing=0
      for f in "${idx[@]}"; do [[ -f "$f" ]] || missing=1; done
      [[ "$missing" -eq 0 ]] && break
    done

    missing=0
    for f in "${idx[@]}"; do [[ -f "$f" ]] || missing=1; done
    [[ "$missing" -eq 0 ]] || die "Timed out waiting for BWA index. Lock: $lock"
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
first_header_line() { "${ZCAT[@]}" "$1" | head -n 1 || true; }

normalize_qname() {
  local hdr="$1"
  local first="${hdr%% *}"
  first="${first#@}"
  first="${first%/1}"
  first="${first%/2}"
  echo "$first"
}

classify_fastq_direction() {
  local fq="$1"
  local hdr firsttok rest secondtok
  hdr="$(first_header_line "$fq")"
  [[ -n "$hdr" ]] || { echo "UNKNOWN"; return; }

  if [[ "$hdr" == *" 1:"* ]]; then echo "R1"; return; fi
  if [[ "$hdr" == *" 2:"* ]]; then echo "R2"; return; fi

  firsttok="${hdr%% *}"
  rest="${hdr#"$firsttok"}"; rest="${rest# }"
  secondtok="${rest%% *}"

  if [[ "$secondtok" == */1 ]]; then echo "R1"; return; fi
  if [[ "$secondtok" == */2 ]]; then echo "R2"; return; fi
  if [[ "$firsttok" == */1 ]]; then echo "R1"; return; fi
  if [[ "$firsttok" == */2 ]]; then echo "R2"; return; fi
  if [[ "$secondtok" == "1" ]]; then echo "R1"; return; fi
  if [[ "$secondtok" == "2" ]]; then echo "R2"; return; fi

  echo "UNKNOWN"
}

check_headers_identical_sra() {
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
# BUILD UNITS SAFELY (SRA unknown handled)
###############################################################################
declare -A R1_BY_KEY
declare -A R2_BY_KEY
declare -A UNK_BY_KEY

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

  dir="$(classify_fastq_direction "$f")"
  key="$(normalize_qname "$hdr")"
  [[ -n "$key" ]] || die "Could not derive header key for: $f"

  if [[ "$dir" == "R1" ]]; then
    [[ -z "${R1_BY_KEY[$key]:-}" ]] || die "Multiple R1 files share key '$key'. Offending: $f and ${R1_BY_KEY[$key]}"
    R1_BY_KEY["$key"]="$f"
  elif [[ "$dir" == "R2" ]]; then
    [[ -z "${R2_BY_KEY[$key]:-}" ]] || die "Multiple R2 files share key '$key'. Offending: $f and ${R2_BY_KEY[$key]}"
    R2_BY_KEY["$key"]="$f"
  else
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

for k in "${!R2_BY_KEY[@]}"; do
  [[ -n "${R1_BY_KEY[$k]:-}" ]] || die "Found R2 without matching R1 for key '$k': ${R2_BY_KEY[$k]}"
done

for k in "${!R1_BY_KEY[@]}"; do
  if [[ -n "${R2_BY_KEY[$k]:-}" ]]; then
    PE_KEYS+=( "$k" )
  else
    SE_FILES+=( "${R1_BY_KEY[$k]}" )
  fi
done

for k in "${!UNK_BY_KEY[@]}"; do
  mapfile -t entries <<< "${UNK_BY_KEY[$k]}"
  if [[ "${#entries[@]}" -eq 1 ]]; then
    f="${entries[0]#*:}"
    SE_FILES+=( "$f" )
  elif [[ "${#entries[@]}" -eq 2 ]]; then
    e1="${entries[0]}"; e2="${entries[1]}"
    i1="${e1%%:*}"; f1="${e1#*:}"
    i2="${e2%%:*}"; f2="${e2#*:}"
    if [[ "$i2" -lt "$i1" ]]; then
      tmpi="$i1"; tmpf="$f1"
      i1="$i2"; f1="$f2"
      i2="$tmpi"; f2="$tmpf"
    fi

    check_headers_identical_sra "$f1" "$f2"
    n1=$(count_reads_fastq_gz "$f1")
    n2=$(count_reads_fastq_gz "$f2")
    [[ "$n1" -eq "$n2" ]] || die "SRA-style PE sanity check failed: read counts differ for key '$k' (R1=$n1 R2=$n2):
$f1
$f2"

    [[ -z "${R1_BY_KEY[$k]:-}" && -z "${R2_BY_KEY[$k]:-}" ]] || die "Internal error: UNKNOWN key '$k' already classified"
    R1_BY_KEY["$k"]="$f1"
    R2_BY_KEY["$k"]="$f2"
    PE_KEYS+=( "$k" )
  else
    die "Ambiguous SRA-style grouping for key '$k': found ${#entries[@]} files with no mate markers."
  fi
done

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
# CREATE PRE-FASTP CHUNKS PER UNIT (PE stays paired) [ckpt: chunk]
###############################################################################
make_chunks() {
  log "STEP: pre-fastp chunking"
  ckpt_clear chunk
  rm -f "$CHUNKS/${SAMPLE}."*.fastq.gz 2>/dev/null || true

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
# FASTP PARTS (RESUMABLE) + POOLING (ATOMIC) [ckpt: fastp_pool]
###############################################################################
POOL_PE_R1="$RAW/${SAMPLE}.pool.PE.R1.fastq.gz"
POOL_PE_R2="$RAW/${SAMPLE}.pool.PE.R2.fastq.gz"
POOL_SE_ALL="$RAW/${SAMPLE}.pool.SE.all.fastq.gz"

RAW_R1_READS=0
RAW_R2_READS=0
RAW_FRAGMENTS=0
TRIMMED_FRAGMENTS=0
MERGED_READS=0
UNPAIRED_READS=0

declare -A PILOT_LEFT_BY_UNIT

fastp_parts_exist_ok_pe() {
  local tag="$1"
  local p1="$PARTS/${SAMPLE}.${tag}.trim.R1.fastq.gz"
  local p2="$PARTS/${SAMPLE}.${tag}.trim.R2.fastq.gz"
  local pj="$REPORTS/${SAMPLE}.${tag}.fastp.json"
  file_nonempty "$p1" && file_nonempty "$p2" && file_nonempty "$pj" && gz_ok "$p1" && gz_ok "$p2"
}

fastp_parts_exist_ok_se() {
  local tag="$1"
  local p1="$PARTS/${SAMPLE}.${tag}.trim.SE.fastq.gz"
  local pj="$REPORTS/${SAMPLE}.${tag}.fastp.json"
  file_nonempty "$p1" && file_nonempty "$pj" && gz_ok "$p1"
}

pool_from_parts_atomic() {
  log "  Pooling from parts (atomic)"

  local tmp_pe1="$POOL_PE_R1.tmp.$$"
  local tmp_pe2="$POOL_PE_R2.tmp.$$"
  local tmp_se="$POOL_SE_ALL.tmp.$$"

  : > "$tmp_pe1"
  : > "$tmp_pe2"
  : > "$tmp_se"

  shopt -s nullglob
  # Deterministic order: sort by filename
  local pe_r1_parts=( "$PARTS/${SAMPLE}."*.trim.R1.fastq.gz )
  local pe_r2_parts=( "$PARTS/${SAMPLE}."*.trim.R2.fastq.gz )
  local se_parts=( "$PARTS/${SAMPLE}."*.trim.SE.fastq.gz )
  local merged_parts=( "$PARTS/${SAMPLE}."*.merged.fastq.gz )
  local up1_parts=( "$PARTS/${SAMPLE}."*.unpaired1.fastq.gz )
  local up2_parts=( "$PARTS/${SAMPLE}."*.unpaired2.fastq.gz )

  # Sort in-place (bash arrays -> printf/sort -> mapfile)
  if [[ ${#pe_r1_parts[@]} -gt 0 ]]; then mapfile -t pe_r1_parts < <(printf "%s\n" "${pe_r1_parts[@]}" | sort); fi
  if [[ ${#pe_r2_parts[@]} -gt 0 ]]; then mapfile -t pe_r2_parts < <(printf "%s\n" "${pe_r2_parts[@]}" | sort); fi
  if [[ ${#se_parts[@]} -gt 0 ]]; then mapfile -t se_parts < <(printf "%s\n" "${se_parts[@]}" | sort); fi
  if [[ ${#merged_parts[@]} -gt 0 ]]; then mapfile -t merged_parts < <(printf "%s\n" "${merged_parts[@]}" | sort); fi
  if [[ ${#up1_parts[@]} -gt 0 ]]; then mapfile -t up1_parts < <(printf "%s\n" "${up1_parts[@]}" | sort); fi
  if [[ ${#up2_parts[@]} -gt 0 ]]; then mapfile -t up2_parts < <(printf "%s\n" "${up2_parts[@]}" | sort); fi

  for f in "${pe_r1_parts[@]}"; do cat "$f" >> "$tmp_pe1"; done
  for f in "${pe_r2_parts[@]}"; do cat "$f" >> "$tmp_pe2"; done

  for f in "${se_parts[@]}"; do cat "$f" >> "$tmp_se"; done
  for f in "${merged_parts[@]}"; do file_nonempty "$f" && cat "$f" >> "$tmp_se" || true; done
  for f in "${up1_parts[@]}"; do file_nonempty "$f" && cat "$f" >> "$tmp_se" || true; done
  for f in "${up2_parts[@]}"; do file_nonempty "$f" && cat "$f" >> "$tmp_se" || true; done

  mv -f "$tmp_pe1" "$POOL_PE_R1"
  mv -f "$tmp_pe2" "$POOL_PE_R2"
  mv -f "$tmp_se"  "$POOL_SE_ALL"

  # Validate pools (may be empty if no data of that type)
  [[ ! -s "$POOL_PE_R1" ]] || gz_ok "$POOL_PE_R1" || die "Pooled PE R1 gzip corrupt: $POOL_PE_R1"
  [[ ! -s "$POOL_PE_R2" ]] || gz_ok "$POOL_PE_R2" || die "Pooled PE R2 gzip corrupt: $POOL_PE_R2"
  [[ ! -s "$POOL_SE_ALL" ]] || gz_ok "$POOL_SE_ALL" || die "Pooled SE gzip corrupt: $POOL_SE_ALL"
}

run_fastp_parts_and_pool() {
  log "STEP: fastp parts (resumable) + pooling"

  # Recompute pilot state deterministically per run (pilot is per-unit)
  declare -A PILOT_LEFT_BY_UNIT
  RAW_R1_READS=0; RAW_R2_READS=0; RAW_FRAGMENTS=0
  TRIMMED_FRAGMENTS=0; MERGED_READS=0; UNPAIRED_READS=0

  shopt -s nullglob

  # PE parts
  local pe_r1_chunks=( "$CHUNKS/${SAMPLE}.PE"*".R1_"*.fastq.gz )
  for r1c in "${pe_r1_chunks[@]}"; do
    local base r2c unit_id unit_left pairs_in_chunk take_pairs
    base=$(basename "$r1c")
    r2c="${r1c/.R1_/.R2_}"
    [[ -f "$r2c" ]] || die "Missing PE mate chunk for: $r1c (expected $r2c)"

    unit_id="${base%.fastq.gz}"
    unit_id="${unit_id%%.R1_*}"

    if [[ "$PILOT_FRAGMENTS" -gt 0 && -z "${PILOT_LEFT_BY_UNIT[$unit_id]:-}" ]]; then
      PILOT_LEFT_BY_UNIT["$unit_id"]="$PILOT_FRAGMENTS"
    fi
    unit_left="${PILOT_LEFT_BY_UNIT[$unit_id]:-0}"
    if [[ "$PILOT_FRAGMENTS" -gt 0 && "$unit_left" -le 0 ]]; then
      continue
    fi

    pairs_in_chunk=$(count_reads_fastq_gz "$r1c")
    take_pairs="$pairs_in_chunk"

    local in_r1="$r1c"
    local in_r2="$r2c"

    if [[ "$PILOT_FRAGMENTS" -gt 0 ]]; then
      take_pairs=$(pilot_take "$pairs_in_chunk" "$unit_left")
      if [[ "$take_pairs" -le 0 ]]; then
        PILOT_LEFT_BY_UNIT["$unit_id"]=0
        continue
      fi
      if [[ "$take_pairs" -lt "$pairs_in_chunk" ]]; then
        log "  Pilot (per-unit): limiting PE unit $unit_id chunk $(basename "$r1c") to first $take_pairs pairs"
      fi
      local lim_r1="$TMP/${SAMPLE}.$(basename "${r1c%.fastq.gz}").pilot.R1.fastq.gz"
      local lim_r2="$TMP/${SAMPLE}.$(basename "${r2c%.fastq.gz}").pilot.R2.fastq.gz"
      limit_fastq_gz "$in_r1" "$lim_r1" "$take_pairs"
      limit_fastq_gz "$in_r2" "$lim_r2" "$take_pairs"
      in_r1="$lim_r1"
      in_r2="$lim_r2"
      unit_left="$(py_sub_nonneg "$unit_left" "$take_pairs")"
      PILOT_LEFT_BY_UNIT["$unit_id"]="$unit_left"
    fi

    RAW_R1_READS=$(py_add "$RAW_R1_READS" "$take_pairs")
    RAW_R2_READS=$(py_add "$RAW_R2_READS" "$take_pairs")
    RAW_FRAGMENTS=$(py_add "$RAW_FRAGMENTS" "$take_pairs")

    local tag
    tag="$(basename "${r1c%.fastq.gz}")"
    tag="${tag//\//_}"

    if [[ $RESUME -eq 1 ]] && fastp_parts_exist_ok_pe "$tag"; then
      log "  fastp PE $tag (skipped; parts exist)"
    else
      log "  fastp PE $tag"
      local fastp_err="$TMP/fastp_${SAMPLE}.${tag}.stderr.log"
      local fp_json="$REPORTS/${SAMPLE}.${tag}.fastp.json"

      local out_trim1="$PARTS/${SAMPLE}.${tag}.trim.R1.fastq.gz"
      local out_trim2="$PARTS/${SAMPLE}.${tag}.trim.R2.fastq.gz"
      local out_merged="$PARTS/${SAMPLE}.${tag}.merged.fastq.gz"
      local out_up1="$PARTS/${SAMPLE}.${tag}.unpaired1.fastq.gz"
      local out_up2="$PARTS/${SAMPLE}.${tag}.unpaired2.fastq.gz"

      # Write to tmp then rename (atomic-ish)
      local t_trim1="${out_trim1}.tmp.$$"
      local t_trim2="${out_trim2}.tmp.$$"
      local t_merged="${out_merged}.tmp.$$"
      local t_up1="${out_up1}.tmp.$$"
      local t_up2="${out_up2}.tmp.$$"
      local t_json="${fp_json}.tmp.$$"

      local extra_pe=()
      if [[ -z "$ADAPTER_R1" && -z "$ADAPTER_R2" ]]; then
        extra_pe+=( --detect_adapter_for_pe )
      fi

      run_fastp "$fastp_err" "${FASTP_ARGS[@]}" "${extra_pe[@]}" \
        -i "$in_r1" -I "$in_r2" \
        --out1 "$t_trim1" --out2 "$t_trim2" \
        -m --merged_out "$t_merged" \
        --unpaired1 "$t_up1" --unpaired2 "$t_up2" \
        -j "$t_json"

      file_nonempty "$t_trim1" || die "fastp produced empty trimmed R1 for $tag"
      file_nonempty "$t_trim2" || die "fastp produced empty trimmed R2 for $tag"
      gz_ok "$t_trim1" || die "Corrupt trimmed R1 gzip for $tag"
      gz_ok "$t_trim2" || die "Corrupt trimmed R2 gzip for $tag"

      mv -f "$t_trim1" "$out_trim1"
      mv -f "$t_trim2" "$out_trim2"
      mv -f "$t_merged" "$out_merged" 2>/dev/null || true
      mv -f "$t_up1" "$out_up1" 2>/dev/null || true
      mv -f "$t_up2" "$out_up2" 2>/dev/null || true
      mv -f "$t_json" "$fp_json"
    fi

    # update counts from json (whether we skipped or ran)
    local fp_json="$REPORTS/${SAMPLE}.${tag}.fastp.json"
    local out1_reads merged_reads
    out1_reads=$(json_get_int_any "$fp_json" "read1_after_filtering.total_reads")
    [[ "$out1_reads" -le 0 ]] && out1_reads=$(count_reads_fastq_gz "$PARTS/${SAMPLE}.${tag}.trim.R1.fastq.gz")
    TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$out1_reads")

    merged_reads=$(json_get_int_any "$fp_json" "merged_and_filtered.total_reads")
    if [[ "$merged_reads" -le 0 && -s "$PARTS/${SAMPLE}.${tag}.merged.fastq.gz" ]]; then
      merged_reads=$(count_reads_fastq_gz "$PARTS/${SAMPLE}.${tag}.merged.fastq.gz")
    fi
    MERGED_READS=$(py_add "$MERGED_READS" "$merged_reads")
    TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$merged_reads")

    local up1n=0 up2n=0
    [[ -s "$PARTS/${SAMPLE}.${tag}.unpaired1.fastq.gz" ]] && up1n=$(count_reads_fastq_gz "$PARTS/${SAMPLE}.${tag}.unpaired1.fastq.gz")
    [[ -s "$PARTS/${SAMPLE}.${tag}.unpaired2.fastq.gz" ]] && up2n=$(count_reads_fastq_gz "$PARTS/${SAMPLE}.${tag}.unpaired2.fastq.gz")
    UNPAIRED_READS=$(py_add "$UNPAIRED_READS" "$up1n")
    UNPAIRED_READS=$(py_add "$UNPAIRED_READS" "$up2n")
    TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$up1n")
    TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$up2n")

    # cleanup pilot-limited inputs
    if [[ "$PILOT_FRAGMENTS" -gt 0 ]]; then
      rm -f "${lim_r1:-}" "${lim_r2:-}" 2>/dev/null || true
    fi
  done

  # SE parts
  local se_chunks=( "$CHUNKS/${SAMPLE}.SE"*".SE_"*.fastq.gz )
  for sec in "${se_chunks[@]}"; do
    local base unit_id unit_left reads_in_chunk take_reads
    base=$(basename "$sec")
    unit_id="${base%.fastq.gz}"
    unit_id="${unit_id%%.SE_*}"

    if [[ "$PILOT_FRAGMENTS" -gt 0 && -z "${PILOT_LEFT_BY_UNIT[$unit_id]:-}" ]]; then
      PILOT_LEFT_BY_UNIT["$unit_id"]="$PILOT_FRAGMENTS"
    fi
    unit_left="${PILOT_LEFT_BY_UNIT[$unit_id]:-0}"
    if [[ "$PILOT_FRAGMENTS" -gt 0 && "$unit_left" -le 0 ]]; then
      continue
    fi

    reads_in_chunk=$(count_reads_fastq_gz "$sec")
    take_reads="$reads_in_chunk"

    local in_se="$sec"
    if [[ "$PILOT_FRAGMENTS" -gt 0 ]]; then
      take_reads=$(pilot_take "$reads_in_chunk" "$unit_left")
      if [[ "$take_reads" -le 0 ]]; then
        PILOT_LEFT_BY_UNIT["$unit_id"]=0
        continue
      fi
      if [[ "$take_reads" -lt "$reads_in_chunk" ]]; then
        log "  Pilot (per-unit): limiting SE unit $unit_id chunk $(basename "$sec") to first $take_reads reads"
      fi
      local lim_se="$TMP/${SAMPLE}.$(basename "${sec%.fastq.gz}").pilot.SE.fastq.gz"
      limit_fastq_gz "$in_se" "$lim_se" "$take_reads"
      in_se="$lim_se"
      unit_left="$(py_sub_nonneg "$unit_left" "$take_reads")"
      PILOT_LEFT_BY_UNIT["$unit_id"]="$unit_left"
    fi

    RAW_R1_READS=$(py_add "$RAW_R1_READS" "$take_reads")
    RAW_FRAGMENTS=$(py_add "$RAW_FRAGMENTS" "$take_reads")

    local tag
    tag="$(basename "${sec%.fastq.gz}")"
    tag="${tag//\//_}"

    if [[ $RESUME -eq 1 ]] && fastp_parts_exist_ok_se "$tag"; then
      log "  fastp SE $tag (skipped; parts exist)"
    else
      log "  fastp SE $tag"
      local fastp_err="$TMP/fastp_${SAMPLE}.${tag}.stderr.log"
      local fp_json="$REPORTS/${SAMPLE}.${tag}.fastp.json"

      local out_trim="$PARTS/${SAMPLE}.${tag}.trim.SE.fastq.gz"
      local t_trim="${out_trim}.tmp.$$"
      local t_json="${fp_json}.tmp.$$"

      run_fastp "$fastp_err" "${FASTP_ARGS[@]}" \
        -i "$in_se" --out1 "$t_trim" \
        -j "$t_json"

      file_nonempty "$t_trim" || die "fastp produced empty trimmed SE for $tag"
      gz_ok "$t_trim" || die "Corrupt trimmed SE gzip for $tag"

      mv -f "$t_trim" "$out_trim"
      mv -f "$t_json" "$fp_json"
    fi

    local fp_json="$REPORTS/${SAMPLE}.${tag}.fastp.json"
    local se_after
    se_after=$(json_get_int_any "$fp_json" "read1_after_filtering.total_reads" "summary.after_filtering.total_reads")
    [[ "$se_after" -le 0 ]] && se_after=$(count_reads_fastq_gz "$PARTS/${SAMPLE}.${tag}.trim.SE.fastq.gz")
    TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$se_after")

    if [[ "$PILOT_FRAGMENTS" -gt 0 ]]; then
      rm -f "${lim_se:-}" 2>/dev/null || true
    fi
  done

  # Pool (atomic)
  pool_from_parts_atomic

  ckpt_mark fastp_pool

  # Space saving: pre-fastp chunks no longer needed after pooling checkpoint
  if [[ $KEEP_INTERMEDIATE -eq 0 ]]; then
    log "  Space: deleting pre-fastp chunks (pooling complete)"
    rm -f "$CHUNKS/${SAMPLE}."*.fastq.gz 2>/dev/null || true
  else
    log "  Keeping intermediates: not deleting pre-fastp chunks"
  fi
}

if ! ckpt_ok fastp_pool "$POOL_PE_R1" "$POOL_PE_R2" "$POOL_SE_ALL"; then
  run_fastp_parts_and_pool
else
  log "STEP: fastp parts + pooling (skipped; checkpoint present)"
  # We still need RAW_* and TRIMMED_* counters for stats. Recompute quickly from parts jsons.
  RAW_R1_READS=0; RAW_R2_READS=0; RAW_FRAGMENTS=0
  TRIMMED_FRAGMENTS=0; MERGED_READS=0; UNPAIRED_READS=0

  # RAW_* cannot be reconstructed from parts reliably without re-reading inputs/pilot;
  # to keep correctness, we recompute raw counts from fastp JSON when possible (fallback to 0).
  # This is acceptable because pilot-aware raw counts are the main reason; for strict accounting,
  # rerun with --no-resume or keep raw counts in a file. Here we persist them in ckpt metadata.
  RAW_META="$CKPT/${SAMPLE}.raw_counts.tsv"
  if [[ -f "$RAW_META" ]]; then
    read -r RAW_R1_READS RAW_R2_READS RAW_FRAGMENTS TRIMMED_FRAGMENTS MERGED_READS UNPAIRED_READS < "$RAW_META" || true
  else
    # Best-effort reconstruct from parts (trimmed counts)
    shopt -s nullglob
    for j in "$REPORTS/${SAMPLE}."*.fastp.json; do
      # if PE json has read1_after_filtering.total_reads, count it
      out1=$(json_get_int_any "$j" "read1_after_filtering.total_reads")
      if [[ "$out1" -gt 0 ]]; then
        TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$out1")
      else
        se_after=$(json_get_int_any "$j" "summary.after_filtering.total_reads" "read1_after_filtering.total_reads")
        TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$se_after")
      fi
      merged=$(json_get_int_any "$j" "merged_and_filtered.total_reads")
      MERGED_READS=$(py_add "$MERGED_READS" "$merged")
    done
    RAW_R1_READS=0; RAW_R2_READS=0; RAW_FRAGMENTS=0
    UNPAIRED_READS=0
  fi
fi

# Persist pilot-aware raw counters for reliable resume of stats (after a fresh run_fastp_parts_and_pool)
RAW_META="$CKPT/${SAMPLE}.raw_counts.tsv"
if [[ ! -f "$RAW_META" || $RESET -eq 1 || $RESUME -eq 0 ]]; then
  echo -e "${RAW_R1_READS}\t${RAW_R2_READS}\t${RAW_FRAGMENTS}\t${TRIMMED_FRAGMENTS}\t${MERGED_READS}\t${UNPAIRED_READS}" > "$RAW_META"
fi

###############################################################################
# TRIM-ONLY MODE
###############################################################################
if [[ $TRIM_ONLY -eq 1 ]]; then
  log "Trim-only mode enabled; writing stats and exiting"
  tmp_stats="${STATS}.tmp.$$"
  {
    printf "sample\tlibrary_type\tseq_mode\tpilot_fragments\tmax_reads_per_chunk\tmapq\traw_R1_reads\traw_R2_reads\traw_fragments\ttrimmed_fragments\tmerged_reads\tunpaired_reads\tmapped_reads_all\tmapped_reads_unique\tmapped_fragments_all\tmapped_fragments_unique\tendog\tduprate\tavg_readlen\tmapped_bp\tavg_depth\tpct_covered\n"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n" \
      "$SAMPLE" "$LIBTYPE" "$SEQ_MODE" "$PILOT_FRAGMENTS" "$MAX_READS_PER_CHUNK" "$MAPQ" \
      "$RAW_R1_READS" "$RAW_R2_READS" "$RAW_FRAGMENTS" \
      "$TRIMMED_FRAGMENTS" "$MERGED_READS" "$UNPAIRED_READS"
  } > "$tmp_stats"
  mv -f "$tmp_stats" "$STATS"

  PIPE_T1=$(date +%s)
  log "STATS written: $STATS"
  log "Total wall time: $((PIPE_T1 - PIPE_T0)) seconds"
  exit 0
fi

###############################################################################
# POST-FASTP MAPPING CHUNKING (OPTIONAL; REUSE MAX_READS_PER_CHUNK)
###############################################################################
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
  [[ "$n1" -eq "$n2" ]] || die "Pooled PE R1/R2 read counts differ (R1=$n1 R2=$n2)."

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
if [[ $RESUME -eq 1 ]]; then
  shopt -s nullglob
  existing_any=( "$MAPCHUNKS/${SAMPLE}."*.fastq.gz )
  if [[ ${#existing_any[@]} -gt 0 ]]; then
    log "Resume enabled: existing mapchunks found; not regenerating mapchunks"
  else
    make_mapchunks_pe "$POOL_PE_R1" "$POOL_PE_R2" "$MAPCHUNKS/${SAMPLE}.PE.R1_" "$MAPCHUNKS/${SAMPLE}.PE.R2_"
    make_mapchunks_se "$POOL_SE_ALL" "$MAPCHUNKS/${SAMPLE}.SE_"
  fi
else
  rm -f "$MAPCHUNKS/${SAMPLE}."*.fastq.gz 2>/dev/null || true
  make_mapchunks_pe "$POOL_PE_R1" "$POOL_PE_R2" "$MAPCHUNKS/${SAMPLE}.PE.R1_" "$MAPCHUNKS/${SAMPLE}.PE.R2_"
  make_mapchunks_se "$POOL_SE_ALL" "$MAPCHUNKS/${SAMPLE}.SE_"
fi

###############################################################################
# MAPPING (CHUNK-RESUMABLE + DELETE MAPCHUNKS AS YOU GO) [ckpt: map]
###############################################################################
map_se_mem_chunk() {
  local fq="$1" outbam="$2" tag="$3"
  local sam="$TMP/${SAMPLE}.${tag}.sam"
  local bam="$TMP/${SAMPLE}.${tag}.bam"
  "$BWA" mem -t "$THREADS" "$REF" "$fq" > "$sam"
  "$SAMTOOLS" view -q "$MAPQ" -u "$sam" > "$bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$outbam" "$bam"
  rm -f "$sam" "$bam"
}

map_pe_mem_chunk() {
  local fq1="$1" fq2="$2" outbam="$3" tag="$4"
  local sam="$TMP/${SAMPLE}.${tag}.sam"
  local bam="$TMP/${SAMPLE}.${tag}.bam"
  "$BWA" mem -t "$THREADS" "$REF" "$fq1" "$fq2" > "$sam"
  "$SAMTOOLS" view -q "$MAPQ" -u "$sam" > "$bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$outbam" "$bam"
  rm -f "$sam" "$bam"
}

map_se_aln_chunk() {
  local fq="$1" outbam="$2" tag="$3"
  local sai="$TMP/${SAMPLE}.${tag}.sai"
  local sam="$TMP/${SAMPLE}.${tag}.sam"
  local bam="$TMP/${SAMPLE}.${tag}.bam.tmp"
  "$BWA" aln -l 999 -n "$MISMATCH" -t "$THREADS" "$REF" "$fq" > "$sai"
  "$BWA" samse "$REF" "$sai" "$fq" > "$sam"
  "$SAMTOOLS" view -q "$MAPQ" -u "$sam" > "$bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$outbam" "$bam"
  rm -f "$sai" "$sam" "$bam"
}

map_pe_aln_chunk() {
  local fq1="$1" fq2="$2" outbam="$3" tag="$4"
  local sai1="$TMP/${SAMPLE}.${tag}.1.sai"
  local sai2="$TMP/${SAMPLE}.${tag}.2.sai"
  local sam="$TMP/${SAMPLE}.${tag}.sam"
  local bam="$TMP/${SAMPLE}.${tag}.bam.tmp"

  local h1 h2 n1x n2x
  h1="$(first_header_line "$fq1")"; h2="$(first_header_line "$fq2")"
  n1x="$(normalize_qname "$h1")"; n2x="$(normalize_qname "$h2")"
  [[ "$n1x" == "$n2x" ]] || die "Refusing sampe: PE chunks first read names differ: $(basename "$fq1")='$n1x' vs $(basename "$fq2")='$n2x'"

  "$BWA" aln -l 999 -n "$MISMATCH" -t "$THREADS" "$REF" "$fq1" > "$sai1"
  "$BWA" aln -l 999 -n "$MISMATCH" -t "$THREADS" "$REF" "$fq2" > "$sai2"
  "$BWA" sampe "$REF" "$sai1" "$sai2" "$fq1" "$fq2" > "$sam"
  "$SAMTOOLS" view -q "$MAPQ" -u "$sam" > "$bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$outbam" "$bam"
  rm -f "$sai1" "$sai2" "$sam" "$bam"
}

expected_pe=0
expected_se=0

if [[ -s "$POOL_PE_R1" && -s "$POOL_PE_R2" ]]; then
  shopt -s nullglob
  pe_r1_chunks=( "$MAPCHUNKS/${SAMPLE}.PE.R1_"*.fastq.gz )
  expected_pe="${#pe_r1_chunks[@]}"
fi
if [[ -s "$POOL_SE_ALL" ]]; then
  shopt -s nullglob
  se_chunks=( "$MAPCHUNKS/${SAMPLE}.SE_"*.fastq.gz )
  expected_se="${#se_chunks[@]}"
fi

do_mapping() {
  log "STEP: mapping (chunk-resume + delete-as-you-go)"
  shopt -s nullglob

  if [[ "$expected_pe" -eq 0 && "$expected_se" -eq 0 ]]; then
    die "No mapchunks found; nothing to map."
  fi

  # PE mapping
  if [[ "$expected_pe" -gt 0 ]]; then
    pe_r1_chunks=( "$MAPCHUNKS/${SAMPLE}.PE.R1_"*.fastq.gz )
    for r1c in "${pe_r1_chunks[@]}"; do
      base=$(basename "$r1c")
      cid="${base#${SAMPLE}.PE.R1_}"
      cid="${cid%.fastq.gz}"
      r2c="$MAPCHUNKS/${SAMPLE}.PE.R2_${cid}.fastq.gz"
      [[ -f "$r2c" ]] || die "Missing PE R2 mapchunk: $r2c"

      outbam="$BAMS/${SAMPLE}.MAP.PE.${cid}.bam"

      if [[ $RESUME -eq 1 && -s "$outbam" ]] && "$SAMTOOLS" quickcheck -q "$outbam"; then
        log "  Mapping PE chunk ${cid} (skipped; valid BAM exists)"
        # If mapchunk still exists, we can delete it to reclaim space.
        if [[ $KEEP_INTERMEDIATE -eq 0 ]]; then
          rm -f "$r1c" "$r2c" 2>/dev/null || true
        fi
        continue
      fi

      if [[ "$LIBTYPE" == "modern" ]]; then
        log "  Mapping PE chunk ${cid} (mem)"
        map_pe_mem_chunk "$r1c" "$r2c" "$outbam" "mem.PE.${cid}"
      elif [[ "$LIBTYPE" == "historical" ]]; then
        log "  Mapping PE chunk ${cid} (aln/sampe)"
        map_pe_aln_chunk "$r1c" "$r2c" "$outbam" "aln.PE.${cid}"
      else
        # ancient ignores unmerged PE
        continue
      fi

      file_nonempty "$outbam" || die "Mapping produced empty BAM: $outbam"
      "$SAMTOOLS" quickcheck -q "$outbam" || die "Mapped chunk BAM failed quickcheck: $outbam"

      # Delete mapchunk FASTQs as soon as BAM is valid (space)
      if [[ $KEEP_INTERMEDIATE -eq 0 ]]; then
        rm -f "$r1c" "$r2c" 2>/dev/null || true
      fi
    done
  fi

  # SE mapping
  if [[ "$expected_se" -gt 0 ]]; then
    se_chunks=( "$MAPCHUNKS/${SAMPLE}.SE_"*.fastq.gz )
    for sec in "${se_chunks[@]}"; do
      base=$(basename "$sec")
      cid="${base#${SAMPLE}.SE_}"
      cid="${cid%.fastq.gz}"
      outbam="$BAMS/${SAMPLE}.MAP.SE.${cid}.bam"

      if [[ $RESUME -eq 1 && -s "$outbam" ]] && "$SAMTOOLS" quickcheck -q "$outbam"; then
        log "  Mapping SE chunk ${cid} (skipped; valid BAM exists)"
        if [[ $KEEP_INTERMEDIATE -eq 0 ]]; then
          rm -f "$sec" 2>/dev/null || true
        fi
        continue
      fi

      if [[ "$LIBTYPE" == "modern" ]]; then
        log "  Mapping SE chunk ${cid} (mem)"
        map_se_mem_chunk "$sec" "$outbam" "mem.SE.${cid}"
      else
        log "  Mapping SE chunk ${cid} (aln/samse)"
        map_se_aln_chunk "$sec" "$outbam" "aln.SE.${cid}"
      fi

      file_nonempty "$outbam" || die "Mapping produced empty BAM: $outbam"
      "$SAMTOOLS" quickcheck -q "$outbam" || die "Mapped chunk BAM failed quickcheck: $outbam"

      if [[ $KEEP_INTERMEDIATE -eq 0 ]]; then
        rm -f "$sec" 2>/dev/null || true
      fi
    done
  fi

  # Sanity: check expected number of BAMs exist and are valid
  local got_pe=0 got_se=0
  local b
  for b in "$BAMS/${SAMPLE}.MAP.PE."*.bam; do
    [[ -e "$b" ]] || continue
    got_pe=$((got_pe+1))
    "$SAMTOOLS" quickcheck -q "$b" || die "Invalid PE chunk BAM after mapping: $b"
  done
  for b in "$BAMS/${SAMPLE}.MAP.SE."*.bam; do
    [[ -e "$b" ]] || continue
    got_se=$((got_se+1))
    "$SAMTOOLS" quickcheck -q "$b" || die "Invalid SE chunk BAM after mapping: $b"
  done

  # For ancient: expected_pe may be >0 but we intentionally skip PE; adjust expectation
  if [[ "$LIBTYPE" == "ancient" ]]; then
    expected_pe=0
  fi

  [[ "$got_pe" -eq "$expected_pe" ]] || die "Mapping incomplete (PE): expected $expected_pe chunk BAMs, got $got_pe"
  [[ "$got_se" -eq "$expected_se" ]] || die "Mapping incomplete (SE): expected $expected_se chunk BAMs, got $got_se"

  ckpt_mark map

  # Space: pooled FASTQs no longer needed once mapping is complete
  if [[ $KEEP_INTERMEDIATE -eq 0 ]]; then
    log "  Space: deleting pooled FASTQs (mapping complete)"
    rm -f "$POOL_PE_R1" "$POOL_PE_R2" "$POOL_SE_ALL" 2>/dev/null || true
  fi
}

if ! ckpt_ok map "$BAMS"; then
  do_mapping
else
  log "STEP: mapping (skipped; checkpoint present)"
fi

###############################################################################
# MERGE (OVERWRITE PARTIAL OUTPUTS WITH -f, ATOMIC) [ckpt: merge]
###############################################################################
MERGED_BAM="$FINAL/${SAMPLE}.merged.bam"

do_merge() {
  log "STEP: merge"
  shopt -s nullglob
  bam_list=( "$BAMS/${SAMPLE}.MAP."*.bam )
  [[ ${#bam_list[@]} -gt 0 ]] || die "No chunk BAMs found to merge."

  local tmp="${MERGED_BAM}.tmp.$$"
  "$SAMTOOLS" merge -f -@ "$SORT_THREADS" "$tmp" "${bam_list[@]}"
  "$SAMTOOLS" quickcheck -q "$tmp" || die "Merged BAM failed quickcheck: $tmp"
  mv -f "$tmp" "$MERGED_BAM"
  file_nonempty "$MERGED_BAM" || die "Merge failed: $MERGED_BAM"

  ckpt_mark merge

  # Space: delete per-chunk BAMs once merge succeeds (unless keeping)
  if [[ $KEEP_INTERMEDIATE -eq 0 ]]; then
    log "  Space: deleting mapped chunk BAMs (merge complete)"
    rm -f "${bam_list[@]}" 2>/dev/null || true
  fi
}

if ! ckpt_ok merge "$MERGED_BAM"; then
  do_merge
else
  log "STEP: merge (skipped; checkpoint present)"
fi

###############################################################################
# DEDUPLICATION (FLAG DUPLICATES, THEN FILTER) [ckpt: dedup]
###############################################################################
MARKDUP_BAM="$FINAL/${SAMPLE}.markdup_flagged.bam"
DEDUP_BAM="$FINAL/${SAMPLE}.dedup.bam"

do_dedup() {
  log "STEP: dedup (mark duplicates, then filter; duprate from duplicate flags)"

  local coordsort="$FINAL/${SAMPLE}.coordsort.bam"
  local namesort="$FINAL/${SAMPLE}.namesort.bam"
  local fixmate="$FINAL/${SAMPLE}.fixmate.bam"

  if [[ ( "$LIBTYPE" == "modern" || "$LIBTYPE" == "historical" ) && "$SEQ_MODE" != "SE" ]]; then
    "$SAMTOOLS" sort -n -@ "$SORT_THREADS" -o "$namesort" "$MERGED_BAM"
    "$SAMTOOLS" fixmate -m -@ "$SORT_THREADS" "$namesort" "$fixmate"
    "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$coordsort" "$fixmate"
    rm -f "$namesort" "$fixmate"
  else
    "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$coordsort" "$MERGED_BAM"
  fi

  local tmp_mark="${MARKDUP_BAM}.tmp.$$"
  "$SAMTOOLS" markdup -s -@ "$SORT_THREADS" "$coordsort" "$tmp_mark"
  rm -f "$coordsort"
  "$SAMTOOLS" quickcheck -q "$tmp_mark" || die "markdup output failed quickcheck: $tmp_mark"
  mv -f "$tmp_mark" "$MARKDUP_BAM"
  file_nonempty "$MARKDUP_BAM" || die "markdup failed: $MARKDUP_BAM"

  local nodup="$TMP/${SAMPLE}.nodup.unsorted.bam"
  "$SAMTOOLS" view -b -F 1024 "$MARKDUP_BAM" > "$nodup"

  local tmp_ded="${DEDUP_BAM}.tmp.$$"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$tmp_ded" "$nodup"
  rm -f "$nodup"
  "$SAMTOOLS" quickcheck -q "$tmp_ded" || die "Dedup BAM failed quickcheck: $tmp_ded"
  mv -f "$tmp_ded" "$DEDUP_BAM"
  file_nonempty "$DEDUP_BAM" || die "Dedup BAM missing/empty: $DEDUP_BAM"

  ckpt_mark dedup
}

if ! ckpt_ok dedup "$MARKDUP_BAM" "$DEDUP_BAM"; then
  do_dedup
else
  log "STEP: dedup (skipped; checkpoint present)"
fi

###############################################################################
# FINAL SORT+INDEX + RG (ATOMIC) [ckpt: final]
###############################################################################
FINAL_BAM="$FINAL/${SAMPLE}.final.bam"
OUT_BAM="$OUT/${SAMPLE}.final.RG.bam"

do_final() {
  log "STEP: final sort+index"
  local tmp_final="${FINAL_BAM}.tmp.$$"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$tmp_final" "$DEDUP_BAM"
  "$SAMTOOLS" index "$tmp_final"
  "$SAMTOOLS" quickcheck -q "$tmp_final" || die "Final BAM failed quickcheck: $tmp_final"
  mv -f "$tmp_final" "$FINAL_BAM"
  mv -f "${tmp_final}.bai" "${FINAL_BAM}.bai" 2>/dev/null || true
  "$SAMTOOLS" index "$FINAL_BAM" 2>/dev/null || true
  file_nonempty "$FINAL_BAM" || die "Final BAM missing/empty"

  log "STEP: read groups"
  local tmp_out="${OUT_BAM}.tmp.$$"
  RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina\tLB:${SAMPLE}"
  "$SAMTOOLS" addreplacerg -r "$RG" -o "$tmp_out" "$FINAL_BAM"
  "$SAMTOOLS" index "$tmp_out"
  "$SAMTOOLS" quickcheck -q "$tmp_out" || die "RG BAM failed quickcheck: $tmp_out"
  mv -f "$tmp_out" "$OUT_BAM"
  mv -f "${tmp_out}.bai" "${OUT_BAM}.bai" 2>/dev/null || true
  "$SAMTOOLS" index "$OUT_BAM" 2>/dev/null || true
  file_nonempty "$OUT_BAM" || die "RG BAM missing/empty"

  ckpt_mark final
}

if ! ckpt_ok final "$OUT_BAM" "${OUT_BAM}.bai"; then
  do_final
else
  log "STEP: final + RG (skipped; checkpoint present)"
fi

###############################################################################
# OPTIONAL mapDamage [ckpt: mapdamage]
###############################################################################
do_mapdamage() {
  log "STEP: mapDamage"
  local outdir="$OUT/${SAMPLE}_mapdamage"
  mkdir -p "$outdir"
  # mapDamage doesn't have a universal "done" file; we checkpoint it.
  "$MAPDAMAGE" -i "$OUT_BAM" --merge-reference-sequences --no-stats -r "$REF" -d "$outdir"
  ckpt_mark mapdamage
}

if [[ "$RUN_MAPDAMAGE" -eq 1 ]]; then
  if ! ckpt_ok mapdamage "$OUT_BAM"; then
    do_mapdamage
  else
    log "STEP: mapDamage (skipped; checkpoint present)"
  fi
fi

###############################################################################
# COVERAGE/DEPTH (ATOMIC) [ckpt: coverage]
###############################################################################
do_coverage() {
  log "STEP: coverage statistics"
  local tmp_cov="${COV_TSV}.tmp.$$"
  "$SAMTOOLS" coverage -q "$MAPQ" "$OUT_BAM" > "$tmp_cov"
  file_nonempty "$tmp_cov" || die "Coverage output empty: $tmp_cov"
  mv -f "$tmp_cov" "$COV_TSV"
  ckpt_mark coverage
}

if ! ckpt_ok coverage "$COV_TSV"; then
  do_coverage
else
  log "STEP: coverage (skipped; checkpoint present)"
fi

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
# SUMMARY STATS (ATOMIC) [ckpt: stats]
###############################################################################
do_stats() {
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

  # Exact mapped bp on final RG BAM (dedup + RG)
  MAPPED_BP=$(
    "$SAMTOOLS" view -F 2308 "$OUT_BAM" \
    | awk '{l=length($10); if(l>0){sum+=l}} END{printf "%.0f", sum+0}'
  )

  ENDOG=$("$PYTHON" - "$MAPPED_FRAGMENTS_UNIQUE" "$RAW_FRAGMENTS" <<'PY'
import sys
n=float(sys.argv[1]); d=float(sys.argv[2])
print("0" if d==0 else f"{(n/d):.6f}")
PY
  )

  local tmp_stats="${STATS}.tmp.$$"
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
  } > "$tmp_stats"
  mv -f "$tmp_stats" "$STATS"
  ckpt_mark stats
}

if ! ckpt_ok stats "$STATS"; then
  do_stats
else
  log "STEP: stats (skipped; checkpoint present)"
fi

log "STATS written: $STATS"
log "Coverage written: $COV_TSV"

###############################################################################
# CLEANUP
###############################################################################
if [[ $KEEP_INTERMEDIATE -eq 0 ]]; then
  log "Cleaning remaining intermediates (work dir only)"
  rm -rf "$WORK" 2>/dev/null || true
  [[ -n "$TMPDIR_USER" ]] && rm -rf "$TMPDIR_USER/plainmap_${SAMPLE}" 2>/dev/null || true
else
  log "Keeping intermediates: work dir retained at $WORK"
fi

PIPE_T1=$(date +%s)
log "STATUS: SUCCESS"
log "Total wall time: $((PIPE_T1 - PIPE_T0)) seconds"
log "Final BAM: $OUT_BAM"
