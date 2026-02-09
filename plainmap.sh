#!/usr/bin/env bash
# PlainMap v0.1
# Transparent, failure-aware mapping pipeline for ancient and modern DNA
#
# Key properties (current dev version):
# - Manifest can contain mixed SE + PE, mixed platforms, mixed naming conventions
# - SE/PE classification is derived from FASTQ headers:
#     * Illumina/CASAVA: header contains ' 1:' / ' 2:'
#     * ENA/older Illumina: token1 or token2 ends with /1 or /2
#     * SRA-export FASTQ (often lacks mate markers): R1/R2 headers are identical;
#       PlainMap pairs files by identical normalized first-read name (key) and
#       assigns mates deterministically by manifest order (first=R1, second=R2),
#       with strict sanity checks.
# - fastp runs per unit (or per unit-chunk), then pools:
#     * pooled unmerged PE (R1+R2) kept separate
#     * pooled SE-like reads (merged + unpaired + SE) kept separate
# - mapping is restartable:
#     * per-chunk resume (won't remap already-produced valid BAMs)
#     * chunk ID lists persisted, so resume works even if mapchunk FASTQs were deleted
# - pigz is used automatically if available (decompress/test/compress)
# - fragment-aware stats + duplication rate derived from duplicate flags (samtools markdup)
# - disk-space friendly: can delete chunk FASTQs and mapchunk FASTQs as soon as safe
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
PILOT_FRAGMENTS=0          # 0 disables pilot; otherwise PER-UNIT cap on fragments BEFORE fastp
MAPQ=20                    # mapping quality threshold (samtools view -q)

ADAPTER_R1=""
ADAPTER_R2=""
TRIM_ONLY=0

# mapDamage toggle (default OFF for all library types)
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

# Disk-space behavior:
#   1 = delete chunk/mapchunk FASTQs as soon as no longer needed
#   0 = keep chunk/mapchunk FASTQs until final cleanup (or forever with --keep-intermediate)
DELETE_AS_YOU_GO=1

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
  --pilot-fragments INT                      (default: 0; disabled) PER-UNIT cap on fragments before fastp
  --adapter-r1 SEQ
  --adapter-r2 SEQ
  --trim-only                                Trim only (fastp); exit before mapping
  --run-mapdamage                            Run mapDamage (any library type; default OFF)
  --no-mapdamage                             Do not run mapDamage (default)
  --tmpdir DIR                               Custom temp dir (e.g. node-local scratch)
  --keep-intermediate                        Keep <outdir>/<prefix>/work
  --resume | --no-resume
  --delete-as-you-go | --no-delete-as-you-go (default: delete-as-you-go)
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
    --run-mapdamage) RUN_MAPDAMAGE=1; shift ;;
    --no-mapdamage) RUN_MAPDAMAGE=0; shift ;;
    --tmpdir) TMPDIR_USER="$2"; shift 2 ;;
    --keep-intermediate) KEEP_INTERMEDIATE=1; shift ;;
    --resume) RESUME=1; shift ;;
    --no-resume) RESUME=0; shift ;;
    --delete-as-you-go) DELETE_AS_YOU_GO=1; shift ;;
    --no-delete-as-you-go) DELETE_AS_YOU_GO=0; shift ;;
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

RAW_META="$CKPT/${SAMPLE}.raw_counts.tsv"
PE_IDLIST="$CKPT/${SAMPLE}.mapchunk_ids.pe.txt"
SE_IDLIST="$CKPT/${SAMPLE}.mapchunk_ids.se.txt"

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
# BASIC FILE/VALIDITY HELPERS
###############################################################################
file_nonempty() { [[ -e "$1" && -s "$1" ]]; }
file_exists() { [[ -e "$1" ]]; }

bam_ok() {
  local b="$1"
  file_nonempty "$b" || return 1
  "$SAMTOOLS" quickcheck -v "$b" >/dev/null 2>&1
}

###############################################################################
# CHECKPOINTS
###############################################################################
ckpt_file() { echo "$CKPT/${SAMPLE}.$1.done"; }
ckpt_mark() { [[ $RESUME -eq 1 ]] && date -Is > "$(ckpt_file "$1")"; }
ckpt_clear() { rm -f "$(ckpt_file "$1")" 2>/dev/null || true; }

# Use with real output files that must be non-empty
ckpt_ok() {
  local step="$1"; shift
  [[ $RESUME -eq 1 ]] || return 1
  [[ -f "$(ckpt_file "$step")" ]] || return 1
  local f
  for f in "$@"; do
    file_nonempty "$f" || return 1
  done
  return 0
}

# Use for metadata/list files that may legitimately be empty (must exist)
ckpt_ok_exists() {
  local step="$1"; shift
  [[ $RESUME -eq 1 ]] || return 1
  [[ -f "$(ckpt_file "$step")" ]] || return 1
  local f
  for f in "$@"; do
    file_exists "$f" || return 1
  done
  return 0
}

if [[ $RESET -eq 1 ]]; then
  log "Reset requested: clearing all checkpoints for sample $SAMPLE"
  rm -f "$CKPT/${SAMPLE}."*.done 2>/dev/null || true
  rm -f "$RAW_META" "$PE_IDLIST" "$SE_IDLIST" 2>/dev/null || true
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
  log "  run_mapdamage:       $RUN_MAPDAMAGE (any library type)"
  log "  delete_as_you_go:    $DELETE_AS_YOU_GO"
  log "Plan: manifest -> build SE/PE units -> pre-fastp chunk (optional) -> pilot per-unit (optional) -> fastp -> pool -> mapchunk -> map (per-chunk resume) -> merge (-f overwrite) -> dedup -> RG -> (mapDamage optional) -> coverage -> stats"
  exit 0
fi

quick_preflight

###############################################################################
# BWA INDEX (PARALLEL SAFE + STALE LOCK RECOVERY)
###############################################################################
ensure_bwa_index() {
  local ref="$1"
  local lock="${ref}.bwa.lock"
  local idx=( "${ref}.bwt" "${ref}.pac" "${ref}.ann" "${ref}.amb" "${ref}.sa" )
  local missing=0
  local f

  missing=0
  for f in "${idx[@]}"; do [[ -f "$f" ]] || missing=1; done
  [[ "$missing" -eq 0 ]] && { log "Reference index: present"; return; }

  # If lock exists and is stale (>1h), remove it.
  if [[ -d "$lock" ]]; then
    local now mtime age
    now=$(date +%s)
    mtime=$(stat -c %Y "$lock" 2>/dev/null || echo 0)
    age=$(( now - mtime ))
    if [[ "$mtime" -gt 0 && "$age" -gt 3600 ]]; then
      log "WARNING: BWA lock appears stale (>1h) and index is missing. Removing stale lock: $lock"
      rm -rf "$lock" || true
    fi
  fi

  if mkdir "$lock" 2>/dev/null; then
    trap 'rm -rf "$lock"' EXIT
    log "Indexing reference: $ref"
    "$BWA" index "$ref"
    rm -rf "$lock"; trap - EXIT
  else
    log "Waiting for BWA index to finish (lock exists: $lock)"
    while :; do
      sleep 15
      missing=0
      for f in "${idx[@]}"; do [[ -f "$f" ]] || missing=1; done
      [[ "$missing" -eq 0 ]] && break

      if [[ -d "$lock" ]]; then
        local now2 mtime2 age2
        now2=$(date +%s)
        mtime2=$(stat -c %Y "$lock" 2>/dev/null || echo 0)
        age2=$(( now2 - mtime2 ))
        if [[ "$mtime2" -gt 0 && "$age2" -gt 3600 ]]; then
          log "WARNING: BWA lock became stale while waiting (>1h). Removing and retrying indexing."
          rm -rf "$lock" || true
          ensure_bwa_index "$ref"
          return
        fi
      fi
    done
    log "Reference index: present (built by another process)"
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
  local status=$?
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
first_header_line() {
  local fq="$1"
  "${ZCAT[@]}" "$fq" | head -n 1 || true
}

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

TOTAL_FILES=0
IDX=0

while read -r f; do
  [[ -z "$f" || "$f" =~ ^# ]] && continue
  [[ "$f" = /* ]] || f="$MANIFEST_DIR/$f"
  f="$(resolve_path "$f")"
  [[ -r "$f" ]] || die "FASTQ not readable: $f"

  TOTAL_FILES=$((TOTAL_FILES+1))
  IDX=$((IDX+1))

  hdr="$(first_header_line "$f")"
  [[ -n "$hdr" ]] || die "Could not read FASTQ header for: $f"

  dir="$(classify_fastq_direction "$f")"
  key="$(normalize_qname "$hdr")"
  [[ -n "$key" ]] || die "Could not derive header key for: $f"

  if [[ "$dir" == "R1" ]]; then
    [[ -z "${R1_BY_KEY[$key]:-}" ]] || die "Multiple R1 files share the same first-read key '$key'. Offending: $f and ${R1_BY_KEY[$key]}"
    R1_BY_KEY["$key"]="$f"
  elif [[ "$dir" == "R2" ]]; then
    [[ -z "${R2_BY_KEY[$key]:-}" ]] || die "Multiple R2 files share the same first-read key '$key'. Offending: $f and ${R2_BY_KEY[$key]}"
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
    [[ "$n1" -eq "$n2" ]] || die "SRA-style PE sanity check failed: read counts differ for key '$k' (R1?=$n1 R2?=$n2):
$f1
$f2"

    [[ -z "${R1_BY_KEY[$k]:-}" && -z "${R2_BY_KEY[$k]:-}" ]] || die "Internal error: UNKNOWN key '$k' already classified"
    R1_BY_KEY["$k"]="$f1"
    R2_BY_KEY["$k"]="$f2"
    PE_KEYS+=( "$k" )
  else
    die "Ambiguous SRA-style grouping for key '$k': found ${#entries[@]} files with no mate markers. PlainMap refuses to guess."
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

chunks_exist() {
  shopt -s nullglob
  local a=( "$CHUNKS/${SAMPLE}."*.fastq.gz )
  [[ ${#a[@]} -gt 0 ]]
}

###############################################################################
# STEP 1: CREATE PRE-FASTP CHUNKS PER UNIT (PE stays paired)
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

# BUGFIX: only skip chunking if checkpoint exists AND chunks still exist.
if [[ $RESUME -eq 1 && -f "$(ckpt_file chunk)" && chunks_exist ]]; then
  log "STEP: pre-fastp chunking (skipped; checkpoint present + chunks exist)"
else
  make_chunks
fi

###############################################################################
# STEP 2: FASTP PER PRE-FASTP CHUNK -> POOL TRIMMED OUTPUTS
###############################################################################
POOL_PE_R1="$RAW/${SAMPLE}.pool.PE.R1.fastq.gz"
POOL_PE_R2="$RAW/${SAMPLE}.pool.PE.R2.fastq.gz"
POOL_SE_ALL="$RAW/${SAMPLE}.pool.SE.all.fastq.gz"

run_fastp_and_pool() {
  log "STEP: fastp per unit-chunk (JSON only) + pooling (unmerged PE separate from SE-like)"
  ckpt_clear fastp_pool

  : > "$POOL_PE_R1"
  : > "$POOL_PE_R2"
  : > "$POOL_SE_ALL"

  local RAW_R1_READS=0
  local RAW_R2_READS=0
  local RAW_FRAGMENTS=0

  local TRIMMED_FRAGMENTS=0
  local MERGED_READS=0
  local UNPAIRED_READS=0

  declare -A PILOT_LEFT_BY_UNIT
  shopt -s nullglob

  # PE chunks
  local pe_r1_chunks=( "$CHUNKS/${SAMPLE}.PE"*".R1_"*.fastq.gz )
  for r1c in "${pe_r1_chunks[@]}"; do
    local base r2c unit_id unit_left in_r1 in_r2 pairs_in_chunk take_pairs
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
      [[ $DELETE_AS_YOU_GO -eq 1 ]] && rm -f "$r1c" "$r2c" 2>/dev/null || true
      continue
    fi

    in_r1="$r1c"
    in_r2="$r2c"

    pairs_in_chunk=$(count_reads_fastq_gz "$in_r1")
    take_pairs="$pairs_in_chunk"

    local lim_r1="" lim_r2=""
    if [[ "$PILOT_FRAGMENTS" -gt 0 ]]; then
      take_pairs=$(pilot_take "$pairs_in_chunk" "$unit_left")
      if [[ "$take_pairs" -le 0 ]]; then
        PILOT_LEFT_BY_UNIT["$unit_id"]=0
        [[ $DELETE_AS_YOU_GO -eq 1 ]] && rm -f "$r1c" "$r2c" 2>/dev/null || true
        continue
      fi
      if [[ "$take_pairs" -lt "$pairs_in_chunk" ]]; then
        log "  Pilot (per-unit): limiting PE unit $unit_id chunk $(basename "$r1c") to first $take_pairs pairs (unit remaining before: $unit_left)"
      fi
      lim_r1="$TMP/${SAMPLE}.$(basename "${r1c%.fastq.gz}").pilot.R1.fastq.gz"
      lim_r2="$TMP/${SAMPLE}.$(basename "${r2c%.fastq.gz}").pilot.R2.fastq.gz"
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

    local tag fastp_err fp_json trim1 trim2 merged up1 up2
    tag="$(basename "${r1c%.fastq.gz}")"
    tag="${tag//\//_}"

    fastp_err="$TMP/fastp_${SAMPLE}.${tag}.stderr.log"
    fp_json="$REPORTS/${SAMPLE}.${tag}.fastp.json"

    trim1="$TMP/${SAMPLE}.${tag}.trim.R1.fastq.gz"
    trim2="$TMP/${SAMPLE}.${tag}.trim.R2.fastq.gz"
    merged="$TMP/${SAMPLE}.${tag}.merged.fastq.gz"
    up1="$TMP/${SAMPLE}.${tag}.unpaired1.fastq.gz"
    up2="$TMP/${SAMPLE}.${tag}.unpaired2.fastq.gz"

    local extra_pe=()
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

    local out1_reads merged_reads up1n up2n
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
    [[ -n "$lim_r1" ]] && rm -f "$lim_r1" "$lim_r2" 2>/dev/null || true

    if [[ $DELETE_AS_YOU_GO -eq 1 ]]; then
      rm -f "$r1c" "$r2c" 2>/dev/null || true
    fi
  done

  # SE chunks
  local se_chunks=( "$CHUNKS/${SAMPLE}.SE"*".SE_"*.fastq.gz )
  for sec in "${se_chunks[@]}"; do
    local base unit_id unit_left in_se reads_in_chunk take_reads
    base=$(basename "$sec")
    unit_id="${base%.fastq.gz}"
    unit_id="${unit_id%%.SE_*}"

    if [[ "$PILOT_FRAGMENTS" -gt 0 && -z "${PILOT_LEFT_BY_UNIT[$unit_id]:-}" ]]; then
      PILOT_LEFT_BY_UNIT["$unit_id"]="$PILOT_FRAGMENTS"
    fi
    unit_left="${PILOT_LEFT_BY_UNIT[$unit_id]:-0}"

    if [[ "$PILOT_FRAGMENTS" -gt 0 && "$unit_left" -le 0 ]]; then
      [[ $DELETE_AS_YOU_GO -eq 1 ]] && rm -f "$sec" 2>/dev/null || true
      continue
    fi

    in_se="$sec"
    reads_in_chunk=$(count_reads_fastq_gz "$in_se")
    take_reads="$reads_in_chunk"

    local lim_se=""
    if [[ "$PILOT_FRAGMENTS" -gt 0 ]]; then
      take_reads=$(pilot_take "$reads_in_chunk" "$unit_left")
      if [[ "$take_reads" -le 0 ]]; then
        PILOT_LEFT_BY_UNIT["$unit_id"]=0
        [[ $DELETE_AS_YOU_GO -eq 1 ]] && rm -f "$sec" 2>/dev/null || true
        continue
      fi
      if [[ "$take_reads" -lt "$reads_in_chunk" ]]; then
        log "  Pilot (per-unit): limiting SE unit $unit_id chunk $(basename "$sec") to first $take_reads reads (unit remaining before: $unit_left)"
      fi
      lim_se="$TMP/${SAMPLE}.$(basename "${sec%.fastq.gz}").pilot.SE.fastq.gz"
      limit_fastq_gz "$in_se" "$lim_se" "$take_reads"
      in_se="$lim_se"
      unit_left="$(py_sub_nonneg "$unit_left" "$take_reads")"
      PILOT_LEFT_BY_UNIT["$unit_id"]="$unit_left"
    fi

    RAW_R1_READS=$(py_add "$RAW_R1_READS" "$take_reads")
    RAW_FRAGMENTS=$(py_add "$RAW_FRAGMENTS" "$take_reads")

    local tag fastp_err fp_json trim_se
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

    local se_after
    se_after=$(json_get_int_any "$fp_json" "read1_after_filtering.total_reads" "summary.after_filtering.total_reads")
    [[ "$se_after" -le 0 ]] && se_after=$(count_reads_fastq_gz "$trim_se")
    TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$se_after")

    rm -f "$trim_se" 2>/dev/null || true
    [[ -n "$lim_se" ]] && rm -f "$lim_se" 2>/dev/null || true

    if [[ $DELETE_AS_YOU_GO -eq 1 ]]; then
      rm -f "$sec" 2>/dev/null || true
    fi
  done

  {
    printf "raw_R1_reads\traw_R2_reads\traw_fragments\ttrimmed_fragments\tmerged_reads\tunpaired_reads\n"
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$RAW_R1_READS" "$RAW_R2_READS" "$RAW_FRAGMENTS" \
      "$TRIMMED_FRAGMENTS" "$MERGED_READS" "$UNPAIRED_READS"
  } > "$RAW_META"

  ckpt_mark fastp_pool
  log "fastp+pool complete. Raw counts written: $RAW_META"
}

need_fastp_files=( "$RAW_META" "$POOL_SE_ALL" )
if [[ "$SEQ_MODE" != "SE" ]]; then
  need_fastp_files+=( "$POOL_PE_R1" "$POOL_PE_R2" )
fi

if ! ckpt_ok fastp_pool "${need_fastp_files[@]}"; then
  run_fastp_and_pool
else
  log "STEP: fastp+pool (skipped; checkpoint present)"
fi

###############################################################################
# TRIM-ONLY MODE
###############################################################################
if [[ $TRIM_ONLY -eq 1 ]]; then
  log "Trim-only mode enabled; writing stats and exiting"

  RAW_R1_READS=$(awk 'NR==2{print $1}' "$RAW_META")
  RAW_R2_READS=$(awk 'NR==2{print $2}' "$RAW_META")
  RAW_FRAGMENTS=$(awk 'NR==2{print $3}' "$RAW_META")
  TRIMMED_FRAGMENTS=$(awk 'NR==2{print $4}' "$RAW_META")
  MERGED_READS=$(awk 'NR==2{print $5}' "$RAW_META")
  UNPAIRED_READS=$(awk 'NR==2{print $6}' "$RAW_META")

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
# STEP 3: POST-FASTP MAPPING CHUNKING (OPTIONAL; REUSE MAX_READS_PER_CHUNK)
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
  [[ "$n1" -eq "$n2" ]] || die "Pooled PE R1/R2 read counts differ (R1=$n1 R2=$n2). This would break sampe."

  if [[ "$MAX_READS_PER_CHUNK" -gt 0 ]]; then
    log "  Mapping chunking PE: splitting pooled PE (max reads/chunk: $MAX_READS_PER_CHUNK)"
    split_fastq_gz "$in1" "$prefix1" "$MAX_READS_PER_CHUNK"
    split_fastq_gz "$in2" "$prefix2" "$MAX_READS_PER_CHUNK"
  else
    ln -sf "$in1" "${prefix1}0000.fastq.gz"
    ln -sf "$in2" "${prefix2}0000.fastq.gz"
  fi
}

persist_mapchunk_ids() {
  # BUGFIX: These lists can be legitimately empty. We must create them regardless.
  : > "$PE_IDLIST"
  : > "$SE_IDLIST"

  shopt -s nullglob
  local f base cid

  if [[ "$LIBTYPE" != "ancient" ]]; then
    for f in "$MAPCHUNKS/${SAMPLE}.PE.R1_"*.fastq.gz; do
      base=$(basename "$f")
      cid="${base#${SAMPLE}.PE.R1_}"
      cid="${cid%.fastq.gz}"
      echo "$cid" >> "$PE_IDLIST"
    done
  fi

  for f in "$MAPCHUNKS/${SAMPLE}.SE_"*.fastq.gz; do
    base=$(basename "$f")
    cid="${base#${SAMPLE}.SE_}"
    cid="${cid%.fastq.gz}"
    echo "$cid" >> "$SE_IDLIST"
  done

  sort -u -o "$PE_IDLIST" "$PE_IDLIST" 2>/dev/null || true
  sort -u -o "$SE_IDLIST" "$SE_IDLIST" 2>/dev/null || true

  ckpt_mark mapchunk
}

log "STEP: post-fastp mapping chunking"

# BUGFIX: mapchunk checkpoint must accept empty idlists; also only require PE list when PE is expected.
need_mapchunk_lists=( "$SE_IDLIST" )
if [[ "$LIBTYPE" != "ancient" && "$SEQ_MODE" != "SE" ]]; then
  need_mapchunk_lists+=( "$PE_IDLIST" )
fi

if ! ckpt_ok_exists mapchunk "${need_mapchunk_lists[@]}"; then
  ckpt_clear mapchunk
  rm -f "$MAPCHUNKS/${SAMPLE}."*.fastq.gz 2>/dev/null || true

  if [[ "$LIBTYPE" != "ancient" && "$SEQ_MODE" != "SE" ]]; then
    make_mapchunks_pe "$POOL_PE_R1" "$POOL_PE_R2" "$MAPCHUNKS/${SAMPLE}.PE.R1_" "$MAPCHUNKS/${SAMPLE}.PE.R2_"
  fi
  make_mapchunks_se "$POOL_SE_ALL" "$MAPCHUNKS/${SAMPLE}.SE_"

  persist_mapchunk_ids
else
  log "  mapchunk (skipped; checkpoint present)"
fi

###############################################################################
# STEP 4: MAPPING (PER-CHUNK RESUME) + DELETE MAPCHUNK FASTQs AS SOON AS MAPPED
###############################################################################
log "STEP: mapping (per-chunk resume)"
mkdir -p "$BAMS"

map_se_mem_chunk() {
  local fq="$1" outbam="$2" tag="$3"
  local tmpbam="${outbam}.tmp"
  "$BWA" mem -t "$THREADS" "$REF" "$fq" > "$TMP/${SAMPLE}.${tag}.sam"
  "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.${tag}.sam" > "$TMP/${SAMPLE}.${tag}.bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$tmpbam" "$TMP/${SAMPLE}.${tag}.bam"
  rm -f "$TMP/${SAMPLE}.${tag}.sam" "$TMP/${SAMPLE}.${tag}.bam"
  bam_ok "$tmpbam" || die "Produced invalid BAM (SE mem): $tmpbam"
  mv -f "$tmpbam" "$outbam"
}

map_pe_mem_chunk() {
  local fq1="$1" fq2="$2" outbam="$3" tag="$4"
  local tmpbam="${outbam}.tmp"
  "$BWA" mem -t "$THREADS" "$REF" "$fq1" "$fq2" > "$TMP/${SAMPLE}.${tag}.sam"
  "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.${tag}.sam" > "$TMP/${SAMPLE}.${tag}.bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$tmpbam" "$TMP/${SAMPLE}.${tag}.bam"
  rm -f "$TMP/${SAMPLE}.${tag}.sam" "$TMP/${SAMPLE}.${tag}.bam"
  bam_ok "$tmpbam" || die "Produced invalid BAM (PE mem): $tmpbam"
  mv -f "$tmpbam" "$outbam"
}

map_se_aln_chunk() {
  local fq="$1" outbam="$2" tag="$3"
  local tmpbam="${outbam}.tmp"
  "$BWA" aln -l 999 -n "$MISMATCH" -t "$THREADS" "$REF" "$fq" > "$TMP/${SAMPLE}.${tag}.sai"
  "$BWA" samse "$REF" "$TMP/${SAMPLE}.${tag}.sai" "$fq" > "$TMP/${SAMPLE}.${tag}.sam"
  "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.${tag}.sam" > "$TMP/${SAMPLE}.${tag}.bam.tmp"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$tmpbam" "$TMP/${SAMPLE}.${tag}.bam.tmp"
  rm -f "$TMP/${SAMPLE}.${tag}.sai" "$TMP/${SAMPLE}.${tag}.sam" "$TMP/${SAMPLE}.${tag}.bam.tmp"
  bam_ok "$tmpbam" || die "Produced invalid BAM (SE aln): $tmpbam"
  mv -f "$tmpbam" "$outbam"
}

map_pe_aln_chunk() {
  local fq1="$1" fq2="$2" outbam="$3" tag="$4"
  local tmpbam="${outbam}.tmp"

  local h1 h2 n1x n2x
  h1="$(first_header_line "$fq1")"; h2="$(first_header_line "$fq2")"
  n1x="$(normalize_qname "$h1")"; n2x="$(normalize_qname "$h2")"
  [[ "$n1x" == "$n2x" ]] || die "Refusing to run sampe: PE chunks first read names differ: $(basename "$fq1")='$n1x' vs $(basename "$fq2")='$n2x'"

  "$BWA" aln -l 999 -n "$MISMATCH" -t "$THREADS" "$REF" "$fq1" > "$TMP/${SAMPLE}.${tag}.1.sai"
  "$BWA" aln -l 999 -n "$MISMATCH" -t "$THREADS" "$REF" "$fq2" > "$TMP/${SAMPLE}.${tag}.2.sai"
  "$BWA" sampe "$REF" "$TMP/${SAMPLE}.${tag}.1.sai" "$TMP/${SAMPLE}.${tag}.2.sai" "$fq1" "$fq2" > "$TMP/${SAMPLE}.${tag}.sam"
  "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.${tag}.sam" > "$TMP/${SAMPLE}.${tag}.bam.tmp"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$tmpbam" "$TMP/${SAMPLE}.${tag}.bam.tmp"
  rm -f "$TMP/${SAMPLE}.${tag}.1.sai" "$TMP/${SAMPLE}.${tag}.2.sai" "$TMP/${SAMPLE}.${tag}.sam" "$TMP/${SAMPLE}.${tag}.bam.tmp"
  bam_ok "$tmpbam" || die "Produced invalid BAM (PE aln): $tmpbam"
  mv -f "$tmpbam" "$outbam"
}

map_one_pe_id() {
  local cid="$1"
  local r1c="$MAPCHUNKS/${SAMPLE}.PE.R1_${cid}.fastq.gz"
  local r2c="$MAPCHUNKS/${SAMPLE}.PE.R2_${cid}.fastq.gz"
  local outbam="$BAMS/${SAMPLE}.MAP.PE.${cid}.bam"

  if bam_ok "$outbam"; then
    log "  PE chunk ${cid}: already mapped (valid BAM exists) -> skipping"
    [[ $DELETE_AS_YOU_GO -eq 1 ]] && rm -f "$r1c" "$r2c" 2>/dev/null || true
    return 0
  fi

  [[ -f "$r1c" && -f "$r2c" ]] || die "PE chunk ${cid} is missing FASTQs needed to map:
  $r1c
  $r2c
If you deleted mapchunk FASTQs, resume requires the BAM to exist and be valid."

  if [[ "$LIBTYPE" == "modern" ]]; then
    log "  Mapping PE chunk ${cid} (mem)"
    map_pe_mem_chunk "$r1c" "$r2c" "$outbam" "mem.PE.${cid}"
  elif [[ "$LIBTYPE" == "historical" ]]; then
    log "  Mapping PE chunk ${cid} (aln/sampe)"
    map_pe_aln_chunk "$r1c" "$r2c" "$outbam" "aln.PE.${cid}"
  else
    die "Internal error: PE mapping called for ancient"
  fi

  bam_ok "$outbam" || die "PE chunk ${cid}: produced BAM is invalid: $outbam"

  if [[ $DELETE_AS_YOU_GO -eq 1 ]]; then
    rm -f "$r1c" "$r2c" 2>/dev/null || true
  fi
}

map_one_se_id() {
  local cid="$1"
  local sec="$MAPCHUNKS/${SAMPLE}.SE_${cid}.fastq.gz"
  local outbam="$BAMS/${SAMPLE}.MAP.SE.${cid}.bam"

  if bam_ok "$outbam"; then
    log "  SE chunk ${cid}: already mapped (valid BAM exists) -> skipping"
    [[ $DELETE_AS_YOU_GO -eq 1 ]] && rm -f "$sec" 2>/dev/null || true
    return 0
  fi

  [[ -f "$sec" ]] || die "SE chunk ${cid} is missing FASTQ needed to map:
  $sec
If you deleted mapchunk FASTQs, resume requires the BAM to exist and be valid."

  if [[ "$LIBTYPE" == "modern" ]]; then
    log "  Mapping SE chunk ${cid} (mem)"
    map_se_mem_chunk "$sec" "$outbam" "mem.SE.${cid}"
  else
    log "  Mapping SE chunk ${cid} (aln/samse)"
    map_se_aln_chunk "$sec" "$outbam" "aln.SE.${cid}"
  fi

  bam_ok "$outbam" || die "SE chunk ${cid}: produced BAM is invalid: $outbam"

  if [[ $DELETE_AS_YOU_GO -eq 1 ]]; then
    rm -f "$sec" 2>/dev/null || true
  fi
}

# Map PE chunk IDs (skip if ancient)
if [[ "$LIBTYPE" != "ancient" && "$SEQ_MODE" != "SE" && -s "$PE_IDLIST" ]]; then
  while read -r cid; do
    [[ -z "$cid" ]] && continue
    map_one_pe_id "$cid"
  done < "$PE_IDLIST"
fi

# Map SE chunk IDs
if [[ -s "$SE_IDLIST" ]]; then
  while read -r cid; do
    [[ -z "$cid" ]] && continue
    map_one_se_id "$cid"
  done < "$SE_IDLIST"
fi

# Verify BAMs exist for every expected chunk ID
missing=0
if [[ "$LIBTYPE" != "ancient" && "$SEQ_MODE" != "SE" && -s "$PE_IDLIST" ]]; then
  while read -r cid; do
    [[ -z "$cid" ]] && continue
    bam_ok "$BAMS/${SAMPLE}.MAP.PE.${cid}.bam" || { log "ERROR: Missing/invalid PE BAM for chunk ${cid}"; missing=1; }
  done < "$PE_IDLIST"
fi
if [[ -s "$SE_IDLIST" ]]; then
  while read -r cid; do
    [[ -z "$cid" ]] && continue
    bam_ok "$BAMS/${SAMPLE}.MAP.SE.${cid}.bam" || { log "ERROR: Missing/invalid SE BAM for chunk ${cid}"; missing=1; }
  done < "$SE_IDLIST"
fi
[[ "$missing" -eq 0 ]] || die "Mapping incomplete: one or more expected chunk BAMs are missing/invalid"

###############################################################################
# STEP 5: MERGE (OVERWRITE SAFE)
###############################################################################
MERGED_BAM="$FINAL/${SAMPLE}.merged.bam"
log "STEP: merge"

bam_list=()

if [[ "$LIBTYPE" != "ancient" && "$SEQ_MODE" != "SE" && -s "$PE_IDLIST" ]]; then
  while read -r cid; do
    [[ -z "$cid" ]] && continue
    bam_list+=( "$BAMS/${SAMPLE}.MAP.PE.${cid}.bam" )
  done < "$PE_IDLIST"
fi
if [[ -s "$SE_IDLIST" ]]; then
  while read -r cid; do
    [[ -z "$cid" ]] && continue
    bam_list+=( "$BAMS/${SAMPLE}.MAP.SE.${cid}.bam" )
  done < "$SE_IDLIST"
fi

[[ ${#bam_list[@]} -gt 0 ]] || die "No BAMs produced by mapping"

if ckpt_ok merge "$MERGED_BAM" && bam_ok "$MERGED_BAM"; then
  log "  merge (skipped; checkpoint present)"
else
  ckpt_clear merge
  mkdir -p "$FINAL"
  "$SAMTOOLS" merge -f -@ "$SORT_THREADS" "$MERGED_BAM" "${bam_list[@]}"
  bam_ok "$MERGED_BAM" || die "Merge failed/invalid: $MERGED_BAM"
  ckpt_mark merge
fi

###############################################################################
# STEP 6: DEDUPLICATION (FLAG DUPLICATES, THEN FILTER)
###############################################################################
MARKDUP_BAM="$FINAL/${SAMPLE}.markdup_flagged.bam"
DEDUP_BAM="$FINAL/${SAMPLE}.dedup.bam"

log "STEP: dedup (mark duplicates, then filter; duprate from duplicate flags)"

if ckpt_ok dedup "$MARKDUP_BAM" "$DEDUP_BAM" && bam_ok "$MARKDUP_BAM" && bam_ok "$DEDUP_BAM"; then
  log "  dedup (skipped; checkpoint present)"
else
  ckpt_clear dedup

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
  bam_ok "$MARKDUP_BAM" || die "markdup failed/invalid: $MARKDUP_BAM"

  "$SAMTOOLS" view -b -F 1024 "$MARKDUP_BAM" > "$TMP/${SAMPLE}.nodup.unsorted.bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$DEDUP_BAM" "$TMP/${SAMPLE}.nodup.unsorted.bam"
  rm -f "$TMP/${SAMPLE}.nodup.unsorted.bam"
  bam_ok "$DEDUP_BAM" || die "Dedup BAM missing/invalid: $DEDUP_BAM"

  ckpt_mark dedup
fi

###############################################################################
# STEP 7: FINAL SORT + RG (RG OUTPUT GUARANTEED SORTED) + INDEX
###############################################################################
FINAL_BAM="$FINAL/${SAMPLE}.final.bam"
OUT_BAM="$OUT/${SAMPLE}.final.RG.bam"

log "STEP: final sort+index + read groups"

if ckpt_ok rg "$OUT_BAM" "${OUT_BAM}.bai" && bam_ok "$OUT_BAM"; then
  log "  final+RG (skipped; checkpoint present)"
else
  ckpt_clear rg

  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$FINAL_BAM" "$DEDUP_BAM"
  "$SAMTOOLS" index "$FINAL_BAM"
  bam_ok "$FINAL_BAM" || die "Final BAM invalid: $FINAL_BAM"

  RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina\tLB:${SAMPLE}"
  tmp_rg="$TMP/${SAMPLE}.rg.tmp.bam"
  "$SAMTOOLS" addreplacerg -r "$RG" -o "$tmp_rg" "$FINAL_BAM"

  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$OUT_BAM" "$tmp_rg"
  rm -f "$tmp_rg"

  "$SAMTOOLS" index "$OUT_BAM"
  bam_ok "$OUT_BAM" || die "RG BAM invalid: $OUT_BAM"

  ckpt_mark rg
fi

###############################################################################
# OPTIONAL: mapDamage
###############################################################################
if [[ "$RUN_MAPDAMAGE" -eq 1 ]]; then
  if [[ $RESUME -eq 1 && -f "$(ckpt_file mapdamage)" ]]; then
    log "STEP: mapDamage (skipped; checkpoint present)"
  else
    log "STEP: mapDamage"
    ckpt_clear mapdamage
    "$MAPDAMAGE" -i "$OUT_BAM" --merge-reference-sequences --no-stats -r "$REF" -d "$OUT/${SAMPLE}_mapdamage"
    ckpt_mark mapdamage
  fi
fi

###############################################################################
# STEP 8: COVERAGE/DEPTH
###############################################################################
log "STEP: coverage statistics"
if ckpt_ok coverage "$COV_TSV"; then
  log "  coverage (skipped; checkpoint present)"
else
  ckpt_clear coverage
  "$SAMTOOLS" coverage -q "$MAPQ" "$OUT_BAM" > "$COV_TSV"
  file_nonempty "$COV_TSV" || die "Coverage TSV missing/empty: $COV_TSV"
  ckpt_mark coverage
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
# STEP 9: SUMMARY STATS
###############################################################################
log "STEP: summary statistics"

RAW_R1_READS=$(awk 'NR==2{print $1}' "$RAW_META")
RAW_R2_READS=$(awk 'NR==2{print $2}' "$RAW_META")
RAW_FRAGMENTS=$(awk 'NR==2{print $3}' "$RAW_META")
TRIMMED_FRAGMENTS=$(awk 'NR==2{print $4}' "$RAW_META")
MERGED_READS=$(awk 'NR==2{print $5}' "$RAW_META")
UNPAIRED_READS=$(awk 'NR==2{print $6}' "$RAW_META")

# Exclude: unmapped(4) + secondary(256) + QC-fail(512) + supplementary(2048) = 2820
EXCL_FLAGS=2820

MAPPED_READS_ALL=$("$SAMTOOLS" view -c -F "$EXCL_FLAGS" "$MERGED_BAM")
MAPPED_READS_UNIQUE=$("$SAMTOOLS" view -c -F "$EXCL_FLAGS" "$OUT_BAM")

MAPPED_FRAGMENTS_ALL=$("$SAMTOOLS" view -F "$EXCL_FLAGS" "$MARKDUP_BAM" \
  | awk 'BEGIN{c=0}
         function hasbit(x,b){return and(x,b)}
         {flag=$2;
          paired=hasbit(flag,1);
          read1=hasbit(flag,64);
          if(paired){ if(read1) c++ } else { c++ }
         }
         END{print c}'
)

DUP_FRAGMENTS=$("$SAMTOOLS" view -F "$EXCL_FLAGS" "$MARKDUP_BAM" \
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
  "$SAMTOOLS" view -F "$EXCL_FLAGS" "$OUT_BAM" \
  | awk '{l=length($10); if(l>0){sum+=l;n++}} END{ if(n>0) printf "%.2f", sum/n; else printf "0.00"}'
)

# Precise mapped_bp computed on the final BAM using samtools view.
MAPPED_BP=$(
  "$SAMTOOLS" view -F "$EXCL_FLAGS" "$OUT_BAM" \
  | awk 'BEGIN{sum=0} {sum+=length($10)} END{printf "%.0f", sum}'
)

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
