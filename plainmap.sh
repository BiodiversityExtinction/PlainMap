#!/usr/bin/env bash
# PlainMap v0.1
# Transparent, failure-aware mapping pipeline for ancient and modern DNA
# Fix in this version:
#   - Resolve OUTDIR/REF/MANIFEST and manifest FASTQ paths to ABSOLUTE paths
#     so symlinks created during "chunking disabled" are not broken.
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
MISMATCH=0.01
LIBTYPE="modern"            # modern|ancient
MAX_READS_PER_CHUNK=0       # 0 disables chunking
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
PlainMap v$VERSION

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
  --adapter-r1 SEQ
  --adapter-r2 SEQ
  --trim-only                     Trim only (fastp); exit before mapping
  --tmpdir DIR                    Custom temp dir (e.g. node-local scratch)
  --keep-intermediate             Keep <outdir>/<prefix>/work
  --resume | --no-resume
  --dry-run                       Print plan only (no gzip -t, no work)
  --validate                      Check tools + gzip -t all manifest FASTQs, then exit
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
# PATH RESOLUTION (CRITICAL FIX)
###############################################################################
resolve_path() {
  # canonical absolute path (does not require realpath)
  "$PYTHON" - "$1" <<'PY'
import os,sys
p=sys.argv[1]
p=os.path.expanduser(p)
print(os.path.realpath(os.path.abspath(p)))
PY
}

# OUT must exist before we can safely canonicalize to a real directory
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
# SAMPLE-SPECIFIC WORK DIRS (PARALLEL-SAFE)
###############################################################################
WORK="$OUT/${SAMPLE}/work"
RAW="$WORK/raw"
CHUNKS="$WORK/chunks"
BAMS="$WORK/chunk_bams"
FINAL="$WORK/final"
CKPT="$WORK/.checkpoints"
REPORTS="$OUT/${SAMPLE}/fastp_reports"

STATS="$OUT/${SAMPLE}.plainmap.stats.tsv"
COV_TSV="$OUT/${SAMPLE}.plainmap.coverage.tsv"

if [[ -n "$TMPDIR_USER" ]]; then
  TMP="$TMPDIR_USER/plainmap_${SAMPLE}"
else
  TMP="$WORK/tmp"
fi

mkdir -p "$RAW" "$CHUNKS" "$BAMS" "$FINAL" "$CKPT" "$TMP" "$REPORTS"

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
# CHECKPOINTS (VALID ONLY IF OUTPUTS EXIST)
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
  require_cmd zcat; require_cmd gzip; require_cmd split
  require_cmd "$FASTP"; require_cmd "$BWA"; require_cmd "$SAMTOOLS"; require_cmd "$PYTHON"
  check_split_filter_support
  [[ "$LIBTYPE" != "ancient" ]] || require_cmd "$MAPDAMAGE"
}

if [[ $VALIDATE_ONLY -eq 1 ]]; then
  log "Validation mode: checking tools + gzip -t for all manifest FASTQs"
  quick_preflight
  n=0
  while read -r f; do
    [[ -z "$f" || "$f" =~ ^# ]] && continue
    n=$((n+1))
    # resolve relative paths relative to manifest directory
    [[ "$f" = /* ]] || f="$MANIFEST_DIR/$f"
    f="$(resolve_path "$f")"
    [[ -r "$f" ]] || die "FASTQ not readable: $f"
    gzip -t "$f" >/dev/null || die "Corrupt/truncated gzip: $f"
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
  log "Planned: concat/link -> chunk -> fastp -> map -> merge -> dedup -> RG -> stats -> cleanup"
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
# FASTP WRAPPER + ARGS
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

FASTP_ARGS="-l $MINLEN -g -w $THREADS"
[[ -n "$ADAPTER_R1" ]] && FASTP_ARGS+=" --adapter_sequence $ADAPTER_R1"
[[ -n "$ADAPTER_R2" ]] && FASTP_ARGS+=" --adapter_sequence_r2 $ADAPTER_R2"

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
pct() { "$PYTHON" - "$1" "$2" <<'PY'
import sys
n=float(sys.argv[1]); d=float(sys.argv[2])
print("0" if d==0 else f"{(100.0*n/d):.3f}")
PY
}
ratio() { "$PYTHON" - "$1" "$2" <<'PY'
import sys
n=float(sys.argv[1]); d=float(sys.argv[2])
print("0" if d==0 else f"{(n/d):.6f}")
PY
}

###############################################################################
# READ MANIFEST + DETECT SE/PE
###############################################################################
R1_LIST=(); R2_LIST=()
while read -r f; do
  [[ -z "$f" || "$f" =~ ^# ]] && continue
  [[ "$f" = /* ]] || f="$MANIFEST_DIR/$f"
  f="$(resolve_path "$f")"
  [[ -r "$f" ]] || die "FASTQ not readable: $f"
  hdr=$(zcat "$f" | head -n1 || true)
  [[ "$hdr" =~ " 1:" ]] && R1_LIST+=("$f")
  [[ "$hdr" =~ " 2:" ]] && R2_LIST+=("$f")
done < "$MANIFEST"

[[ ${#R1_LIST[@]} -gt 0 ]] || die "No R1 reads found in manifest"

HAS_R2=0
[[ ${#R2_LIST[@]} -gt 0 ]] && HAS_R2=1
log "Detected sequencing mode: $([[ $HAS_R2 -eq 1 ]] && echo paired-end || echo single-end)"

###############################################################################
# CONCATENATE / LINK
###############################################################################
out_r1="$RAW/${SAMPLE}_R1.fastq.gz"
out_r2="$RAW/${SAMPLE}_R2.fastq.gz"

concat_gz_list() {
  local out="$1"; shift
  local arr=( "$@" )
  : > "$out"
  for f in "${arr[@]}"; do
    cat "$f" >> "$out"
  done
  file_nonempty "$out" || die "Concatenation produced empty output: $out"
}

if [[ $HAS_R2 -eq 1 ]]; then
  if ! ckpt_ok concat "$out_r1" "$out_r2"; then
    log "STEP: concatenate/link"
    ckpt_clear concat

    if [[ ${#R1_LIST[@]} -eq 1 ]]; then
      log "Linking single R1 -> $out_r1"
      ln -sf "${R1_LIST[0]}" "$out_r1"
    else
      log "Concatenating ${#R1_LIST[@]} R1 FASTQs -> $out_r1"
      concat_gz_list "$out_r1" "${R1_LIST[@]}"
    fi
    file_nonempty "$out_r1" || die "R1 output missing/empty: $out_r1"

    if [[ ${#R2_LIST[@]} -eq 1 ]]; then
      log "Linking single R2 -> $out_r2"
      ln -sf "${R2_LIST[0]}" "$out_r2"
    else
      log "Concatenating ${#R2_LIST[@]} R2 FASTQs -> $out_r2"
      concat_gz_list "$out_r2" "${R2_LIST[@]}"
    fi
    file_nonempty "$out_r2" || die "R2 output missing/empty: $out_r2"

    ckpt_mark concat
  else
    log "STEP: concatenate/link (skipped; outputs present)"
  fi
else
  if ! ckpt_ok concat "$out_r1"; then
    log "STEP: concatenate/link"
    ckpt_clear concat
    if [[ ${#R1_LIST[@]} -eq 1 ]]; then
      log "Linking single R1 -> $out_r1"
      ln -sf "${R1_LIST[0]}" "$out_r1"
    else
      log "Concatenating ${#R1_LIST[@]} R1 FASTQs -> $out_r1"
      concat_gz_list "$out_r1" "${R1_LIST[@]}"
    fi
    file_nonempty "$out_r1" || die "R1 output missing/empty: $out_r1"
    ckpt_mark concat
  else
    log "STEP: concatenate/link (skipped; outputs present)"
  fi
fi

###############################################################################
# COUNT RAW READS
###############################################################################
log "STEP: counting raw reads"
RAW_R1_READS=$(count_reads_fastq_gz "$out_r1")
RAW_R2_READS=0
[[ $HAS_R2 -eq 1 ]] && RAW_R2_READS=$(count_reads_fastq_gz "$out_r2")
RAW_PAIRS="$RAW_R1_READS"

###############################################################################
# CHUNKING
###############################################################################
split_fastq() {
  local fq="$1" prefix="$2" reads="$3"
  local lines=$((reads * 4))
  zcat "$fq" | split -l "$lines" -d -a 4 --filter='gzip > $FILE.fastq.gz' - "$prefix"
}

log "STEP: chunking"
first_r1_chunk="$CHUNKS/${SAMPLE}.R1_0000.fastq.gz"

if ! ckpt_ok chunk "$first_r1_chunk"; then
  ckpt_clear chunk
  rm -f "$CHUNKS/${SAMPLE}.R1_"*.fastq.gz 2>/dev/null || true
  rm -f "$CHUNKS/${SAMPLE}.R2_"*.fastq.gz 2>/dev/null || true

  if [[ "$MAX_READS_PER_CHUNK" -gt 0 ]]; then
    log "Splitting FASTQs (max reads per chunk: $MAX_READS_PER_CHUNK)"
    split_fastq "$out_r1" "$CHUNKS/${SAMPLE}.R1_" "$MAX_READS_PER_CHUNK"
    [[ $HAS_R2 -eq 1 ]] && split_fastq "$out_r2" "$CHUNKS/${SAMPLE}.R2_" "$MAX_READS_PER_CHUNK"
  else
    log "Chunking disabled; linking as single chunk"
    # IMPORTANT: out_r1/out_r2 are ABSOLUTE => symlink will not break
    ln -sf "$out_r1" "$first_r1_chunk"
    [[ $HAS_R2 -eq 1 ]] && ln -sf "$out_r2" "$CHUNKS/${SAMPLE}.R2_0000.fastq.gz"
  fi

  file_nonempty "$first_r1_chunk" || die "Chunking failed: chunk file missing/empty (symlink broken?)"
  ckpt_mark chunk
else
  log "Chunking (skipped; outputs present)"
fi

###############################################################################
# PER-CHUNK FASTP + MAPPING
###############################################################################
log "STEP: per-chunk fastp (and mapping unless --trim-only)"
shopt -s nullglob
R1_CHUNKS=( "$CHUNKS"/"${SAMPLE}.R1_"*.fastq.gz )
[[ ${#R1_CHUNKS[@]} -gt 0 ]] || die "No R1 chunks found in $CHUNKS"

TOTAL_R1_AFTER=0
TOTAL_R2_AFTER=0
TOTAL_MERGED=0

for r1 in "${R1_CHUNKS[@]}"; do
  base=$(basename "$r1")
  id="${base#${SAMPLE}.R1_}"
  id="${id%.fastq.gz}"
  r2="$CHUNKS/${SAMPLE}.R2_${id}.fastq.gz"

  fastp_err="$TMP/fastp_${SAMPLE}_${id}.stderr.log"
  fp_json="$REPORTS/${SAMPLE}.${id}.fastp.json"
  fp_html="$REPORTS/${SAMPLE}.${id}.fastp.html"
  outbam="$BAMS/${SAMPLE}.${id}.bam"

  log "Processing chunk $id"

  if [[ $HAS_R2 -eq 0 ]]; then
    trim_se="$TMP/${SAMPLE}.${id}.SE.trim.fastq"
    run_fastp "$fastp_err" $FASTP_ARGS -i "$r1" --out1 "$trim_se" -j "$fp_json" -h "$fp_html"
    file_nonempty "$trim_se" || die "fastp produced empty trimmed SE FASTQ"

    r1_after=$(json_get_int "$fp_json" "read1_after_filtering.total_reads")
    TOTAL_R1_AFTER=$("$PYTHON" - <<PY
print(int($TOTAL_R1_AFTER) + int($r1_after))
PY
)

    if [[ $TRIM_ONLY -eq 0 ]]; then
      if [[ "$LIBTYPE" == "ancient" ]]; then
        "$BWA" aln -l 999 -n "$MISMATCH" -t "$THREADS" "$REF" "$trim_se" > "$TMP/${SAMPLE}.${id}.sai"
        "$BWA" samse "$REF" "$TMP/${SAMPLE}.${id}.sai" "$trim_se" > "$TMP/${SAMPLE}.${id}.sam"
        "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.${id}.sam" > "$TMP/${SAMPLE}.${id}.bam.tmp"
        "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$outbam" "$TMP/${SAMPLE}.${id}.bam.tmp"
        rm -f "$TMP/${SAMPLE}.${id}.sai" "$TMP/${SAMPLE}.${id}.sam" "$TMP/${SAMPLE}.${id}.bam.tmp"
      else
        "$BWA" mem -t "$THREADS" "$REF" "$trim_se" > "$TMP/${SAMPLE}.${id}.SE.sam"
        "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.${id}.SE.sam" > "$TMP/${SAMPLE}.${id}.SE.bam"
        "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$outbam" "$TMP/${SAMPLE}.${id}.SE.bam"
        rm -f "$TMP/${SAMPLE}.${id}.SE.sam" "$TMP/${SAMPLE}.${id}.SE.bam"
      fi
    fi
    rm -f "$trim_se"

  else
    [[ -f "$r2" ]] || die "Missing chunk R2: $r2"

    if [[ "$LIBTYPE" == "ancient" ]]; then
      merged="$TMP/${SAMPLE}.${id}.merged.fastq"
      run_fastp "$fastp_err" $FASTP_ARGS -i "$r1" -I "$r2" -m --merged_out "$merged" -j "$fp_json" -h "$fp_html"
      file_nonempty "$merged" || die "fastp produced empty merged FASTQ"

      r1_after=$(json_get_int "$fp_json" "read1_after_filtering.total_reads")
      merged_reads=$(json_get_int "$fp_json" "merged_and_filtered.total_reads")
      TOTAL_R1_AFTER=$("$PYTHON" - <<PY
print(int($TOTAL_R1_AFTER) + int($r1_after))
PY
)
      TOTAL_MERGED=$("$PYTHON" - <<PY
print(int($TOTAL_MERGED) + int($merged_reads))
PY
)

      if [[ $TRIM_ONLY -eq 0 ]]; then
        "$BWA" aln -l 999 -n "$MISMATCH" -t "$THREADS" "$REF" "$merged" > "$TMP/${SAMPLE}.${id}.sai"
        "$BWA" samse "$REF" "$TMP/${SAMPLE}.${id}.sai" "$merged" > "$TMP/${SAMPLE}.${id}.sam"
        "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.${id}.sam" > "$TMP/${SAMPLE}.${id}.bam.tmp"
        "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$outbam" "$TMP/${SAMPLE}.${id}.bam.tmp"
        rm -f "$TMP/${SAMPLE}.${id}.sai" "$TMP/${SAMPLE}.${id}.sam" "$TMP/${SAMPLE}.${id}.bam.tmp"
      fi
      rm -f "$merged"

    else
      trim1="$TMP/${SAMPLE}.${id}.R1.trim.fastq"
      trim2="$TMP/${SAMPLE}.${id}.R2.trim.fastq"
      merged="$TMP/${SAMPLE}.${id}.merged.fastq"
      up1="$TMP/${SAMPLE}.${id}.unpaired1.fastq"
      up2="$TMP/${SAMPLE}.${id}.unpaired2.fastq"
      secat="$TMP/${SAMPLE}.${id}.SE.concat.fastq"

      run_fastp "$fastp_err" $FASTP_ARGS \
        -i "$r1" -I "$r2" \
        --out1 "$trim1" --out2 "$trim2" \
        -m --merged_out "$merged" \
        --unpaired1 "$up1" --unpaired2 "$up2" \
        -j "$fp_json" -h "$fp_html"

      file_nonempty "$trim1" || die "fastp produced empty trimmed R1"
      file_nonempty "$trim2" || die "fastp produced empty trimmed R2"

      r1_after=$(json_get_int "$fp_json" "read1_after_filtering.total_reads")
      r2_after=$(json_get_int "$fp_json" "read2_after_filtering.total_reads")
      merged_reads=$(json_get_int "$fp_json" "merged_and_filtered.total_reads")
      TOTAL_R1_AFTER=$("$PYTHON" - <<PY
print(int($TOTAL_R1_AFTER) + int($r1_after))
PY
)
      TOTAL_R2_AFTER=$("$PYTHON" - <<PY
print(int($TOTAL_R2_AFTER) + int($r2_after))
PY
)
      TOTAL_MERGED=$("$PYTHON" - <<PY
print(int($TOTAL_MERGED) + int($merged_reads))
PY
)

      if [[ $TRIM_ONLY -eq 0 ]]; then
        "$BWA" mem -t "$THREADS" "$REF" "$trim1" "$trim2" > "$TMP/${SAMPLE}.${id}.PE.sam"
        "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.${id}.PE.sam" > "$TMP/${SAMPLE}.${id}.PE.bam"
        "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$TMP/${SAMPLE}.${id}.PE.sort.bam" "$TMP/${SAMPLE}.${id}.PE.bam"
        rm -f "$TMP/${SAMPLE}.${id}.PE.sam" "$TMP/${SAMPLE}.${id}.PE.bam"

        : > "$secat"
        [[ -s "$merged" ]] && cat "$merged" >> "$secat"
        [[ -s "$up1" ]] && cat "$up1" >> "$secat"
        [[ -s "$up2" ]] && cat "$up2" >> "$secat"

        if [[ -s "$secat" ]]; then
          "$BWA" mem -t "$THREADS" "$REF" "$secat" > "$TMP/${SAMPLE}.${id}.SE.sam"
          "$SAMTOOLS" view -q "$MAPQ" -u "$TMP/${SAMPLE}.${id}.SE.sam" > "$TMP/${SAMPLE}.${id}.SE.bam"
          "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$TMP/${SAMPLE}.${id}.SE.sort.bam" "$TMP/${SAMPLE}.${id}.SE.bam"
          rm -f "$TMP/${SAMPLE}.${id}.SE.sam" "$TMP/${SAMPLE}.${id}.SE.bam"
          "$SAMTOOLS" merge -@ "$SORT_THREADS" "$outbam" "$TMP/${SAMPLE}.${id}.PE.sort.bam" "$TMP/${SAMPLE}.${id}.SE.sort.bam"
          rm -f "$TMP/${SAMPLE}.${id}.PE.sort.bam" "$TMP/${SAMPLE}.${id}.SE.sort.bam"
        else
          cp "$TMP/${SAMPLE}.${id}.PE.sort.bam" "$outbam"
          rm -f "$TMP/${SAMPLE}.${id}.PE.sort.bam"
        fi
      fi

      rm -f "$trim1" "$trim2" "$merged" "$up1" "$up2" "$secat"
    fi
  fi
done

###############################################################################
# TRIMMING STATS (CHUNK-SAFE)
###############################################################################
if [[ $HAS_R2 -eq 1 ]]; then
  TOTAL_TRIMMED_PAIRS=$("$PYTHON" - <<PY
print(min(int($TOTAL_R1_AFTER), int($TOTAL_R2_AFTER)))
PY
)
  TOTAL_UNPAIRED=$("$PYTHON" - <<PY
tp=min(int($TOTAL_R1_AFTER), int($TOTAL_R2_AFTER))
print(max(0, (int($TOTAL_R1_AFTER)+int($TOTAL_R2_AFTER)) - 2*tp))
PY
)
else
  TOTAL_TRIMMED_PAIRS="$TOTAL_R1_AFTER"
  TOTAL_UNPAIRED=0
fi

if [[ $TRIM_ONLY -eq 1 ]]; then
  log "Trim-only mode enabled; exiting"
  PCT_REMAIN=$(pct "$TOTAL_TRIMMED_PAIRS" "$RAW_PAIRS")
  RATIO_TRIM=$(ratio "$TOTAL_TRIMMED_PAIRS" "$RAW_PAIRS")
  {
    printf "sample\tlibrary_type\tseq_mode\traw_R1_reads\traw_R2_reads\traw_pairs\ttrimmed_pairs\tmerged_reads\tunpaired_reads\tpct_remaining_after_trim\ttrimmed_over_raw\n"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$SAMPLE" "$LIBTYPE" "$([[ $HAS_R2 -eq 1 ]] && echo PE || echo SE)" \
      "$RAW_R1_READS" "$RAW_R2_READS" "$RAW_PAIRS" \
      "$TOTAL_TRIMMED_PAIRS" "$TOTAL_MERGED" "$TOTAL_UNPAIRED" \
      "$PCT_REMAIN" "$RATIO_TRIM"
  } > "$STATS"
  PIPE_T1=$(date +%s)
  log "STATS written: $STATS"
  log "Total wall time: $((PIPE_T1 - PIPE_T0)) seconds"
  exit 0
fi

###############################################################################
# MERGE + DEDUP + RG
###############################################################################
MERGED_BAM="$FINAL/${SAMPLE}.merged.bam"
DEDUP_BAM="$FINAL/${SAMPLE}.dedup.bam"
FINAL_BAM="$FINAL/${SAMPLE}.final.bam"
OUT_BAM="$OUT/${SAMPLE}.final.RG.bam"

log "STEP: merge"
"$SAMTOOLS" merge -@ "$SORT_THREADS" "$MERGED_BAM" "$BAMS"/*.bam
file_nonempty "$MERGED_BAM" || die "Merge failed: $MERGED_BAM"

log "STEP: dedup"
if [[ $HAS_R2 -eq 0 || "$LIBTYPE" == "ancient" ]]; then
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$FINAL/${SAMPLE}.coordsort.bam" "$MERGED_BAM"
  "$SAMTOOLS" markdup -r -s -@ "$SORT_THREADS" "$FINAL/${SAMPLE}.coordsort.bam" "$DEDUP_BAM"
  rm -f "$FINAL/${SAMPLE}.coordsort.bam"
else
  "$SAMTOOLS" sort -n -@ "$SORT_THREADS" -o "$FINAL/${SAMPLE}.namesort.bam" "$MERGED_BAM"
  "$SAMTOOLS" fixmate -m -@ "$SORT_THREADS" "$FINAL/${SAMPLE}.namesort.bam" "$FINAL/${SAMPLE}.fixmate.bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$FINAL/${SAMPLE}.coordsort.bam" "$FINAL/${SAMPLE}.fixmate.bam"
  "$SAMTOOLS" markdup -r -@ "$SORT_THREADS" "$FINAL/${SAMPLE}.coordsort.bam" "$DEDUP_BAM"
  rm -f "$FINAL/${SAMPLE}.namesort.bam" "$FINAL/${SAMPLE}.fixmate.bam" "$FINAL/${SAMPLE}.coordsort.bam"
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

log "STEP: coverage statistics"
"$SAMTOOLS" coverage -q "$MAPQ" "$OUT_BAM" > "$COV_TSV"

# genome-wide mean coverage (length-weighted), fast
AVG_COV=$(
  awk 'BEGIN{sum=0;len=0}
       NR==1{next}                    # skip header
       {sum += $7*$6; len += $6}      # mean_depth * ref_len
       END{if(len>0) printf "%.6f", sum/len; else print 0}' \
    "$COV_TSV"
)

log "STEP: summary statistics"
MAPPED_ALL=$("$SAMTOOLS" view -c -F 4 "$MERGED_BAM")
MAPPED_UNIQUE=$("$SAMTOOLS" view -c -F 4 "$DEDUP_BAM")

AVG_READLEN=$(
  "$SAMTOOLS" view -F 4 "$DEDUP_BAM" \
  | awk '{l=length($10); if(l>0){sum+=l;n++}} END{ if(n>0) printf "%.2f", sum/n; else printf "0.00"}'
)

MAPPED_BP=$(awk -v n="$MAPPED_ALL" -v l="$AVG_READLEN" 'BEGIN{printf "%.0f", n*l}')

PCT_REMAIN=$(pct "$TOTAL_TRIMMED_PAIRS" "$RAW_PAIRS")
RATIO_TRIM=$(ratio "$TOTAL_TRIMMED_PAIRS" "$RAW_PAIRS")
ENDOG=$(ratio "$MAPPED_UNIQUE" "$RAW_PAIRS")
DUPRATE=$(
  echo "1 - $(ratio "$MAPPED_UNIQUE" "$MAPPED_ALL")" \
  | bc -l \
  | awk '{printf "%.6f", $0}'
)

{
  printf "sample\tlibrary_type\tseq_mode\traw_R1_reads\traw_R2_reads\traw_pairs\ttrimmed_pairs\tmerged_reads\tunpaired_reads\tpct_pairs_remaining\tratio_trim\tmapped_all\tmapped_unique\tendog\tduprate\tavg_readlen\tmapped_bp\tavg_cov\n"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$SAMPLE" "$LIBTYPE" "$([[ $HAS_R2 -eq 1 ]] && echo PE || echo SE)" \
    "$RAW_R1_READS" "$RAW_R2_READS" "$RAW_PAIRS" \
    "$TOTAL_TRIMMED_PAIRS" "$TOTAL_MERGED" "$TOTAL_UNPAIRED" \
    "$PCT_REMAIN" "$RATIO_TRIM" \
    "$MAPPED_ALL" "$MAPPED_UNIQUE" "$ENDOG" "$DUPRATE" \
    "$AVG_READLEN" "$MAPPED_BP" "$AVG_COV"
} > "$STATS"


if [[ $KEEP_INTERMEDIATE -eq 0 ]]; then
  log "Cleaning intermediates (work dir only)"
  rm -rf "$WORK"
  [[ -n "$TMPDIR_USER" ]] && rm -rf "$TMPDIR_USER/plainmap_${SAMPLE}" || true
fi

PIPE_T1=$(date +%s)
log "STATUS: SUCCESS"
log "Total wall time: $((PIPE_T1 - PIPE_T0)) seconds"
log "Final BAM: $OUT_BAM"
