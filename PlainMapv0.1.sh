#!/usr/bin/env bash
# PlainMap v0.1
# Transparent, failure-aware mapping pipeline for ancient and modern DNA

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
LIBTYPE="modern"
MAX_READS_PER_CHUNK=0

MAPQ=20

DRYRUN=0
RESUME=1
VALIDATE_ONLY=0
KEEP_INTERMEDIATE=0

TMPDIR_USER=""

FASTP=${FASTP:-fastp}
BWA=${BWA:-bwa}
SAMTOOLS=${SAMTOOLS:-samtools}

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
  -library-type ancient|modern    (default: modern)
  -t, --threads INT               (default: 1)
  -minlength INT                  (default: 30)
  -mismatch FLOAT                 (ancient only; default: 0.01)
  -max-reads-per-chunk INT        (default: 0; disabled)
  --tmpdir DIR                    Custom temporary directory
  --keep-intermediate             Keep all intermediate files

Execution control:
  --resume | --no-resume
  --dry-run
  --validate

Tool overrides:
  --fastp CMD
  --bwa CMD
  --samtools CMD
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
    --tmpdir) TMPDIR_USER="$2"; shift 2 ;;
    --keep-intermediate) KEEP_INTERMEDIATE=1; shift ;;
    --resume) RESUME=1; shift ;;
    --no-resume) RESUME=0; shift ;;
    --dry-run) DRYRUN=1; shift ;;
    --validate) VALIDATE_ONLY=1; shift ;;
    --fastp) FASTP="$2"; shift 2 ;;
    --bwa) BWA="$2"; shift 2 ;;
    --samtools) SAMTOOLS="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown option $1"; usage ;;
  esac
done

[[ -z "$MANIFEST" || -z "$SAMPLE" || -z "$REF" || -z "$OUT" ]] && usage

###############################################################################
# LOGGING
###############################################################################
mkdir -p "$OUT"
LOG="$OUT/${SAMPLE}_plainmap.log"
exec > >(tee -a "$LOG") 2>&1

log() { echo "[$(date '+%F %T')] $*"; }
die() { log "ERROR: $*"; exit 1; }

run() {
  log "$*"
  [[ $DRYRUN -eq 0 ]] && eval "$*"
}

PIPE_T0=$(date +%s)

###############################################################################
# THREAD SPLIT
###############################################################################
BWA_THREADS="$THREADS"
SORT_THREADS=$(( THREADS > 2 ? THREADS / 2 : 1 ))
VIEW_THREADS=1

###############################################################################
# DIRECTORIES
###############################################################################
RAW="$OUT/raw"
CHUNKS="$OUT/chunks"
BAMS="$OUT/chunk_bams"
FINAL="$OUT/final"
CKPT="$OUT/.checkpoints"

if [[ -n "$TMPDIR_USER" ]]; then
  TMP="$TMPDIR_USER/plainmap_${SAMPLE}"
else
  TMP="$OUT/tmp"
fi

mkdir -p "$RAW" "$CHUNKS" "$BAMS" "$FINAL" "$CKPT" "$TMP"

###############################################################################
# CHECKPOINTS
###############################################################################
ckpt() { echo "$CKPT/${SAMPLE}.$1.done"; }
ckpt_ok() { [[ $RESUME -eq 1 && -f "$(ckpt "$1")" ]]; }
ckpt_mark() { [[ $DRYRUN -eq 0 && $RESUME -eq 1 ]] && date -Is > "$(ckpt "$1")"; }

###############################################################################
# FASTP WRAPPER (ROBUST ERROR HANDLING)
###############################################################################
run_fastp() {
  local stderr_file="$1"; shift

  if [[ $DRYRUN -eq 1 ]]; then
    log "[DRY-RUN] $FASTP $*"
    return 0
  fi

  : > "$stderr_file"

  set +e
  "$FASTP" "$@" 2> >(tee -a "$stderr_file" >&2)
  status=$?
  set -e

  if [[ $status -ne 0 ]]; then
    die "fastp exited with non-zero status ($status)"
  fi

  if grep -qiE "igzip|invalid gzip|premature|truncated|error" "$stderr_file"; then
    die "fastp reported gzip/igzip error (see $stderr_file)"
  fi
}

###############################################################################
# VALIDATE
###############################################################################
if [[ $VALIDATE_ONLY -eq 1 ]]; then
  log "Validation mode"
  command -v "$FASTP" >/dev/null || die "fastp not found"
  command -v "$BWA" >/dev/null || die "bwa not found"
  command -v "$SAMTOOLS" >/dev/null || die "samtools not found"
  while read -r f; do
    [[ -z "$f" || "$f" =~ ^# ]] && continue
    gzip -t "$f" >/dev/null || die "Corrupt gzip: $f"
  done < "$MANIFEST"
  log "Validation OK"
  exit 0
fi

###############################################################################
# BWA INDEX (SAFE)
###############################################################################
ensure_bwa_index() {
  local ref="$1"
  local lock="${ref}.bwa.lock"
  local idx=( "${ref}.bwt" "${ref}.pac" "${ref}.ann" "${ref}.amb" "${ref}.sa" )

  local missing=0
  for f in "${idx[@]}"; do [[ -f "$f" ]] || missing=1; done
  [[ "$missing" -eq 0 ]] && return

  if mkdir "$lock" 2>/dev/null; then
    trap 'rm -rf "$lock"' EXIT
    log "Indexing reference"
    "$BWA" index "$ref"
    rm -rf "$lock"
    trap - EXIT
  else
    log "Waiting for BWA index"
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
# CONCATENATE FASTQS
###############################################################################
if ! ckpt_ok concat; then
  log "Concatenating FASTQs"
  R1=()
  R2=()
  while read -r f; do
    [[ -z "$f" || "$f" =~ ^# ]] && continue
    hdr=$(zcat "$f" | head -n1)
    [[ "$hdr" =~ " 1:" ]] && R1+=("$f")
    [[ "$hdr" =~ " 2:" ]] && R2+=("$f")
  done < "$MANIFEST"

  [[ "${#R1[@]}" -gt 0 ]] || die "No R1 reads found"

  cat "${R1[@]}" > "$RAW/${SAMPLE}_R1.fastq.gz"
  [[ "${#R2[@]}" -gt 0 ]] && cat "${R2[@]}" > "$RAW/${SAMPLE}_R2.fastq.gz"
  ckpt_mark concat
fi

###############################################################################
# SPLIT INTO CHUNKS
###############################################################################
split_fastq() {
  local fq="$1"
  local prefix="$2"
  local reads="$3"
  local lines=$((reads * 4))
  zcat "$fq" | split -l "$lines" -d -a 4 --filter='gzip > $FILE.fastq.gz' - "$prefix"
}

if [[ "$MAX_READS_PER_CHUNK" -gt 0 ]]; then
  log "Splitting FASTQs"
  split_fastq "$RAW/${SAMPLE}_R1.fastq.gz" "$CHUNKS/R1_" "$MAX_READS_PER_CHUNK"
  [[ -f "$RAW/${SAMPLE}_R2.fastq.gz" ]] && \
    split_fastq "$RAW/${SAMPLE}_R2.fastq.gz" "$CHUNKS/R2_" "$MAX_READS_PER_CHUNK"
else
  ln -sf "$RAW/${SAMPLE}_R1.fastq.gz" "$CHUNKS/R1_0000.fastq.gz"
  [[ -f "$RAW/${SAMPLE}_R2.fastq.gz" ]] && \
    ln -sf "$RAW/${SAMPLE}_R2.fastq.gz" "$CHUNKS/R2_0000.fastq.gz"
fi

###############################################################################
# PER-CHUNK PROCESSING
###############################################################################
for r1 in "$CHUNKS"/R1_*.fastq.gz; do
  id=$(basename "$r1" | sed 's/R1_//;s/.fastq.gz//')
  r2="$CHUNKS/R2_${id}.fastq.gz"
  bam="$BAMS/${SAMPLE}.${id}.bam"
  fastp_err="$TMP/fastp_${id}.stderr.log"

  log "Processing chunk $id"

  if [[ "$LIBTYPE" == "ancient" ]]; then
    merged="$TMP/${SAMPLE}.${id}.merged.fastq"
    run_fastp "$fastp_err" -g -l "$MINLEN" -i "$r1" -I "$r2" --merged_out "$merged"
    [[ -s "$merged" ]] || die "fastp produced empty merged FASTQ"
    $BWA aln -l 999 -n "$MISMATCH" -t "$BWA_THREADS" "$REF" "$merged" |
      $BWA samse "$REF" - "$merged" |
      $SAMTOOLS view -q "$MAPQ" -u |
      $SAMTOOLS sort -@ "$SORT_THREADS" -o "$bam"
    rm -f "$merged"
  else
    trim1="$TMP/${SAMPLE}.${id}.R1.fastq"
    trim2="$TMP/${SAMPLE}.${id}.R2.fastq"
    run_fastp "$fastp_err" -g -l "$MINLEN" -i "$r1" -I "$r2" --out1 "$trim1" --out2 "$trim2"
    [[ -s "$trim1" && -s "$trim2" ]] || die "fastp produced empty trimmed FASTQ"
    $BWA mem -t "$BWA_THREADS" "$REF" "$trim1" "$trim2" |
      $SAMTOOLS view -q "$MAPQ" -u |
      $SAMTOOLS sort -@ "$SORT_THREADS" -o "$bam"
    rm -f "$trim1" "$trim2"
  fi
done

###############################################################################
# MERGE + DEDUP
###############################################################################
log "Merging and removing duplicates"
$SAMTOOLS merge -@ "$SORT_THREADS" "$FINAL/${SAMPLE}.merged.bam" "$BAMS"/*.bam
$SAMTOOLS sort -n "$FINAL/${SAMPLE}.merged.bam" |
  $SAMTOOLS markdup -r -s - "$FINAL/${SAMPLE}.dedup.bam"
$SAMTOOLS sort "$FINAL/${SAMPLE}.dedup.bam" -o "$FINAL/${SAMPLE}.final.bam"
$SAMTOOLS index "$FINAL/${SAMPLE}.final.bam"

###############################################################################
# READ GROUPS
###############################################################################
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina\tLB:${SAMPLE}"
$SAMTOOLS addreplacerg -r "$RG" \
  -o "$OUT/${SAMPLE}.final.RG.bam" \
  "$FINAL/${SAMPLE}.final.bam"
$SAMTOOLS index "$OUT/${SAMPLE}.final.RG.bam"

###############################################################################
# CLEANUP
###############################################################################
if [[ "$KEEP_INTERMEDIATE" -eq 0 ]]; then
  log "Cleaning intermediate files"
  rm -rf "$RAW" "$CHUNKS" "$BAMS" "$TMP" "$FINAL"
fi

PIPE_T1=$(date +%s)
log "STATUS: SUCCESS"
log "Total wall time: $((PIPE_T1 - PIPE_T0)) seconds"
log "Final BAM: $OUT/${SAMPLE}.final.RG.bam"
