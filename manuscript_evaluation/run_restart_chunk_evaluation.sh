#!/usr/bin/env bash
set -euo pipefail

: "${PLAINMAP:?Set PLAINMAP to plainmap.sh}"
: "${DATA_DIR:?Set DATA_DIR to the evaluation FASTQ directory}"
: "${REF:?Set REF to GCF_009762305.2_mZalCal1.pri.v2_genomic.fna}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EVAL_ROOT="${EVAL_ROOT:-$SCRIPT_DIR/evaluation_runs/restart_chunk}"
FASTP="${FASTP:-fastp}"
BWA="${BWA:-bwa}"
SAMTOOLS="${SAMTOOLS:-samtools}"
PYTHON="${PYTHON:-python3}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
PREFIX=Z3_modern_pilot1M_validation
MAX_READS=100000
MANIFEST="$EVAL_ROOT/Z3_manifest_test.resolved.txt"
UNCHUNKED_OUT="$EVAL_ROOT/results/unchunked"
CHUNKED_OUT="$EVAL_ROOT/results/chunked"
INTERRUPTED_OUT="$EVAL_ROOT/results/interrupted"
LOG_DIR="$EVAL_ROOT/logs"
VALIDATION_DIR="$EVAL_ROOT/validation"

[[ ! -e "$EVAL_ROOT" ]] || {
  echo "ERROR: evaluation output already exists: $EVAL_ROOT" >&2
  echo "Choose a new EVAL_ROOT or remove the previous evaluation explicitly." >&2
  exit 1
}

for cmd in setsid /usr/bin/time sha256sum sort; do
  command -v "$cmd" >/dev/null 2>&1 || {
    echo "ERROR: required command not found: $cmd" >&2
    exit 1
  }
done

mkdir -p "$UNCHUNKED_OUT" "$CHUNKED_OUT" "$INTERRUPTED_OUT" "$LOG_DIR" "$VALIDATION_DIR"
DATA_DIR="$DATA_DIR" "$SCRIPT_DIR/resolve_manifest.sh" "$MANIFEST" >/dev/null

COMMON_ARGS=(
  --manifest "$MANIFEST"
  --prefix "$PREFIX"
  --ref "$REF"
  --threads "$THREADS"
  --library-type modern
  --pilot-fragments 1000000
  --fastp "$FASTP"
  --bwa "$BWA"
  --samtools "$SAMTOOLS"
  --python "$PYTHON"
)

{
  printf 'plainmap_version\t%s\n' "$("$PLAINMAP" --help | awk 'NR==1{print $2}')"
  printf 'plainmap_sha256\t%s\n' "$(sha256sum "$PLAINMAP" | awk '{print $1}')"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$MANIFEST" | awk '{print $1}')"
  printf 'reference\t%s\n' "$REF"
  printf 'threads\t%s\n' "$THREADS"
  printf 'pilot_fragments_per_unit\t1000000\n'
  printf 'forced_chunk_size\t%s\n' "$MAX_READS"
} > "$VALIDATION_DIR/provenance.tsv"

run_timed() {
  local label="$1"
  local outdir="$2"
  shift 2

  /usr/bin/time \
    -f $'wall_seconds\t%e\nmax_rss_kb\t%M\nexit_status\t%x' \
    -o "$VALIDATION_DIR/${label}.time.tsv" \
    "$PLAINMAP" "${COMMON_ARGS[@]}" --outdir "$outdir" "$@" \
    > "$LOG_DIR/${label}.out" 2>&1
}

alignment_checksum() {
  "$SAMTOOLS" view "$1" | LC_ALL=C sort | sha256sum | awk '{print $1}'
}

validate_bam() {
  local bam="$1"
  "$SAMTOOLS" quickcheck "$bam"
  test -s "${bam}.bai"
  [[ "$("$SAMTOOLS" view -c -f 4 "$bam")" -eq 0 ]]
  [[ "$("$SAMTOOLS" view -c -f 1024 "$bam")" -eq 0 ]]
  "$SAMTOOLS" view -H "$bam" | grep -F 'SO:coordinate' >/dev/null
  "$SAMTOOLS" view -H "$bam" | grep -F $'@RG\tID:'"$PREFIX"$'\t' >/dev/null
}

count_valid_chunk_bams() {
  local count=0
  local bam
  shopt -s nullglob
  for bam in "$INTERRUPTED_OUT/$PREFIX/work/bams/"*.bam; do
    if "$SAMTOOLS" quickcheck "$bam" >/dev/null 2>&1; then
      count=$((count + 1))
    fi
  done
  shopt -u nullglob
  printf '%s\n' "$count"
}

echo "Running unchunked baseline"
run_timed unchunked "$UNCHUNKED_OUT"
UNCHUNKED_BAM="$UNCHUNKED_OUT/${PREFIX}.final.RG.bam"
validate_bam "$UNCHUNKED_BAM"

echo "Running uninterrupted forced-chunk evaluation"
run_timed chunked "$CHUNKED_OUT" --max-reads-per-chunk "$MAX_READS"
CHUNKED_BAM="$CHUNKED_OUT/${PREFIX}.final.RG.bam"
CHUNKED_STATS="$CHUNKED_OUT/${PREFIX}.plainmap.stats.tsv"
validate_bam "$CHUNKED_BAM"

echo "Starting run that will be terminated during partial mapping"
INTERRUPTED_ARGS=(
  "${COMMON_ARGS[@]}"
  --outdir "$INTERRUPTED_OUT"
  --max-reads-per-chunk "$MAX_READS"
  --keep-intermediate
)

setsid "$PLAINMAP" "${INTERRUPTED_ARGS[@]}" \
  > "$LOG_DIR/interrupted_attempt.out" 2>&1 &
ATTEMPT_PID=$!

PE_IDS="$INTERRUPTED_OUT/$PREFIX/work/.checkpoints/${PREFIX}.mapchunk_ids.pe.txt"
SE_IDS="$INTERRUPTED_OUT/$PREFIX/work/.checkpoints/${PREFIX}.mapchunk_ids.se.txt"
EXPECTED_CHUNKS=0
COMPLETED_AT_TERMINATION=0
TERMINATION_REQUESTED=0

# Poll for up to 40 minutes. Terminate only after at least one valid chunk is
# complete and at least one mapping chunk remains unfinished.
for ((poll=1; poll<=2400; poll++)); do
  if ! kill -0 "$ATTEMPT_PID" 2>/dev/null; then
    break
  fi

  if [[ -e "$PE_IDS" || -e "$SE_IDS" ]]; then
    EXPECTED_CHUNKS=$(cat "$PE_IDS" "$SE_IDS" 2>/dev/null | awk 'NF{n++} END{print n+0}')
  fi
  COMPLETED_AT_TERMINATION="$(count_valid_chunk_bams)"

  if [[ "$EXPECTED_CHUNKS" -gt 1 && "$COMPLETED_AT_TERMINATION" -gt 0 && "$COMPLETED_AT_TERMINATION" -lt "$EXPECTED_CHUNKS" ]]; then
    kill -TERM -- "-$ATTEMPT_PID"
    TERMINATION_REQUESTED=1
    break
  fi
  sleep 1
done

[[ "$TERMINATION_REQUESTED" -eq 1 ]] || {
  echo "ERROR: failed to terminate during partial mapping" >&2
  kill -TERM -- "-$ATTEMPT_PID" 2>/dev/null || true
  exit 1
}

for _ in {1..10}; do
  kill -0 -- "-$ATTEMPT_PID" 2>/dev/null || break
  sleep 1
done
kill -KILL -- "-$ATTEMPT_PID" 2>/dev/null || true

set +e
wait "$ATTEMPT_PID"
ATTEMPT_STATUS=$?
set -e

COMPLETED_AT_TERMINATION="$(count_valid_chunk_bams)"
[[ "$ATTEMPT_STATUS" -eq 143 ]]
[[ "$COMPLETED_AT_TERMINATION" -gt 0 ]]
[[ "$COMPLETED_AT_TERMINATION" -lt "$EXPECTED_CHUNKS" ]]

echo "Resuming with the identical command"
/usr/bin/time \
  -f $'wall_seconds\t%e\nmax_rss_kb\t%M\nexit_status\t%x' \
  -o "$VALIDATION_DIR/interrupted_resume.time.tsv" \
  "$PLAINMAP" "${INTERRUPTED_ARGS[@]}" \
  > "$LOG_DIR/interrupted_resume.out" 2>&1

INTERRUPTED_BAM="$INTERRUPTED_OUT/${PREFIX}.final.RG.bam"
INTERRUPTED_STATS="$INTERRUPTED_OUT/${PREFIX}.plainmap.stats.tsv"
validate_bam "$INTERRUPTED_BAM"

CHUNKED_HASH="$(alignment_checksum "$CHUNKED_BAM")"
INTERRUPTED_HASH="$(alignment_checksum "$INTERRUPTED_BAM")"
UNCHUNKED_HASH="$(alignment_checksum "$UNCHUNKED_BAM")"
[[ "$CHUNKED_HASH" == "$INTERRUPTED_HASH" ]]
cmp -s "$CHUNKED_STATS" "$INTERRUPTED_STATS"

REUSED_CHUNKS=$(grep -c 'already mapped (valid BAM exists) -> skipping' \
  "$INTERRUPTED_OUT/${PREFIX}_plainmap.log" || true)
[[ "$REUSED_CHUNKS" -ge "$COMPLETED_AT_TERMINATION" ]]

{
  printf 'check\tvalue\n'
  printf 'forced_chunk_size\t%s\n' "$MAX_READS"
  printf 'expected_mapping_chunks\t%s\n' "$EXPECTED_CHUNKS"
  printf 'chunks_completed_before_termination\t%s\n' "$COMPLETED_AT_TERMINATION"
  printf 'chunks_reused_on_resume\t%s\n' "$REUSED_CHUNKS"
  printf 'interrupted_attempt_exit_status\t%s\n' "$ATTEMPT_STATUS"
  printf 'unchunked_alignment_checksum\t%s\n' "$UNCHUNKED_HASH"
  printf 'chunked_alignment_checksum\t%s\n' "$CHUNKED_HASH"
  printf 'resumed_alignment_checksum\t%s\n' "$INTERRUPTED_HASH"
  printf 'resumed_vs_uninterrupted_alignments\tPASS\n'
  printf 'resumed_vs_uninterrupted_statistics\tPASS\n'
  printf 'final_bam_quickcheck\tPASS\n'
  printf 'final_bam_unmapped_records\t0\n'
  printf 'final_bam_duplicate_records\t0\n'
} > "$VALIDATION_DIR/restart_checks.tsv"

echo "Restart/chunk evaluation passed"
