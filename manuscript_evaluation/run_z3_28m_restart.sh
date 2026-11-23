#!/usr/bin/env bash
set -euo pipefail

: "${PLAINMAP:?Set PLAINMAP to the exact v0.3.0 script}"
: "${MANIFEST:?Set MANIFEST to the 28M-fragment Z3 manifest}"
: "${REF:?Set REF to the California sea lion reference}"
: "${DEDUP:?Set DEDUP to the pinned DeDup executable}"
: "${BASELINE_OUT:?Set BASELINE_OUT to the completed 1M-chunk modern run}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${EVAL_ROOT:-$SCRIPT_DIR/evaluation_runs/z3_v030_28m_restart}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
FASTP="${FASTP:-fastp}"
BWA="${BWA:-bwa}"
SAMTOOLS="${SAMTOOLS:-samtools}"
PYTHON="${PYTHON:-python3}"
PREFIX=Z3_large_modern
CHUNK_SIZE=1000000
PILOT_FRAGMENTS=10000000
RUN_OUT="$OUT/resumed"
BASELINE_BAM="$BASELINE_OUT/$PREFIX.final.RG.bam"
BASELINE_STATS="$BASELINE_OUT/$PREFIX.plainmap.stats.tsv"
BASELINE_PROVENANCE="$BASELINE_OUT/../provenance.tsv"
RESUMED_BAM="$RUN_OUT/$PREFIX.final.RG.bam"
RESUMED_STATS="$RUN_OUT/$PREFIX.plainmap.stats.tsv"
LOG="$RUN_OUT/${PREFIX}_plainmap.log"

[[ ! -e "$OUT" ]] || { echo "ERROR: output exists: $OUT" >&2; exit 1; }
for file in "$PLAINMAP" "$MANIFEST" "$REF" "$REF.bwt" "$DEDUP" "$BASELINE_BAM" "$BASELINE_STATS" "$BASELINE_PROVENANCE"; do
  [[ -s "$file" ]] || { echo "ERROR: required file is missing or empty: $file" >&2; exit 1; }
done
[[ "$(sha256sum "$PLAINMAP" | cut -d ' ' -f 1)" == \
   adbeb8f15c8100938344ac012f68989b62095f406e381050d9dbec649d9dd492 ]] || {
  echo "ERROR: PLAINMAP is not the exact v0.3.0 release" >&2; exit 1;
}
baseline_value() {
  awk -F '\t' -v key="$1" '$1 == key { print $2; exit }' "$BASELINE_PROVENANCE"
}
baseline_script_hash="$(baseline_value plainmap_sha256)"
[[ ( "$baseline_script_hash" == \
      2723d9eed4ac2d39fe33af11c0cd2f3046362b2963f5632b1d2662daf91fedb6 ||
      "$baseline_script_hash" == \
      adbeb8f15c8100938344ac012f68989b62095f406e381050d9dbec649d9dd492 ) &&
   "$(baseline_value manifest_sha256)" == "$(sha256sum "$MANIFEST" | cut -d ' ' -f 1)" &&
   "$(baseline_value reference)" == "$REF" &&
   "$(baseline_value threads)" == "$THREADS" &&
   "$(baseline_value pilot_fragments_per_unit)" == "$PILOT_FRAGMENTS" &&
   "$(baseline_value chunk_size)" == "$CHUNK_SIZE" &&
   "$(baseline_value dedup_method)" == hybrid ]] || {
  echo "ERROR: baseline provenance does not match this evaluation" >&2; exit 1;
}
mkdir -p "$RUN_OUT"

{
  printf 'plainmap_sha256\t%s\n' "$(sha256sum "$PLAINMAP" | cut -d ' ' -f 1)"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$MANIFEST" | cut -d ' ' -f 1)"
  printf 'reference\t%s\n' "$REF"
  printf 'baseline\t%s\n' "$BASELINE_OUT"
  printf 'mode\tmodern\nthreads\t%s\npilot_fragments_per_unit\t%s\n' "$THREADS" "$PILOT_FRAGMENTS"
  printf 'chunk_size\t%s\ndedup_method\thybrid\n' "$CHUNK_SIZE"
  printf 'fastp_version\t%s\n' "$("$FASTP" --version 2>&1)"
  printf 'bwa_version\t%s\n' "$("$BWA" 2>&1 | awk '/Version:/ {print $2; exit}' || true)"
  printf 'samtools_version\t%s\n' "$("$SAMTOOLS" --version | sed -n '1p')"
  printf 'dedup_version\t%s\n' "$("$DEDUP" --version 2>&1 | sed -n '1p')"
  printf 'baseline_plainmap_sha256\t%s\n' "$(baseline_value plainmap_sha256)"
  printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-none}"
} > "$OUT/provenance.tsv"

ARGS=(
  --manifest "$MANIFEST" --prefix "$PREFIX" --ref "$REF" --outdir "$RUN_OUT"
  --threads "$THREADS" --library-type modern --pilot-fragments "$PILOT_FRAGMENTS"
  --max-reads-per-chunk "$CHUNK_SIZE" --dedup-method hybrid --dedup "$DEDUP"
  --keep-intermediate
  --fastp "$FASTP" --bwa "$BWA" --samtools "$SAMTOOLS" --python "$PYTHON"
)

count_valid_chunks() {
  local bam count=0
  shopt -s nullglob
  for bam in "$RUN_OUT/$PREFIX/work/bams/"*.bam; do
    "$SAMTOOLS" quickcheck "$bam" >/dev/null 2>&1 && count=$((count + 1))
  done
  shopt -u nullglob
  printf '%s\n' "$count"
}

echo "Starting exact v0.3.0 run; will interrupt after a completed mapping chunk"
setsid bash "$PLAINMAP" "${ARGS[@]}" > "$OUT/first_attempt.log" 2>&1 &
PID=$!
TERMINATED=0
EXPECTED=0
COMPLETED=0
PE_IDS="$RUN_OUT/$PREFIX/work/.checkpoints/$PREFIX.mapchunk_ids.pe.txt"
SE_IDS="$RUN_OUT/$PREFIX/work/.checkpoints/$PREFIX.mapchunk_ids.se.txt"
for ((poll=0; poll<3600; poll++)); do
  kill -0 "$PID" 2>/dev/null || break
  if [[ -e "$PE_IDS" || -e "$SE_IDS" ]]; then
    id_lists=()
    [[ ! -f "$PE_IDS" ]] || id_lists+=("$PE_IDS")
    [[ ! -f "$SE_IDS" ]] || id_lists+=("$SE_IDS")
    EXPECTED=$(awk 'NF { n++ } END { print n+0 }' "${id_lists[@]}")
  fi
  COMPLETED="$(count_valid_chunks)"
  if ((EXPECTED > 1 && COMPLETED > 0 && COMPLETED < EXPECTED)); then
    kill -TERM -- "-$PID"
    TERMINATED=1
    break
  fi
  sleep 1
done

if ((TERMINATED == 0)); then
  kill -TERM -- "-$PID" 2>/dev/null || true
  wait "$PID" 2>/dev/null || true
  echo "ERROR: could not interrupt during partial mapping" >&2
  exit 1
fi
for _ in {1..10}; do
  kill -0 -- "-$PID" 2>/dev/null || break
  sleep 1
done
kill -KILL -- "-$PID" 2>/dev/null || true
set +e
wait "$PID"
FIRST_STATUS=$?
set -e
COMPLETED="$(count_valid_chunks)"
((FIRST_STATUS == 143 && COMPLETED > 0 && COMPLETED < EXPECTED)) || {
  echo "ERROR: interrupted run did not retain valid partial mapping chunks" >&2
  exit 1
}

echo "Resuming with identical command: $COMPLETED/$EXPECTED valid chunks already exist"
/usr/bin/time -f $'wall_seconds\t%e\nuser_cpu_seconds\t%U\nsystem_cpu_seconds\t%S\nmax_process_rss_kb\t%M\nexit_status\t%x' \
  -o "$OUT/resume_time.tsv" \
  bash "$PLAINMAP" "${ARGS[@]}" > "$OUT/resume.log" 2>&1

"$SAMTOOLS" quickcheck "$RESUMED_BAM"
[[ -s "$RESUMED_BAM.bai" && -s "$RESUMED_STATS" ]]
[[ "$("$SAMTOOLS" view -c -f 4 "$RESUMED_BAM")" -eq 0 ]]
[[ "$("$SAMTOOLS" view -c -f 1024 "$RESUMED_BAM")" -eq 0 ]]
REUSED=$(grep -c 'already mapped (valid BAM exists) -> skipping' "$LOG" || true)
((REUSED >= COMPLETED)) || { echo "ERROR: completed chunks were not reused" >&2; exit 1; }

alignment_hash() {
  "$SAMTOOLS" view "$1" | LC_ALL=C sort | sha256sum | cut -d ' ' -f 1
}
BASELINE_HASH="$(alignment_hash "$BASELINE_BAM")"
RESUMED_HASH="$(alignment_hash "$RESUMED_BAM")"
[[ "$BASELINE_HASH" == "$RESUMED_HASH" ]] || {
  echo "ERROR: resumed alignments differ from uninterrupted baseline" >&2; exit 1;
}
cmp -s "$BASELINE_STATS" "$RESUMED_STATS" || {
  echo "ERROR: resumed statistics differ from uninterrupted baseline" >&2; exit 1;
}

{
  printf 'check\tvalue\n'
  printf 'raw_fragments\t27943970\n'
  printf 'mapping_chunks_expected\t%s\n' "$EXPECTED"
  printf 'chunks_valid_at_interruption\t%s\n' "$COMPLETED"
  printf 'chunks_reused_on_resume\t%s\n' "$REUSED"
  printf 'interrupted_attempt_exit_status\t%s\n' "$FIRST_STATUS"
  printf 'baseline_alignment_checksum\t%s\n' "$BASELINE_HASH"
  printf 'resumed_alignment_checksum\t%s\n' "$RESUMED_HASH"
  printf 'resumed_vs_uninterrupted_alignments\tPASS\n'
  printf 'resumed_vs_uninterrupted_statistics\tPASS\n'
  printf 'final_bam_quickcheck\tPASS\n'
  printf 'final_bam_unmapped_records\t0\n'
  printf 'final_bam_duplicate_records\t0\n'
} > "$OUT/restart_checks.tsv"
echo "Exact v0.3.0 large-workload restart validation passed"
