#!/usr/bin/env bash
set -euo pipefail

: "${PLAINMAP:?Set PLAINMAP to the frozen PlainMap script}"
: "${MANIFEST:?Set MANIFEST to the validated Z3 manifest}"
: "${REF:?Set REF to the indexed California sea lion reference}"
: "${DEDUP:?Set DEDUP to the pinned DeDup wrapper}"
: "${SLURM_ARRAY_TASK_ID:?Submit as a Slurm array}"

[[ "$SLURM_ARRAY_TASK_ID" =~ ^[0-9]+$ ]] || { echo "Invalid array task ID" >&2; exit 1; }
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EVAL_BASE="${EVAL_BASE:-$SCRIPT_DIR/evaluation_runs/z3_unchunked_repeats}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
PILOT_FRAGMENTS="${PILOT_FRAGMENTS:-10000000}"
FASTP="${FASTP:-fastp}"
BWA="${BWA:-bwa}"
SAMTOOLS="${SAMTOOLS:-samtools}"
PYTHON="${PYTHON:-python3}"
out="$EVAL_BASE/rep_$SLURM_ARRAY_TASK_ID"
[[ ! -e "$out" ]] || { echo "ERROR: output exists: $out" >&2; exit 1; }
mkdir -p "$out"

{
  printf 'plainmap_sha256\t%s\n' "$(sha256sum "$PLAINMAP" | cut -d ' ' -f 1)"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$MANIFEST" | cut -d ' ' -f 1)"
  printf 'reference\t%s\n' "$REF"
  printf 'threads\t%s\n' "$THREADS"
  printf 'pilot_fragments_per_unit\t%s\n' "$PILOT_FRAGMENTS"
  printf 'chunk_size\t0\n'
  printf 'library_type\tmodern\n'
  printf 'dedup_method\thybrid\n'
  printf 'dedup_version\t%s\n' "$("$DEDUP" --version 2>&1 | sed -n '1p')"
  printf 'fastp_version\t%s\n' "$("$FASTP" --version 2>&1)"
  printf 'bwa_version\t%s\n' "$("$BWA" 2>&1 | awk '/Version:/ {print $2; exit}' || true)"
  printf 'samtools_version\t%s\n' "$("$SAMTOOLS" --version | sed -n '1p')"
} > "$out/provenance.tsv"

/usr/bin/time -f $'wall_seconds\t%e\nuser_cpu_seconds\t%U\nsystem_cpu_seconds\t%S\nmax_rss_kb\t%M\nexit_status\t%x' \
  -o "$out/time.tsv" \
  "$PLAINMAP" \
    --manifest "$MANIFEST" \
    --prefix Z3_large_modern \
    --ref "$REF" \
    --outdir "$out" \
    --threads "$THREADS" \
    --library-type modern \
    --dedup-method hybrid --dedup "$DEDUP" \
    --pilot-fragments "$PILOT_FRAGMENTS" \
    --fastp "$FASTP" --bwa "$BWA" --samtools "$SAMTOOLS" --python "$PYTHON" \
    > "$out/driver.log" 2>&1

bam="$out/Z3_large_modern.final.RG.bam"
"$SAMTOOLS" quickcheck "$bam"
[[ "$("$SAMTOOLS" view -c -f 4 "$bam")" -eq 0 ]]
[[ "$("$SAMTOOLS" view -c -f 1024 "$bam")" -eq 0 ]]
"$SAMTOOLS" view "$bam" | LC_ALL=C sort | sha256sum | cut -d ' ' -f 1 \
  > "$out/alignment_records.sha256"
echo "Completed and validated repeat $SLURM_ARRAY_TASK_ID: $out"
