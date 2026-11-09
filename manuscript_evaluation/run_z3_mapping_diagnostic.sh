#!/usr/bin/env bash
set -euo pipefail

: "${PLAINMAP:?Set PLAINMAP to the frozen PlainMap script}"
: "${MANIFEST:?Set MANIFEST to the validated Z3 manifest}"
: "${REF:?Set REF to the indexed California sea lion reference}"
: "${DEDUP:?Set DEDUP to the pinned DeDup wrapper}"
: "${SLURM_ARRAY_TASK_ID:?Submit as a 0-1 Slurm array}"
[[ "$SLURM_ARRAY_TASK_ID" =~ ^[0-1]$ ]] || { echo "Invalid array task ID" >&2; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EVAL_BASE="${EVAL_BASE:-$SCRIPT_DIR/evaluation_runs/z3_mapping_diagnostic}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
PILOT_FRAGMENTS="${PILOT_FRAGMENTS:-10000000}"
FASTP="${FASTP:-fastp}"
BWA="${BWA:-bwa}"
SAMTOOLS="${SAMTOOLS:-samtools}"
PYTHON="${PYTHON:-python3}"
label=(unchunked chunk_10m)
chunk_size=(0 10000000)
out="$EVAL_BASE/${label[SLURM_ARRAY_TASK_ID]}"
[[ ! -e "$out" ]] || { echo "ERROR: output exists: $out" >&2; exit 1; }
mkdir -p "$out"

{
  printf 'plainmap_sha256\t%s\n' "$(sha256sum "$PLAINMAP" | cut -d ' ' -f 1)"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$MANIFEST" | cut -d ' ' -f 1)"
  printf 'reference\t%s\n' "$REF"
  printf 'threads\t%s\n' "$THREADS"
  printf 'pilot_fragments_per_unit\t%s\n' "$PILOT_FRAGMENTS"
  printf 'chunk_size\t%s\n' "${chunk_size[SLURM_ARRAY_TASK_ID]}"
  printf 'library_type\tmodern\n'
  printf 'dedup_method\thybrid\n'
  printf 'dedup_version\t%s\n' "$("$DEDUP" --version 2>&1 | sed -n '1p')"
  printf 'fastp_version\t%s\n' "$("$FASTP" --version 2>&1)"
  printf 'bwa_version\t%s\n' "$("$BWA" 2>&1 | awk '/Version:/ {print $2; exit}' || true)"
  printf 'samtools_version\t%s\n' "$("$SAMTOOLS" --version | sed -n '1p')"
} > "$out/provenance.tsv"

chunk_arg=()
if [[ "${chunk_size[SLURM_ARRAY_TASK_ID]}" -gt 0 ]]; then
  chunk_arg=(--max-reads-per-chunk "${chunk_size[SLURM_ARRAY_TASK_ID]}")
fi

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
    --keep-intermediate --no-delete-as-you-go \
    "${chunk_arg[@]}" \
    --fastp "$FASTP" --bwa "$BWA" --samtools "$SAMTOOLS" --python "$PYTHON" \
    > "$out/driver.log" 2>&1

work="$out/Z3_large_modern/work"
pooled="$work/raw/Z3_large_modern.pool"
merged="$work/final/Z3_large_modern.merged.bam"
final="$out/Z3_large_modern.final.RG.bam"
"$SAMTOOLS" quickcheck "$merged" "$final"

{
  printf 'file\tsha256_uncompressed_or_sorted_records\n'
  for type in PE.R1 PE.R2 SE.all; do
    path="$pooled.$type.fastq.gz"
    [[ -s "$path" ]] || { echo "Missing pooled FASTQ: $path" >&2; exit 1; }
    printf '%s\t%s\n' "$type" "$(gzip -dc "$path" | sha256sum | cut -d ' ' -f 1)"
  done
  for type in merged final; do
    if [[ "$type" == merged ]]; then path="$merged"; else path="$final"; fi
    printf '%s\t%s\n' "$type" "$("$SAMTOOLS" view "$path" | LC_ALL=C sort | sha256sum | cut -d ' ' -f 1)"
  done
} > "$out/stage_checksums.tsv"

echo "Completed mapping diagnostic ${label[SLURM_ARRAY_TASK_ID]}: $out"
