#!/usr/bin/env bash
set -euo pipefail

: "${PLAINMAP:?Set PLAINMAP to the frozen PlainMap script}"
: "${MANIFEST:?Set MANIFEST to the validated Z3 manifest}"
: "${REF:?Set REF to the indexed California sea lion reference}"
: "${DEDUP:?Set DEDUP to the pinned DeDup wrapper}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EVAL_BASE="${EVAL_BASE:-$SCRIPT_DIR/evaluation_runs/z3_chunk_size_sweep}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
PILOT_FRAGMENTS="${PILOT_FRAGMENTS:-10000000}"
FASTP="${FASTP:-fastp}"
BWA="${BWA:-bwa}"
SAMTOOLS="${SAMTOOLS:-samtools}"
PYTHON="${PYTHON:-python3}"

sizes=(5000000 10000000)
if [[ -z "${CHUNK_SIZE:-}" ]]; then
  : "${SLURM_ARRAY_TASK_ID:?Set CHUNK_SIZE or submit as a 0-1 Slurm array}"
  [[ "$SLURM_ARRAY_TASK_ID" =~ ^[0-1]$ ]] || { echo "Invalid array task ID" >&2; exit 1; }
  CHUNK_SIZE="${sizes[SLURM_ARRAY_TASK_ID]}"
fi
[[ "$CHUNK_SIZE" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid CHUNK_SIZE: $CHUNK_SIZE" >&2; exit 1; }

out="$EVAL_BASE/size_$CHUNK_SIZE"
[[ ! -e "$out" ]] || { echo "ERROR: output exists: $out" >&2; exit 1; }
mkdir -p "$out"

{
  printf 'plainmap\t%s\n' "$PLAINMAP"
  printf 'plainmap_sha256\t%s\n' "$(sha256sum "$PLAINMAP" | cut -d ' ' -f 1)"
  printf 'manifest\t%s\n' "$MANIFEST"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$MANIFEST" | cut -d ' ' -f 1)"
  printf 'reference\t%s\n' "$REF"
  printf 'threads\t%s\n' "$THREADS"
  printf 'pilot_fragments_per_unit\t%s\n' "$PILOT_FRAGMENTS"
  printf 'chunk_size\t%s\n' "$CHUNK_SIZE"
  printf 'library_type\tmodern\n'
  printf 'dedup_method\thybrid\n'
  printf 'dedup_executable\t%s\n' "$DEDUP"
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
    --max-reads-per-chunk "$CHUNK_SIZE" \
    --fastp "$FASTP" --bwa "$BWA" --samtools "$SAMTOOLS" --python "$PYTHON" \
    > "$out/driver.log" 2>&1

bam="$out/Z3_large_modern.final.RG.bam"
"$SAMTOOLS" quickcheck "$bam"
[[ "$("$SAMTOOLS" view -c -f 4 "$bam")" -eq 0 ]]
[[ "$("$SAMTOOLS" view -c -f 1024 "$bam")" -eq 0 ]]

"$PYTHON" - "$out" "$CHUNK_SIZE" <<'PY'
import csv
import pathlib
import sys

out = pathlib.Path(sys.argv[1])
with (out / "Z3_large_modern.plainmap.stats.tsv").open() as stream:
    stats = next(csv.DictReader(stream, delimiter="\t"))
times = dict(line.rstrip().split("\t", 1) for line in (out / "time.tsv").open())
with (out / "summary.tsv").open("w") as stream:
    writer = csv.writer(stream, delimiter="\t")
    columns = ("chunk_size", "wall_seconds", "user_cpu_seconds", "system_cpu_seconds",
               "max_rss_kb", "raw_fragments", "trimmed_fragments", "merged_reads",
               "unpaired_reads", "mapped_fragments_unique")
    writer.writerow(columns)
    writer.writerow((sys.argv[2], *(times[key] if key in times else stats[key]
                                    for key in columns[1:])))
PY

echo "Completed and validated chunk size $CHUNK_SIZE: $out"
