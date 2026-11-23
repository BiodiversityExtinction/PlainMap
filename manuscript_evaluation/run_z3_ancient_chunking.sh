#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HERE/evaluation_runs"
PLAINMAP="$BASE/plainmap_integrated_v030.sh"
MANIFEST="$BASE/z3_10m_inputs_sra/manifest.txt"
REF="${REF:-/home/people/micwe/Biodiversity_Extinction/Japanese_sea_lion/Reference_genomes/GCF_009762305.2_mZalCal1.pri.v2_genomic.fna}"
DEDUP="$BASE/tools/dedup-release"
OUT="${EVAL_ROOT:-$BASE/z3_ancient_chunking_v030}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

[[ ! -e "$OUT" ]] || { echo "Output already exists: $OUT" >&2; exit 1; }
for file in "$PLAINMAP" "$MANIFEST" "$REF" "$REF.bwt" "$DEDUP"; do
  [[ -s "$file" ]] || { echo "Missing input: $file" >&2; exit 1; }
done
[[ "$(sha256sum "$PLAINMAP" | cut -d ' ' -f 1)" == \
   "$(cat "$BASE/plainmap_integrated_v030.sha256")" ]] || {
  echo "Frozen PlainMap script has changed" >&2; exit 1;
}
[[ "$(sha256sum "$MANIFEST" | cut -d ' ' -f 1)" == \
   c80a82b78bebe1cf55dff94ce84c4f54f20469bc4fbe282b467a16585ed3cfab ]] || {
  echo "Validated SRA-derived manifest has changed" >&2; exit 1;
}
mkdir -p "$OUT"
{
  printf 'field\tvalue\n'
  printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-none}"
  printf 'node\t%s\n' "$(hostname)"
  printf 'plainmap_sha256\t%s\n' "$(sha256sum "$PLAINMAP" | cut -d ' ' -f 1)"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$MANIFEST" | cut -d ' ' -f 1)"
  printf 'reference\t%s\n' "$REF"
  printf 'threads\t%s\n' "$THREADS"
  printf 'pilot_fragments_per_unit\t10000000\n'
  printf 'chunk_sizes\t1000000,10000000\n'
  printf 'minlength\t30\nmapq\t20\ndedup_method\thybrid\n'
  printf 'fastp_version\t%s\n' "$(fastp --version 2>&1)"
  printf 'bwa_version\t%s\n' "$(bwa 2>&1 | awk '/Version:/ {print $2; exit}' || true)"
  printf 'samtools_version\t%s\n' "$(samtools --version | sed -n '1p')"
  printf 'dedup_version\t%s\n' "$($DEDUP --version 2>&1 | head -1 || true)"
} > "$OUT/provenance.tsv"

for mode in ancient-fast ancient-all; do
  for size in 1000000 10000000; do
    label="${mode}_${size}"
    prefix="Z3_${mode//-/_}"
    dir="$OUT/$label"
    mkdir -p "$dir"
    echo "Starting $label on $(hostname)"
    /usr/bin/time -f $'wall_seconds\t%e\nuser_cpu_seconds\t%U\nsystem_cpu_seconds\t%S\nmax_rss_kb\t%M\nexit_status\t%x' \
      -o "$dir/time.tsv" \
      "$PLAINMAP" --manifest "$MANIFEST" --prefix "$prefix" --ref "$REF" \
        --outdir "$dir" --threads "$THREADS" --library-type "$mode" \
        --pilot-fragments 10000000 --max-reads-per-chunk "$size" \
        --minlength 30 --mapq 20 --dedup-method hybrid --dedup "$DEDUP" \
        --fastp fastp --bwa bwa --samtools samtools --python python3 \
        > "$dir/driver.log" 2>&1
    bam="$dir/$prefix.final.RG.bam"
    samtools quickcheck "$bam"
    [[ "$(samtools view -c -f 4 "$bam")" -eq 0 ]]
    [[ "$(samtools view -c -f 1024 "$bam")" -eq 0 ]]
    [[ -s "$dir/$prefix.plainmap.stats.tsv" ]]
    [[ -s "$dir/$prefix.plainmap.adapters.tsv" ]]
    echo "Validated $label"
  done
done

python3 - "$OUT" <<'PY'
import csv
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
fields = (
    "raw_fragments", "trimmed_fragments", "merged_reads",
    "pe_fragments_all", "pe_fragments_unique",
    "se_like_fragments_all", "se_like_fragments_unique",
    "mapped_fragments_unique",
)
with (root / "results.tsv").open("w") as output:
    writer = csv.writer(output, delimiter="\t")
    writer.writerow(("mode", "chunk_size", *fields, "wall_seconds", "max_rss_kb"))
    for mode in ("ancient-fast", "ancient-all"):
        for size in (1000000, 10000000):
            case = root / f"{mode}_{size}"
            prefix = f"Z3_{mode.replace('-', '_')}"
            with (case / f"{prefix}.plainmap.stats.tsv").open() as stream:
                stats = next(csv.DictReader(stream, delimiter="\t"))
            times = dict(line.rstrip().split("\t", 1) for line in (case / "time.tsv").open())
            writer.writerow((mode, size, *(stats[field] for field in fields),
                             times["wall_seconds"], times["max_rss_kb"]))
PY
cat "$OUT/results.tsv"
echo "Ancient-mode chunking comparison passed: $OUT"
