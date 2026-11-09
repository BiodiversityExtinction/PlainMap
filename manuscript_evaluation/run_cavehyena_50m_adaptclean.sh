#!/usr/bin/env bash
set -euo pipefail

case "${1:-${PLAINMAP_CONDITION:-}}" in
  off) label=off; adaptclean_enabled=0 ;;
  on) label=on; adaptclean_enabled=1 ;;
  *) echo "Usage: $0 off|on" >&2; exit 2 ;;
esac

REPO_DIR="${PLAINMAP_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
SCRIPT_DIR="$REPO_DIR/manuscript_evaluation"
BASE="$SCRIPT_DIR/evaluation_runs"
PLAINMAP="${PLAINMAP:-$BASE/plainmap_current_20261105.sh}"
FASTQ="${FASTQ:-/home/projects2/Biodiversity_Extinction/Cave_hyena/Raw/Ancient/Ccsp043_R1.fastq.gz}"
REF="${REF:-/home/projects2/Biodiversity_Extinction/Cave_hyena/Striped_map/GCF_003009895.1_ASM300989v1_genomic.fna}"
DEDUP="${DEDUP:-$BASE/tools/dedup-release}"
ADAPTCLEAN="${ADAPTCLEAN:-$BASE/tools/AdaptClean/bin/AdaptClean}"
EVAL_ROOT="${EVAL_ROOT:-$BASE/cavehyena_50m_adaptclean}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
FASTP="${FASTP:-fastp}"
BWA="${BWA:-bwa}"
SAMTOOLS="${SAMTOOLS:-samtools}"
PYTHON="${PYTHON:-python3}"
PILOT_FRAGMENTS="${PILOT_FRAGMENTS:-50000000}"
PREFIX=Ccsp043_50m

out="$EVAL_ROOT/$label"
[[ ! -e "$out" ]] || { echo "ERROR: output already exists: $out" >&2; exit 1; }
for file in "$PLAINMAP" "$FASTQ" "$REF" "$REF.bwt" "$DEDUP" "$ADAPTCLEAN"; do
  [[ -s "$file" ]] || { echo "ERROR: required file is missing or empty: $file" >&2; exit 1; }
done
[[ -x "$DEDUP" && -x "$ADAPTCLEAN" ]] || {
  echo "ERROR: DeDup or AdaptClean is not executable" >&2; exit 1;
}
mkdir -p "$out"
printf '%s\n' "$FASTQ" > "$out/manifest.txt"

{
  printf 'sample\tCcsp043\n'
  printf 'source_fastq\t%s\n' "$FASTQ"
  printf 'source_fastq_bytes\t%s\n' "$(stat -c %s "$FASTQ")"
  printf 'reference\t%s\n' "$REF"
  printf 'reference_assembly\tGCF_003009895.1\n'
  printf 'plainmap\t%s\n' "$PLAINMAP"
  printf 'plainmap_sha256\t%s\n' "$(sha256sum "$PLAINMAP" | cut -d ' ' -f 1)"
  printf 'adaptclean_enabled\t%s\n' "$adaptclean_enabled"
  printf 'adaptclean_sha256\t%s\n' "$(sha256sum "$ADAPTCLEAN" | cut -d ' ' -f 1)"
  printf 'adaptclean_source_commit\t3fffb7b60e2032f4b2ed4203c320cfed1ccea868\n'
  printf 'pilot_fragments\t%s\n' "$PILOT_FRAGMENTS"
  printf 'library_type\tancient-fast\n'
  printf 'dedup_method\thybrid\n'
  printf 'threads\t%s\n' "$THREADS"
  printf 'fastp_version\t%s\n' "$("$FASTP" --version 2>&1)"
  printf 'bwa_version\t%s\n' "$("$BWA" 2>&1 | awk '/Version:/ {print $2; exit}' || true)"
  printf 'samtools_version\t%s\n' "$("$SAMTOOLS" --version | sed -n '1p')"
} > "$out/provenance.tsv"

adaptclean_args=()
if [[ "$adaptclean_enabled" -eq 1 ]]; then
  adaptclean_args=(--run-adaptclean --adaptclean "$ADAPTCLEAN")
fi

/usr/bin/time -f $'wall_seconds\t%e\nuser_cpu_seconds\t%U\nsystem_cpu_seconds\t%S\nmax_rss_kb\t%M\nexit_status\t%x' \
  -o "$out/time.tsv" \
  bash "$PLAINMAP" \
    --manifest "$out/manifest.txt" \
    --prefix "$PREFIX" --ref "$REF" --outdir "$out" \
    --threads "$THREADS" --library-type ancient-fast \
    --dedup-method hybrid --dedup "$DEDUP" \
    --pilot-fragments "$PILOT_FRAGMENTS" \
    --keep-intermediate \
    "${adaptclean_args[@]}" \
    --fastp "$FASTP" --bwa "$BWA" --samtools "$SAMTOOLS" --python "$PYTHON" \
    > "$out/driver.log" 2>&1

bam="$out/$PREFIX.final.RG.bam"
"$SAMTOOLS" quickcheck "$bam"
[[ "$("$SAMTOOLS" view -c -f 4 "$bam")" -eq 0 ]]
[[ "$("$SAMTOOLS" view -c -f 1024 "$bam")" -eq 0 ]]

pool="$out/$PREFIX/work/raw/$PREFIX.pool.SE.all.fastq.gz"
[[ -s "$pool" ]] || { echo "ERROR: missing trimmed SE pool: $pool" >&2; exit 1; }
{
  printf 'length_bp\treads\n'
  gzip -dc "$pool" | awk 'NR % 4 == 2 {counts[length($0)]++} END {for (len in counts) print len "\t" counts[len]}' | sort -n
} > "$out/trimmed_length_histogram.tsv"

echo "Completed and validated cave hyena 50M AdaptClean $label: $out"
