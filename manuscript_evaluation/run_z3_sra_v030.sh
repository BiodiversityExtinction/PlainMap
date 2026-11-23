#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$SCRIPT_DIR/evaluation_runs"
INPUT_ROOT="${INPUT_ROOT:-$BASE/z3_10m_inputs_sra}"
MANIFEST="$INPUT_ROOT/manifest.txt"
PLAINMAP="${PLAINMAP:-$BASE/z3_v030_28m_modes/plainmap_v0.3.0.sh}"
REF="${REF:-/home/people/micwe/Biodiversity_Extinction/Japanese_sea_lion/Reference_genomes/GCF_009762305.2_mZalCal1.pri.v2_genomic.fna}"
DEDUP="${DEDUP:-$BASE/tools/dedup-release}"
EVAL_ROOT="${EVAL_ROOT:-$BASE/z3_sra_v030}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
export PLAINMAP MANIFEST REF DEDUP THREADS

[[ ! -e "$EVAL_ROOT" ]] || {
  echo "ERROR: output exists: $EVAL_ROOT" >&2; exit 1;
}
for file in "$MANIFEST" "$PLAINMAP" "$REF" "$REF.bwt" "$DEDUP" "$INPUT_ROOT/provenance.tsv"; do
  [[ -s "$file" ]] || { echo "ERROR: required input missing or empty: $file" >&2; exit 1; }
done
[[ "$(sha256sum "$PLAINMAP" | cut -d ' ' -f 1)" == \
   adbeb8f15c8100938344ac012f68989b62095f406e381050d9dbec649d9dd492 ]] || {
  echo "ERROR: PlainMap script is not exact v0.3.0" >&2; exit 1;
}
grep -F $'mate_headers_identical\tPASS (all pairs)' "$INPUT_ROOT/provenance.tsv" >/dev/null || {
  echo "ERROR: SRA mate-header validation is missing" >&2; exit 1;
}
input_value() {
  awk -F '\t' -v key="$1" '$1 == key { print $2; exit }' "$INPUT_ROOT/provenance.tsv"
}
[[ "$(input_value manifest_sha256)" == "$(sha256sum "$MANIFEST" | cut -d ' ' -f 1)" ]] || {
  echo "ERROR: validated SRA manifest changed" >&2; exit 1;
}
for mate in 1 2; do
  file="$INPUT_ROOT/bulk/SRR25607793_${mate}.fastq.gz"
  expected="$(input_value "SRR25607793_${mate}.fastq.gz_sha256")"
  [[ -s "$file" && -n "$expected" && "$(sha256sum "$file" | cut -d ' ' -f 1)" == "$expected" ]] || {
    echo "ERROR: validated SRA mate $mate changed" >&2; exit 1;
  }
done
mkdir -p "$EVAL_ROOT"
{
  printf 'field\tvalue\n'
  printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-none}"
  printf 'node\t%s\n' "$(hostname)"
  printf 'plainmap_sha256\t%s\n' "$(sha256sum "$PLAINMAP" | cut -d ' ' -f 1)"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$MANIFEST" | cut -d ' ' -f 1)"
  printf 'sra_accession\tSRR25607793\n'
  printf 'sra_pair_count\t10000000\n'
  printf 'reference\t%s\n' "$REF"
  printf 'threads\t%s\n' "$THREADS"
  printf 'pilot_fragments_per_unit\t10000000\n'
  printf 'minlength\t30\nmapq\t20\ndedup_method\thybrid\n'
} > "$EVAL_ROOT/provenance.tsv"

echo "Running exact v0.3.0 mapping-mode comparison on SRA-derived pairs"
for mode in modern ancient-fast ancient-all; do
  EVAL_ROOT="$EVAL_ROOT/modes" PILOT_FRAGMENTS=10000000 \
    bash "$SCRIPT_DIR/run_z3_28m_modes.sh" "$mode"
done
python3 "$SCRIPT_DIR/summarize_z3_28m_modes.py" "$EVAL_ROOT/modes"

echo "Running modern unchunked and 1M-chunk conditions sequentially"
EVAL_ROOT="$EVAL_ROOT/chunking_baseline" PILOT_FRAGMENTS=10000000 \
  CHUNK_SIZE=1000000 DEDUP_METHOD=hybrid \
  bash "$SCRIPT_DIR/run_large_z3_benchmark.sh"

for size in 5000000 10000000; do
  echo "Running modern ${size}-fragment chunk condition"
  EVAL_BASE="$EVAL_ROOT/chunk_size_sweep" PILOT_FRAGMENTS=10000000 \
    CHUNK_SIZE="$size" bash "$SCRIPT_DIR/run_z3_chunk_size_sweep.sh"
done

echo "Running interrupted restart against the SRA-derived 1M-chunk baseline"
BASELINE_OUT="$EVAL_ROOT/chunking_baseline/chunked" \
  EVAL_ROOT="$EVAL_ROOT/restart" \
  bash "$SCRIPT_DIR/run_z3_28m_restart.sh"

python3 "$SCRIPT_DIR/summarize_z3_sra_v030.py" "$EVAL_ROOT"
echo "All SRA-derived Z3 v0.3.0 evaluations passed: $EVAL_ROOT"
