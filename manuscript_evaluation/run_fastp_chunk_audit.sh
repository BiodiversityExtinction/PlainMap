#!/usr/bin/env bash
set -euo pipefail

: "${R1:?Set R1 to the SRA paired-end R1 FASTQ}"
: "${R2:?Set R2 to the corresponding R2 FASTQ}"
: "${REF:?Set REF to the indexed reference FASTA}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
EVAL_ROOT="${EVAL_ROOT:-$SCRIPT_DIR/evaluation_runs/fastp_chunk_audit}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
CHUNK_SIZE="${CHUNK_SIZE:-100000}"
FASTP="${FASTP:-fastp}"
BWA="${BWA:-bwa}"
SAMTOOLS="${SAMTOOLS:-samtools}"
PYTHON="${PYTHON:-python3}"
ADAPTER_R1="${ADAPTER_R1:-GATCGGAAGAGCACACGTCTGAACTCCAGTCAC}"
ADAPTER_R2="${ADAPTER_R2:-AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT}"

[[ ! -e "$EVAL_ROOT" ]] || { echo "ERROR: output exists: $EVAL_ROOT" >&2; exit 1; }
mkdir -p "$EVAL_ROOT"

# Reproduce the v0.1.0 fastp path used for the manuscript's original table.
git -C "$REPO_ROOT" show v0.1.0:plainmap.sh > "$EVAL_ROOT/plainmap_v0.1.0.sh"
chmod +x "$EVAL_ROOT/plainmap_v0.1.0.sh"
printf '%s\n%s\n' "$R1" "$R2" > "$EVAL_ROOT/manifest.txt"

{
  printf 'plainmap_tag\tv0.1.0\n'
  printf 'plainmap_sha256\t%s\n' "$(sha256sum "$EVAL_ROOT/plainmap_v0.1.0.sh" | cut -d ' ' -f 1)"
  printf 'r1\t%s\n' "$R1"
  printf 'r2\t%s\n' "$R2"
  printf 'reference\t%s\n' "$REF"
  printf 'fastp_version\t%s\n' "$("$FASTP" --version 2>&1)"
  printf 'threads\t%s\n' "$THREADS"
  printf 'chunk_size\t%s\n' "$CHUNK_SIZE"
  printf 'fixed_adapter_r1\t%s\n' "$ADAPTER_R1"
  printf 'fixed_adapter_r2\t%s\n' "$ADAPTER_R2"
} > "$EVAL_ROOT/provenance.tsv"

run_case() {
  local label="$1"; shift
  local out="$EVAL_ROOT/$label"
  mkdir -p "$out"
  /usr/bin/time -f $'wall_seconds\t%e\nmax_rss_kb\t%M' -o "$out/time.tsv" \
    "$EVAL_ROOT/plainmap_v0.1.0.sh" \
      --manifest "$EVAL_ROOT/manifest.txt" \
      --prefix Z3_sra_chunk_audit \
      --ref "$REF" \
      --outdir "$out" \
      --threads "$THREADS" \
      --library-type modern \
      --trim-only \
      --keep-intermediate \
      --fastp "$FASTP" --bwa "$BWA" --samtools "$SAMTOOLS" --python "$PYTHON" \
      "$@" > "$out/driver.log" 2>&1
  echo "Completed $label"
}

run_case auto_unchunked
run_case auto_chunked --max-reads-per-chunk "$CHUNK_SIZE"
run_case fixed_unchunked --adapter-r1 "$ADAPTER_R1" --adapter-r2 "$ADAPTER_R2"
run_case fixed_chunked --max-reads-per-chunk "$CHUNK_SIZE" \
  --adapter-r1 "$ADAPTER_R1" --adapter-r2 "$ADAPTER_R2"

"$PYTHON" "$SCRIPT_DIR/summarize_fastp_chunk_audit.py" "$EVAL_ROOT"
