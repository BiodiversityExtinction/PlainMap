#!/usr/bin/env bash
set -euo pipefail

: "${PLAINMAP:?Set PLAINMAP to the PlainMap version being evaluated}"
: "${MANIFEST:?Set MANIFEST to the complete Z3 FASTQ manifest}"
: "${REF:?Set REF to the indexed reference FASTA}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EVAL_ROOT="${EVAL_ROOT:-$SCRIPT_DIR/evaluation_runs/large_z3_benchmark}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
PILOT_FRAGMENTS="${PILOT_FRAGMENTS:-10000000}"
CHUNK_SIZE="${CHUNK_SIZE:-1000000}"
DEDUP_METHOD="${DEDUP_METHOD:-}"
FASTP="${FASTP:-fastp}"
BWA="${BWA:-bwa}"
SAMTOOLS="${SAMTOOLS:-samtools}"
PYTHON="${PYTHON:-python3}"

[[ ! -e "$EVAL_ROOT" ]] || { echo "ERROR: output exists: $EVAL_ROOT" >&2; exit 1; }
mkdir -p "$EVAL_ROOT"

{
  printf 'plainmap\t%s\n' "$PLAINMAP"
  printf 'plainmap_sha256\t%s\n' "$(sha256sum "$PLAINMAP" | cut -d ' ' -f 1)"
  printf 'manifest\t%s\n' "$MANIFEST"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$MANIFEST" | cut -d ' ' -f 1)"
  printf 'reference\t%s\n' "$REF"
  printf 'threads\t%s\n' "$THREADS"
  printf 'pilot_fragments_per_unit\t%s\n' "$PILOT_FRAGMENTS"
  printf 'chunk_size\t%s\n' "$CHUNK_SIZE"
  printf 'dedup_method\t%s\n' "${DEDUP_METHOD:-script_default}"
  printf 'dedup_executable\t%s\n' "${DEDUP:-script_default}"
  if [[ -n "${DEDUP:-}" ]]; then
    printf 'dedup_version\t%s\n' "$("$DEDUP" --version 2>&1 | head -1)"
  fi
  printf 'fastp_version\t%s\n' "$("$FASTP" --version 2>&1)"
  printf 'bwa_version\t%s\n' "$("$BWA" 2>&1 | awk '/Version:/ {print $2; exit}' || true)"
  printf 'samtools_version\t%s\n' "$("$SAMTOOLS" --version | sed -n '1p')"
} > "$EVAL_ROOT/provenance.tsv"

run_case() {
  local label="$1"; shift
  local out="$EVAL_ROOT/$label"
  local dedup_args=()
  [[ -z "$DEDUP_METHOD" ]] || dedup_args=( --dedup-method "$DEDUP_METHOD" )
  mkdir -p "$out"
  /usr/bin/time -f $'wall_seconds\t%e\nuser_cpu_seconds\t%U\nsystem_cpu_seconds\t%S\nmax_rss_kb\t%M\nexit_status\t%x' \
    -o "$out/time.tsv" \
    "$PLAINMAP" \
      --manifest "$MANIFEST" \
      --prefix Z3_large_modern \
      --ref "$REF" \
      --outdir "$out" \
      --threads "$THREADS" \
      --library-type modern \
      "${dedup_args[@]}" \
      --pilot-fragments "$PILOT_FRAGMENTS" \
      --fastp "$FASTP" --bwa "$BWA" --samtools "$SAMTOOLS" --python "$PYTHON" \
      "$@" > "$out/driver.log" 2>&1
  local bam="$out/Z3_large_modern.final.RG.bam"
  "$SAMTOOLS" quickcheck "$bam"
  [[ "$("$SAMTOOLS" view -c -f 4 "$bam")" -eq 0 ]]
  [[ "$("$SAMTOOLS" view -c -f 1024 "$bam")" -eq 0 ]]
  echo "Completed and validated $label"
}

run_case unchunked
run_case chunked --max-reads-per-chunk "$CHUNK_SIZE"

"$PYTHON" - "$EVAL_ROOT" <<'PY'
import csv
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
with (root / "summary.tsv").open("w") as output:
    writer = csv.writer(output, delimiter="\t")
    writer.writerow(("case", "wall_seconds", "user_cpu_seconds", "system_cpu_seconds",
                     "max_rss_kb", "raw_fragments", "trimmed_fragments", "mapped_fragments_unique"))
    for label in ("unchunked", "chunked"):
        out = root / label
        times = dict(line.rstrip().split("\t", 1) for line in (out / "time.tsv").open())
        with (out / "Z3_large_modern.plainmap.stats.tsv").open() as stream:
            stats = next(csv.DictReader(stream, delimiter="\t"))
        writer.writerow((label, times["wall_seconds"], times["user_cpu_seconds"],
                         times["system_cpu_seconds"], times["max_rss_kb"],
                         stats["raw_fragments"], stats["trimmed_fragments"],
                         stats["mapped_fragments_unique"]))
PY
