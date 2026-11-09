#!/usr/bin/env bash
set -euo pipefail

MODE="${1:?Usage: run_z3_28m_modes.sh modern|ancient-fast|ancient-all}"
case "$MODE" in
  modern|ancient-fast|ancient-all) ;;
  *) echo "ERROR: unsupported mode: $MODE" >&2; exit 2 ;;
esac

: "${PLAINMAP:?Set PLAINMAP to the frozen v0.3.0 script}"
: "${MANIFEST:?Set MANIFEST to the validated Z3 28M-fragment manifest}"
: "${REF:?Set REF to the indexed California sea lion reference}"
: "${DEDUP:?Set DEDUP to the pinned DeDup executable}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EVAL_ROOT="${EVAL_ROOT:-$SCRIPT_DIR/evaluation_runs/z3_v030_28m_modes}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
PILOT_FRAGMENTS="${PILOT_FRAGMENTS:-10000000}"
FASTP="${FASTP:-fastp}"
BWA="${BWA:-bwa}"
SAMTOOLS="${SAMTOOLS:-samtools}"
PYTHON="${PYTHON:-python3}"
MODE_TAG="${MODE//-/_}"
PREFIX="Z3_28m_${MODE_TAG}"
OUT="$EVAL_ROOT/$MODE"

[[ ! -e "$OUT" ]] || { echo "ERROR: output already exists: $OUT" >&2; exit 1; }
for file in "$PLAINMAP" "$MANIFEST" "$REF" "$REF.bwt" "$DEDUP"; do
  [[ -s "$file" ]] || { echo "ERROR: required file is missing or empty: $file" >&2; exit 1; }
done
mkdir -p "$OUT"

{
  printf 'plainmap\t%s\n' "$PLAINMAP"
  printf 'plainmap_sha256\t%s\n' "$(sha256sum "$PLAINMAP" | cut -d ' ' -f 1)"
  printf 'manifest\t%s\n' "$MANIFEST"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$MANIFEST" | cut -d ' ' -f 1)"
  printf 'reference\t%s\n' "$REF"
  printf 'mode\t%s\n' "$MODE"
  printf 'threads\t%s\n' "$THREADS"
  printf 'pilot_fragments_per_unit\t%s\n' "$PILOT_FRAGMENTS"
  printf 'minlength\t30\nmapq\t20\nmax_reads_per_chunk\t0\n'
  printf 'dedup_method\thybrid\n'
  printf 'ancient_bwa_aln\t-l 999 -n 0.01 -o 2 -t %s\n' "$THREADS"
  printf 'fastp_version\t%s\n' "$("$FASTP" --version 2>&1)"
  printf 'bwa_version\t%s\n' "$("$BWA" 2>&1 | awk '/Version:/ {print $2; exit}' || true)"
  printf 'samtools_version\t%s\n' "$("$SAMTOOLS" --version | sed -n '1p')"
  printf 'dedup_version\t%s\n' "$("$DEDUP" --version 2>&1 | sed -n '1p')"
} > "$OUT/provenance.tsv"

/usr/bin/time -f $'wall_seconds\t%e\nuser_cpu_seconds\t%U\nsystem_cpu_seconds\t%S\nmax_process_rss_kb\t%M\nexit_status\t%x' \
  -o "$OUT/time.tsv" \
  bash "$PLAINMAP" \
    --manifest "$MANIFEST" --prefix "$PREFIX" --ref "$REF" --outdir "$OUT" \
    --threads "$THREADS" --library-type "$MODE" \
    --pilot-fragments "$PILOT_FRAGMENTS" \
    --dedup-method hybrid --dedup "$DEDUP" \
    --fastp "$FASTP" --bwa "$BWA" --samtools "$SAMTOOLS" --python "$PYTHON" \
    > "$OUT/driver.log" 2>&1

BAM="$OUT/$PREFIX.final.RG.bam"
"$SAMTOOLS" quickcheck "$BAM"
[[ "$("$SAMTOOLS" view -c -f 4 "$BAM")" -eq 0 ]]
[[ "$("$SAMTOOLS" view -c -f 1024 "$BAM")" -eq 0 ]]
[[ -s "$OUT/$PREFIX.plainmap.stats.tsv" ]]
printf 'Completed and validated %s: %s\n' "$MODE" "$OUT"
