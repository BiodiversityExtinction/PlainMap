#!/usr/bin/env bash
set -euo pipefail

MODE="${1:?Usage: run_mapping_mode.sh modern|ancient-fast|ancient-all}"
case "$MODE" in
  modern|ancient-fast|ancient-all) ;;
  *) echo "ERROR: unsupported mode: $MODE" >&2; exit 2 ;;
esac

: "${PLAINMAP:?Set PLAINMAP to plainmap.sh}"
: "${DATA_DIR:?Set DATA_DIR to the evaluation FASTQ directory}"
: "${REF:?Set REF to GCF_009762305.2_mZalCal1.pri.v2_genomic.fna}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EVAL_ROOT="${EVAL_ROOT:-$SCRIPT_DIR/evaluation_runs/mapping_modes}"
FASTP="${FASTP:-fastp}"
BWA="${BWA:-bwa}"
SAMTOOLS="${SAMTOOLS:-samtools}"
PYTHON="${PYTHON:-python3}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
MODE_TAG="${MODE//-/_}"
PREFIX="Z3_${MODE_TAG}_pilot1M"
OUTDIR="$EVAL_ROOT/results/$MODE"
MANIFEST="$EVAL_ROOT/Z3_manifest_test.${MODE_TAG}.resolved.txt"

mkdir -p "$OUTDIR" "$EVAL_ROOT/logs" "$EVAL_ROOT/validation"
DATA_DIR="$DATA_DIR" "$SCRIPT_DIR/resolve_manifest.sh" "$MANIFEST" >/dev/null

{
  printf 'plainmap_version\t%s\n' "$("$PLAINMAP" --help | awk 'NR==1{print $2}')"
  printf 'plainmap_sha256\t%s\n' "$(sha256sum "$PLAINMAP" | awk '{print $1}')"
  printf 'mode\t%s\n' "$MODE"
  printf 'pilot_fragments_per_unit\t1000000\n'
  printf 'threads\t%s\n' "$THREADS"
  printf 'manifest\t%s\n' "$MANIFEST"
  printf 'reference\t%s\n' "$REF"
} > "$EVAL_ROOT/validation/${PREFIX}.provenance.tsv"

/usr/bin/time \
  -f $'wall_seconds\t%e\nmax_rss_kb\t%M\nexit_status\t%x' \
  -o "$EVAL_ROOT/validation/${PREFIX}.time.tsv" \
  "$PLAINMAP" \
    --manifest "$MANIFEST" \
    --prefix "$PREFIX" \
    --ref "$REF" \
    --outdir "$OUTDIR" \
    --threads "$THREADS" \
    --library-type "$MODE" \
    --pilot-fragments 1000000 \
    --keep-intermediate \
    --fastp "$FASTP" \
    --bwa "$BWA" \
    --samtools "$SAMTOOLS" \
    --python "$PYTHON" \
    > "$EVAL_ROOT/logs/${PREFIX}.out" 2>&1

BAM="$OUTDIR/${PREFIX}.final.RG.bam"
"$SAMTOOLS" quickcheck "$BAM"
[[ "$("$SAMTOOLS" view -c -f 4 "$BAM")" -eq 0 ]]
[[ "$("$SAMTOOLS" view -c -f 1024 "$BAM")" -eq 0 ]]

printf 'Completed and validated %s\n' "$MODE"
