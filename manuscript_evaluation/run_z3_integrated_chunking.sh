#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HERE/evaluation_runs"
PLAINMAP="${PLAINMAP:-$BASE/plainmap_integrated_v030.sh}"
MANIFEST="$BASE/z3_10m_inputs_sra/manifest.txt"
REF="${REF:-/home/people/micwe/Biodiversity_Extinction/Japanese_sea_lion/Reference_genomes/GCF_009762305.2_mZalCal1.pri.v2_genomic.fna}"
DEDUP="$BASE/tools/dedup-release"
OUT="${EVAL_ROOT:-$BASE/z3_integrated_chunking_v030}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
PREFIX=Z3_large_modern

[[ ! -e "$OUT" ]] || { echo "Output already exists: $OUT" >&2; exit 1; }
for file in "$PLAINMAP" "$MANIFEST" "$REF" "$REF.bwt" "$DEDUP"; do
  [[ -s "$file" ]] || { echo "Missing input: $file" >&2; exit 1; }
done
[[ -s "$BASE/plainmap_integrated_v030.sha256" ]] || {
  echo "Missing frozen PlainMap checksum" >&2; exit 1;
}
[[ "$(sha256sum "$PLAINMAP" | cut -d ' ' -f 1)" == \
  "$(cat "$BASE/plainmap_integrated_v030.sha256")" ]] || {
  echo "Frozen PlainMap script changed" >&2; exit 1;
}
[[ "$(sha256sum "$MANIFEST" | cut -d ' ' -f 1)" == \
  c80a82b78bebe1cf55dff94ce84c4f54f20469bc4fbe282b467a16585ed3cfab ]] || {
  echo "Validated SRA-derived manifest has changed" >&2; exit 1;
}
grep -F $'mate_headers_identical\tPASS (all pairs)' \
  "$BASE/z3_10m_inputs_sra/provenance.tsv" >/dev/null
mkdir -p "$OUT"
{
  printf 'field\tvalue\n'
  printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-none}"
  printf 'node\t%s\n' "$(hostname)"
  printf 'plainmap_version\t%s\n' '0.3.0 (integrated direct-chunk update)'
  printf 'plainmap_sha256\t%s\n' "$(sha256sum "$PLAINMAP" | cut -d ' ' -f 1)"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$MANIFEST" | cut -d ' ' -f 1)"
  printf 'reference\t%s\n' "$REF"
  printf 'threads\t%s\n' "$THREADS"
  printf 'pilot_fragments_per_unit\t10000000\n'
  printf 'minlength\t30\nmapq\t20\ndedup_method\thybrid\n'
  printf 'fastp_version\t%s\n' "$(fastp --version 2>&1)"
  printf 'bwa_version\t%s\n' "$(bwa 2>&1 | awk '/Version:/ {print $2; exit}' || true)"
  printf 'samtools_version\t%s\n' "$(samtools --version | sed -n '1p')"
  printf 'dedup_version\t%s\n' "$($DEDUP --version 2>&1 | head -1 || true)"
} > "$OUT/provenance.tsv"

common_args=( --manifest "$MANIFEST" --prefix "$PREFIX" --ref "$REF"
  --threads "$THREADS" --library-type modern --pilot-fragments 10000000
  --dedup-method hybrid --dedup "$DEDUP" --fastp fastp --bwa bwa
  --samtools samtools --python python3 )

run_case() {
  local label="$1" chunk_size="$2" dir="$OUT/$1"
  mkdir -p "$dir"
  /usr/bin/time -f $'wall_seconds\t%e\nuser_cpu_seconds\t%U\nsystem_cpu_seconds\t%S\nmax_rss_kb\t%M\nexit_status\t%x' \
    -o "$dir/time.tsv" \
    "$PLAINMAP" "${common_args[@]}" --outdir "$dir" \
      --max-reads-per-chunk "$chunk_size" > "$dir/driver.log" 2>&1
  local bam="$dir/$PREFIX.final.RG.bam"
  samtools quickcheck "$bam"
  [[ "$(samtools view -c -f 4 "$bam")" -eq 0 ]]
  [[ "$(samtools view -c -f 1024 "$bam")" -eq 0 ]]
  [[ -s "$dir/$PREFIX.plainmap.stats.tsv" ]]
  [[ -s "$dir/$PREFIX.plainmap.adapters.tsv" ]]
  echo "Validated $label"
}

# Sequential cases share the same Slurm allocation, avoiding cross-node timing comparisons.
for size in 1000000 5000000 10000000; do
  run_case "chunk_${size}" "$size"
done

BASELINE="$OUT/chunk_1000000"
RESTART="$OUT/restart"
mkdir -p "$RESTART"
restart_args=( "${common_args[@]}" --outdir "$RESTART" --max-reads-per-chunk 1000000 --keep-intermediate )
BAMS="$RESTART/$PREFIX/work/bams"
PE_IDS="$RESTART/$PREFIX/work/.checkpoints/$PREFIX.mapchunk_ids.pe.txt"
SE_IDS="$RESTART/$PREFIX/work/.checkpoints/$PREFIX.mapchunk_ids.se.txt"

count_valid_chunks() {
  local bam count=0
  shopt -s nullglob
  for bam in "$BAMS"/*.bam; do
    samtools quickcheck "$bam" >/dev/null 2>&1 && count=$((count + 1))
  done
  shopt -u nullglob
  printf '%s\n' "$count"
}

echo "Starting restart case; interrupting after at least three valid mapping chunks"
setsid "$PLAINMAP" "${restart_args[@]}" > "$RESTART/first_attempt.log" 2>&1 &
pid=$!
interrupted=0 expected=0 completed=0
for ((poll=0; poll<3600; poll++)); do
  kill -0 "$pid" 2>/dev/null || break
  if [[ -f "$PE_IDS" && -f "$SE_IDS" ]]; then
    expected=$(awk 'NF {n++} END {print n+0}' "$PE_IDS" "$SE_IDS")
  fi
  completed=$(count_valid_chunks)
  if (( expected > 3 && completed >= 3 && completed < expected )); then
    kill -TERM -- "-$pid"
    interrupted=1
    break
  fi
  sleep 2
done
if (( interrupted == 0 )); then
  kill -TERM -- "-$pid" 2>/dev/null || true
  wait "$pid" 2>/dev/null || true
  echo "Could not interrupt a partially completed mapping run" >&2
  exit 1
fi
for _ in {1..10}; do
  kill -0 -- "-$pid" 2>/dev/null || break
  sleep 1
done
kill -KILL -- "-$pid" 2>/dev/null || true
set +e
wait "$pid"
first_status=$?
set -e
completed=$(count_valid_chunks)
(( first_status != 0 && completed >= 3 && completed < expected )) || {
  echo "Interrupted run did not retain valid partial mapping" >&2; exit 1;
}
printf 'bam\tsha256\n' > "$RESTART/completed_before.tsv"
shopt -s nullglob
for bam in "$BAMS"/*.bam; do
  if samtools quickcheck "$bam" >/dev/null 2>&1; then
    printf '%s\t%s\n' "$(basename "$bam")" "$(sha256sum "$bam" | cut -d ' ' -f 1)" \
      >> "$RESTART/completed_before.tsv"
  fi
done
shopt -u nullglob

/usr/bin/time -f $'wall_seconds\t%e\nuser_cpu_seconds\t%U\nsystem_cpu_seconds\t%S\nmax_rss_kb\t%M\nexit_status\t%x' \
  -o "$RESTART/resume_time.tsv" \
  "$PLAINMAP" "${restart_args[@]}" > "$RESTART/resume.log" 2>&1
reused=$(grep -c 'already mapped (valid BAM exists) -> skipping' "$RESTART/resume.log" || true)
(( reused >= completed )) || { echo "Previously mapped chunks were not reused" >&2; exit 1; }
while IFS=$'\t' read -r bam expected_hash; do
  [[ "$bam" == bam ]] && continue
  [[ "$(sha256sum "$BAMS/$bam" | cut -d ' ' -f 1)" == "$expected_hash" ]] || {
    echo "Previously mapped chunk changed: $bam" >&2; exit 1;
  }
done < "$RESTART/completed_before.tsv"
baseline_hash=$(samtools view "$BASELINE/$PREFIX.final.RG.bam" | sha256sum | cut -d ' ' -f 1)
resumed_bam="$RESTART/$PREFIX.final.RG.bam"
samtools quickcheck "$resumed_bam"
[[ "$(samtools view -c -f 4 "$resumed_bam")" -eq 0 ]]
[[ "$(samtools view -c -f 1024 "$resumed_bam")" -eq 0 ]]
resumed_hash=$(samtools view "$resumed_bam" | sha256sum | cut -d ' ' -f 1)
[[ "$baseline_hash" == "$resumed_hash" ]] || { echo "Restart BAM records differ" >&2; exit 1; }
cmp -s "$BASELINE/$PREFIX.plainmap.stats.tsv" "$RESTART/$PREFIX.plainmap.stats.tsv" || {
  echo "Restart statistics differ" >&2; exit 1;
}
{
  printf 'check\tvalue\n'
  printf 'mapping_chunks_expected\t%s\n' "$expected"
  printf 'valid_chunks_at_interruption\t%s\n' "$completed"
  printf 'reused_chunks_on_resume\t%s\n' "$reused"
  printf 'first_attempt_exit_status\t%s\n' "$first_status"
  printf 'baseline_bam_records_sha256\t%s\n' "$baseline_hash"
  printf 'resumed_bam_records_sha256\t%s\n' "$resumed_hash"
  printf 'completed_chunk_checksums_preserved\tPASS\n'
  printf 'final_bam_records_identical\tPASS\n'
  printf 'statistics_identical\tPASS\n'
} > "$RESTART/restart_checks.tsv"

python3 - "$OUT" "$PREFIX" <<'PY'
import csv
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
prefix = sys.argv[2]
with (root / "results.tsv").open("w") as output:
    writer = csv.writer(output, delimiter="\t")
    writer.writerow(("case", "chunk_size", "raw_fragments", "trimmed_fragments",
                     "mapped_fragments_unique", "wall_seconds", "user_cpu_seconds",
                     "system_cpu_seconds", "max_rss_kb"))
    for size in (1000000, 5000000, 10000000):
        case = root / f"chunk_{size}"
        with (case / f"{prefix}.plainmap.stats.tsv").open() as stream:
            stats = next(csv.DictReader(stream, delimiter="\t"))
        times = dict(line.rstrip().split("\t", 1) for line in (case / "time.tsv").open())
        writer.writerow((case.name, size, stats["raw_fragments"],
                         stats["trimmed_fragments"], stats["mapped_fragments_unique"],
                         times["wall_seconds"], times["user_cpu_seconds"],
                         times["system_cpu_seconds"], times["max_rss_kb"]))
PY
cat "$OUT/results.tsv"
echo "All integrated v0.3.0 chunking and restart evaluations passed: $OUT"
