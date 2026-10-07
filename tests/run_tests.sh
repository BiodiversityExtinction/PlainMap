#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PLAINMAP="$REPO_ROOT/plainmap.sh"
TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/plainmap-tests.XXXXXX")"
trap 'rm -rf "$TMP_ROOT"' EXIT

DEDUP_CMD="${DEDUP:-dedup}"

for cmd in fastp bwa samtools "$DEDUP_CMD" python3 gzip sha256sum; do
  command -v "$cmd" >/dev/null 2>&1 || {
    echo "ERROR: required test command not found: $cmd" >&2
    exit 1
  }
done

bash -n "$PLAINMAP"
gzip -t "$REPO_ROOT/example/reads_R1.fastq.gz" "$REPO_ROOT/example/reads_R2.fastq.gz"

cp "$REPO_ROOT/example/reference.fa" "$TMP_ROOT/reference.fa"

# Add a second, unrelated contig and a matching SE read so the mapDamage test
# exercises the real parallel split/merge path rather than its serial fallback.
chr2_seq="TTGACCATGATCCGTTACGAGGCTAACCTGATGCGTTCAGATCCAGTACGTTGACCTAGGATCCGTAGACCTTGCATGGTACCTGAGTTCGATACCGTTCAGGATGCTAACGTAGCTTGACGGTACAGTCCGATGATCCTGAC"
printf '>toy_chr2\n%s\n' "$chr2_seq" >> "$TMP_ROOT/reference.fa"

# SRA paired files commonly carry the same read name in both mates. Retain the
# toy sequences while replacing both mate headers with identical SRA-style names.
for mate in R1 R2; do
  gzip -dc "$REPO_ROOT/example/reads_${mate}.fastq.gz" \
    | awk 'NR % 4 == 1 {n++; print "@ERRTEST." n " " n " length=52"; next} {print}' \
    | gzip -c > "$TMP_ROOT/sra_${mate}.fastq.gz"
done

gzip -dc "$REPO_ROOT/example/reads_R1.fastq.gz" \
  | awk 'NR <= 4 {line[NR]=$0}
         END {
           for (copy=1; copy<=2; copy++) {
             print "@SETEST." copy " length=52"
             print line[2]
             print "+"
             print line[4]
           }
         }' \
  | gzip -c > "$TMP_ROOT/single_end.fastq.gz"

quality52=$(printf '%*s' 52 '' | tr ' ' I)
printf '@SETEST.2 length=52\n%s\n+\n%s\n' "${chr2_seq:0:52}" "$quality52" \
  | gzip -c > "$TMP_ROOT/single_end_chr2.fastq.gz"

cat > "$TMP_ROOT/manifest.txt" <<'EOF'
sra_R1.fastq.gz
sra_R2.fastq.gz
single_end.fastq.gz
single_end_chr2.fastq.gz
EOF

cd "$TMP_ROOT"

run_plainmap() {
  local prefix="$1"
  local mode="$2"
  local outdir="$3"
  shift 3
  local driver_log="$TMP_ROOT/${prefix}.${mode}.$(basename "$outdir").driver.log"

  if ! bash "$PLAINMAP" \
      --manifest "$TMP_ROOT/manifest.txt" \
      --prefix "$prefix" \
      --ref "$TMP_ROOT/reference.fa" \
      --outdir "$outdir" \
      --threads 1 \
      --library-type "$mode" \
      "$@" > "$driver_log" 2>&1; then
    cat "$driver_log" >&2
    return 1
  fi
}

stats_value() {
  local stats="$1"
  local field="$2"
  awk -F '\t' -v field="$field" '
    NR == 1 {for (i=1; i<=NF; i++) column[$i]=i; next}
    NR == 2 {print $column[field]}
  ' "$stats"
}

alignment_checksum() {
  samtools view "$1" | LC_ALL=C sort | sha256sum | awk '{print $1}'
}

validate_outputs() {
  local prefix="$1"
  local mode="$2"
  local outdir="$3"
  local dedup_method="${4:-hybrid}"
  local bam="$outdir/${prefix}.final.RG.bam"
  local stats="$outdir/${prefix}.plainmap.stats.tsv"

  samtools quickcheck "$bam"
  test -s "${bam}.bai"
  test -s "$outdir/${prefix}.plainmap.coverage.tsv"
  test -s "$outdir/${prefix}.plainmap.run-config.tsv"
  samtools idxstats "$bam" >/dev/null

  [[ "$(samtools view -c -f 4 "$bam")" -eq 0 ]]
  [[ "$(samtools view -c -f 1024 "$bam")" -eq 0 ]]
  samtools view -H "$bam" | grep -F $'@HD\tVN:' | grep -F $'SO:coordinate' >/dev/null
  samtools view -H "$bam" | grep -F $'@RG\tID:'"$prefix"$'\t' >/dev/null

  local mapped_reads
  mapped_reads="$(samtools view -c -F 2820 "$bam")"
  [[ "$mapped_reads" -gt 0 ]]
  [[ "$mapped_reads" -eq "$(stats_value "$stats" mapped_reads_unique)" ]]
  [[ "$(stats_value "$stats" library_type)" == "$mode" ]]
  [[ "$(stats_value "$stats" seq_mode)" == "MIX" ]]

  if samtools view "$bam" | awk -v tag="RG:Z:${prefix}" '
      {found=0; for (i=12; i<=NF; i++) if ($i == tag) found=1; if (!found) exit 1}
    '; then
    :
  else
    echo "ERROR: one or more alignments lack the expected read-group tag" >&2
    return 1
  fi

  python3 - "$stats" "$dedup_method" <<'PY'
import csv
import math
import sys

with open(sys.argv[1], newline="") as handle:
    row = next(csv.DictReader(handle, delimiter="\t"))
expected_method = sys.argv[2]

raw = int(row["raw_fragments"])
mapped_all = int(row["mapped_fragments_all"])
mapped_unique = int(row["mapped_fragments_unique"])
duplicates = int(row["duplicate_fragments"])

assert mapped_unique == mapped_all - duplicates
expected_endog = 0.0 if raw == 0 else mapped_unique / raw
expected_duprate = 0.0 if mapped_all == 0 else duplicates / mapped_all
assert math.isclose(float(row["endog"]), expected_endog, abs_tol=5.1e-7)
assert math.isclose(float(row["duprate"]), expected_duprate, abs_tol=5.1e-7)
assert float(row["avg_depth"]) >= 0
assert 0 <= float(row["pct_covered"]) <= 100

assert row["pe_dedup_method"] == "samtools-markdup"
assert row["dedup_method"] == expected_method
expected_se_method = "DeDup" if expected_method == "hybrid" else "samtools-markdup"
assert row["se_like_dedup_method"] == expected_se_method

for prefix in ("pe", "se_like"):
    all_fragments = int(row[f"{prefix}_fragments_all"])
    unique_fragments = int(row[f"{prefix}_fragments_unique"])
    duplicate_fragments = int(row[f"{prefix}_duplicate_fragments"])
    duplicate_rate = float(row[f"{prefix}_duprate"])
    assert unique_fragments == all_fragments - duplicate_fragments
    expected_rate = 0.0 if all_fragments == 0 else duplicate_fragments / all_fragments
    assert math.isclose(duplicate_rate, expected_rate, abs_tol=5.1e-7)

assert int(row["mapped_fragments_all"]) == (
    int(row["pe_fragments_all"]) + int(row["se_like_fragments_all"])
)
assert int(row["duplicate_fragments"]) == (
    int(row["pe_duplicate_fragments"]) + int(row["se_like_duplicate_fragments"])
)
assert int(row["se_like_duplicate_fragments"]) >= 1
if row["library_type"] == "ancient-fast":
    assert int(row["pe_fragments_all"]) == 0
else:
    assert int(row["pe_fragments_all"]) > 0
PY
}

echo "Testing all mapping modes with mixed SRA-style PE and SE input"
for mode in modern ancient-fast ancient-all; do
  outdir="$TMP_ROOT/modes/$mode"
  run_plainmap "$mode" "$mode" "$outdir"
  validate_outputs "$mode" "$mode" "$outdir"
done

echo "Testing v0.1-compatible samtools-only deduplication in all mapping modes"
for mode in modern ancient-fast ancient-all; do
  outdir="$TMP_ROOT/samtools-modes/$mode"
  run_plainmap "samtools_${mode}" "$mode" "$outdir" \
    --dedup-method samtools --dedup /definitely/not/installed/dedup
  validate_outputs "samtools_${mode}" "$mode" "$outdir" samtools
  grep -F $'dedup_method\tsamtools' \
    "$outdir/samtools_${mode}.plainmap.run-config.tsv" >/dev/null
  if find "$outdir/samtools_${mode}/dedup_reports" -type f -print -quit | grep -q .; then
    echo "ERROR: samtools-only mode unexpectedly produced DeDup reports" >&2
    exit 1
  fi
done

if git -C "$REPO_ROOT" rev-parse --verify --quiet refs/tags/v0.1.0 >/dev/null; then
  echo "Testing samtools-only output against the immutable v0.1.0 release"
  old_plainmap="$TMP_ROOT/plainmap-v0.1.0.sh"
  git -C "$REPO_ROOT" show v0.1.0:plainmap.sh > "$old_plainmap"
  chmod +x "$old_plainmap"
  old_out="$TMP_ROOT/v0.1-compatibility"
  if ! bash "$old_plainmap" \
      --manifest "$TMP_ROOT/manifest.txt" \
      --prefix samtools_modern \
      --ref "$TMP_ROOT/reference.fa" \
      --outdir "$old_out" \
      --threads 1 \
      --library-type modern > "$TMP_ROOT/v0.1-compatibility.driver.log" 2>&1; then
    cat "$TMP_ROOT/v0.1-compatibility.driver.log" >&2
    exit 1
  fi

  current_out="$TMP_ROOT/samtools-modes/modern"
  [[ "$(alignment_checksum "$old_out/samtools_modern.final.RG.bam")" == \
     "$(alignment_checksum "$current_out/samtools_modern.final.RG.bam")" ]]

  python3 - \
    "$old_out/samtools_modern.plainmap.stats.tsv" \
    "$current_out/samtools_modern.plainmap.stats.tsv" <<'PY'
import csv
import sys

def read(path):
    with open(path, newline="") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))

old = read(sys.argv[1])
current = read(sys.argv[2])
shared = old.keys() & current.keys()
differences = {key: (old[key], current[key]) for key in shared if old[key] != current[key]}
assert not differences, differences
assert len(shared) == 23, shared
PY
else
  echo "Skipping v0.1.0 compatibility test (tag not available in this checkout)"
fi

echo "Testing strict resume configuration and reset"
config_out="$TMP_ROOT/config-output"
run_plainmap config modern "$config_out" --keep-intermediate
run_plainmap config modern "$config_out" --keep-intermediate
grep -F "Run configuration matches" "$config_out/config_plainmap.log" >/dev/null

if bash "$PLAINMAP" \
    --manifest "$TMP_ROOT/manifest.txt" \
    --prefix config \
    --ref "$TMP_ROOT/reference.fa" \
    --outdir "$config_out" \
    --threads 1 \
    --library-type modern \
    --mapq 25 \
    --keep-intermediate > "$TMP_ROOT/config-mismatch.log" 2>&1; then
  echo "ERROR: changed resume configuration unexpectedly succeeded" >&2
  exit 1
fi
grep -F "Resume configuration mismatch" "$TMP_ROOT/config-mismatch.log" >/dev/null
grep -F "mapq: previous='20'; current='25'" "$TMP_ROOT/config-mismatch.log" >/dev/null

if bash "$PLAINMAP" \
    --manifest "$TMP_ROOT/manifest.txt" \
    --prefix config \
    --ref "$TMP_ROOT/reference.fa" \
    --outdir "$config_out" \
    --threads 1 \
    --library-type modern \
    --dedup-method samtools \
    --keep-intermediate > "$TMP_ROOT/dedup-method-mismatch.log" 2>&1; then
  echo "ERROR: changed dedup method unexpectedly resumed" >&2
  exit 1
fi
grep -F "dedup_method: previous='hybrid'; current='samtools'" \
  "$TMP_ROOT/dedup-method-mismatch.log" >/dev/null

run_plainmap config modern "$config_out" --mapq 25 --keep-intermediate --reset
grep -F $'mapq\t25' "$config_out/config.plainmap.run-config.tsv" >/dev/null

echo "Testing forced chunking equivalence"
base_out="$TMP_ROOT/equivalence/base"
chunk_out="$TMP_ROOT/equivalence/chunked"
run_plainmap equivalence modern "$base_out" --keep-intermediate
run_plainmap equivalence modern "$chunk_out" --max-reads-per-chunk 1 --keep-intermediate
validate_outputs equivalence modern "$base_out"
validate_outputs equivalence modern "$chunk_out"

base_checksum="$(alignment_checksum "$base_out/equivalence.final.RG.bam")"
chunk_checksum="$(alignment_checksum "$chunk_out/equivalence.final.RG.bam")"
[[ "$base_checksum" == "$chunk_checksum" ]]

python3 - \
  "$base_out/equivalence.plainmap.stats.tsv" \
  "$chunk_out/equivalence.plainmap.stats.tsv" <<'PY'
import csv
import sys

def read(path):
    with open(path, newline="") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))

base = read(sys.argv[1])
chunked = read(sys.argv[2])
base.pop("max_reads_per_chunk")
chunked.pop("max_reads_per_chunk")
assert base == chunked, (base, chunked)
PY

echo "Testing recovery of a missing persistent DeDup report"
dedup_json="$chunk_out/equivalence/dedup_reports/equivalence.se_like.coordsort.dedup.json"
rm -f "$dedup_json"
run_plainmap equivalence modern "$chunk_out" --max-reads-per-chunk 1 --keep-intermediate
test -s "$dedup_json"
if grep -F "class-specific dedup (skipped; checkpoint present)" \
    "$chunk_out/equivalence_plainmap.log" >/dev/null; then
  echo "ERROR: dedup checkpoint was reused despite a missing DeDup report" >&2
  exit 1
fi

echo "Testing reuse of completed mapping chunks during downstream recovery"
checkpoint_dir="$chunk_out/equivalence/work/.checkpoints"
rm -f \
  "$checkpoint_dir/equivalence.merge.done" \
  "$checkpoint_dir/equivalence.dedup.done" \
  "$checkpoint_dir/equivalence.rg.done" \
  "$checkpoint_dir/equivalence.coverage.done" \
  "$chunk_out/equivalence.final.RG.bam" \
  "$chunk_out/equivalence.final.RG.bam.bai" \
  "$chunk_out/equivalence.plainmap.stats.tsv" \
  "$chunk_out/equivalence.plainmap.coverage.tsv"
rm -f "$chunk_out/equivalence/work/final/"*

run_plainmap equivalence modern "$chunk_out" --max-reads-per-chunk 1 --keep-intermediate
grep -F "already mapped (valid BAM exists) -> skipping" \
  "$chunk_out/equivalence_plainmap.log" >/dev/null
[[ "$chunk_checksum" == "$(alignment_checksum "$chunk_out/equivalence.final.RG.bam")" ]]
validate_outputs equivalence modern "$chunk_out"

if command -v mapDamage >/dev/null 2>&1; then
  echo "Testing optional mapDamage integration (serial and parallel results agree)"
  md_serial_out="$TMP_ROOT/mapdamage/serial"
  md_auto_out="$TMP_ROOT/mapdamage/auto"

  # Shared references are often read-only. An existing index must be reused,
  # not rebuilt by the parallel implementation.
  samtools faidx "$TMP_ROOT/reference.fa"
  chmod a-w "$TMP_ROOT/reference.fa" "$TMP_ROOT/reference.fa.fai"

  run_plainmap md modern "$md_serial_out" --run-mapdamage --mapdamage-serial --keep-intermediate
  run_plainmap md modern "$md_auto_out" --run-mapdamage --threads 2 --keep-intermediate
  grep -F "splitting 2 contig(s) into 2 parallel chunk(s)" "$md_auto_out/md_plainmap.log" >/dev/null

  for out in "$md_serial_out" "$md_auto_out"; do
    for f in misincorporation.txt dnacomp.txt lgdistribution.txt; do
      test -s "$out/md_mapdamage/$f"
    done
    if command -v Rscript >/dev/null 2>&1; then
      for f in Fragmisincorporation_plot.pdf Length_plot.pdf 5pCtoT_freq.txt 3pGtoA_freq.txt; do
        test -s "$out/md_mapdamage/$f"
      done
    fi
  done

  # Comments and row order may differ after merging, so compare table schemas
  # and complete row multisets rather than requiring byte-identical files.
  python3 - "$md_serial_out/md_mapdamage" "$md_auto_out/md_mapdamage" <<'PY'
import sys

def read_table(path):
    with open(path) as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip() and not line.startswith("#")]
    return lines[0], sorted(lines[1:])

serial_dir, parallel_dir = sys.argv[1:]
for name in ("misincorporation.txt", "dnacomp.txt", "lgdistribution.txt"):
    serial = read_table(f"{serial_dir}/{name}")
    parallel = read_table(f"{parallel_dir}/{name}")
    assert serial == parallel, name
PY

  # A stale top-level checkpoint must not hide missing final output. Complete
  # chunk outputs should be reused to reconstruct only the merged result.
  mv "$md_auto_out/md_mapdamage/dnacomp.txt" "$md_auto_out/md_mapdamage/dnacomp.txt.missing"
  run_plainmap md modern "$md_auto_out" --run-mapdamage --threads 2 --keep-intermediate
  grep -F "checkpoint is present but outputs are incomplete" "$md_auto_out/md_plainmap.log" >/dev/null
  test -s "$md_auto_out/md_mapdamage/dnacomp.txt"
else
  echo "Skipping mapDamage tests (mapDamage not found in PATH)"
fi

echo "Testing corrupt FASTQ validation failure"
cp "$TMP_ROOT/sra_R1.fastq.gz" "$TMP_ROOT/corrupt.fastq.gz"
truncate -s -8 "$TMP_ROOT/corrupt.fastq.gz"
printf 'corrupt.fastq.gz\n' > "$TMP_ROOT/corrupt_manifest.txt"
if bash "$PLAINMAP" \
    --manifest "$TMP_ROOT/corrupt_manifest.txt" \
    --prefix corrupt \
    --ref "$TMP_ROOT/reference.fa" \
    --outdir "$TMP_ROOT/corrupt-output" \
    --validate > "$TMP_ROOT/corrupt.log" 2>&1; then
  echo "ERROR: corrupt FASTQ unexpectedly passed validation" >&2
  exit 1
fi
grep -F "Corrupt/truncated gzip" "$TMP_ROOT/corrupt.log" >/dev/null

echo "All PlainMap tests passed"
