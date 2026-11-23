#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$SCRIPT_DIR/evaluation_runs"
OUT="${OUT:-$BASE/z3_10m_inputs_sra}"
FASTQ_DUMP="${FASTQ_DUMP:-$BASE/tools/sratoolkit/sratoolkit.3.4.1-ubuntu64/bin/fastq-dump}"
SEQKIT="${SEQKIT:-seqkit}"
ACCESSION=SRR25607793
TARGET=10000000
SE="$BASE/z3_10m_inputs/SRR25607782_10M.fastq.gz"
PE2_R1=/home/people/micwe/Biodiversity_Extinction/Japanese_sea_lion/Raw/SRR25607804_1.fastq.gz
PE2_R2=/home/people/micwe/Biodiversity_Extinction/Japanese_sea_lion/Raw/SRR25607804_2.fastq.gz
R1="$OUT/bulk/${ACCESSION}_1.fastq.gz"
R2="$OUT/bulk/${ACCESSION}_2.fastq.gz"

for file in "$FASTQ_DUMP" "$SE" "$PE2_R1" "$PE2_R2"; do
  [[ -s "$file" ]] || { echo "ERROR: missing input: $file" >&2; exit 1; }
done
[[ ! -e "$R1" && ! -e "$R2" && ! -e "$OUT/manifest.txt" ]] || {
  echo "ERROR: SRA output already exists; use a fresh OUT" >&2; exit 1;
}
mkdir -p "$OUT/bulk"

"$FASTQ_DUMP" --split-files --skip-technical --gzip \
  -N 1 -X "$TARGET" --outdir "$OUT/bulk" "$ACCESSION"

for file in "$R1" "$R2"; do
  gzip -t "$file"
  count=$("$SEQKIT" stats -T "$file" | awk -F '\t' 'NR == 2 { gsub(/,/, "", $4); print $4 }')
  [[ "$count" == "$TARGET" ]] || {
    echo "ERROR: expected $TARGET reads in $file, found $count" >&2; exit 1;
  }
done

# Compare every defline, not just the first pair: this is the SRA-header test.
paste <(gzip -dc "$R1" | awk 'NR % 4 == 1') \
      <(gzip -dc "$R2" | awk 'NR % 4 == 1') \
  | awk -F '\t' '$1 != $2 { print "ERROR: mate headers differ at pair " NR > "/dev/stderr"; exit 1 }'

printf '%s\n' "$SE" "$R1" "$R2" "$PE2_R1" "$PE2_R2" > "$OUT/manifest.txt"
{
  printf 'field\tvalue\n'
  printf 'source_accession\t%s\n' "$ACCESSION"
  printf 'source_archive\tNCBI SRA (fastq-dump direct accession access)\n'
  printf 'fastq_dump_version\t%s\n' "$("$FASTQ_DUMP" --version | awk 'NF { print; exit }')"
  printf 'fastq_dump_options\t--split-files --skip-technical --gzip -N 1 -X %s\n' "$TARGET"
  printf 'pairs\t%s\n' "$TARGET"
  printf 'mate_headers_identical\tPASS (all pairs)\n'
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$OUT/manifest.txt" | cut -d ' ' -f 1)"
  for file in "$R1" "$R2"; do
    printf '%s_sha256\t%s\n' "$(basename "$file")" "$(sha256sum "$file" | cut -d ' ' -f 1)"
  done
} > "$OUT/provenance.tsv"
echo "Validated $TARGET SRA pairs with identical mate headers: $OUT/manifest.txt"
