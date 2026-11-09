#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${OUT:-$SCRIPT_DIR/evaluation_runs/z3_10m_inputs}"
TARGET="${TARGET:-10000000}"
RANGE_BYTES="${RANGE_BYTES:-65000000}"
RANGE_PARTS="${RANGE_PARTS:-10}"
PARALLEL="${PARALLEL:-4}"
SE_SOURCE="${SE_SOURCE:-/home/people/micwe/Biodiversity_Extinction/Japanese_sea_lion/Raw/SRR25607782.fastq.gz}"
ENA_R1="${ENA_R1:-/home/people/micwe/Biodiversity_Extinction/Japanese_sea_lion/Raw/SRR25607804_1.fastq.gz}"
ENA_R2="${ENA_R2:-/home/people/micwe/Biodiversity_Extinction/Japanese_sea_lion/Raw/SRR25607804_2.fastq.gz}"
SRA_URL="${SRA_URL:-https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/093/SRR25607793/SRR25607793}"

for cmd in curl seqkit gzip sha256sum awk stat; do
  command -v "$cmd" >/dev/null || { echo "Missing command: $cmd" >&2; exit 1; }
done
[[ "$TARGET" =~ ^[1-9][0-9]*$ && "$RANGE_BYTES" =~ ^[1-9][0-9]*$ && "$RANGE_PARTS" =~ ^[1-9][0-9]*$ && "$PARALLEL" =~ ^[1-9][0-9]*$ ]] || {
  echo "TARGET, RANGE_BYTES, RANGE_PARTS, and PARALLEL must be positive integers" >&2
  exit 1
}

mkdir -p "$OUT/parts"

count_reads() {
  seqkit stats -T "$1" | awk -F '\t' 'NR == 2 {gsub(/,/, "", $4); print $4}'
}

validate_fastq() {
  local file="$1" expected="$2" actual
  gzip -t "$file"
  actual=$(count_reads "$file")
  [[ "$actual" == "$expected" ]] || {
    echo "Wrong read count in $file: expected $expected, found $actual" >&2
    return 1
  }
}

fetch_part() {
  local mate="$1" index="$2" start end expected part tmp tail status have
  start=$((index * RANGE_BYTES))
  end=$((start + RANGE_BYTES - 1))
  expected="$RANGE_BYTES"
  part="$OUT/parts/SRR25607793_${mate}.$(printf '%03d' "$index")"
  tmp="${part}.partial"
  if [[ -f "$part" && "$(stat -c %s "$part")" -eq "$expected" ]]; then
    return
  fi
  have=0
  [[ ! -f "$tmp" ]] || have=$(stat -c %s "$tmp")
  if (( have > expected )); then
    echo "Oversized partial ENA range: $tmp" >&2
    return 1
  fi
  if (( have == expected )); then
    mv "$tmp" "$part"
    return
  fi
  tail="${tmp}.tail"
  status=$(curl --silent --show-error --fail --location --retry 5 --retry-delay 5 \
    --connect-timeout 30 --max-time 1800 --range "$((start + have))-$end" \
    --output "$tail" --write-out '%{http_code}' "${SRA_URL}_${mate}.fastq.gz") || return 1
  if [[ "$status" != 206 || "$(stat -c %s "$tail")" -ne "$((expected - have))" ]]; then
    echo "Incomplete ENA range $mate:$((start + have))-$end (HTTP $status)" >&2
    return 1
  fi
  cat "$tail" >> "$tmp"
  rm -f "$tail"
  mv "$tmp" "$part"
  echo "Downloaded SRA mate $mate range $start-$end"
}

fetch_mate() {
  local mate="$1" index pid
  local pids=()
  for ((index=0; index<RANGE_PARTS; index++)); do
    fetch_part "$mate" "$index" &
    pids+=("$!")
    if ((${#pids[@]} >= PARALLEL)); then
      for pid in "${pids[@]}"; do wait "$pid" || return 1; done
      pids=()
    fi
  done
  for pid in "${pids[@]}"; do wait "$pid" || return 1; done
}

make_sra_subset() {
  local mate="$1" output partial tmp
  output="$OUT/SRR25607793_10M_${mate}.fastq.gz"
  if [[ -f "$output" ]] && validate_fastq "$output" "$TARGET"; then
    echo "Validated existing $output"
    return
  fi
  fetch_mate "$mate"
  partial="$OUT/parts/SRR25607793_${mate}.prefix.fastq.gz"
  cat "$OUT/parts/SRR25607793_${mate}."[0-9][0-9][0-9] > "$partial"
  tmp="${output%.gz}.partial.gz"
  seqkit head -n "$TARGET" -o "$tmp" "$partial"
  validate_fastq "$tmp" "$TARGET"
  mv "$tmp" "$output"
  echo "Prepared $output"
}

se_output="$OUT/SRR25607782_10M.fastq.gz"
if [[ ! -f "$se_output" ]] || ! validate_fastq "$se_output" "$TARGET"; then
  seqkit head -n "$TARGET" -o "${se_output}.partial.gz" "$SE_SOURCE"
  validate_fastq "${se_output}.partial.gz" "$TARGET"
  mv "${se_output}.partial.gz" "$se_output"
fi
echo "Prepared $se_output"

make_sra_subset 1
make_sra_subset 2

mate1_ids=$(seqkit seq -n -i "$OUT/SRR25607793_10M_1.fastq.gz" | sha256sum | awk '{print $1}')
mate2_ids=$(seqkit seq -n -i "$OUT/SRR25607793_10M_2.fastq.gz" | sha256sum | awk '{print $1}')
[[ "$mate1_ids" == "$mate2_ids" ]] || {
  echo "SRA paired FASTQ IDs do not match in order" >&2
  exit 1
}

ena_count=$(count_reads "$ENA_R1")
validate_fastq "$ENA_R1" "$ena_count"
validate_fastq "$ENA_R2" "$ena_count"

printf '%s\n' "$se_output" "$OUT/SRR25607793_10M_1.fastq.gz" \
  "$OUT/SRR25607793_10M_2.fastq.gz" "$ENA_R1" "$ENA_R2" > "$OUT/manifest.txt"

{
  printf 'source\tselection\treads_or_pairs\tsha256\tsource_url_or_path\n'
  printf 'SRR25607782\tfirst SE reads\t%s\t%s\t%s\n' \
    "$TARGET" "$(sha256sum "$se_output" | awk '{print $1}')" "$SE_SOURCE"
  for mate in 1 2; do
    file="$OUT/SRR25607793_10M_${mate}.fastq.gz"
    printf 'SRR25607793_%s\tfirst PE mate reads from ENA FASTQ\t%s\t%s\t%s\n' \
      "$mate" "$TARGET" "$(sha256sum "$file" | awk '{print $1}')" "${SRA_URL}_${mate}.fastq.gz"
  done
  for file in "$ENA_R1" "$ENA_R2"; do
    printf 'SRR25607804_%s\tall PE mate reads\t%s\t%s\t%s\n' \
      "$(basename "$file" | sed -E 's/.*_([12])\.fastq\.gz/\1/')" "$ena_count" "$(sha256sum "$file" | awk '{print $1}')" "$file"
  done
} > "$OUT/provenance.tsv"

echo "Input manifest: $OUT/manifest.txt"
