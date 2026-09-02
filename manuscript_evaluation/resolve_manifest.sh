#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: resolve_manifest.sh OUTPUT_MANIFEST" >&2
  exit 2
fi
: "${DATA_DIR:?Set DATA_DIR to the directory containing the evaluation FASTQs}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT="$1"
mkdir -p "$(dirname "$OUTPUT")"
: > "$OUTPUT"

while IFS= read -r name; do
  [[ -z "$name" || "$name" == \#* ]] && continue
  file="$DATA_DIR/$name"
  [[ -r "$file" ]] || {
    echo "ERROR: evaluation FASTQ is not readable: $file" >&2
    exit 1
  }
  printf '%s\n' "$file" >> "$OUTPUT"
done < "$SCRIPT_DIR/Z3_manifest_test.txt"

echo "$OUTPUT"
