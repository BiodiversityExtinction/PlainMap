#!/usr/bin/env bash
set -euo pipefail

repo="${PLAINMAP_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
root="${EVAL_ROOT:-$repo/manuscript_evaluation/evaluation_runs/cavehyena_50m_adaptclean}"
prefix=Ccsp043_50m
outfile="$root/comparison.tsv"

{
  printf 'condition\traw_fragments\ttrimmed_fragments\tadaptclean_shortened\tmapped_fragments_all\tmapped_fragments_unique\tduprate\tmean_trimmed_length_bp\twall_seconds\tuser_cpu_seconds\tsystem_cpu_seconds\tmax_process_rss_kb\n'
  for condition in off on; do
    dir="$root/$condition"
    stats="$dir/$prefix.plainmap.stats.tsv"
    hist="$dir/trimmed_length_histogram.tsv"
    timefile="$dir/time.tsv"
    [[ -s "$stats" && -s "$hist" && -s "$timefile" ]] || {
      echo "ERROR: incomplete $condition run" >&2; exit 1;
    }
    values=$(awk -F '\t' 'NR==1 {for(i=1;i<=NF;i++) col[$i]=i; next} NR==2 {print $(col["raw_fragments"]), $(col["trimmed_fragments"]), $(col["mapped_fragments_all"]), $(col["mapped_fragments_unique"]), $(col["duprate"])}' OFS='\t' "$stats")
    mean_len=$(awk -F '\t' 'NR>1 {total += $1*$2; reads += $2} END {if(reads == 0) exit 1; printf "%.2f", total/reads}' "$hist")
    shortened=0
    if [[ "$condition" == on ]]; then
      report="$dir/$prefix.adaptclean.tsv"
      [[ -s "$report" ]] || { echo "ERROR: missing AdaptClean report" >&2; exit 1; }
      shortened=$(awk -F '\t' 'NR>1 {total += $6} END {print total+0}' "$report")
    fi
    timing=$(awk -F '\t' '$1=="wall_seconds" || $1=="user_cpu_seconds" || $1=="system_cpu_seconds" || $1=="max_rss_kb" {value[$1]=$2} END {print value["wall_seconds"], value["user_cpu_seconds"], value["system_cpu_seconds"], value["max_rss_kb"]}' OFS='\t' "$timefile")
    IFS=$'\t' read -r raw trimmed mapped unique duprate <<< "$values"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$condition" "$raw" "$trimmed" "$shortened" "$mapped" "$unique" "$duprate" "$mean_len" "$timing"
  done
} > "$outfile.tmp"
mv "$outfile.tmp" "$outfile"
echo "Wrote $outfile"
