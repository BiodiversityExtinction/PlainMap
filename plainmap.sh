#!/usr/bin/env bash
# PlainMap v0.2.1
# Lightweight, restartable mapping pipeline for ancient-fast, ancient-all, and modern DNA
#
# Key properties (current dev version):
# - Manifest can contain mixed SE + PE, mixed platforms, mixed naming conventions
# - SE/PE classification is derived from FASTQ headers:
#     * Illumina/CASAVA: header contains ' 1:' / ' 2:'
#     * ENA/older Illumina: token1 or token2 ends with /1 or /2
#     * SRA-export FASTQ (often lacks mate markers): R1/R2 headers are identical;
#       PlainMap pairs files by identical normalized first-read name (key) and
#       assigns mates deterministically by manifest order (first=R1, second=R2),
#       with strict sanity checks.
# - fastp runs per unit (or per unit-chunk), then pools:
#     * pooled unmerged PE (R1+R2) kept separate
#     * pooled SE-like reads (merged + unpaired + SE) kept separate
# - mapping is restartable:
#     * per-chunk resume (won't remap already-produced valid BAMs)
#     * chunk ID lists persisted, so resume works even if mapchunk FASTQs were deleted
# - pigz is used automatically if available (decompress/test/compress)
# - selectable duplicate removal: hybrid DeDup/samtools or v0.1-compatible samtools-only
# - disk-space friendly: can delete chunk FASTQs and mapchunk FASTQs as soon as safe
#
# Library types:
#   modern       : bwa mem; maps unmerged PE as PE and SE-like as SE
#   ancient-fast : bwa aln/samse; maps SE-like only (merged + unpaired + SE)
#   ancient-all  : bwa aln; maps unmerged PE with aln+sampe AND SE-like with aln+samse
#
set -euo pipefail

###############################################################################
# REQUIRE BASH
###############################################################################
if [ -z "${BASH_VERSION:-}" ]; then
  echo "ERROR: PlainMap must be run with bash (not sh)."
  exit 1
fi

###############################################################################
# DEFAULTS
###############################################################################
VERSION="0.2.1"

THREADS=1
MINLEN=30
MISMATCH=0.01              # bwa aln -n (ancient-fast/ancient-all)
LIBTYPE="modern"           # modern|ancient-fast|ancient-all
MAX_READS_PER_CHUNK=0      # 0 disables pre-fastp chunking safety valve
PILOT_FRAGMENTS=0          # 0 disables pilot; otherwise PER-UNIT cap on fragments BEFORE fastp
MAPQ=20                    # mapping quality threshold (samtools view -q)
DEDUP_METHOD="hybrid"      # hybrid|samtools
EMIT_UNMAPPED_FASTQ=0      # 1 exports unmapped reads to FASTQ without changing final BAM/stats behavior

ADAPTER_R1=""
ADAPTER_R2=""
TRIM_ONLY=0

# mapDamage toggle (default OFF for all library types)
RUN_MAPDAMAGE=0
# When RUN_MAPDAMAGE=1 and --threads > 1, mapDamage is split by whole contig and
# run in parallel, then merged; set to 1 to force the old single-process behavior.
MAPDAMAGE_SERIAL=0

DRYRUN=0
VALIDATE_ONLY=0
RESUME=1
KEEP_INTERMEDIATE=0
TMPDIR_USER=""
RESET=0

FASTP=${FASTP:-fastp}
BWA=${BWA:-bwa}
SAMTOOLS=${SAMTOOLS:-samtools}
DEDUP=${DEDUP:-dedup}
MAPDAMAGE=${MAPDAMAGE:-mapDamage}
PYTHON=${PYTHON:-python3}

# Disk-space behavior:
#   1 = delete chunk/mapchunk FASTQs as soon as no longer needed
#   0 = keep chunk/mapchunk FASTQs until final cleanup (or forever with --keep-intermediate)
DELETE_AS_YOU_GO=1

###############################################################################
# HELP
###############################################################################
usage() {
cat <<EOF
PlainMap v$VERSION

Required:
  --manifest FILE
  --prefix STRING
  --ref FILE
  --outdir DIR

Optional:
  --library-type modern|ancient-fast|ancient-all
  --threads INT                              (default: 1)
  --minlength INT                            (default: 30)
  --mismatch FLOAT                           (ancient-fast/ancient-all; bwa aln -n; default: 0.01)
  --mapq INT                                 Mapping quality threshold (default: 20)
  --dedup-method hybrid|samtools             (default: hybrid; samtools matches v0.1 behavior)
  --emit-unmapped-fastq                      Export unmapped reads to FASTQ
  --max-reads-per-chunk INT                  (default: 0; disabled) pre-fastp chunking safety valve
  --pilot-fragments INT                      (default: 0; disabled) PER-UNIT cap on fragments before fastp
  --adapter-r1 SEQ
  --adapter-r2 SEQ
  --trim-only                                Trim only (fastp); exit before mapping
  --run-mapdamage                            Run mapDamage (any library type; default OFF)
  --mapdamage-serial                         Disable automatic contig-split mapDamage parallelization
  --tmpdir DIR                               Custom temp dir (e.g. node-local scratch)
  --keep-intermediate                        Keep <outdir>/<prefix>/work
  --no-resume                                Disable resume (default is resume ON)
  --no-delete-as-you-go                      Keep chunk/mapchunk FASTQs until cleanup
  --dry-run                                  Print plan only
  --validate                                 Check tools + gzip/pigz -t all manifest FASTQs, then exit
  --reset                                    Restart from scratch and accept the current configuration

Tool overrides:
  --fastp CMD
  --bwa CMD
  --samtools CMD
  --dedup CMD
  --mapdamage CMD
  --python CMD

Resume safety:
  PlainMap records the effective run configuration for each prefix. Any changed
  parameter, input/output path, or tool command causes an error before processing.
  Use --reset to intentionally restart with a different configuration.
EOF
exit 0
}

###############################################################################
# PARSE ARGS (double-dash only)
###############################################################################
MANIFEST=""
SAMPLE=""
REF=""
OUT=""
ORIG_ARGS=( "$@" )

need_arg() {
  local opt="$1"
  local val="${2-}"
  if [[ -z "$val" || "$val" == --* ]]; then
    echo "ERROR: $opt requires a value"
    exit 1
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) need_arg "$1" "${2-}"; MANIFEST="$2"; shift 2 ;;
    --prefix) need_arg "$1" "${2-}"; SAMPLE="$2"; shift 2 ;;
    --ref) need_arg "$1" "${2-}"; REF="$2"; shift 2 ;;
    --outdir) need_arg "$1" "${2-}"; OUT="$2"; shift 2 ;;
    --library-type) need_arg "$1" "${2-}"; LIBTYPE="$2"; shift 2 ;;
    --threads) need_arg "$1" "${2-}"; THREADS="$2"; shift 2 ;;
    --minlength) need_arg "$1" "${2-}"; MINLEN="$2"; shift 2 ;;
    --mismatch) need_arg "$1" "${2-}"; MISMATCH="$2"; shift 2 ;;
    --mapq) need_arg "$1" "${2-}"; MAPQ="$2"; shift 2 ;;
    --dedup-method) need_arg "$1" "${2-}"; DEDUP_METHOD="$2"; shift 2 ;;
    --emit-unmapped-fastq) EMIT_UNMAPPED_FASTQ=1; shift ;;
    --no-emit-unmapped-fastq) EMIT_UNMAPPED_FASTQ=0; shift ;;  # backward-compatible no-op
    --max-reads-per-chunk) need_arg "$1" "${2-}"; MAX_READS_PER_CHUNK="$2"; shift 2 ;;
    --pilot-fragments) need_arg "$1" "${2-}"; PILOT_FRAGMENTS="$2"; shift 2 ;;
    --adapter-r1) need_arg "$1" "${2-}"; ADAPTER_R1="$2"; shift 2 ;;
    --adapter-r2) need_arg "$1" "${2-}"; ADAPTER_R2="$2"; shift 2 ;;
    --trim-only) TRIM_ONLY=1; shift ;;
    --run-mapdamage) RUN_MAPDAMAGE=1; shift ;;
    --no-mapdamage) RUN_MAPDAMAGE=0; shift ;;                  # backward-compatible no-op
    --mapdamage-serial) MAPDAMAGE_SERIAL=1; shift ;;
    --tmpdir) need_arg "$1" "${2-}"; TMPDIR_USER="$2"; shift 2 ;;
    --keep-intermediate) KEEP_INTERMEDIATE=1; shift ;;
    --resume) RESUME=1; shift ;;
    --no-resume) RESUME=0; shift ;;
    --delete-as-you-go) DELETE_AS_YOU_GO=1; shift ;;
    --no-delete-as-you-go) DELETE_AS_YOU_GO=0; shift ;;
    --dry-run) DRYRUN=1; shift ;;
    --validate) VALIDATE_ONLY=1; shift ;;
    --reset) RESET=1; shift ;;
    --fastp) need_arg "$1" "${2-}"; FASTP="$2"; shift 2 ;;
    --bwa) need_arg "$1" "${2-}"; BWA="$2"; shift 2 ;;
    --samtools) need_arg "$1" "${2-}"; SAMTOOLS="$2"; shift 2 ;;
    --dedup) need_arg "$1" "${2-}"; DEDUP="$2"; shift 2 ;;
    --mapdamage) need_arg "$1" "${2-}"; MAPDAMAGE="$2"; shift 2 ;;
    --python) need_arg "$1" "${2-}"; PYTHON="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown option $1"; usage ;;
  esac
done

[[ -z "$MANIFEST" || -z "$SAMPLE" || -z "$REF" || -z "$OUT" ]] && usage
[[ "$LIBTYPE" == "modern" || "$LIBTYPE" == "ancient-fast" || "$LIBTYPE" == "ancient-all" ]] || {
  echo "ERROR: --library-type must be modern|ancient-fast|ancient-all"
  exit 1
}
[[ "$DEDUP_METHOD" == "hybrid" || "$DEDUP_METHOD" == "samtools" ]] || {
  echo "ERROR: --dedup-method must be hybrid|samtools"
  exit 1
}
[[ "$MAPQ" =~ ^[0-9]+$ ]] || { echo "ERROR: --mapq must be an integer"; exit 1; }
[[ "$THREADS" =~ ^[0-9]+$ ]] || { echo "ERROR: --threads must be an integer"; exit 1; }
[[ "$MAX_READS_PER_CHUNK" =~ ^[0-9]+$ ]] || { echo "ERROR: --max-reads-per-chunk must be an integer"; exit 1; }
[[ "$PILOT_FRAGMENTS" =~ ^[0-9]+$ ]] || { echo "ERROR: --pilot-fragments must be an integer"; exit 1; }

###############################################################################
# PATH RESOLUTION
###############################################################################
resolve_path() {
  "$PYTHON" - "$1" <<'PY'
import os,sys
p=sys.argv[1]
p=os.path.expanduser(p)
print(os.path.realpath(os.path.abspath(p)))
PY
}

resolve_abs_path() {
  "$PYTHON" - "$1" <<'PY'
import os,sys
p=sys.argv[1]
p=os.path.expanduser(p)
print(os.path.abspath(p))
PY
}

mkdir -p "$OUT"
MANIFEST="$(resolve_path "$MANIFEST")"
REF="$(resolve_path "$REF")"
OUT="$(resolve_abs_path "$OUT")"
if [[ -n "$TMPDIR_USER" ]]; then
  TMPDIR_USER="$(resolve_abs_path "$TMPDIR_USER")"
fi
MANIFEST_DIR="$(dirname "$MANIFEST")"

###############################################################################
# LOGGING
###############################################################################
LOG="$OUT/${SAMPLE}_plainmap.log"
if [[ -f "$LOG" ]]; then
  ts=$(date '+%Y%m%d-%H%M%S')
  mv "$LOG" "$OUT/${SAMPLE}_plainmap.${ts}.log"
fi
exec > >(tee "$LOG") 2>&1

log() { echo "[$(date '+%F %T')] $*"; }
die() { log "ERROR: $*"; exit 1; }

PIPE_T0=$(date +%s)
SORT_THREADS=$(( THREADS > 2 ? THREADS / 2 : 1 ))

###############################################################################
# WORK DIRS
###############################################################################
WORK="$OUT/${SAMPLE}/work"
RAW="$WORK/raw"
CHUNKS="$WORK/chunks"          # pre-fastp chunks
TMPBASE="$WORK/tmp"
MAPCHUNKS="$WORK/mapchunks"    # post-fastp mapping chunks
BAMS="$WORK/bams"
MDCHUNKS="$WORK/mapdamage_chunks"
UNMAPPED="$WORK/unmapped"
FINAL="$WORK/final"
CKPT="$WORK/.checkpoints"
REPORTS="$OUT/${SAMPLE}/fastp_reports"
DEDUP_REPORTS="$OUT/${SAMPLE}/dedup_reports"

STATS="$OUT/${SAMPLE}.plainmap.stats.tsv"
COV_TSV="$OUT/${SAMPLE}.plainmap.coverage.tsv"
OUT_BAM="$OUT/${SAMPLE}.final.RG.bam"
RUN_CONFIG="$OUT/${SAMPLE}.plainmap.run-config.tsv"
UNMAPPED_R1_FQ="$OUT/${SAMPLE}.unmapped.R1.fastq.gz"
UNMAPPED_R2_FQ="$OUT/${SAMPLE}.unmapped.R2.fastq.gz"
UNMAPPED_SE_FQ="$OUT/${SAMPLE}.unmapped.SE.fastq.gz"

RAW_META="$CKPT/${SAMPLE}.raw_counts.tsv"
PE_IDLIST="$CKPT/${SAMPLE}.mapchunk_ids.pe.txt"
SE_IDLIST="$CKPT/${SAMPLE}.mapchunk_ids.se.txt"

if [[ -n "$TMPDIR_USER" ]]; then
  TMP="$TMPDIR_USER/plainmap_${SAMPLE}"
else
  TMP="$TMPBASE"
fi

mkdir -p "$RAW" "$CHUNKS" "$MAPCHUNKS" "$BAMS" "$MDCHUNKS" "$UNMAPPED" "$FINAL" "$CKPT" "$TMP" "$REPORTS" "$DEDUP_REPORTS"

###############################################################################
# COMPRESSION HELPERS (pigz if available)
###############################################################################
HAVE_PIGZ=0
if command -v pigz >/dev/null 2>&1; then
  HAVE_PIGZ=1
fi

if [[ $HAVE_PIGZ -eq 1 ]]; then
  GZ_TEST=( pigz -t )
  ZCAT=( pigz -dc )
  GZIP=( pigz -p "$THREADS" )
else
  GZ_TEST=( gzip -t )
  ZCAT=( zcat )
  GZIP=( gzip )
fi

###############################################################################
# BASIC FILE/VALIDITY HELPERS
###############################################################################
file_nonempty() { [[ -e "$1" && -s "$1" ]]; }
file_exists() { [[ -e "$1" ]]; }

bam_ok() {
  local b="$1"
  file_nonempty "$b" || return 1
  "$SAMTOOLS" quickcheck -v "$b" >/dev/null 2>&1
}

bam_ok_with_retries() {
  local b="$1"
  local attempts="${2:-1}"
  local sleep_s="${3:-0}"
  local i=1

  while [[ "$i" -le "$attempts" ]]; do
    if bam_ok "$b"; then
      return 0
    fi
    if [[ "$i" -lt "$attempts" ]]; then
      log "  BAM validation not yet passing for $(basename "$b") (attempt $i/$attempts); retrying after ${sleep_s}s"
      sleep "$sleep_s"
    fi
    i=$((i+1))
  done
  return 1
}

###############################################################################
# CHECKPOINTS
###############################################################################
ckpt_file() { echo "$CKPT/${SAMPLE}.$1.done"; }
ckpt_mark() { [[ $RESUME -eq 1 ]] && date -Is > "$(ckpt_file "$1")"; }
ckpt_clear() { rm -f "$(ckpt_file "$1")" 2>/dev/null || true; }

# Use with real output files that must be non-empty
ckpt_ok() {
  local step="$1"; shift
  [[ $RESUME -eq 1 ]] || return 1
  [[ -f "$(ckpt_file "$step")" ]] || return 1
  local f
  for f in "$@"; do
    file_nonempty "$f" || return 1
  done
  return 0
}

# Use for metadata/list files that may legitimately be empty (must exist)
ckpt_ok_exists() {
  local step="$1"; shift
  [[ $RESUME -eq 1 ]] || return 1
  [[ -f "$(ckpt_file "$step")" ]] || return 1
  local f
  for f in "$@"; do
    file_exists "$f" || return 1
  done
  return 0
}

if [[ $RESET -eq 1 ]]; then
  log "Reset requested: restarting from scratch for sample $SAMPLE"
  rm -rf "$WORK" "$REPORTS" "$DEDUP_REPORTS" "$OUT/${SAMPLE}_mapdamage" 2>/dev/null || true
  rm -f "$OUT/${SAMPLE}.final.RG.bam" "$OUT/${SAMPLE}.final.RG.bam.bai" \
        "$STATS" "$COV_TSV" "$RUN_CONFIG" \
        "$UNMAPPED_R1_FQ" "$UNMAPPED_R2_FQ" "$UNMAPPED_SE_FQ" 2>/dev/null || true
  rm -f "$CKPT/${SAMPLE}."*.done "$RAW_META" "$PE_IDLIST" "$SE_IDLIST" 2>/dev/null || true
  mkdir -p "$RAW" "$CHUNKS" "$MAPCHUNKS" "$BAMS" "$MDCHUNKS" "$UNMAPPED" "$FINAL" "$CKPT" "$TMP" "$REPORTS" "$DEDUP_REPORTS"
fi

###############################################################################
# PREFLIGHT
###############################################################################
require_cmd() { command -v "$1" >/dev/null 2>&1 || die "Tool not found in PATH: $1"; }
check_split_filter_support() { split --help 2>/dev/null | grep -q -- '--filter' || die "Need GNU split with --filter support"; }

quick_preflight() {
  [[ -r "$MANIFEST" ]] || die "Manifest not readable: $MANIFEST"
  [[ -r "$REF" ]] || die "Reference not readable: $REF"
  require_cmd split
  require_cmd "$FASTP"; require_cmd "$BWA"; require_cmd "$SAMTOOLS"; require_cmd "$PYTHON"
  [[ "$DEDUP_METHOD" == "samtools" ]] || require_cmd "$DEDUP"
  check_split_filter_support
  if [[ $HAVE_PIGZ -eq 1 ]]; then
    require_cmd pigz
  else
    require_cmd gzip
    require_cmd zcat
  fi
  if [[ "$RUN_MAPDAMAGE" -eq 1 ]]; then
    require_cmd "$MAPDAMAGE"
  fi
}

###############################################################################
# VALIDATE / DRYRUN
###############################################################################
if [[ $VALIDATE_ONLY -eq 1 ]]; then
  log "Validation mode: checking tools + gzip/pigz -t for all manifest FASTQs"
  quick_preflight
  n=0
  while read -r f; do
    f="${f#"${f%%[![:space:]]*}"}"
    f="${f%"${f##*[![:space:]]}"}"
    [[ -z "$f" || "${f:0:1}" == "#" ]] && continue
    n=$((n+1))
    [[ "$f" = /* ]] || f="$MANIFEST_DIR/$f"
    f="$(resolve_path "$f")"
    [[ -r "$f" ]] || die "FASTQ not readable: $f"
    "${GZ_TEST[@]}" "$f" >/dev/null || die "Corrupt/truncated gzip: $f"
  done < "$MANIFEST"
  [[ $n -gt 0 ]] || die "Manifest empty: $MANIFEST"
  log "Validation OK"
  exit 0
fi

if [[ $DRYRUN -eq 1 ]]; then
  log "Dry-run mode: printing plan only"
  quick_preflight
  log "Resolved paths:"
  log "  manifest: $MANIFEST"
  log "  ref:      $REF"
  log "  outdir:   $OUT"
  log "Settings:"
  log "  library_type:        $LIBTYPE"
  log "  threads:             $THREADS"
  log "  minlength:           $MINLEN"
  log "  mismatch (aln):      $MISMATCH"
  log "  mapq:                $MAPQ"
  log "  emit_unmapped_fastq: $EMIT_UNMAPPED_FASTQ"
  log "  pilot_fragments:     $PILOT_FRAGMENTS (PER-UNIT cap)"
  log "  max_reads_per_chunk: $MAX_READS_PER_CHUNK"
  log "  run_mapdamage:       $RUN_MAPDAMAGE (any library type)"
  log "  dedup_method:        $DEDUP_METHOD"
  log "  delete_as_you_go:    $DELETE_AS_YOU_GO"
  log "Plan: manifest -> build SE/PE units -> pre-fastp chunk (optional) -> pilot per-unit (optional) -> fastp -> pool -> mapchunk -> map (per-chunk resume) -> merge (-f overwrite) -> dedup ($DEDUP_METHOD) -> RG -> (mapDamage optional) -> coverage -> stats"
  exit 0
fi

###############################################################################
# RESUME CONFIGURATION SAFETY
###############################################################################
file_sha256() {
  "$PYTHON" - "$1" <<'PY'
import hashlib
import sys

h = hashlib.sha256()
with open(sys.argv[1], "rb") as handle:
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        h.update(block)
print(h.hexdigest())
PY
}

file_identity() {
  "$PYTHON" - "$1" <<'PY'
import os
import sys

st = os.stat(sys.argv[1])
print(f"{st.st_size}:{st.st_mtime_ns}")
PY
}

write_run_config() {
  local destination="$1"
  local script_path manifest_hash script_hash ref_identity
  script_path="$(resolve_path "$0")"
  manifest_hash="$(file_sha256 "$MANIFEST")"
  script_hash="$(file_sha256 "$script_path")"
  ref_identity="$(file_identity "$REF")"

  {
    printf 'config_schema\t3\n'
    printf 'plainmap_version\t%s\n' "$VERSION"
    printf 'plainmap_script\t%s\n' "$script_path"
    printf 'plainmap_script_sha256\t%s\n' "$script_hash"
    printf 'manifest\t%s\n' "$MANIFEST"
    printf 'manifest_sha256\t%s\n' "$manifest_hash"
    printf 'prefix\t%s\n' "$SAMPLE"
    printf 'reference\t%s\n' "$REF"
    printf 'reference_size_mtime_ns\t%s\n' "$ref_identity"
    printf 'outdir\t%s\n' "$OUT"
    printf 'library_type\t%s\n' "$LIBTYPE"
    printf 'threads\t%s\n' "$THREADS"
    printf 'minlength\t%s\n' "$MINLEN"
    printf 'mismatch\t%s\n' "$MISMATCH"
    printf 'mapq\t%s\n' "$MAPQ"
    printf 'dedup_method\t%s\n' "$DEDUP_METHOD"
    printf 'emit_unmapped_fastq\t%s\n' "$EMIT_UNMAPPED_FASTQ"
    printf 'max_reads_per_chunk\t%s\n' "$MAX_READS_PER_CHUNK"
    printf 'pilot_fragments\t%s\n' "$PILOT_FRAGMENTS"
    printf 'adapter_r1\t%s\n' "${ADAPTER_R1:-<auto>}"
    printf 'adapter_r2\t%s\n' "${ADAPTER_R2:-<auto>}"
    printf 'trim_only\t%s\n' "$TRIM_ONLY"
    printf 'run_mapdamage\t%s\n' "$RUN_MAPDAMAGE"
    printf 'mapdamage_serial\t%s\n' "$MAPDAMAGE_SERIAL"
    printf 'tmpdir\t%s\n' "${TMPDIR_USER:-<default>}"
    printf 'keep_intermediate\t%s\n' "$KEEP_INTERMEDIATE"
    printf 'resume\t%s\n' "$RESUME"
    printf 'delete_as_you_go\t%s\n' "$DELETE_AS_YOU_GO"
    printf 'fastp\t%s\n' "$FASTP"
    printf 'bwa\t%s\n' "$BWA"
    printf 'samtools\t%s\n' "$SAMTOOLS"
    printf 'dedup\t%s\n' "$DEDUP"
    printf 'mapdamage\t%s\n' "$MAPDAMAGE"
    printf 'python\t%s\n' "$PYTHON"
  } > "$destination"
}

prior_run_artifacts_exist() {
  [[ -e "$OUT_BAM" || -e "${OUT_BAM}.bai" || -e "$STATS" || -e "$COV_TSV" ]] && return 0
  [[ -n "$(find "$WORK" "$REPORTS" "$DEDUP_REPORTS" -type f -print -quit 2>/dev/null)" ]] && return 0
  return 1
}

report_config_changes() {
  "$PYTHON" - "$RUN_CONFIG" "$1" <<'PY'
import sys

def read_config(path):
    values = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            key, separator, value = line.rstrip("\n").partition("\t")
            if separator:
                values[key] = value
    return values

previous = read_config(sys.argv[1])
current = read_config(sys.argv[2])
for key in sorted(previous.keys() | current.keys()):
    if previous.get(key) != current.get(key):
        print(f"  {key}: previous={previous.get(key, '<missing>')!r}; current={current.get(key, '<missing>')!r}")
PY
}

enforce_run_config() {
  local current_config="${RUN_CONFIG}.current.$$"
  write_run_config "$current_config"

  if [[ -e "$RUN_CONFIG" ]]; then
    if ! cmp -s "$RUN_CONFIG" "$current_config"; then
      log "Resume configuration mismatch for sample $SAMPLE:"
      report_config_changes "$current_config"
      rm -f "$current_config"
      die "Refusing to reuse existing outputs or checkpoints with changed settings. Re-run with the original configuration, or use --reset to restart from scratch."
    fi
    rm -f "$current_config"
    log "Run configuration matches: $RUN_CONFIG"
    return
  fi

  if prior_run_artifacts_exist; then
    rm -f "$current_config"
    die "Existing outputs or checkpoints were found, but no run-configuration file exists. PlainMap cannot verify that this legacy run used the same settings. Use --reset to restart from scratch, or use a new --prefix."
  fi

  mv "$current_config" "$RUN_CONFIG"
  log "Run configuration recorded: $RUN_CONFIG"
}

enforce_run_config
quick_preflight

log_run_configuration() {
  local cmdline_q=""
  if [[ ${#ORIG_ARGS[@]} -gt 0 ]]; then
    cmdline_q="$(printf '%q ' "${ORIG_ARGS[@]}")"
  fi

  log "Run configuration:"
  log "  sample:              $SAMPLE"
  log "  library_type:        $LIBTYPE"
  log "  manifest:            $MANIFEST"
  log "  ref:                 $REF"
  log "  outdir:              $OUT"
  log "  threads:             $THREADS"
  log "  minlength:           $MINLEN"
  log "  mismatch (aln):      $MISMATCH"
  log "  mapq:                $MAPQ"
  log "  dedup_method:        $DEDUP_METHOD"
  log "  pilot_fragments:     $PILOT_FRAGMENTS"
  log "  max_reads_per_chunk: $MAX_READS_PER_CHUNK"
  log "  emit_unmapped_fastq: $EMIT_UNMAPPED_FASTQ"
  log "  run_mapdamage:       $RUN_MAPDAMAGE"
  log "  mapdamage_serial:    $MAPDAMAGE_SERIAL"
  log "  resume:              $RESUME"
  log "  delete_as_you_go:    $DELETE_AS_YOU_GO"
  log "  keep_intermediate:   $KEEP_INTERMEDIATE"
  log "  reset:               $RESET"
  log "  trim_only:           $TRIM_ONLY"
  log "  fastp:               $FASTP"
  log "  bwa:                 $BWA"
  log "  samtools:            $SAMTOOLS"
  [[ "$DEDUP_METHOD" == "hybrid" ]] && log "  dedup:               $DEDUP" || log "  dedup:               (not used)"
  log "  mapdamage:           $MAPDAMAGE"
  log "  python:              $PYTHON"
  [[ -n "$TMPDIR_USER" ]] && log "  tmpdir override:     $TMPDIR_USER" || log "  tmpdir override:     (none)"
  [[ -n "$ADAPTER_R1" ]] && log "  adapter_r1:          set" || log "  adapter_r1:          auto"
  [[ -n "$ADAPTER_R2" ]] && log "  adapter_r2:          set" || log "  adapter_r2:          auto"
  [[ -n "$cmdline_q" ]] && log "  cli_args:            $cmdline_q"
}

log_run_configuration

###############################################################################
# BWA INDEX (PARALLEL SAFE + STALE LOCK RECOVERY)
###############################################################################
ensure_bwa_index() {
  local ref="$1"
  local lock="${ref}.bwa.lock"
  local idx=( "${ref}.bwt" "${ref}.pac" "${ref}.ann" "${ref}.amb" "${ref}.sa" )
  local missing=0
  local f

  missing=0
  for f in "${idx[@]}"; do [[ -f "$f" ]] || missing=1; done
  [[ "$missing" -eq 0 ]] && { log "Reference index: present"; return; }

  # If lock exists and is stale (>1h), remove it.
  if [[ -d "$lock" ]]; then
    local now mtime age
    now=$(date +%s)
    mtime=$(stat -c %Y "$lock" 2>/dev/null || echo 0)
    age=$(( now - mtime ))
    if [[ "$mtime" -gt 0 && "$age" -gt 3600 ]]; then
      log "WARNING: BWA lock appears stale (>1h) and index is missing. Removing stale lock: $lock"
      rm -rf "$lock" || true
    fi
  fi

  if mkdir "$lock" 2>/dev/null; then
    trap 'rm -rf "$lock"' EXIT
    log "Indexing reference: $ref"
    "$BWA" index "$ref"
    rm -rf "$lock"; trap - EXIT
  else
    log "Waiting for BWA index to finish (lock exists: $lock)"
    while :; do
      sleep 15
      missing=0
      for f in "${idx[@]}"; do [[ -f "$f" ]] || missing=1; done
      [[ "$missing" -eq 0 ]] && break

      if [[ -d "$lock" ]]; then
        local now2 mtime2 age2
        now2=$(date +%s)
        mtime2=$(stat -c %Y "$lock" 2>/dev/null || echo 0)
        age2=$(( now2 - mtime2 ))
        if [[ "$mtime2" -gt 0 && "$age2" -gt 3600 ]]; then
          log "WARNING: BWA lock became stale while waiting (>1h). Removing and retrying indexing."
          rm -rf "$lock" || true
          ensure_bwa_index "$ref"
          return
        fi
      fi
    done
    log "Reference index: present (built by another process)"
  fi
}
ensure_bwa_index "$REF"

###############################################################################
# FASTP WRAPPER + ARGS (JSON only)
###############################################################################
run_fastp() {
  local stderr_file="$1"; shift
  log "$FASTP $*"
  : > "$stderr_file"
  set +e
  "$FASTP" "$@" 2> >(tee -a "$stderr_file" >&2)
  local status=$?
  set -e
  [[ $status -eq 0 ]] || die "fastp exited non-zero ($status)"
  if grep -qiE "igzip|invalid gzip|premature|truncated|error" "$stderr_file"; then
    die "fastp reported gzip/igzip error (see $stderr_file)"
  fi
}

FASTP_ARGS=( -l "$MINLEN" -g -w "$THREADS" )
[[ -n "$ADAPTER_R1" ]] && FASTP_ARGS+=( --adapter_sequence "$ADAPTER_R1" )
[[ -n "$ADAPTER_R2" ]] && FASTP_ARGS+=( --adapter_sequence_r2 "$ADAPTER_R2" )

###############################################################################
# JSON HELPERS (robust across fastp schemas)
###############################################################################
json_get_int_any() {
  "$PYTHON" - "$@" <<'PY'
import json,sys
path=sys.argv[1]
keys_list=sys.argv[2:]
with open(path) as f:
    d=json.load(f)

def get_path(obj, path):
    cur=obj
    for k in path.split("."):
        if isinstance(cur, dict) and k in cur:
            cur=cur[k]
        else:
            return None
    return cur

for k in keys_list:
    v=get_path(d,k)
    if v is None:
        continue
    try:
        iv=int(v)
    except Exception:
        continue
    if iv != 0:
        print(iv)
        sys.exit(0)

for k in keys_list:
    v=get_path(d,k)
    if v is None:
        continue
    try:
        print(int(v))
        sys.exit(0)
    except Exception:
        pass

print(0)
PY
}

###############################################################################
# STATS HELPERS
###############################################################################
count_reads_fastq_gz() { "${ZCAT[@]}" "$1" | awk 'END{print NR/4}'; }

py_add() { "$PYTHON" - "$1" "$2" <<'PY'
import sys
print(int(sys.argv[1]) + int(sys.argv[2]))
PY
}

py_sub_nonneg() {
  "$PYTHON" - "$1" "$2" <<'PY'
import sys
a=int(sys.argv[1]); b=int(sys.argv[2])
x=a-b
print(x if x>=0 else 0)
PY
}

###############################################################################
# FASTQ HEADER PARSING / DIRECTION / PAIR KEYS (HEADER-ONLY, SRA-SAFE)
###############################################################################
first_header_line() {
  local fq="$1"
  "${ZCAT[@]}" "$fq" | head -n 1 || true
}

sample_header_lines() {
  local fq="$1" n="${2:-3}"
  "${ZCAT[@]}" "$fq" | awk -v max="$n" 'NR % 4 == 1 { print; c++; if (c >= max) exit }' || true
}

normalize_qname() {
  local hdr="$1"
  local first="${hdr%% *}"
  first="${first#@}"
  first="${first%/1}"
  first="${first%/2}"
  echo "$first"
}

classify_fastq_direction() {
  local fq="$1"
  local hdr firsttok rest secondtok
  hdr="$(first_header_line "$fq")"
  [[ -n "$hdr" ]] || { echo "UNKNOWN"; return; }

  if [[ "$hdr" == *" 1:"* ]]; then echo "R1"; return; fi
  if [[ "$hdr" == *" 2:"* ]]; then echo "R2"; return; fi

  firsttok="${hdr%% *}"
  rest="${hdr#"$firsttok"}"; rest="${rest# }"
  secondtok="${rest%% *}"

  if [[ "$secondtok" == */1 ]]; then echo "R1"; return; fi
  if [[ "$secondtok" == */2 ]]; then echo "R2"; return; fi

  if [[ "$firsttok" == */1 ]]; then echo "R1"; return; fi
  if [[ "$firsttok" == */2 ]]; then echo "R2"; return; fi

  if [[ "$secondtok" == "1" || "$secondtok" == "2" ]]; then
    local sampled consistent value s_hdr s_first s_rest s_second
    sampled=0
    consistent=1
    value="$secondtok"
    while IFS= read -r s_hdr; do
      [[ -n "$s_hdr" ]] || continue
      s_first="${s_hdr%% *}"
      s_rest="${s_hdr#"$s_first"}"; s_rest="${s_rest# }"
      s_second="${s_rest%% *}"
      sampled=$((sampled+1))
      if [[ "$s_second" != "$value" ]]; then
        consistent=0
        break
      fi
    done < <(sample_header_lines "$fq" 3)

    if [[ "$sampled" -gt 0 && "$consistent" -eq 1 ]]; then
      if [[ "$value" == "1" ]]; then echo "R1"; return; fi
      if [[ "$value" == "2" ]]; then echo "R2"; return; fi
    fi
  fi

  echo "UNKNOWN"
}

check_headers_identical_sra() {
  local f1="$1" f2="$2"
  local h1 h2 k1 k2
  h1="$(first_header_line "$f1")"
  h2="$(first_header_line "$f2")"
  k1="$(normalize_qname "$h1")"
  k2="$(normalize_qname "$h2")"
  [[ "$k1" == "$k2" ]] || die "SRA/unknown-PE sanity check failed: normalized first read names differ:
$f1 -> $k1
$f2 -> $k2"
}

###############################################################################
# BUILD UNITS SAFELY (SRA unknown handled)
###############################################################################
declare -A R1_BY_KEY
declare -A R2_BY_KEY
declare -A UNK_BY_KEY

TOTAL_FILES=0
IDX=0

while read -r f; do
  f="${f#"${f%%[![:space:]]*}"}"
  f="${f%"${f##*[![:space:]]}"}"
  [[ -z "$f" || "${f:0:1}" == "#" ]] && continue
  [[ "$f" = /* ]] || f="$MANIFEST_DIR/$f"
  f="$(resolve_path "$f")"
  [[ -r "$f" ]] || die "FASTQ not readable: $f"

  TOTAL_FILES=$((TOTAL_FILES+1))
  IDX=$((IDX+1))

  hdr="$(first_header_line "$f")"
  [[ -n "$hdr" ]] || die "Could not read FASTQ header for: $f"

  dir="$(classify_fastq_direction "$f")"
  key="$(normalize_qname "$hdr")"
  [[ -n "$key" ]] || die "Could not derive header key for: $f"

  if [[ "$dir" == "R1" ]]; then
    [[ -z "${R1_BY_KEY[$key]:-}" ]] || die "Multiple R1 files share the same first-read key '$key'. Offending: $f and ${R1_BY_KEY[$key]}"
    R1_BY_KEY["$key"]="$f"
  elif [[ "$dir" == "R2" ]]; then
    [[ -z "${R2_BY_KEY[$key]:-}" ]] || die "Multiple R2 files share the same first-read key '$key'. Offending: $f and ${R2_BY_KEY[$key]}"
    R2_BY_KEY["$key"]="$f"
  else
    if [[ -z "${UNK_BY_KEY[$key]:-}" ]]; then
      UNK_BY_KEY["$key"]="${IDX}:${f}"
    else
      UNK_BY_KEY["$key"]+=$'\n'"${IDX}:${f}"
    fi
  fi
done < "$MANIFEST"

[[ $TOTAL_FILES -gt 0 ]] || die "Manifest empty: $MANIFEST"

PE_KEYS=()
SE_FILES=()

for k in "${!R2_BY_KEY[@]}"; do
  [[ -n "${R1_BY_KEY[$k]:-}" ]] || die "Found R2 without matching R1 for key '$k': ${R2_BY_KEY[$k]}"
done

for k in "${!R1_BY_KEY[@]}"; do
  if [[ -n "${R2_BY_KEY[$k]:-}" ]]; then
    PE_KEYS+=( "$k" )
  else
    SE_FILES+=( "${R1_BY_KEY[$k]}" )
  fi
done

for k in "${!UNK_BY_KEY[@]}"; do
  mapfile -t entries <<< "${UNK_BY_KEY[$k]}"
  if [[ "${#entries[@]}" -eq 1 ]]; then
    f="${entries[0]#*:}"
    log "SRA/unknown-header fallback: treating key '$k' as SE ($f)"
    SE_FILES+=( "$f" )
  elif [[ "${#entries[@]}" -eq 2 ]]; then
    e1="${entries[0]}"; e2="${entries[1]}"
    i1="${e1%%:*}"; f1="${e1#*:}"
    i2="${e2%%:*}"; f2="${e2#*:}"
    if [[ "$i2" -lt "$i1" ]]; then
      tmpi="$i1"; tmpf="$f1"
      i1="$i2"; f1="$f2"
      i2="$tmpi"; f2="$tmpf"
    fi

    check_headers_identical_sra "$f1" "$f2"
    n1=$(count_reads_fastq_gz "$f1")
    n2=$(count_reads_fastq_gz "$f2")
    [[ "$n1" -eq "$n2" ]] || die "SRA-style PE sanity check failed: read counts differ for key '$k' (R1?=$n1 R2?=$n2):
$f1
$f2"

    log "SRA/unknown-header fallback: assigning PE by manifest order for key '$k'"
    log "  R1: $f1"
    log "  R2: $f2"

    [[ -z "${R1_BY_KEY[$k]:-}" && -z "${R2_BY_KEY[$k]:-}" ]] || die "Internal error: UNKNOWN key '$k' already classified"
    R1_BY_KEY["$k"]="$f1"
    R2_BY_KEY["$k"]="$f2"
    PE_KEYS+=( "$k" )
  else
    die "Ambiguous SRA-style grouping for key '$k': found ${#entries[@]} files with no mate markers:
${UNK_BY_KEY[$k]}
PlainMap refuses to guess."
  fi
done

SEQ_MODE="SE"
if [[ ${#PE_KEYS[@]} -gt 0 && ${#SE_FILES[@]} -gt 0 ]]; then
  SEQ_MODE="MIX"
elif [[ ${#PE_KEYS[@]} -gt 0 ]]; then
  SEQ_MODE="PE"
fi

log "Detected sequencing mode: $SEQ_MODE (PE keys: ${#PE_KEYS[@]}, SE files: ${#SE_FILES[@]})"

###############################################################################
# PRE-FASTP CHUNKING HELPERS
###############################################################################
split_fastq_gz() {
  local in="$1" outprefix="$2" reads="$3"
  local lines=$((reads * 4))
  local filter_cmd
  if [[ $HAVE_PIGZ -eq 1 ]]; then
    filter_cmd="pigz -p ${THREADS} > \$FILE.fastq.gz"
  else
    filter_cmd="gzip > \$FILE.fastq.gz"
  fi
  "${ZCAT[@]}" "$in" | split -l "$lines" -d -a 4 --filter="$filter_cmd" - "$outprefix"
}

limit_fastq_gz() {
  local in="$1" out="$2" reads="$3"
  local lines=$((reads * 4))
  set +o pipefail
  "${ZCAT[@]}" "$in" | head -n "$lines" | "${GZIP[@]}" > "$out"
  set -o pipefail
  file_nonempty "$out" || die "Pilot limiting produced empty FASTQ: $out"
}

pilot_take() {
  "$PYTHON" - "$1" "$2" <<'PY'
import sys
c=int(sys.argv[1]); p=int(sys.argv[2])
if p<=0: print(0)
else: print(min(c,p))
PY
}

chunks_exist() {
  shopt -s nullglob
  local a=( "$CHUNKS/${SAMPLE}."*.fastq.gz )
  [[ ${#a[@]} -gt 0 ]]
}

###############################################################################
# STEP 1: CREATE PRE-FASTP CHUNKS PER UNIT (PE stays paired)
###############################################################################
make_chunks() {
  log "STEP: pre-fastp chunking"
  ckpt_clear chunk
  rm -f "$CHUNKS/${SAMPLE}."*.fastq.gz 2>/dev/null || true

  local idx=0
  for k in "${PE_KEYS[@]}"; do
    local r1="${R1_BY_KEY[$k]}"
    local r2="${R2_BY_KEY[$k]}"
    [[ -n "$r1" && -n "$r2" ]] || die "Internal error: PE key '$k' missing mate"

    local unit="PE$(printf '%04d' $idx)"
    local prefix="$CHUNKS/${SAMPLE}.${unit}."

    if [[ "$MAX_READS_PER_CHUNK" -gt 0 ]]; then
      log "  PE unit $unit (key=$k): splitting into chunks (max reads/chunk: $MAX_READS_PER_CHUNK)"
      split_fastq_gz "$r1" "${prefix}R1_" "$MAX_READS_PER_CHUNK"
      split_fastq_gz "$r2" "${prefix}R2_" "$MAX_READS_PER_CHUNK"
    else
      log "  PE unit $unit (key=$k): chunking disabled; linking as single chunk"
      ln -sf "$r1" "${prefix}R1_0000.fastq.gz"
      ln -sf "$r2" "${prefix}R2_0000.fastq.gz"
    fi

    file_nonempty "${prefix}R1_0000.fastq.gz" || die "Missing PE R1 chunk for $unit"
    file_nonempty "${prefix}R2_0000.fastq.gz" || die "Missing PE R2 chunk for $unit"
    idx=$((idx+1))
  done

  local sidx=0
  for fq in "${SE_FILES[@]}"; do
    local unit="SE$(printf '%04d' $sidx)"
    local prefix="$CHUNKS/${SAMPLE}.${unit}."
    if [[ "$MAX_READS_PER_CHUNK" -gt 0 ]]; then
      log "  SE unit $unit: splitting into chunks (max reads/chunk: $MAX_READS_PER_CHUNK)"
      split_fastq_gz "$fq" "${prefix}SE_" "$MAX_READS_PER_CHUNK"
    else
      log "  SE unit $unit: chunking disabled; linking as single chunk"
      ln -sf "$fq" "${prefix}SE_0000.fastq.gz"
    fi
    file_nonempty "${prefix}SE_0000.fastq.gz" || die "Missing SE chunk for $unit"
    sidx=$((sidx+1))
  done

  ckpt_mark chunk
}

# BUGFIX: only skip chunking if checkpoint exists AND chunks still exist.
if [[ $RESUME -eq 1 && -f "$(ckpt_file chunk)" && chunks_exist ]]; then
  log "STEP: pre-fastp chunking (skipped; checkpoint present + chunks exist)"
else
  make_chunks
fi

###############################################################################
# STEP 2: FASTP PER PRE-FASTP CHUNK -> POOL TRIMMED OUTPUTS
###############################################################################
POOL_PE_R1="$RAW/${SAMPLE}.pool.PE.R1.fastq.gz"
POOL_PE_R2="$RAW/${SAMPLE}.pool.PE.R2.fastq.gz"
POOL_SE_ALL="$RAW/${SAMPLE}.pool.SE.all.fastq.gz"

run_fastp_and_pool() {
  log "STEP: fastp per unit-chunk (JSON only) + pooling (unmerged PE separate from SE-like)"
  ckpt_clear fastp_pool

  : > "$POOL_PE_R1"
  : > "$POOL_PE_R2"
  : > "$POOL_SE_ALL"

  local RAW_R1_READS=0
  local RAW_R2_READS=0
  local RAW_FRAGMENTS=0

  local TRIMMED_FRAGMENTS=0
  local MERGED_READS=0
  local UNPAIRED_READS=0

  declare -A PILOT_LEFT_BY_UNIT
  shopt -s nullglob

  # PE chunks
  local pe_r1_chunks=( "$CHUNKS/${SAMPLE}.PE"*".R1_"*.fastq.gz )
  for r1c in "${pe_r1_chunks[@]}"; do
    local base r2c unit_id unit_left in_r1 in_r2 pairs_in_chunk take_pairs
    base=$(basename "$r1c")
    r2c="${r1c/.R1_/.R2_}"
    [[ -f "$r2c" ]] || die "Missing PE mate chunk for: $r1c (expected $r2c)"

    unit_id="${base%.fastq.gz}"
    unit_id="${unit_id%%.R1_*}"

    if [[ "$PILOT_FRAGMENTS" -gt 0 && -z "${PILOT_LEFT_BY_UNIT[$unit_id]:-}" ]]; then
      PILOT_LEFT_BY_UNIT["$unit_id"]="$PILOT_FRAGMENTS"
    fi
    unit_left="${PILOT_LEFT_BY_UNIT[$unit_id]:-0}"

    if [[ "$PILOT_FRAGMENTS" -gt 0 && "$unit_left" -le 0 ]]; then
      [[ $DELETE_AS_YOU_GO -eq 1 ]] && rm -f "$r1c" "$r2c" 2>/dev/null || true
      continue
    fi

    in_r1="$r1c"
    in_r2="$r2c"

    pairs_in_chunk=$(count_reads_fastq_gz "$in_r1")
    take_pairs="$pairs_in_chunk"

    local lim_r1="" lim_r2=""
    if [[ "$PILOT_FRAGMENTS" -gt 0 ]]; then
      take_pairs=$(pilot_take "$pairs_in_chunk" "$unit_left")
      if [[ "$take_pairs" -le 0 ]]; then
        PILOT_LEFT_BY_UNIT["$unit_id"]=0
        [[ $DELETE_AS_YOU_GO -eq 1 ]] && rm -f "$r1c" "$r2c" 2>/dev/null || true
        continue
      fi
      if [[ "$take_pairs" -lt "$pairs_in_chunk" ]]; then
        log "  Pilot (per-unit): limiting PE unit $unit_id chunk $(basename "$r1c") to first $take_pairs pairs (unit remaining before: $unit_left)"
      fi
      lim_r1="$TMP/${SAMPLE}.$(basename "${r1c%.fastq.gz}").pilot.R1.fastq.gz"
      lim_r2="$TMP/${SAMPLE}.$(basename "${r2c%.fastq.gz}").pilot.R2.fastq.gz"
      limit_fastq_gz "$in_r1" "$lim_r1" "$take_pairs"
      limit_fastq_gz "$in_r2" "$lim_r2" "$take_pairs"
      in_r1="$lim_r1"
      in_r2="$lim_r2"
      unit_left="$(py_sub_nonneg "$unit_left" "$take_pairs")"
      PILOT_LEFT_BY_UNIT["$unit_id"]="$unit_left"
    fi

    RAW_R1_READS=$(py_add "$RAW_R1_READS" "$take_pairs")
    RAW_R2_READS=$(py_add "$RAW_R2_READS" "$take_pairs")
    RAW_FRAGMENTS=$(py_add "$RAW_FRAGMENTS" "$take_pairs")

    local tag fastp_err fp_json fp_html trim1 trim2 merged up1 up2
    tag="$(basename "${r1c%.fastq.gz}")"
    tag="${tag//\//_}"

    fastp_err="$TMP/fastp_${SAMPLE}.${tag}.stderr.log"
    fp_json="$REPORTS/${SAMPLE}.${tag}.fastp.json"
    fp_html="$REPORTS/${SAMPLE}.${tag}.fastp.html"

    trim1="$TMP/${SAMPLE}.${tag}.trim.R1.fastq.gz"
    trim2="$TMP/${SAMPLE}.${tag}.trim.R2.fastq.gz"
    merged="$TMP/${SAMPLE}.${tag}.merged.fastq.gz"
    up1="$TMP/${SAMPLE}.${tag}.unpaired1.fastq.gz"
    up2="$TMP/${SAMPLE}.${tag}.unpaired2.fastq.gz"

    local extra_pe=()
    if [[ -z "$ADAPTER_R1" && -z "$ADAPTER_R2" ]]; then
      extra_pe+=( --detect_adapter_for_pe )
    fi

    run_fastp "$fastp_err" "${FASTP_ARGS[@]}" "${extra_pe[@]}" \
      -i "$in_r1" -I "$in_r2" \
      --out1 "$trim1" --out2 "$trim2" \
      -m --merged_out "$merged" \
      --unpaired1 "$up1" --unpaired2 "$up2" \
      -j "$fp_json" --html "$fp_html"

    file_nonempty "$trim1" || die "fastp produced empty trimmed R1 for $tag"
    file_nonempty "$trim2" || die "fastp produced empty trimmed R2 for $tag"

    cat "$trim1" >> "$POOL_PE_R1"
    cat "$trim2" >> "$POOL_PE_R2"

    file_nonempty "$merged" && cat "$merged" >> "$POOL_SE_ALL"
    file_nonempty "$up1"   && cat "$up1"   >> "$POOL_SE_ALL"
    file_nonempty "$up2"   && cat "$up2"   >> "$POOL_SE_ALL"

    local out1_reads merged_reads up1n up2n
    out1_reads=$(json_get_int_any "$fp_json" "read1_after_filtering.total_reads")
    [[ "$out1_reads" -le 0 ]] && out1_reads=$(count_reads_fastq_gz "$trim1")
    TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$out1_reads")

    merged_reads=$(json_get_int_any "$fp_json" "merged_and_filtered.total_reads")
    if [[ "$merged_reads" -le 0 && -s "$merged" ]]; then
      merged_reads=$(count_reads_fastq_gz "$merged")
    fi
    MERGED_READS=$(py_add "$MERGED_READS" "$merged_reads")
    TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$merged_reads")

    up1n=0; up2n=0
    [[ -s "$up1" ]] && up1n=$(count_reads_fastq_gz "$up1")
    [[ -s "$up2" ]] && up2n=$(count_reads_fastq_gz "$up2")
    UNPAIRED_READS=$(py_add "$UNPAIRED_READS" "$up1n")
    UNPAIRED_READS=$(py_add "$UNPAIRED_READS" "$up2n")
    TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$up1n")
    TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$up2n")

    rm -f "$trim1" "$trim2" "$merged" "$up1" "$up2" 2>/dev/null || true
    [[ -n "$lim_r1" ]] && rm -f "$lim_r1" "$lim_r2" 2>/dev/null || true

    if [[ $DELETE_AS_YOU_GO -eq 1 ]]; then
      rm -f "$r1c" "$r2c" 2>/dev/null || true
    fi
  done

  # SE chunks
  local se_chunks=( "$CHUNKS/${SAMPLE}.SE"*".SE_"*.fastq.gz )
  for sec in "${se_chunks[@]}"; do
    local base unit_id unit_left in_se reads_in_chunk take_reads
    base=$(basename "$sec")
    unit_id="${base%.fastq.gz}"
    unit_id="${unit_id%%.SE_*}"

    if [[ "$PILOT_FRAGMENTS" -gt 0 && -z "${PILOT_LEFT_BY_UNIT[$unit_id]:-}" ]]; then
      PILOT_LEFT_BY_UNIT["$unit_id"]="$PILOT_FRAGMENTS"
    fi
    unit_left="${PILOT_LEFT_BY_UNIT[$unit_id]:-0}"

    if [[ "$PILOT_FRAGMENTS" -gt 0 && "$unit_left" -le 0 ]]; then
      [[ $DELETE_AS_YOU_GO -eq 1 ]] && rm -f "$sec" 2>/dev/null || true
      continue
    fi

    in_se="$sec"
    reads_in_chunk=$(count_reads_fastq_gz "$in_se")
    take_reads="$reads_in_chunk"

    local lim_se=""
    if [[ "$PILOT_FRAGMENTS" -gt 0 ]]; then
      take_reads=$(pilot_take "$reads_in_chunk" "$unit_left")
      if [[ "$take_reads" -le 0 ]]; then
        PILOT_LEFT_BY_UNIT["$unit_id"]=0
        [[ $DELETE_AS_YOU_GO -eq 1 ]] && rm -f "$sec" 2>/dev/null || true
        continue
      fi
      if [[ "$take_reads" -lt "$reads_in_chunk" ]]; then
        log "  Pilot (per-unit): limiting SE unit $unit_id chunk $(basename "$sec") to first $take_reads reads (unit remaining before: $unit_left)"
      fi
      lim_se="$TMP/${SAMPLE}.$(basename "${sec%.fastq.gz}").pilot.SE.fastq.gz"
      limit_fastq_gz "$in_se" "$lim_se" "$take_reads"
      in_se="$lim_se"
      unit_left="$(py_sub_nonneg "$unit_left" "$take_reads")"
      PILOT_LEFT_BY_UNIT["$unit_id"]="$unit_left"
    fi

    RAW_R1_READS=$(py_add "$RAW_R1_READS" "$take_reads")
    RAW_FRAGMENTS=$(py_add "$RAW_FRAGMENTS" "$take_reads")

    local tag fastp_err fp_json fp_html trim_se
    tag="$(basename "${sec%.fastq.gz}")"
    tag="${tag//\//_}"

    fastp_err="$TMP/fastp_${SAMPLE}.${tag}.stderr.log"
    fp_json="$REPORTS/${SAMPLE}.${tag}.fastp.json"
    fp_html="$REPORTS/${SAMPLE}.${tag}.fastp.html"
    trim_se="$TMP/${SAMPLE}.${tag}.trim.SE.fastq.gz"

    run_fastp "$fastp_err" "${FASTP_ARGS[@]}" \
      -i "$in_se" --out1 "$trim_se" \
      -j "$fp_json" --html "$fp_html"

    file_nonempty "$trim_se" || die "fastp produced empty trimmed SE for $tag"
    cat "$trim_se" >> "$POOL_SE_ALL"

    local se_after
    se_after=$(json_get_int_any "$fp_json" "read1_after_filtering.total_reads" "summary.after_filtering.total_reads")
    [[ "$se_after" -le 0 ]] && se_after=$(count_reads_fastq_gz "$trim_se")
    TRIMMED_FRAGMENTS=$(py_add "$TRIMMED_FRAGMENTS" "$se_after")

    rm -f "$trim_se" 2>/dev/null || true
    [[ -n "$lim_se" ]] && rm -f "$lim_se" 2>/dev/null || true

    if [[ $DELETE_AS_YOU_GO -eq 1 ]]; then
      rm -f "$sec" 2>/dev/null || true
    fi
  done

  {
    printf "raw_R1_reads\traw_R2_reads\traw_fragments\ttrimmed_fragments\tmerged_reads\tunpaired_reads\n"
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$RAW_R1_READS" "$RAW_R2_READS" "$RAW_FRAGMENTS" \
      "$TRIMMED_FRAGMENTS" "$MERGED_READS" "$UNPAIRED_READS"
  } > "$RAW_META"

  ckpt_mark fastp_pool
  log "fastp+pool complete. Raw counts written: $RAW_META"
}

need_fastp_files=( "$RAW_META" "$POOL_SE_ALL" )
if [[ "$SEQ_MODE" != "SE" ]]; then
  need_fastp_files+=( "$POOL_PE_R1" "$POOL_PE_R2" )
fi

if ! ckpt_ok fastp_pool "${need_fastp_files[@]}"; then
  run_fastp_and_pool
else
  log "STEP: fastp+pool (skipped; checkpoint present)"
fi

###############################################################################
# TRIM-ONLY MODE
###############################################################################
if [[ $TRIM_ONLY -eq 1 ]]; then
  log "Trim-only mode enabled; writing stats and exiting"

  RAW_R1_READS=$(awk 'NR==2{print $1}' "$RAW_META")
  RAW_R2_READS=$(awk 'NR==2{print $2}' "$RAW_META")
  RAW_FRAGMENTS=$(awk 'NR==2{print $3}' "$RAW_META")
  TRIMMED_FRAGMENTS=$(awk 'NR==2{print $4}' "$RAW_META")
  MERGED_READS=$(awk 'NR==2{print $5}' "$RAW_META")
  UNPAIRED_READS=$(awk 'NR==2{print $6}' "$RAW_META")

  {
    printf "sample\tlibrary_type\tseq_mode\tpilot_fragments\tmax_reads_per_chunk\tmapq\traw_R1_reads\traw_R2_reads\traw_fragments\ttrimmed_fragments\tmerged_reads\tunpaired_reads\tmapped_reads_all\tmapped_reads_unique\tmapped_fragments_all\tmapped_fragments_unique\tendog\tduprate\tavg_readlen\tmapped_bp\tavg_depth\tpct_covered\n"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n" \
      "$SAMPLE" "$LIBTYPE" "$SEQ_MODE" "$PILOT_FRAGMENTS" "$MAX_READS_PER_CHUNK" "$MAPQ" \
      "$RAW_R1_READS" "$RAW_R2_READS" "$RAW_FRAGMENTS" \
      "$TRIMMED_FRAGMENTS" "$MERGED_READS" "$UNPAIRED_READS"
  } > "$STATS"

  PIPE_T1=$(date +%s)
  log "STATS written: $STATS"
  log "Total wall time: $((PIPE_T1 - PIPE_T0)) seconds"
  exit 0
fi

###############################################################################
# STEP 3: POST-FASTP MAPPING CHUNKING (OPTIONAL; REUSE MAX_READS_PER_CHUNK)
###############################################################################
make_mapchunks_se() {
  local in="$1" prefix="$2"
  [[ -s "$in" ]] || return 0
  if [[ "$MAX_READS_PER_CHUNK" -gt 0 ]]; then
    log "  Mapping chunking SE: splitting pooled SE-like (max reads/chunk: $MAX_READS_PER_CHUNK)"
    split_fastq_gz "$in" "$prefix" "$MAX_READS_PER_CHUNK"
  else
    ln -sf "$in" "${prefix}0000.fastq.gz"
  fi
}

make_mapchunks_pe() {
  local in1="$1" in2="$2" prefix1="$3" prefix2="$4"
  [[ -s "$in1" && -s "$in2" ]] || return 0

  local n1 n2
  n1=$(count_reads_fastq_gz "$in1")
  n2=$(count_reads_fastq_gz "$in2")
  [[ "$n1" -eq "$n2" ]] || die "Pooled PE R1/R2 read counts differ (R1=$n1 R2=$n2). This would break sampe."

  if [[ "$MAX_READS_PER_CHUNK" -gt 0 ]]; then
    log "  Mapping chunking PE: splitting pooled PE (max reads/chunk: $MAX_READS_PER_CHUNK)"
    split_fastq_gz "$in1" "$prefix1" "$MAX_READS_PER_CHUNK"
    split_fastq_gz "$in2" "$prefix2" "$MAX_READS_PER_CHUNK"
  else
    ln -sf "$in1" "${prefix1}0000.fastq.gz"
    ln -sf "$in2" "${prefix2}0000.fastq.gz"
  fi
}

persist_mapchunk_ids() {
  # BUGFIX: These lists can be legitimately empty. We must create them regardless.
  : > "$PE_IDLIST"
  : > "$SE_IDLIST"

  shopt -s nullglob
  local f base cid

  if [[ "$LIBTYPE" != "ancient-fast" ]]; then
    for f in "$MAPCHUNKS/${SAMPLE}.PE.R1_"*.fastq.gz; do
      base=$(basename "$f")
      cid="${base#${SAMPLE}.PE.R1_}"
      cid="${cid%.fastq.gz}"
      echo "$cid" >> "$PE_IDLIST"
    done
  fi

  for f in "$MAPCHUNKS/${SAMPLE}.SE_"*.fastq.gz; do
    base=$(basename "$f")
    cid="${base#${SAMPLE}.SE_}"
    cid="${cid%.fastq.gz}"
    echo "$cid" >> "$SE_IDLIST"
  done

  sort -u -o "$PE_IDLIST" "$PE_IDLIST" 2>/dev/null || true
  sort -u -o "$SE_IDLIST" "$SE_IDLIST" 2>/dev/null || true

  ckpt_mark mapchunk
}

log "STEP: post-fastp mapping chunking"

# BUGFIX: mapchunk checkpoint must accept empty idlists; also only require PE list when PE is expected.
need_mapchunk_lists=( "$SE_IDLIST" )
if [[ "$LIBTYPE" != "ancient-fast" && "$SEQ_MODE" != "SE" ]]; then
  need_mapchunk_lists+=( "$PE_IDLIST" )
fi

if ! ckpt_ok_exists mapchunk "${need_mapchunk_lists[@]}"; then
  ckpt_clear mapchunk
  rm -f "$MAPCHUNKS/${SAMPLE}."*.fastq.gz 2>/dev/null || true

  if [[ "$LIBTYPE" != "ancient-fast" && "$SEQ_MODE" != "SE" ]]; then
    make_mapchunks_pe "$POOL_PE_R1" "$POOL_PE_R2" "$MAPCHUNKS/${SAMPLE}.PE.R1_" "$MAPCHUNKS/${SAMPLE}.PE.R2_"
  fi
  make_mapchunks_se "$POOL_SE_ALL" "$MAPCHUNKS/${SAMPLE}.SE_"

  persist_mapchunk_ids
else
  log "  mapchunk (skipped; checkpoint present)"
fi

###############################################################################
# STEP 4: MAPPING (PER-CHUNK RESUME) + DELETE MAPCHUNK FASTQs AS SOON AS MAPPED
###############################################################################
log "STEP: mapping (per-chunk resume)"
mkdir -p "$BAMS"

sam_to_chunk_bam() {
  local sam_in="$1" bam_out="$2"
  "$SAMTOOLS" view -F 4 -q "$MAPQ" -u "$sam_in" > "$bam_out"
}

capture_unmapped_from_sam() {
  local sam_in="$1" out_prefix="$2"
  [[ "$EMIT_UNMAPPED_FASTQ" -eq 1 ]] || return 0

  local comp_cmd=( gzip )
  if [[ $HAVE_PIGZ -eq 1 ]]; then
    comp_cmd=( pigz -p "$THREADS" )
  fi

  "$SAMTOOLS" fastq -f 4 -n \
    -1 >("${comp_cmd[@]}" > "${out_prefix}.R1.fastq.gz") \
    -2 >("${comp_cmd[@]}" > "${out_prefix}.R2.fastq.gz") \
    -s >("${comp_cmd[@]}" > "${out_prefix}.SE.fastq.gz") \
    -0 /dev/null \
    "$sam_in" >/dev/null
}

map_se_mem_chunk() {
  local fq="$1" outbam="$2" tag="$3" unmapped_prefix="$4"
  local tmpbam="${outbam}.tmp"
  "$BWA" mem -t "$THREADS" "$REF" "$fq" > "$TMP/${SAMPLE}.${tag}.sam"
  capture_unmapped_from_sam "$TMP/${SAMPLE}.${tag}.sam" "$unmapped_prefix"
  sam_to_chunk_bam "$TMP/${SAMPLE}.${tag}.sam" "$TMP/${SAMPLE}.${tag}.bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$tmpbam" "$TMP/${SAMPLE}.${tag}.bam"
  rm -f "$TMP/${SAMPLE}.${tag}.sam" "$TMP/${SAMPLE}.${tag}.bam"
  bam_ok "$tmpbam" || die "Produced invalid BAM (SE mem): $tmpbam"
  mv -f "$tmpbam" "$outbam"
}

map_pe_mem_chunk() {
  local fq1="$1" fq2="$2" outbam="$3" tag="$4" unmapped_prefix="$5"
  local tmpbam="${outbam}.tmp"
  "$BWA" mem -t "$THREADS" "$REF" "$fq1" "$fq2" > "$TMP/${SAMPLE}.${tag}.sam"
  capture_unmapped_from_sam "$TMP/${SAMPLE}.${tag}.sam" "$unmapped_prefix"
  sam_to_chunk_bam "$TMP/${SAMPLE}.${tag}.sam" "$TMP/${SAMPLE}.${tag}.bam"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$tmpbam" "$TMP/${SAMPLE}.${tag}.bam"
  rm -f "$TMP/${SAMPLE}.${tag}.sam" "$TMP/${SAMPLE}.${tag}.bam"
  bam_ok "$tmpbam" || die "Produced invalid BAM (PE mem): $tmpbam"
  mv -f "$tmpbam" "$outbam"
}

map_se_aln_chunk() {
  local fq="$1" outbam="$2" tag="$3" unmapped_prefix="$4"
  local tmpbam="${outbam}.tmp"
  "$BWA" aln -l 999 -n "$MISMATCH" -t "$THREADS" "$REF" "$fq" > "$TMP/${SAMPLE}.${tag}.sai"
  "$BWA" samse "$REF" "$TMP/${SAMPLE}.${tag}.sai" "$fq" > "$TMP/${SAMPLE}.${tag}.sam"
  capture_unmapped_from_sam "$TMP/${SAMPLE}.${tag}.sam" "$unmapped_prefix"
  sam_to_chunk_bam "$TMP/${SAMPLE}.${tag}.sam" "$TMP/${SAMPLE}.${tag}.bam.tmp"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$tmpbam" "$TMP/${SAMPLE}.${tag}.bam.tmp"
  rm -f "$TMP/${SAMPLE}.${tag}.sai" "$TMP/${SAMPLE}.${tag}.sam" "$TMP/${SAMPLE}.${tag}.bam.tmp"
  bam_ok "$tmpbam" || die "Produced invalid BAM (SE aln): $tmpbam"
  mv -f "$tmpbam" "$outbam"
}

map_pe_aln_chunk() {
  local fq1="$1" fq2="$2" outbam="$3" tag="$4" unmapped_prefix="$5"
  local tmpbam="${outbam}.tmp"

  local h1 h2 n1x n2x
  h1="$(first_header_line "$fq1")"; h2="$(first_header_line "$fq2")"
  n1x="$(normalize_qname "$h1")"; n2x="$(normalize_qname "$h2")"
  [[ "$n1x" == "$n2x" ]] || die "Refusing to run sampe: PE chunks first read names differ: $(basename "$fq1")='$n1x' vs $(basename "$fq2")='$n2x'"

  "$BWA" aln -l 999 -n "$MISMATCH" -t "$THREADS" "$REF" "$fq1" > "$TMP/${SAMPLE}.${tag}.1.sai"
  "$BWA" aln -l 999 -n "$MISMATCH" -t "$THREADS" "$REF" "$fq2" > "$TMP/${SAMPLE}.${tag}.2.sai"
  "$BWA" sampe "$REF" "$TMP/${SAMPLE}.${tag}.1.sai" "$TMP/${SAMPLE}.${tag}.2.sai" "$fq1" "$fq2" > "$TMP/${SAMPLE}.${tag}.sam"
  capture_unmapped_from_sam "$TMP/${SAMPLE}.${tag}.sam" "$unmapped_prefix"
  sam_to_chunk_bam "$TMP/${SAMPLE}.${tag}.sam" "$TMP/${SAMPLE}.${tag}.bam.tmp"
  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$tmpbam" "$TMP/${SAMPLE}.${tag}.bam.tmp"
  rm -f "$TMP/${SAMPLE}.${tag}.1.sai" "$TMP/${SAMPLE}.${tag}.2.sai" "$TMP/${SAMPLE}.${tag}.sam" "$TMP/${SAMPLE}.${tag}.bam.tmp"
  bam_ok "$tmpbam" || die "Produced invalid BAM (PE aln): $tmpbam"
  mv -f "$tmpbam" "$outbam"
}

map_one_pe_id() {
  local cid="$1"
  local r1c="$MAPCHUNKS/${SAMPLE}.PE.R1_${cid}.fastq.gz"
  local r2c="$MAPCHUNKS/${SAMPLE}.PE.R2_${cid}.fastq.gz"
  local outbam="$BAMS/${SAMPLE}.MAP.PE.${cid}.bam"
  local unmapped_prefix="$UNMAPPED/${SAMPLE}.PE.${cid}.unmapped"
  local ur1="${unmapped_prefix}.R1.fastq.gz"
  local ur2="${unmapped_prefix}.R2.fastq.gz"
  local use="${unmapped_prefix}.SE.fastq.gz"
  local have_unmapped=1
  [[ -f "$ur1" && -f "$ur2" && -f "$use" ]] || have_unmapped=0

  if bam_ok "$outbam"; then
    if [[ "$EMIT_UNMAPPED_FASTQ" -eq 1 && "$have_unmapped" -eq 0 ]]; then
      if [[ ! -f "$r1c" || ! -f "$r2c" ]]; then
        die "Cannot reproduce missing unmapped FASTQ outputs for PE chunk ${cid}.
A valid chunk BAM exists, but unmapped chunk outputs and one or both source mapchunk FASTQs are missing.
The recorded run configuration is unchanged, so required intermediate files were likely removed.
Restore the missing files or use --reset to restart from scratch."
      fi
      log "  PE chunk ${cid}: BAM exists but unmapped FASTQ chunk outputs are missing -> remapping required"
    else
      log "  PE chunk ${cid}: already mapped (valid BAM exists) -> skipping"
      [[ $DELETE_AS_YOU_GO -eq 1 ]] && rm -f "$r1c" "$r2c" 2>/dev/null || true
      return 0
    fi
  fi

  [[ -f "$r1c" && -f "$r2c" ]] || die "PE chunk ${cid} is missing FASTQs needed to map:
  $r1c
  $r2c
If you deleted mapchunk FASTQs, resume requires the BAM to exist and be valid."

  if [[ "$LIBTYPE" == "modern" ]]; then
    log "  Mapping PE chunk ${cid} (mem)"
    map_pe_mem_chunk "$r1c" "$r2c" "$outbam" "mem.PE.${cid}" "$unmapped_prefix"
  elif [[ "$LIBTYPE" == "ancient-all" ]]; then
    log "  Mapping PE chunk ${cid} (aln/sampe)"
    map_pe_aln_chunk "$r1c" "$r2c" "$outbam" "aln.PE.${cid}" "$unmapped_prefix"
  else
    die "Internal error: PE mapping called for ancient-fast"
  fi

  bam_ok "$outbam" || die "PE chunk ${cid}: produced BAM is invalid: $outbam"

  if [[ $DELETE_AS_YOU_GO -eq 1 ]]; then
    rm -f "$r1c" "$r2c" 2>/dev/null || true
  fi
}

map_one_se_id() {
  local cid="$1"
  local sec="$MAPCHUNKS/${SAMPLE}.SE_${cid}.fastq.gz"
  local outbam="$BAMS/${SAMPLE}.MAP.SE.${cid}.bam"
  local unmapped_prefix="$UNMAPPED/${SAMPLE}.SE.${cid}.unmapped"
  local ur1="${unmapped_prefix}.R1.fastq.gz"
  local ur2="${unmapped_prefix}.R2.fastq.gz"
  local use="${unmapped_prefix}.SE.fastq.gz"
  local have_unmapped=1
  [[ -f "$ur1" && -f "$ur2" && -f "$use" ]] || have_unmapped=0

  if bam_ok "$outbam"; then
    if [[ "$EMIT_UNMAPPED_FASTQ" -eq 1 && "$have_unmapped" -eq 0 ]]; then
      if [[ ! -f "$sec" ]]; then
        die "Cannot reproduce missing unmapped FASTQ outputs for SE chunk ${cid}.
A valid chunk BAM exists, but unmapped chunk outputs and the source mapchunk FASTQ are missing.
The recorded run configuration is unchanged, so required intermediate files were likely removed.
Restore the missing files or use --reset to restart from scratch."
      fi
      log "  SE chunk ${cid}: BAM exists but unmapped FASTQ chunk outputs are missing -> remapping required"
    else
      log "  SE chunk ${cid}: already mapped (valid BAM exists) -> skipping"
      [[ $DELETE_AS_YOU_GO -eq 1 ]] && rm -f "$sec" 2>/dev/null || true
      return 0
    fi
  fi

  [[ -f "$sec" ]] || die "SE chunk ${cid} is missing FASTQ needed to map:
  $sec
If you deleted mapchunk FASTQs, resume requires the BAM to exist and be valid."

  if [[ "$LIBTYPE" == "modern" ]]; then
    log "  Mapping SE chunk ${cid} (mem)"
    map_se_mem_chunk "$sec" "$outbam" "mem.SE.${cid}" "$unmapped_prefix"
  else
    log "  Mapping SE chunk ${cid} (aln/samse)"
    map_se_aln_chunk "$sec" "$outbam" "aln.SE.${cid}" "$unmapped_prefix"
  fi

  bam_ok "$outbam" || die "SE chunk ${cid}: produced BAM is invalid: $outbam"

  if [[ $DELETE_AS_YOU_GO -eq 1 ]]; then
    rm -f "$sec" 2>/dev/null || true
  fi
}

# Map PE chunk IDs (skip if ancient-fast)
if [[ "$LIBTYPE" != "ancient-fast" && "$SEQ_MODE" != "SE" && -s "$PE_IDLIST" ]]; then
  while read -r cid; do
    [[ -z "$cid" ]] && continue
    map_one_pe_id "$cid"
  done < "$PE_IDLIST"
fi

# Map SE chunk IDs
if [[ -s "$SE_IDLIST" ]]; then
  while read -r cid; do
    [[ -z "$cid" ]] && continue
    map_one_se_id "$cid"
  done < "$SE_IDLIST"
fi

# Verify BAMs exist for every expected chunk ID
missing=0
if [[ "$LIBTYPE" != "ancient-fast" && "$SEQ_MODE" != "SE" && -s "$PE_IDLIST" ]]; then
  while read -r cid; do
    [[ -z "$cid" ]] && continue
    bam_ok "$BAMS/${SAMPLE}.MAP.PE.${cid}.bam" || { log "ERROR: Missing/invalid PE BAM for chunk ${cid}"; missing=1; }
  done < "$PE_IDLIST"
fi
if [[ -s "$SE_IDLIST" ]]; then
  while read -r cid; do
    [[ -z "$cid" ]] && continue
    bam_ok "$BAMS/${SAMPLE}.MAP.SE.${cid}.bam" || { log "ERROR: Missing/invalid SE BAM for chunk ${cid}"; missing=1; }
  done < "$SE_IDLIST"
fi
[[ "$missing" -eq 0 ]] || die "Mapping incomplete: one or more expected chunk BAMs are missing/invalid"

###############################################################################
# OPTIONAL: UNMAPPED FASTQ EXPORT
###############################################################################
make_empty_gz() {
  local out="$1"
  if [[ $HAVE_PIGZ -eq 1 ]]; then
    printf '' | pigz -p "$THREADS" > "$out"
  else
    printf '' | gzip > "$out"
  fi
}

export_unmapped_fastq() {
  ckpt_clear unmapped_fastq
  log "STEP: unmapped FASTQ export"

  rm -f "$UNMAPPED_R1_FQ" "$UNMAPPED_R2_FQ" "$UNMAPPED_SE_FQ" 2>/dev/null || true
  make_empty_gz "$UNMAPPED_R1_FQ"
  make_empty_gz "$UNMAPPED_R2_FQ"
  make_empty_gz "$UNMAPPED_SE_FQ"

  local cid prefix ur1 ur2 use
  if [[ "$LIBTYPE" != "ancient-fast" && "$SEQ_MODE" != "SE" && -s "$PE_IDLIST" ]]; then
    while read -r cid; do
      [[ -z "$cid" ]] && continue
      prefix="$UNMAPPED/${SAMPLE}.PE.${cid}.unmapped"
      ur1="${prefix}.R1.fastq.gz"; ur2="${prefix}.R2.fastq.gz"; use="${prefix}.SE.fastq.gz"
      [[ -f "$ur1" && -f "$ur2" && -f "$use" ]] || die "Missing unmapped FASTQ chunk outputs for PE chunk ${cid}; rerun with --reset if needed."
      cat "$ur1" >> "$UNMAPPED_R1_FQ"
      cat "$ur2" >> "$UNMAPPED_R2_FQ"
      cat "$use" >> "$UNMAPPED_SE_FQ"
    done < "$PE_IDLIST"
  fi

  if [[ -s "$SE_IDLIST" ]]; then
    while read -r cid; do
      [[ -z "$cid" ]] && continue
      prefix="$UNMAPPED/${SAMPLE}.SE.${cid}.unmapped"
      ur1="${prefix}.R1.fastq.gz"; ur2="${prefix}.R2.fastq.gz"; use="${prefix}.SE.fastq.gz"
      [[ -f "$ur1" && -f "$ur2" && -f "$use" ]] || die "Missing unmapped FASTQ chunk outputs for SE chunk ${cid}; rerun with --reset if needed."
      cat "$ur1" >> "$UNMAPPED_R1_FQ"
      cat "$ur2" >> "$UNMAPPED_R2_FQ"
      cat "$use" >> "$UNMAPPED_SE_FQ"
    done < "$SE_IDLIST"
  fi

  ckpt_mark unmapped_fastq
}

if [[ "$EMIT_UNMAPPED_FASTQ" -eq 1 ]]; then
  if ckpt_ok_exists unmapped_fastq "$UNMAPPED_R1_FQ" "$UNMAPPED_R2_FQ" "$UNMAPPED_SE_FQ"; then
    log "STEP: unmapped FASTQ export (skipped; checkpoint present)"
  else
    export_unmapped_fastq
  fi
fi

###############################################################################
# STEP 5: MERGE (OVERWRITE SAFE)
###############################################################################
MERGED_BAM="$FINAL/${SAMPLE}.merged.bam"
log "STEP: merge"

bam_list=()

if [[ "$LIBTYPE" != "ancient-fast" && "$SEQ_MODE" != "SE" && -s "$PE_IDLIST" ]]; then
  while read -r cid; do
    [[ -z "$cid" ]] && continue
    bam_list+=( "$BAMS/${SAMPLE}.MAP.PE.${cid}.bam" )
  done < "$PE_IDLIST"
fi
if [[ -s "$SE_IDLIST" ]]; then
  while read -r cid; do
    [[ -z "$cid" ]] && continue
    bam_list+=( "$BAMS/${SAMPLE}.MAP.SE.${cid}.bam" )
  done < "$SE_IDLIST"
fi

[[ ${#bam_list[@]} -gt 0 ]] || die "No BAMs produced by mapping"

if ckpt_ok merge "$MERGED_BAM" && bam_ok "$MERGED_BAM"; then
  log "  merge (skipped; checkpoint present)"
else
  ckpt_clear merge
  mkdir -p "$FINAL"
  "$SAMTOOLS" merge -f -@ "$SORT_THREADS" "$MERGED_BAM" "${bam_list[@]}"
  bam_ok_with_retries "$MERGED_BAM" 5 15 || die "Merge failed/invalid after retries: $MERGED_BAM"
  ckpt_mark merge
fi

###############################################################################
# STEP 6: DUPLICATE REMOVAL
###############################################################################
PE_INPUT_BAM="$FINAL/${SAMPLE}.pe.input.bam"
SE_INPUT_BAM="$FINAL/${SAMPLE}.se_like.input.bam"
PE_MARKDUP_BAM="$FINAL/${SAMPLE}.pe.markdup_flagged.bam"
PE_DEDUP_BAM="$FINAL/${SAMPLE}.pe.dedup.bam"
SE_DEDUP_BAM="$FINAL/${SAMPLE}.se_like.dedup.bam"
LEGACY_MARKDUP_BAM="$FINAL/${SAMPLE}.markdup_flagged.bam"
DEDUP_BAM="$FINAL/${SAMPLE}.dedup.bam"
DEDUP_META="$FINAL/${SAMPLE}.dedup.stats.tsv"
DEDUP_EXCL_FLAGS=2820
DEDUP_REPORT_BASE="$DEDUP_REPORTS/${SAMPLE}.se_like.coordsort"

count_fragments_in_bam() {
  local bam="$1"
  local include_flags="${2:-0}"
  local exclude_flags="${3:-$DEDUP_EXCL_FLAGS}"
  local view_args=( -F "$exclude_flags" )
  [[ "$include_flags" -eq 0 ]] || view_args+=( -f "$include_flags" )
  "$SAMTOOLS" view "${view_args[@]}" "$bam" \
    | awk 'BEGIN{c=0}
           function hasbit(x,b){return and(x,b)}
           {flag=$2; if(hasbit(flag,1)){if(hasbit(flag,64)) c++} else c++}
           END{print c}'
}

count_duplicate_fragments_in_bam() {
  local bam="$1"
  local include_flags="${2:-0}"
  local exclude_flags="${3:-$DEDUP_EXCL_FLAGS}"
  local view_args=( -F "$exclude_flags" )
  [[ "$include_flags" -eq 0 ]] || view_args+=( -f "$include_flags" )
  "$SAMTOOLS" view "${view_args[@]}" "$bam" \
    | awk 'BEGIN{c=0}
           function hasbit(x,b){return and(x,b)}
           {flag=$2;
            if(hasbit(flag,1024)){
              if(hasbit(flag,1)){if(hasbit(flag,64)) c++} else c++
            }
           }
           END{print c}'
}

dedup_checkpoint_ok() {
  ckpt_ok dedup "$DEDUP_BAM" "$DEDUP_META" || return 1
  bam_ok "$DEDUP_BAM" || return 1

  if [[ "$DEDUP_METHOD" == "samtools" ]]; then
    bam_ok "$LEGACY_MARKDUP_BAM"
    return
  fi

  local se_fragments
  se_fragments=$(awk -F '\t' '$1 == "se_like_fragments_all" {print $2}' "$DEDUP_META")
  [[ "$se_fragments" =~ ^[0-9]+$ ]] || return 1
  if [[ "$se_fragments" -gt 0 ]]; then
    file_nonempty "${DEDUP_REPORT_BASE}.log" || return 1
    file_exists "${DEDUP_REPORT_BASE}.hist" || return 1
    file_nonempty "${DEDUP_REPORT_BASE}.dedup.json" || return 1
    file_nonempty "${DEDUP_REPORT_BASE}.stdout.log" || return 1
  fi
}

if [[ "$DEDUP_METHOD" == "samtools" ]]; then
  log "STEP: dedup (samtools-only v0.1-compatible path on combined BAM)"
else
  log "STEP: dedup (hybrid: DeDup for SE-like; samtools markdup for unmerged PE)"
fi

if dedup_checkpoint_ok; then
  log "  dedup (skipped; checkpoint present; method=$DEDUP_METHOD)"
else
  ckpt_clear dedup
  rm -f "$PE_INPUT_BAM" "$SE_INPUT_BAM" "$PE_MARKDUP_BAM" "$LEGACY_MARKDUP_BAM" \
        "$PE_DEDUP_BAM" "$SE_DEDUP_BAM" "$DEDUP_BAM" "$DEDUP_META" \
        "$FINAL/${SAMPLE}.pe.namesort.bam" "$FINAL/${SAMPLE}.pe.fixmate.bam" \
        "$FINAL/${SAMPLE}.pe.coordsort.bam" "$FINAL/${SAMPLE}.se_like.coordsort.bam" \
        "$FINAL/${SAMPLE}.namesort.bam" "$FINAL/${SAMPLE}.fixmate.bam" \
        "$FINAL/${SAMPLE}.coordsort.bam"

  if [[ "$DEDUP_METHOD" == "samtools" ]]; then
    if [[ ( "$LIBTYPE" == "modern" || "$LIBTYPE" == "ancient-all" ) && "$SEQ_MODE" != "SE" ]]; then
      "$SAMTOOLS" sort -n -@ "$SORT_THREADS" -o "$FINAL/${SAMPLE}.namesort.bam" "$MERGED_BAM"
      "$SAMTOOLS" fixmate -m -@ "$SORT_THREADS" "$FINAL/${SAMPLE}.namesort.bam" "$FINAL/${SAMPLE}.fixmate.bam"
      "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$FINAL/${SAMPLE}.coordsort.bam" "$FINAL/${SAMPLE}.fixmate.bam"
      rm -f "$FINAL/${SAMPLE}.namesort.bam" "$FINAL/${SAMPLE}.fixmate.bam"
    else
      "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$FINAL/${SAMPLE}.coordsort.bam" "$MERGED_BAM"
    fi

    "$SAMTOOLS" markdup -s -@ "$SORT_THREADS" "$FINAL/${SAMPLE}.coordsort.bam" "$LEGACY_MARKDUP_BAM"
    rm -f "$FINAL/${SAMPLE}.coordsort.bam"
    bam_ok "$LEGACY_MARKDUP_BAM" || die "samtools markdup failed/invalid: $LEGACY_MARKDUP_BAM"

    "$SAMTOOLS" view -b -F 1028 "$LEGACY_MARKDUP_BAM" > "$TMP/${SAMPLE}.nodup.unsorted.bam"
    "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$DEDUP_BAM" "$TMP/${SAMPLE}.nodup.unsorted.bam"
    rm -f "$TMP/${SAMPLE}.nodup.unsorted.bam"
    bam_ok "$DEDUP_BAM" || die "Deduplicated BAM missing/invalid: $DEDUP_BAM"

    PE_FRAGMENTS_ALL=$(count_fragments_in_bam "$LEGACY_MARKDUP_BAM" 1 "$DEDUP_EXCL_FLAGS")
    PE_DUP_FRAGMENTS=$(count_duplicate_fragments_in_bam "$LEGACY_MARKDUP_BAM" 1 "$DEDUP_EXCL_FLAGS")
    PE_FRAGMENTS_UNIQUE=$(py_sub_nonneg "$PE_FRAGMENTS_ALL" "$PE_DUP_FRAGMENTS")
    # Excluding flag 0x1 in addition to the standard exclusions selects SE-like records.
    SE_FRAGMENTS_ALL=$(count_fragments_in_bam "$LEGACY_MARKDUP_BAM" 0 2821)
    SE_DUP_FRAGMENTS=$(count_duplicate_fragments_in_bam "$LEGACY_MARKDUP_BAM" 0 2821)
    SE_FRAGMENTS_UNIQUE=$(py_sub_nonneg "$SE_FRAGMENTS_ALL" "$SE_DUP_FRAGMENTS")
    PE_DEDUP_METHOD="samtools-markdup"
    SE_DEDUP_METHOD="samtools-markdup"
  else

    # BAM flag 0x1 distinguishes alignments originating from unmerged PE mapping.
    "$SAMTOOLS" view -@ "$SORT_THREADS" -b -f 1 -o "$PE_INPUT_BAM" "$MERGED_BAM"
    "$SAMTOOLS" view -@ "$SORT_THREADS" -b -F 1 -o "$SE_INPUT_BAM" "$MERGED_BAM"
    bam_ok "$PE_INPUT_BAM" || die "PE-class BAM missing/invalid: $PE_INPUT_BAM"
    bam_ok "$SE_INPUT_BAM" || die "SE-like-class BAM missing/invalid: $SE_INPUT_BAM"

    PE_RECORDS=$("$SAMTOOLS" view -c "$PE_INPUT_BAM")
    SE_RECORDS=$("$SAMTOOLS" view -c "$SE_INPUT_BAM")
    dedup_outputs=()

    PE_FRAGMENTS_ALL=0
    PE_FRAGMENTS_UNIQUE=0
    PE_DUP_FRAGMENTS=0
    if [[ "$PE_RECORDS" -gt 0 ]]; then
      log "  Unmerged PE: samtools fixmate + markdup ($PE_RECORDS alignment records)"
      "$SAMTOOLS" sort -n -@ "$SORT_THREADS" -o "$FINAL/${SAMPLE}.pe.namesort.bam" "$PE_INPUT_BAM"
      "$SAMTOOLS" fixmate -m -@ "$SORT_THREADS" "$FINAL/${SAMPLE}.pe.namesort.bam" "$FINAL/${SAMPLE}.pe.fixmate.bam"
      "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$FINAL/${SAMPLE}.pe.coordsort.bam" "$FINAL/${SAMPLE}.pe.fixmate.bam"
      "$SAMTOOLS" markdup -s -@ "$SORT_THREADS" "$FINAL/${SAMPLE}.pe.coordsort.bam" "$PE_MARKDUP_BAM"
      bam_ok "$PE_MARKDUP_BAM" || die "PE samtools markdup failed/invalid: $PE_MARKDUP_BAM"

      "$SAMTOOLS" view -@ "$SORT_THREADS" -b -F 1028 -o "$TMP/${SAMPLE}.pe.nodup.unsorted.bam" "$PE_MARKDUP_BAM"
      "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$PE_DEDUP_BAM" "$TMP/${SAMPLE}.pe.nodup.unsorted.bam"
      bam_ok "$PE_DEDUP_BAM" || die "PE deduplicated BAM missing/invalid: $PE_DEDUP_BAM"

      PE_FRAGMENTS_ALL=$(count_fragments_in_bam "$PE_MARKDUP_BAM")
      PE_DUP_FRAGMENTS=$(count_duplicate_fragments_in_bam "$PE_MARKDUP_BAM")
      PE_FRAGMENTS_UNIQUE=$(py_sub_nonneg "$PE_FRAGMENTS_ALL" "$PE_DUP_FRAGMENTS")
      dedup_outputs+=( "$PE_DEDUP_BAM" )
    else
      log "  Unmerged PE: no alignments; samtools markdup not required"
    fi

    SE_FRAGMENTS_ALL=0
    SE_FRAGMENTS_UNIQUE=0
    SE_DUP_FRAGMENTS=0
    if [[ "$SE_RECORDS" -gt 0 ]]; then
      log "  SE-like: DeDup --merged ($SE_RECORDS alignment records)"
      "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$FINAL/${SAMPLE}.se_like.coordsort.bam" "$SE_INPUT_BAM"

      dedup_generated_bam="${DEDUP_REPORT_BASE}_rmdup.bam"
      rm -f "$dedup_generated_bam" \
            "${DEDUP_REPORT_BASE}.log" \
            "${DEDUP_REPORT_BASE}.hist" \
            "${DEDUP_REPORT_BASE}.dedup.json" \
            "${DEDUP_REPORT_BASE}.stdout.log"

      "$DEDUP" --merged \
        --input "$FINAL/${SAMPLE}.se_like.coordsort.bam" \
        --output "$DEDUP_REPORTS" \
        > "${DEDUP_REPORT_BASE}.stdout.log"
      bam_ok "$dedup_generated_bam" || die "DeDup failed or produced an invalid SE-like BAM: $dedup_generated_bam"
      file_nonempty "${DEDUP_REPORT_BASE}.log" || die "DeDup log missing/empty: ${DEDUP_REPORT_BASE}.log"
      file_exists "${DEDUP_REPORT_BASE}.hist" || die "DeDup histogram missing: ${DEDUP_REPORT_BASE}.hist"
      file_nonempty "${DEDUP_REPORT_BASE}.dedup.json" || die "DeDup JSON report missing/empty: ${DEDUP_REPORT_BASE}.dedup.json"
      file_nonempty "${DEDUP_REPORT_BASE}.stdout.log" || die "DeDup stdout report missing/empty: ${DEDUP_REPORT_BASE}.stdout.log"

      "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$SE_DEDUP_BAM" "$dedup_generated_bam"
      rm -f "$dedup_generated_bam"
      bam_ok "$SE_DEDUP_BAM" || die "SE-like deduplicated BAM missing/invalid: $SE_DEDUP_BAM"

      SE_FRAGMENTS_ALL=$(count_fragments_in_bam "$SE_INPUT_BAM")
      SE_FRAGMENTS_UNIQUE=$(count_fragments_in_bam "$SE_DEDUP_BAM")
      [[ "$SE_FRAGMENTS_UNIQUE" -le "$SE_FRAGMENTS_ALL" ]] || \
        die "DeDup output has more SE-like fragments than its input ($SE_FRAGMENTS_UNIQUE > $SE_FRAGMENTS_ALL)"
      SE_DUP_FRAGMENTS=$(py_sub_nonneg "$SE_FRAGMENTS_ALL" "$SE_FRAGMENTS_UNIQUE")
      dedup_outputs+=( "$SE_DEDUP_BAM" )
    else
      log "  SE-like: no alignments; DeDup not required"
    fi

    if [[ ${#dedup_outputs[@]} -eq 0 ]]; then
      log "  No mapped alignments remain; producing an empty final BAM"
      cp "$MERGED_BAM" "$DEDUP_BAM"
    elif [[ ${#dedup_outputs[@]} -eq 1 ]]; then
      cp "${dedup_outputs[0]}" "$DEDUP_BAM"
    else
      "$SAMTOOLS" merge -f -@ "$SORT_THREADS" "$DEDUP_BAM" "${dedup_outputs[@]}"
    fi
    bam_ok "$DEDUP_BAM" || die "Combined deduplicated BAM missing/invalid: $DEDUP_BAM"
    PE_DEDUP_METHOD="samtools-markdup"
    SE_DEDUP_METHOD="DeDup"
  fi

  TOTAL_FRAGMENTS_ALL=$(py_add "$PE_FRAGMENTS_ALL" "$SE_FRAGMENTS_ALL")
  TOTAL_FRAGMENTS_UNIQUE=$(py_add "$PE_FRAGMENTS_UNIQUE" "$SE_FRAGMENTS_UNIQUE")
  TOTAL_DUP_FRAGMENTS=$(py_add "$PE_DUP_FRAGMENTS" "$SE_DUP_FRAGMENTS")
  {
    printf 'metric\tvalue\n'
    printf 'pe_fragments_all\t%s\n' "$PE_FRAGMENTS_ALL"
    printf 'pe_fragments_unique\t%s\n' "$PE_FRAGMENTS_UNIQUE"
    printf 'pe_duplicate_fragments\t%s\n' "$PE_DUP_FRAGMENTS"
    printf 'pe_dedup_method\t%s\n' "$PE_DEDUP_METHOD"
    printf 'se_like_fragments_all\t%s\n' "$SE_FRAGMENTS_ALL"
    printf 'se_like_fragments_unique\t%s\n' "$SE_FRAGMENTS_UNIQUE"
    printf 'se_like_duplicate_fragments\t%s\n' "$SE_DUP_FRAGMENTS"
    printf 'se_like_dedup_method\t%s\n' "$SE_DEDUP_METHOD"
    printf 'total_fragments_all\t%s\n' "$TOTAL_FRAGMENTS_ALL"
    printf 'total_fragments_unique\t%s\n' "$TOTAL_FRAGMENTS_UNIQUE"
    printf 'total_duplicate_fragments\t%s\n' "$TOTAL_DUP_FRAGMENTS"
  } > "$DEDUP_META"

  rm -f "$PE_INPUT_BAM" "$SE_INPUT_BAM" \
        "$FINAL/${SAMPLE}.pe.namesort.bam" "$FINAL/${SAMPLE}.pe.fixmate.bam" \
        "$FINAL/${SAMPLE}.pe.coordsort.bam" "$FINAL/${SAMPLE}.se_like.coordsort.bam" \
        "$TMP/${SAMPLE}.pe.nodup.unsorted.bam" "$TMP/${SAMPLE}.nodup.unsorted.bam" \
        "$FINAL/${SAMPLE}.namesort.bam" "$FINAL/${SAMPLE}.fixmate.bam" \
        "$FINAL/${SAMPLE}.coordsort.bam"
  ckpt_mark dedup
fi

###############################################################################
# STEP 7: FINAL SORT + RG (RG OUTPUT GUARANTEED SORTED) + INDEX
###############################################################################
FINAL_BAM="$FINAL/${SAMPLE}.final.bam"

log "STEP: final sort+index + read groups"

if ckpt_ok rg "$OUT_BAM" "${OUT_BAM}.bai" && bam_ok "$OUT_BAM"; then
  log "  final+RG (skipped; checkpoint present)"
else
  ckpt_clear rg

  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$FINAL_BAM" "$DEDUP_BAM"
  "$SAMTOOLS" index "$FINAL_BAM"
  bam_ok "$FINAL_BAM" || die "Final BAM invalid: $FINAL_BAM"

  RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina\tLB:${SAMPLE}"
  tmp_rg="$TMP/${SAMPLE}.rg.tmp.bam"
  "$SAMTOOLS" addreplacerg -r "$RG" -o "$tmp_rg" "$FINAL_BAM"

  "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$OUT_BAM" "$tmp_rg"
  rm -f "$tmp_rg"

  "$SAMTOOLS" index "$OUT_BAM"
  bam_ok "$OUT_BAM" || die "RG BAM invalid: $OUT_BAM"

  ckpt_mark rg
fi

###############################################################################
# OPTIONAL: mapDamage (parallelized by whole-contig BAM splitting)
###############################################################################
# mapDamage2 itself is single-threaded and, unmodified, scans the whole BAM in
# one pass. Its scan produces plain per-contig count tables (misincorporation,
# base composition, fragment-length) with no cross-read state, and (run with
# --no-stats, as PlainMap already does) its plotting stage never rescans the
# BAM -- it only reads those table files back in. That means the scan can be
# split by whole contig (never by sub-contig region, which risks double-
# counting a read that spans a split boundary), run in parallel, and the
# resulting tables merged by summing matching rows before a final, unmodified
# mapDamage --plot-only call regenerates the usual output. mapDamage's own
# code is never modified; only how PlainMap invokes it changes.
MAPDAMAGE_OUT="$OUT/${SAMPLE}_mapdamage"

mapdamage_tables_ok() {
  local dir="$1" table
  for table in misincorporation.txt dnacomp.txt lgdistribution.txt; do
    file_nonempty "$dir/$table" || return 1
  done
}

mapdamage_has_length_data() {
  local table="$1"
  awk '
    /^#/ || NF == 0 {next}
    !header_seen {header_seen=1; next}
    {found=1; exit}
    END {exit(found ? 0 : 1)}
  ' "$table"
}

mapdamage_plots_ok() {
  local dir="$1"
  file_nonempty "$dir/Fragmisincorporation_plot.pdf" || return 1
  file_nonempty "$dir/5pCtoT_freq.txt" || return 1
  file_nonempty "$dir/3pGtoA_freq.txt" || return 1
  if mapdamage_has_length_data "$dir/lgdistribution.txt"; then
    file_nonempty "$dir/Length_plot.pdf" || return 1
  fi
}

mapdamage_results_ok() {
  local dir="$1"
  mapdamage_tables_ok "$dir" || return 1
  command -v Rscript >/dev/null 2>&1 || return 0
  mapdamage_plots_ok "$dir"
}

validate_mapdamage_results() {
  local dir="$1"
  mapdamage_tables_ok "$dir" || die "mapDamage result tables are missing or empty: $dir"
  if command -v Rscript >/dev/null 2>&1; then
    mapdamage_plots_ok "$dir" || die "mapDamage plots are missing or empty despite Rscript being available: $dir"
  fi
}

run_mapdamage_single() {
  rm -rf "$MAPDAMAGE_OUT"
  "$MAPDAMAGE" -i "$OUT_BAM" --merge-reference-sequences --no-stats -r "$REF" -d "$MAPDAMAGE_OUT"
  validate_mapdamage_results "$MAPDAMAGE_OUT"
}

mapdamage_chunk_worker() {
  local idx="$1" regions_file="$2"
  local step="mapdamage_chunk_${idx}"
  local chunk_bam="$MDCHUNKS/${SAMPLE}.chunk${idx}.bam"
  local chunk_dir="$MDCHUNKS/${SAMPLE}.chunk${idx}.results"

  if ckpt_ok "$step" "$chunk_dir/misincorporation.txt" "$chunk_dir/dnacomp.txt" "$chunk_dir/lgdistribution.txt"; then
    log "  mapDamage chunk ${idx}: already complete (checkpoint present) -> skipping"
    return 0
  fi
  ckpt_clear "$step"

  local nregions
  nregions=$(awk 'END {print NR + 0}' "$regions_file")
  [[ "$nregions" -gt 0 ]] || die "mapDamage chunk ${idx}: empty region list"

  log "  mapDamage chunk ${idx}: extracting $nregions contig(s)"
  "$SAMTOOLS" view -b -M -L "$regions_file" -o "${chunk_bam}.tmp" "$OUT_BAM"
  bam_ok "${chunk_bam}.tmp" || die "mapDamage chunk ${idx}: extracted BAM invalid"
  mv -f "${chunk_bam}.tmp" "$chunk_bam"

  rm -rf "$chunk_dir"
  log "  mapDamage chunk ${idx}: scanning"
  "$MAPDAMAGE" -i "$chunk_bam" --merge-reference-sequences --no-stats --no-plot -r "$REF" -d "$chunk_dir"

  file_nonempty "$chunk_dir/misincorporation.txt" || die "mapDamage chunk ${idx}: missing misincorporation.txt"
  file_nonempty "$chunk_dir/dnacomp.txt"          || die "mapDamage chunk ${idx}: missing dnacomp.txt"
  file_nonempty "$chunk_dir/lgdistribution.txt"   || die "mapDamage chunk ${idx}: missing lgdistribution.txt"

  if [[ $DELETE_AS_YOU_GO -eq 1 ]]; then
    rm -f "$chunk_bam"
  fi

  ckpt_mark "$step" || true
}

# Sums matching rows across chunk result tables (each chunk already ran with
# --merge-reference-sequences, so misincorporation.txt/dnacomp.txt rows are all
# keyed "Chr=*"; lgdistribution.txt is keyed by strand+length). Column layout
# matches what mapdamage/tables.py writes.
merge_mapdamage_chunk_tables() {
  local outdir="$1"; shift
  mkdir -p "$outdir"
  "$PYTHON" - "$outdir" "$@" <<'PY'
import sys, os

def read_data(path):
    header, rows = None, []
    with open(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            if header is None:
                header = line
                continue
            rows.append((lineno, line.split("\t")))
    return header, rows

def merge(files, key_idxs, sum_start):
    header, sums, order, ncols = None, {}, [], None
    for f in files:
        h, rows = read_data(f)
        if h is None:
            raise ValueError("missing table header in %s" % f)
        if header is None:
            header = h
            ncols = len(h.split("\t"))
        elif h != header:
            raise ValueError("table header differs in %s" % f)
        for lineno, cols in rows:
            if len(cols) != ncols:
                raise ValueError("expected %d columns in %s:%d, found %d" %
                                 (ncols, f, lineno, len(cols)))
            key = tuple(cols[i] for i in key_idxs)
            if key not in sums:
                sums[key] = [0] * (ncols - sum_start)
                order.append(key)
            for j in range(sum_start, ncols):
                try:
                    sums[key][j - sum_start] += int(cols[j])
                except ValueError as exc:
                    raise ValueError("non-integer count in %s:%d" % (f, lineno)) from exc
    return header, sums, order, ncols

def sort_key(k):
    return [(0, int(v)) if v.lstrip("-").isdigit() else (1, v) for v in k]

def write_out(path, comments, header, key_idxs, const0, ncols, sums, order):
    with open(path, "w") as out:
        for c in comments:
            out.write(c + "\n")
        out.write(header + "\n")
        for key in sorted(order, key=sort_key):
            row = [""] * ncols
            if const0 is not None:
                row[0] = const0
            for pos, idx in enumerate(key_idxs):
                row[idx] = key[pos]
            sum_start = max(key_idxs) + 1 if const0 is None else 4
            for j in range(sum_start, ncols):
                row[j] = str(sums[key][j - sum_start])
            out.write("\t".join(row) + "\n")

outdir, chunk_dirs = sys.argv[1], sys.argv[2:]

for fname in ("misincorporation.txt", "dnacomp.txt"):
    files = [os.path.join(d, fname) for d in chunk_dirs]
    missing = [f for f in files if not os.path.isfile(f)]
    if missing:
        sys.exit("missing chunk table(s): %s" % ", ".join(missing))
    header, sums, order, ncols = merge(files, key_idxs=[1, 2, 3], sum_start=4)
    write_out(os.path.join(outdir, fname),
              ["# table produced by PlainMap (merged from %d parallel mapDamage chunks)" % len(files),
               "# Chr: reference from sam/bam header, End: from which termini of DNA sequences, Std: strand of reads"],
              header, [1, 2, 3], "*", ncols, sums, order)

fname = "lgdistribution.txt"
files = [os.path.join(d, fname) for d in chunk_dirs]
missing = [f for f in files if not os.path.isfile(f)]
if missing:
    sys.exit("missing chunk table(s): %s" % ", ".join(missing))
header, sums, order, ncols = merge(files, key_idxs=[0, 1], sum_start=2)
write_out(os.path.join(outdir, fname),
          ["# table produced by PlainMap (merged from %d parallel mapDamage chunks)" % len(files)],
          header, [0, 1], None, ncols, sums, order)
PY
}

run_mapdamage_parallel() {
  mkdir -p "$MDCHUNKS"

  # Avoid a race between workers creating a missing index, but do not rewrite an
  # existing index: references are commonly shared from read-only directories.
  if ! file_nonempty "${REF}.fai"; then
    log "  mapDamage: creating missing FASTA index"
    "$SAMTOOLS" faidx "$REF"
  fi

  local contig_stats="$MDCHUNKS/${SAMPLE}.contig_stats.tsv"
  "$SAMTOOLS" idxstats "$OUT_BAM" \
    | awk -F'\t' '$1!="*" && $3>0 {print $1"\t"$2"\t"$3}' > "$contig_stats"

  local ncontigs
  ncontigs=$(awk 'END {print NR + 0}' "$contig_stats")

  if [[ "$MAPDAMAGE_SERIAL" -eq 1 || "$THREADS" -le 1 || "$ncontigs" -le 1 ]]; then
    log "  mapDamage: serial mode (threads=$THREADS, contigs-with-reads=$ncontigs, serial-flag=$MAPDAMAGE_SERIAL)"
    run_mapdamage_single
    return
  fi

  local n=$(( THREADS < ncontigs ? THREADS : ncontigs ))
  log "  mapDamage: splitting $ncontigs contig(s) into $n parallel chunk(s)"

  local bin_prefix="$MDCHUNKS/${SAMPLE}.bin"
  rm -f "${bin_prefix}."*.bed

  # This reads the statistics from a file, not piped stdin: "python3 -"
  # already uses stdin to receive this heredoc as the script's own source, so
  # piping data into the same invocation would silently starve sys.stdin.
  "$PYTHON" - "$n" "$bin_prefix" "$contig_stats" <<'PY'
import sys
n, prefix, contig_stats = int(sys.argv[1]), sys.argv[2], sys.argv[3]
contigs = []
with open(contig_stats) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        name, length, reads = line.split("\t")
        contigs.append((name, int(length), int(reads)))
contigs.sort(key=lambda x: x[2], reverse=True)
bins = [[] for _ in range(n)]
loads = [0] * n
for name, length, reads in contigs:
    i = min(range(n), key=lambda j: loads[j])
    bins[i].append((name, length))
    loads[i] += reads
for i, entries in enumerate(bins):
    with open("%s.%d.bed" % (prefix, i), "w") as fh:
        for name, length in entries:
            fh.write("%s\t0\t%d\n" % (name, length))
PY

  local pids=() idx
  for (( idx=0; idx<n; idx++ )); do
    local regions="${bin_prefix}.${idx}.bed"
    [[ -s "$regions" ]] || die "mapDamage: expected region bin file missing/empty: $regions"
    ( mapdamage_chunk_worker "$idx" "$regions" ) &
    pids+=( "$!" )
  done

  local fail=0 p
  for p in "${pids[@]}"; do
    wait "$p" || fail=1
  done
  [[ "$fail" -eq 0 ]] || die "mapDamage: one or more parallel chunks failed (see log above for which)"

  local chunk_dirs=()
  for (( idx=0; idx<n; idx++ )); do
    chunk_dirs+=( "$MDCHUNKS/${SAMPLE}.chunk${idx}.results" )
  done

  log "  mapDamage: merging $n chunk result table(s)"
  rm -rf "$MAPDAMAGE_OUT"
  merge_mapdamage_chunk_tables "$MAPDAMAGE_OUT" "${chunk_dirs[@]}"

  if command -v Rscript >/dev/null 2>&1; then
    log "  mapDamage: generating final plots from merged tables"
    "$MAPDAMAGE" --plot-only -d "$MAPDAMAGE_OUT" --title "${SAMPLE}.final.RG"
  else
    log "WARNING: Rscript not found; mapDamage tables were generated but plotting was skipped"
  fi
  validate_mapdamage_results "$MAPDAMAGE_OUT"
}

if [[ "$RUN_MAPDAMAGE" -eq 1 ]]; then
  if [[ $RESUME -eq 1 && -f "$(ckpt_file mapdamage)" ]] && mapdamage_results_ok "$MAPDAMAGE_OUT"; then
    log "STEP: mapDamage (skipped; checkpoint present)"
  else
    if [[ $RESUME -eq 1 && -f "$(ckpt_file mapdamage)" ]]; then
      log "STEP: mapDamage checkpoint is present but outputs are incomplete; rerunning"
    fi
    log "STEP: mapDamage"
    ckpt_clear mapdamage
    run_mapdamage_parallel
    ckpt_mark mapdamage
  fi
fi

###############################################################################
# STEP 8: COVERAGE/DEPTH
###############################################################################
log "STEP: coverage statistics"
if ckpt_ok coverage "$COV_TSV"; then
  log "  coverage (skipped; checkpoint present)"
else
  ckpt_clear coverage
  "$SAMTOOLS" coverage -q "$MAPQ" "$OUT_BAM" > "$COV_TSV"
  file_nonempty "$COV_TSV" || die "Coverage TSV missing/empty: $COV_TSV"
  ckpt_mark coverage
fi

AVG_DEPTH=$(
  awk 'BEGIN{wsum=0;len=0}
       NR==1{next}
       {ref_len=($3-$2+1); wsum += $7*ref_len; len += ref_len}
       END{if(len>0) printf "%.6f", wsum/len; else print 0}' \
    "$COV_TSV"
)

PCT_COVERED=$(
  awk 'BEGIN{wsum=0;len=0}
       NR==1{next}
       {ref_len=($3-$2+1); wsum += $6*ref_len; len += ref_len}
       END{if(len>0) printf "%.6f", wsum/len; else print 0}' \
    "$COV_TSV"
)

###############################################################################
# STEP 9: SUMMARY STATS
###############################################################################
log "STEP: summary statistics"

RAW_R1_READS=$(awk 'NR==2{print $1}' "$RAW_META")
RAW_R2_READS=$(awk 'NR==2{print $2}' "$RAW_META")
RAW_FRAGMENTS=$(awk 'NR==2{print $3}' "$RAW_META")
TRIMMED_FRAGMENTS=$(awk 'NR==2{print $4}' "$RAW_META")
MERGED_READS=$(awk 'NR==2{print $5}' "$RAW_META")
UNPAIRED_READS=$(awk 'NR==2{print $6}' "$RAW_META")

# Exclude: unmapped(4) + secondary(256) + QC-fail(512) + supplementary(2048) = 2820
EXCL_FLAGS=2820

MAPPED_READS_ALL=$("$SAMTOOLS" view -c -F "$EXCL_FLAGS" "$MERGED_BAM")
MAPPED_READS_UNIQUE=$("$SAMTOOLS" view -c -F "$EXCL_FLAGS" "$OUT_BAM")

dedup_stat() {
  local metric="$1"
  awk -F '\t' -v metric="$metric" '$1 == metric {print $2; found=1} END{if(!found) exit 1}' "$DEDUP_META"
}

PE_FRAGMENTS_ALL=$(dedup_stat pe_fragments_all)
PE_FRAGMENTS_UNIQUE=$(dedup_stat pe_fragments_unique)
PE_DUP_FRAGMENTS=$(dedup_stat pe_duplicate_fragments)
PE_DEDUP_METHOD=$(dedup_stat pe_dedup_method)
SE_FRAGMENTS_ALL=$(dedup_stat se_like_fragments_all)
SE_FRAGMENTS_UNIQUE=$(dedup_stat se_like_fragments_unique)
SE_DUP_FRAGMENTS=$(dedup_stat se_like_duplicate_fragments)
SE_DEDUP_METHOD=$(dedup_stat se_like_dedup_method)
MAPPED_FRAGMENTS_ALL=$(dedup_stat total_fragments_all)
MAPPED_FRAGMENTS_UNIQUE=$(dedup_stat total_fragments_unique)
DUP_FRAGMENTS=$(dedup_stat total_duplicate_fragments)

DUPRATE=$("$PYTHON" - "$DUP_FRAGMENTS" "$MAPPED_FRAGMENTS_ALL" <<'PY'
import sys
d=float(sys.argv[1]); a=float(sys.argv[2])
print("0.000000" if a<=0 else f"{(d/a):.6f}")
PY
)

PE_DUPRATE=$("$PYTHON" - "$PE_DUP_FRAGMENTS" "$PE_FRAGMENTS_ALL" <<'PY'
import sys
d=float(sys.argv[1]); a=float(sys.argv[2])
print("0.000000" if a<=0 else f"{(d/a):.6f}")
PY
)

SE_DUPRATE=$("$PYTHON" - "$SE_DUP_FRAGMENTS" "$SE_FRAGMENTS_ALL" <<'PY'
import sys
d=float(sys.argv[1]); a=float(sys.argv[2])
print("0.000000" if a<=0 else f"{(d/a):.6f}")
PY
)

AVG_READLEN=$(
  "$SAMTOOLS" view -F "$EXCL_FLAGS" "$OUT_BAM" \
  | awk '{l=length($10); if(l>0){sum+=l;n++}} END{ if(n>0) printf "%.2f", sum/n; else printf "0.00"}'
)

# Precise mapped_bp computed on the final BAM using samtools view.
MAPPED_BP=$(
  "$SAMTOOLS" view -F "$EXCL_FLAGS" "$OUT_BAM" \
  | awk 'BEGIN{sum=0} {sum+=length($10)} END{printf "%.0f", sum}'
)

ENDOG=$("$PYTHON" - "$MAPPED_FRAGMENTS_UNIQUE" "$RAW_FRAGMENTS" <<'PY'
import sys
n=float(sys.argv[1]); d=float(sys.argv[2])
print("0" if d==0 else f"{(n/d):.6f}")
PY
)

{
  printf "sample\tlibrary_type\tseq_mode\tpilot_fragments\tmax_reads_per_chunk\tmapq\traw_R1_reads\traw_R2_reads\traw_fragments\ttrimmed_fragments\tmerged_reads\tunpaired_reads\tmapped_reads_all\tmapped_reads_unique\tmapped_fragments_all\tmapped_fragments_unique\tendog\tduprate\tavg_readlen\tmapped_bp\tavg_depth\tpct_covered\tduplicate_fragments\tpe_dedup_method\tpe_fragments_all\tpe_fragments_unique\tpe_duplicate_fragments\tpe_duprate\tse_like_dedup_method\tse_like_fragments_all\tse_like_fragments_unique\tse_like_duplicate_fragments\tse_like_duprate\tdedup_method\n"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$SAMPLE" "$LIBTYPE" "$SEQ_MODE" "$PILOT_FRAGMENTS" "$MAX_READS_PER_CHUNK" "$MAPQ" \
    "$RAW_R1_READS" "$RAW_R2_READS" "$RAW_FRAGMENTS" \
    "$TRIMMED_FRAGMENTS" "$MERGED_READS" "$UNPAIRED_READS" \
    "$MAPPED_READS_ALL" "$MAPPED_READS_UNIQUE" \
    "$MAPPED_FRAGMENTS_ALL" "$MAPPED_FRAGMENTS_UNIQUE" \
    "$ENDOG" "$DUPRATE" \
    "$AVG_READLEN" "$MAPPED_BP" \
    "$AVG_DEPTH" "$PCT_COVERED" "$DUP_FRAGMENTS" \
    "$PE_DEDUP_METHOD" "$PE_FRAGMENTS_ALL" "$PE_FRAGMENTS_UNIQUE" "$PE_DUP_FRAGMENTS" "$PE_DUPRATE" \
    "$SE_DEDUP_METHOD" "$SE_FRAGMENTS_ALL" "$SE_FRAGMENTS_UNIQUE" "$SE_DUP_FRAGMENTS" "$SE_DUPRATE" \
    "$DEDUP_METHOD"
} > "$STATS"

log "STATS written: $STATS"
log "Coverage written: $COV_TSV"

###############################################################################
# CLEANUP
###############################################################################
if [[ $KEEP_INTERMEDIATE -eq 0 ]]; then
  log "Cleaning intermediates (work dir only)"
  rm -rf "$WORK"
  [[ -n "$TMPDIR_USER" ]] && rm -rf "$TMPDIR_USER/plainmap_${SAMPLE}" || true
fi

PIPE_T1=$(date +%s)
log "STATUS: SUCCESS"
log "Total wall time: $((PIPE_T1 - PIPE_T0)) seconds"
log "Final BAM: $OUT_BAM"
