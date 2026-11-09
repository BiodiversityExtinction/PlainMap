# Manuscript evaluation

This directory records the manifests, commands, parameters, and results for the
original evaluation and the later scale, chunking, and AdaptClean tests. The
original [`evaluation_results.tsv`](evaluation_results.tsv) and
[`restart_checks.tsv`](restart_checks.tsv) are committed. Later run outputs
under `evaluation_runs/` are ignored by Git; the large sequence files,
reference genomes, BAMs, and newer result tables are not available in this
repository. The scripts and instructions below document how to regenerate
those results with the required input data.

The original manuscript table used PlainMap v0.1.0 and its samtools-only
duplicate-removal policy. Check out the immutable `v0.1.0` tag for the exact
PlainMap script version. The later hybrid-policy Z3 and cave hyena
tests used a frozen script snapshot from November 5, 2026 (SHA-256
`2723d9eed4ac2d39fe33af11c0cd2f3046362b2963f5632b1d2662daf91fedb6`).
Its executable logic matches v0.3.0; only the version label and comments
differ. It uses hybrid DeDup/samtools duplicate removal by default, so its
duplicate statistics cannot be directly compared with v0.1.0. To retain the
v0.1 mapping and duplicate-removal methodology with v0.3.0, specify
`--gap-opens 1 --dedup-method samtools`.

## Evaluation overview

| Evaluation | Data and question | Code and public materials |
|---|---|---|
| Original mapping modes and restart | Three million pilot fragments from a heterogeneous Japanese sea lion (Z3) manifest containing SE and PE FASTQs with SRA- and ENA-style headers; compare `modern`, `ancient-fast`, and `ancient-all`, then interrupt and resume a chunked `modern` run. | v0.1.0; [manifest](Z3_manifest_test.txt), [mapping-mode runner](run_mapping_mode.sh), [restart runner](run_restart_chunk_evaluation.sh), [reported results](evaluation_results.tsv), and [restart checks](restart_checks.tsv). |
| Larger Z3 workload | 27,943,970 input fragments across the same three manifest units; compare unchunked processing with one-million-fragment chunks under `modern` mapping. | v0.1.0 and, separately, the v0.3.0-equivalent snapshot using hybrid duplicate removal; [input preparation](prepare_z3_10m_inputs.sh) and [benchmark runner](run_large_z3_benchmark.sh). Results from the two duplicate policies are not directly interchangeable. |
| Chunking diagnostics | Compare five- and ten-million-fragment chunks with unchunked runs; repeat unchunked runs and inspect stage checksums to locate count differences. A controlled fastp audit compares automatic with fixed adapters. | Snapshot for the larger-run diagnostics; v0.1.0 for the [fastp audit](run_fastp_chunk_audit.sh). See the [chunk-size sweep](run_z3_chunk_size_sweep.sh), [mapping diagnostic](run_z3_mapping_diagnostic.sh), and [repeat runner](run_z3_unchunked_repeat.sh). |
| Short-read ancient sample | Fifty million pilot fragments from a single-end cave hyena library mapped in `ancient-fast` mode to striped hyena, with AdaptClean off versus on after fastp. | Snapshot; [runner](run_cavehyena_50m_adaptclean.sh) and [summary script](summarize_cavehyena_50m.sh). AdaptClean is not included in the Conda environment. |

These tests check input handling, execution, chunking, and restart behavior.
Mapped-fragment counts do not establish alignment accuracy or biological
endogeneity. Chunking can change fastp's automatic adapter inference, so
small differences in downstream counts are possible; the controlled
fixed-adapter audit is described under Results below.

## Dataset

The heterogeneous Z3 manifest contains:

- one single-end unit: `SRR25607782`
- one SRA-style paired-end unit whose mate files have identical read names:
  `SRR25607793`
- one paired-end unit with explicit `/1` and `/2` header suffixes: `SRR25607804`

The evaluated file order and names are stored in
[`Z3_manifest_test.txt`](Z3_manifest_test.txt). The `SRR25607793_1M_1.fastq.gz` and
`SRR25607793_1M_2.fastq.gz` files contain the first 1,000,000 mate pairs. Using complete
mate files produces the same evaluated subset because PlainMap applies
`--pilot-fragments 1000000` independently to each unit, provided record order is
preserved.

Reads were mapped to the California sea lion (*Zalophus californianus*) mZalCal1.pri.v2
assembly (NCBI RefSeq accession `GCF_009762305.2`).

## Requirements

The reported runs used:

| Software | Version |
|---|---|
| PlainMap | `0.1.0` |
| fastp | `1.0.1` |
| BWA | `0.7.17-r1188` |
| SAMtools/HTSlib | `1.13` / `1.13` |
| Python | `3.13.5` |
| Bash | `5.3.9` |

Set the following variables before running the scripts:

```bash
export PLAINMAP=/path/to/PlainMap/plainmap.sh
export DATA_DIR=/path/to/Z3_fastq_files
export REF=/path/to/GCF_009762305.2_mZalCal1.pri.v2_genomic.fna

# Optional command overrides; PATH commands are used by default.
export FASTP=fastp
export BWA=bwa
export SAMTOOLS=samtools
export PYTHON=python3
export THREADS=16
```

`DATA_DIR` must contain the five files named in `Z3_manifest_test.txt`. The scripts create
an absolute, resolved manifest in their result directory; local institutional paths are
therefore not embedded in this repository.

## Mapping-mode evaluation

Run directly:

```bash
for mode in modern ancient-fast ancient-all; do
    manuscript_evaluation/run_mapping_mode.sh "$mode"
done
```

The output root defaults to `manuscript_evaluation/evaluation_runs/mapping_modes/` and can
be changed with `EVAL_ROOT`.

Equivalent Slurm submissions used 16 CPUs and 48 GB of memory:

```bash
mkdir -p logs
for mode in modern ancient-fast ancient-all; do
    sbatch \
        --cpus-per-task=16 \
        --mem-per-cpu=3G \
        --time=03:00:00 \
        --job-name="PM_${mode}" \
        --output="logs/PM_${mode}-%j.out" \
        -- manuscript_evaluation/run_mapping_mode.sh "$mode"
done
```

Each mode used the following PlainMap settings:

```text
--threads 16
--pilot-fragments 1000000
--library-type MODE
--keep-intermediate
```

All other analysis parameters retained PlainMap defaults: minimum read length 30 bp,
MAPQ 20, `bwa aln -n 0.01 -o 1` for ancient modes, automatic adapter detection, resume enabled,
delete-as-you-go enabled, mapDamage disabled, and unmapped FASTQ export disabled.

## Chunking and interrupted restart

The complete operational evaluation is implemented in
[`run_restart_chunk_evaluation.sh`](run_restart_chunk_evaluation.sh). It performs:

1. An unchunked `modern` baseline.
2. An uninterrupted `modern` run with `--max-reads-per-chunk 100000`.
3. A second chunked run that is terminated with `SIGTERM` after at least one but not all
   mapping chunks have completed.
4. An identical restart of the interrupted command.
5. BAM validation and comparison of the uninterrupted and resumed alignment-record
   checksums and statistics.

Run directly:

```bash
manuscript_evaluation/run_restart_chunk_evaluation.sh
```

Or submit to Slurm:

```bash
sbatch \
    --cpus-per-task=16 \
    --mem-per-cpu=3G \
    --time=03:00:00 \
    --job-name=PM_restart \
    --output=logs/PM_restart-%j.out \
    -- manuscript_evaluation/run_restart_chunk_evaluation.sh
```

The output root defaults to `manuscript_evaluation/evaluation_runs/restart_chunk/` and can
be changed with `EVAL_ROOT`. Forced pre-trimming chunking temporarily required
approximately 12 GB for this dataset, so the output filesystem should have sufficient
free space.

## Results

- [`evaluation_results.tsv`](evaluation_results.tsv) contains the values reported in the
  combined manuscript table.
- [`restart_checks.tsv`](restart_checks.tsv) records the interruption and checksum checks.

The mapping modes intentionally process different read categories and should not be
interpreted as competing estimates of alignment accuracy. Adaptive chunking is a
memory-safety option: in this evaluation it reduced peak resident memory but increased I/O
and produced a small change in fastp output because adapter inference occurred separately
for each pre-trimming chunk.

In a separate controlled fastp audit, automatic adapter detection produced
1,124,015 trimmed fragments without chunking and 1,123,944 with chunking.
Using the same fixed adapter sequences in both conditions produced 1,124,015
trimmed fragments in each, supporting adapter inference as the source of that
specific count difference. This does not establish that every downstream
alignment difference has the same cause. The `max_rss_kb` field in the larger
benchmark is GNU `time`'s maximum process RSS, not the aggregate peak memory of
all concurrent processes on the node.

## Larger Z3 scale evaluation

`prepare_z3_10m_inputs.sh` prepares the first 10,000,000 reads from the Z3 single-end
run and the first 10,000,000 mate pairs from `SRR25607793`. It uses the complete
`SRR25607804` pair, which has 7,943,970 pairs. Thus the workload contains
27,943,970 input fragments across three units, about 9.3 times the original table's
3,000,000. The two `SRR25607793` FASTQs are fetched as bounded ranges from ENA
because the local full copies are truncated. The preparer validates gzip integrity,
exact subset counts, and paired read IDs, and writes source URLs and subset SHA-256
checksums to `evaluation_runs/z3_10m_inputs/provenance.tsv`.

The large run uses the immutable v0.1.0 script, modern mode, 16 threads, and both
unchunked and 1,000,000-read pre-fastp chunked conditions. The two runs share the
same input manifest and reference. Run preparation and evaluation with:

```bash
bash manuscript_evaluation/prepare_z3_10m_inputs.sh
PLAINMAP="$PWD/manuscript_evaluation/evaluation_runs/plainmap_v0.1.0.sh" \
MANIFEST="$PWD/manuscript_evaluation/evaluation_runs/z3_10m_inputs/manifest.txt" \
REF=/path/to/GCF_009762305.2_mZalCal1.pri.v2_genomic.fna \
EVAL_ROOT="$PWD/manuscript_evaluation/evaluation_runs/large_z3_10m" \
PILOT_FRAGMENTS=10000000 CHUNK_SIZE=1000000 THREADS=16 \
  bash manuscript_evaluation/run_large_z3_benchmark.sh
```

The same workload was also evaluated with a snapshot of the current PlainMap
script (`plainmap_current_20261105.sh`, SHA-256
`2723d9eed4ac2d39fe33af11c0cd2f3046362b2963f5632b1d2662daf91fedb6`).
This run explicitly uses `--dedup-method hybrid` with DeDup v0.12.9 for SE-like
reads; AdaptClean remains off. It is a separate evaluation, not a direct
version-to-version mapping-count comparison, because duplicate handling changed.
The `provenance.tsv` written in each output root records the script, tools, and
settings actually used. On the cluster, the current-code run was submitted as:

```bash
root="$PWD/manuscript_evaluation/evaluation_runs"
PLAINMAP="$root/plainmap_current_20261105.sh" \
MANIFEST="$root/z3_10m_inputs/manifest.txt" \
REF=/home/people/micwe/Biodiversity_Extinction/Japanese_sea_lion/Reference_genomes/GCF_009762305.2_mZalCal1.pri.v2_genomic.fna \
EVAL_ROOT="$root/large_z3_current" \
PILOT_FRAGMENTS=10000000 CHUNK_SIZE=1000000 THREADS=16 \
DEDUP_METHOD=hybrid DEDUP="$root/tools/dedup-release" \
  sbatch -c 16 --mem=48G --time=24:00:00 -J PM_Z3_current \
    --output="$root/PM_Z3_current-%j.out" --export=ALL \
    -- manuscript_evaluation/run_large_z3_benchmark.sh
```

The pinned DeDup wrapper invokes the official DeDup v0.12.9 JAR with Java 17.
These locally staged tool binaries and FASTQs are not committed to Git; the
runner and input-preparation scripts are included for reproducibility.

The ENA-generated `SRR25607793` FASTQs preserve the sequence reads but have `/1`
and `/2` mate markers in the header's second field. The original 1M-pair evaluation
used SRA-derived FASTQs with identical full headers across mates, so that earlier
test remains the direct check of SRA-style header fallback. A comparison of the
first 1,000 R1 sequences found them identical, but their quality strings differ
between the ENA-generated and SRA-derived files. Accordingly, the large workload
tests scaling and chunking within its own two conditions; its mapping counts should
not be interpreted as a direct extension of the original 3M-fragment table.

## Short-read cave hyena and AdaptClean evaluation

`run_cavehyena_50m_adaptclean.sh` compares 50,000,000 pilot fragments from the
single-end ancient cave hyena sample `Ccsp043` against the indexed striped hyena
assembly `GCF_003009895.1`. The two conditions use the same raw FASTQ, PlainMap
snapshot, `ancient-fast` mapping mode, hybrid deduplication, and 16 threads.
The `off` run uses fastp alone; the `on` run additionally applies AdaptClean to
fastp's true SE output.
Pre-trimming chunking is disabled in both conditions. The first 10,000 raw reads
sampled were all 76 bp. In a separate 100,000-read preflight, fastp retained
59,595 reads, of which AdaptClean shortened 2,171; this confirms that the
comparison has a measurable adapter-cleaning intervention before scaling up.

Run from the repository root after staging the pinned DeDup and AdaptClean tools
under `manuscript_evaluation/evaluation_runs/tools/`:

```bash
root="$PWD/manuscript_evaluation/evaluation_runs"
mkdir -p "$root"
for condition in off on; do
  sbatch -c 16 --mem=48G --time=48:00:00 \
    -J "PM_Ccsp043_$condition" --output="$root/PM_Ccsp043_$condition-%j.out" \
    -- manuscript_evaluation/run_cavehyena_50m_adaptclean.sh "$condition"
done
```

Outputs default to `evaluation_runs/cavehyena_50m_adaptclean/{off,on}/`. Each
condition records its input paths, software hashes, parameters, wall/CPU time,
PlainMap statistics, AdaptClean report when enabled, and a histogram of
post-trimming read lengths. GNU `time`'s `max_rss_kb` is maximum process RSS,
not aggregate node memory. The full input FASTQ and reference are not tracked
in Git. Once both runs complete, `summarize_cavehyena_50m.sh` writes a compact
`comparison.tsv` from their validated statistics, timing files, and length
histograms.
