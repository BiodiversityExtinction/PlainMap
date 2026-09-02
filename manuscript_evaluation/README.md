# Manuscript evaluation

This directory records the manifest, commands, parameters, and expected results used for
the PlainMap manuscript evaluation. Sequence data and the 2.44 GB reference FASTA are not
committed to the repository.

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
MAPQ 20, `bwa aln -n 0.01` for ancient modes, automatic adapter detection, resume enabled,
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
