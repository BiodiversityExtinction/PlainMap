# Manuscript evaluation

This directory records manifests, commands, parameters, and results for the
original evaluation and later scale, chunking, restart, and AdaptClean tests.
The [combined manuscript table](manuscript_table_v030_sra.md) uses SRA-derived
Z3 reads and the cave hyena dataset. Its unchunked Z3 mapping modes used the
immutable `v0.3.0` tag; the newer direct-chunk and restart runs used commit
`f637f70`. Both scripts print `0.3.0`, but they are different code revisions.
The input, script, and tool provenance for each result group is linked below.
Earlier v0.1.0 and ENA-derived evaluations remain archived, but their counts
should not be mixed with the SRA-derived evaluation. Full run outputs under
`evaluation_runs/` are ignored by Git; large FASTQs, reference genomes, and
BAMs are not in this repository.

The original manuscript table used PlainMap v0.1.0 and its samtools-only
duplicate-removal policy. Check out the immutable `v0.1.0` tag for the exact
PlainMap script version. Some later hybrid-policy Z3 and cave hyena
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
| Current SRA-derived Z3 mapping modes | 27,943,970 input fragments, including 10 million `SRR25607793` SRA pairs with identical mate headers; compare unchunked `modern`, `ancient-fast`, and `ancient-all`. | Tagged `v0.3.0`; [SRA preparation](prepare_z3_10m_sra.sh), [input provenance](z3_sra_10m_input_provenance.tsv), [evaluation runner](run_z3_sra_v030.sh), [results](z3_sra_v030_results.tsv), and [run provenance](z3_sra_v030_provenance.tsv). The same results file also archives the older, superseded chunking series. |
| Current modern direct chunking and restart | Same SRA-derived Z3 input; compare 1M, 5M, and 10M chunk sizes and resume after three of 46 mapping chunks. | Commit `f637f70`; [runner](run_z3_integrated_chunking.sh), [results](z3_modern_direct_chunk_results.tsv), [provenance](z3_modern_direct_chunk_provenance.tsv), [adapter calls](z3_direct_chunk_adapters.tsv), [restart checks](z3_modern_direct_chunk_restart_checks.tsv), and [resume time](z3_modern_direct_chunk_restart_time.tsv). |
| Current ancient-mode direct chunking | Same SRA-derived Z3 input; compare 1M and 10M chunks in `ancient-fast` and `ancient-all`. | Commit `f637f70`; [runner](run_z3_ancient_chunking.sh), [results](z3_ancient_direct_chunk_results.tsv), and [provenance](z3_ancient_direct_chunk_provenance.tsv). |
| Original mapping modes and restart | Three million pilot fragments from a heterogeneous Japanese sea lion (Z3) manifest containing SE and PE FASTQs with SRA- and ENA-style headers; compare `modern`, `ancient-fast`, and `ancient-all`, then interrupt and resume a chunked `modern` run. | v0.1.0; [manifest](Z3_manifest_test.txt), [mapping-mode runner](run_mapping_mode.sh), [restart runner](run_restart_chunk_evaluation.sh), [reported results](evaluation_results.tsv), and [restart checks](restart_checks.tsv). |
| Archived larger ENA-derived Z3 workload | 27,943,970 input fragments across the same three manifest units, but with ENA-generated `SRR25607793` mates; compare unchunked processing with one-million-fragment chunks under `modern` mapping. | v0.1.0 and, separately, the v0.3.0-equivalent snapshot; [input preparation](prepare_z3_10m_inputs.sh) and [benchmark runner](run_large_z3_benchmark.sh). Do not combine counts from these sources. |
| Archived ENA-derived released-version comparison | The same ENA-derived 27,943,970 Z3 fragments, processed without chunking in `modern`, `ancient-fast`, and `ancient-all` modes. | Exact v0.3.0; [mode runner](run_z3_28m_modes.sh), [summary collector](summarize_z3_28m_modes.py), [results](z3_v030_28m_results.tsv), [input checksums](z3_v030_28m_input_checksums.tsv), and [run provenance](z3_v030_28m_provenance.tsv). |
| Archived ENA-derived interrupted restart | Interrupt an ENA-derived 27,943,970-fragment `modern` run after one mapping chunk, then resume and compare with its uninterrupted baseline. | Exact v0.3.0 restart; [runner](run_z3_28m_restart.sh), [checks](z3_v030_28m_restart_checks.tsv), [provenance](z3_v030_28m_restart_provenance.tsv), and [resume resource usage](z3_v030_28m_restart_time.tsv). |
| Chunking diagnostics | Compare five- and ten-million-fragment chunks with unchunked runs; repeat unchunked runs and inspect stage checksums to locate count differences. A controlled fastp audit compares automatic with fixed adapters. | Snapshot for the larger-run diagnostics; v0.1.0 for the [fastp audit](run_fastp_chunk_audit.sh). See the [chunk-size sweep](run_z3_chunk_size_sweep.sh), [mapping diagnostic](run_z3_mapping_diagnostic.sh), and [repeat runner](run_z3_unchunked_repeat.sh). |
| Short-read ancient sample | Fifty million pilot fragments from a single-end cave hyena library mapped in `ancient-fast` mode to striped hyena, with AdaptClean off versus on after fastp. | Snapshot; [runner](run_cavehyena_50m_adaptclean.sh) and [summary script](summarize_cavehyena_50m.sh). AdaptClean is not included in the Conda environment. |

These tests check input handling, execution, chunking, and restart behavior.
Mapped-fragment counts do not establish alignment accuracy or biological
endogeneity. The older per-chunk adapter audit is described under Results;
the current direct-chunk runs inferred adapters once per input unit and still
showed some paired-end alignment differences between chunk sizes.

## Current SRA-derived Z3 evaluation

The primary Z3 analysis uses the first 10,000,000 single-end reads from
`SRR25607782`, the first 10,000,000 paired SRA spots from `SRR25607793`, and
all 7,943,970 pairs from `SRR25607804`. The SRA Toolkit v3.4.1 extracted the
`SRR25607793` mates with `fastq-dump --split-files --skip-technical --gzip -N 1
-X 10000000`. Both files passed gzip and read-count checks; all 10 million
corresponding FASTQ headers were identical. PlainMap logged its SRA-header
fallback and inferred the unit as paired-end in each mapping mode. The source
and compressed FASTQ hashes are in the [input provenance](z3_sra_10m_input_provenance.tsv).

The tagged `v0.3.0` script processed the same 27,943,970 raw fragments without
chunking in `modern`, `ancient-fast`, and `ancient-all` modes. Each run used 16
threads, minimum trimmed length 30 bp, MAPQ 20, and hybrid duplicate removal.
Unique mapped fragments numbered 70,348, 39,983, and 40,837, respectively;
ancient-all exceeded ancient-fast by 854. These runs were performed on
`compute02`. The [SRA evaluation results](z3_sra_v030_results.tsv) also retain
the original tagged-release chunking and restart tests as historical records;
those earlier chunked counts are **not** the current manuscript results.

The current direct-chunk analysis used PlainMap commit `f637f70` with the same
manifest and reference. It probed up to 10,000 raw reads per input unit to
infer [unit-specific adapters](z3_direct_chunk_adapters.tsv), then mapped
fastp output chunks directly. All seven chunked conditions retained
27,928,803 trimmed fragments and 15,036,758 merged reads. Modern-mode unique
mapped fragments numbered 71,006, 70,341, and 70,328 at 1M, 5M, and 10M
chunks, respectively. Ancient-fast retained 39,983 at both 1M and 10M;
the final alignment records were identical between those two runs. Ancient-all
retained 40,863 and 40,837, respectively; its 26-fragment difference came from
paired-end alignments and was present before duplicate removal. The
[modern results](z3_modern_direct_chunk_results.tsv) ran sequentially on
`compute04`; the [ancient results](z3_ancient_direct_chunk_results.tsv) ran
sequentially on `node12`. The groups cannot support a controlled cross-node
timing comparison. Maximum process RSS is not total job memory, and mapped
counts do not establish alignment accuracy or biological endogeneity.

For the current restart validation, the one-million-fragment modern run was
interrupted after three of 46 mapping chunks had completed. PlainMap reused
all three valid BAMs, and the resumed alignment records and complete summary
statistics matched the uninterrupted run exactly. Its 8.32-minute wall time
and 4.85-GiB peak process RSS cover only the remaining work after interruption.
The [checks](z3_modern_direct_chunk_restart_checks.tsv),
[resume time](z3_modern_direct_chunk_restart_time.tsv), and
[modern run provenance](z3_modern_direct_chunk_provenance.tsv) record this
result. All current final Z3 BAMs passed `samtools quickcheck` and contained no
unmapped or duplicate-flagged records.

To regenerate these analyses, stage the indexed California sea lion reference,
DeDup v0.12.9, the three Z3 FASTQ units, and SRA Toolkit v3.4.1. See the
[SRA input preparer](prepare_z3_10m_sra.sh) and its
[input provenance](z3_sra_10m_input_provenance.tsv). The tagged-release runner
requires the exact `v0.3.0` script; the two direct-chunk runners require the
exact `f637f70` script. For example, from the repository root after staging the
inputs and tools under `manuscript_evaluation/evaluation_runs/`:

```bash
root="$PWD/manuscript_evaluation/evaluation_runs"
mkdir -p "$root/z3_v030_28m_modes"
git show v0.3.0:plainmap.sh > "$root/z3_v030_28m_modes/plainmap_v0.3.0.sh"
git show f637f70:plainmap.sh > "$root/plainmap_integrated_v030.sh"
sha256sum "$root/plainmap_integrated_v030.sh" | cut -d ' ' -f 1 \
  > "$root/plainmap_integrated_v030.sha256"
chmod +x "$root/plainmap_integrated_v030.sh"
```

The [tagged-release runner](run_z3_sra_v030.sh),
[modern direct-chunk runner](run_z3_integrated_chunking.sh), and
[ancient direct-chunk runner](run_z3_ancient_chunking.sh) record parameters and
script hashes. Their defaults contain local cluster paths, and the runners
pin the original absolute-path manifest checksum. Exact reproduction therefore
requires the original input paths; on another system, adjust those paths and
the manifest checksum check after verifying the FASTQ hashes in the input
provenance. Set `REF` to the indexed reference location. Full FASTQs,
reference sequences, and BAMs are not committed to Git.

## Archived original pilot dataset

The heterogeneous Z3 manifest contains:

- one single-end unit: `SRR25607782`
- one SRA-style paired-end unit whose mate files have identical read names:
  `SRR25607793`
- one paired-end unit with explicit `/1` and `/2` header suffixes: `SRR25607804`

The evaluated file order and names are stored in
[`Z3_manifest_test.txt`](Z3_manifest_test.txt). Despite their `1M` filenames,
`SRR25607793_1M_1.fastq.gz` and `SRR25607793_1M_2.fastq.gz` each contain
1,125,000 reads. PlainMap uses their first 1,000,000 pairs because it applies
`--pilot-fragments 1000000` independently to each unit, provided record order is
preserved.

Reads were mapped to the California sea lion (*Zalophus californianus*) mZalCal1.pri.v2
assembly (NCBI RefSeq accession `GCF_009762305.2`).

## Archived original v0.1.0 requirements

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
        --wrap="bash '$PWD/manuscript_evaluation/run_mapping_mode.sh' '$mode'"
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
    --wrap="bash '$PWD/manuscript_evaluation/run_restart_chunk_evaluation.sh'"
```

The output root defaults to `manuscript_evaluation/evaluation_runs/restart_chunk/` and can
be changed with `EVAL_ROOT`. Forced pre-trimming chunking temporarily required
approximately 12 GB for this dataset, so the output filesystem should have sufficient
free space.

## Results

- [`evaluation_results.tsv`](evaluation_results.tsv) contains the original
  v0.1.0 evaluation values.
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

The same workload was also evaluated with a November 5 pre-release PlainMap
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

## Released v0.3.0 mapping-mode comparison on larger Z3 input

The [mode runner](run_z3_28m_modes.sh) uses the same validated 27,943,970-fragment
manifest and reference as the larger Z3 scale test, with 16 threads,
`--pilot-fragments 10000000`, no chunking, and hybrid duplicate removal in all
three modes. Ancient modes use the v0.3.0 default `bwa aln -o 2`. The runner
records the script and manifest hashes, tool versions, parameters, elapsed
time, CPU time, and maximum process RSS; it checks each final BAM for
integrity and for unmapped or duplicate-flagged records.

To reproduce from the repository root after preparing `z3_10m_inputs/` and
activating the Conda environment, set `REF` to the indexed California sea lion
reference and run:

```bash
root="$PWD/manuscript_evaluation/evaluation_runs"
mkdir -p "$root/z3_v030_28m_modes"
git show v0.3.0:plainmap.sh > "$root/z3_v030_28m_modes/plainmap_v0.3.0.sh"
export PLAINMAP="$root/z3_v030_28m_modes/plainmap_v0.3.0.sh"
export MANIFEST="$root/z3_10m_inputs/manifest.txt"
export REF=/path/to/GCF_009762305.2_mZalCal1.pri.v2_genomic.fna
export DEDUP="$(command -v dedup)"
export EVAL_ROOT="$root/z3_v030_28m_modes"
export THREADS=16
for mode in modern ancient-fast ancient-all; do
  sbatch -c 16 --mem=48G --time=24:00:00 \
    --job-name="PM28_${mode}" \
    --output="$EVAL_ROOT/PM28_${mode}-%j.out" \
    --export=ALL \
    --wrap="bash '$PWD/manuscript_evaluation/run_z3_28m_modes.sh' '$mode'"
done
```

After all three jobs complete, run
`python3 manuscript_evaluation/summarize_z3_28m_modes.py "$EVAL_ROOT"` to write
`summary.tsv`. The validated released-version run produced 71,089 unique
mapped fragments in `modern`, 39,855 in `ancient-fast`, and 40,763 in
`ancient-all`. Raw, trimmed, merged, and unpaired counts were identical across
the three runs. All final BAMs passed `samtools quickcheck` and contained zero
unmapped or duplicate-flagged records. [Results](z3_v030_28m_results.tsv)
include wall time, CPU time, and maximum **process** RSS, which is not total
node memory. These are single runs on different nodes, so runtime differences
should not be treated as a controlled performance comparison between modes.
The `endog` column is a fraction of raw fragments retained as unique mapped
fragments, not an independent estimate of biological endogeneity.
The input FASTQs, reference, and full run outputs are not committed to Git.
The [input checksums](z3_v030_28m_input_checksums.tsv) are SHA-256 hashes of
the compressed FASTQ files used in these runs. Tool versions, script hash,
parameters, and Slurm job IDs are in the
[run provenance](z3_v030_28m_provenance.tsv).

## Released v0.3.0 interrupted restart on larger Z3 input

The [restart runner](run_z3_28m_restart.sh) uses the same 27,943,970-fragment
manifest and reference, 16 threads, `modern` mode, hybrid duplicate removal,
`--pilot-fragments 10000000`, and one-million-fragment preprocessing chunks.
It launches the exact v0.3.0 script, interrupts it after at least one but not
all mapping-chunk BAMs are valid, then reruns the identical command. The runner
checks that completed chunks were reused and compares the resumed alignment
records and statistics with the previously completed, uninterrupted large
chunked run. The comparator used a frozen script whose executable logic is
identical to v0.3.0; only comments and the version label differ.

From the repository root, after preparing the inputs and completed chunked
baseline described above:

```bash
root="$PWD/manuscript_evaluation/evaluation_runs"
export PLAINMAP="$root/z3_v030_28m_modes/plainmap_v0.3.0.sh"
export MANIFEST="$root/z3_10m_inputs/manifest.txt"
export REF=/path/to/GCF_009762305.2_mZalCal1.pri.v2_genomic.fna
export DEDUP="$(command -v dedup)"
export BASELINE_OUT="$root/large_z3_current/chunked"
export EVAL_ROOT="$root/z3_v030_28m_restart"
sbatch -c 16 --mem=48G --time=08:00:00 \
  --job-name=PM_Z3_v030_restart \
  --output="$root/PM_Z3_v030_restart-%j.out" \
  --export=ALL \
  --wrap="bash '$PWD/manuscript_evaluation/run_z3_28m_restart.sh'"
```

The runner records its checks, script and manifest hashes, tool versions,
parameters, and resume resource usage in its ignored `EVAL_ROOT` output.
The completed job (`231879`, exit 0) stopped after one of 29 mapping chunks
was valid. On resumption, PlainMap reused that chunk and finished successfully.
The resumed alignment records and complete statistics TSV matched the
uninterrupted chunked baseline exactly; `samtools quickcheck` passed, and the
final BAM contained no unmapped or duplicate-flagged records. The resume
attempt took 399.75 seconds; this is remaining-run time, not a measure of
restart overhead or a like-for-like runtime comparison with an uninterrupted
run. The archived [checks](z3_v030_28m_restart_checks.tsv),
[provenance](z3_v030_28m_restart_provenance.tsv), and
[timing](z3_v030_28m_restart_time.tsv) record the result.

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
