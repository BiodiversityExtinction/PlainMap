# PlainMap

PlainMap is a lightweight, restartable mapping pipeline for modern and degraded DNA sequencing data, designed to process heterogeneous FASTQ datasets with minimal configuration.

PlainMap is implemented as a **single Bash script** and focuses on preprocessing, mapping, duplicate handling, BAM generation, and summary statistics without requiring a workflow engine or complex configuration system.

Current version: **v0.3.0**

---

## Design Goals

- Robust to real-world sequencing messiness
- No reliance on filename conventions
- Safe handling of mixed SE/PE datasets
- Restartable and scalable to very large inputs
- Fragment-aware statistics suitable for ancient DNA
- Minimal configuration, minimal hidden state

---

## Key Features

- Manifest may contain **mixed SE and PE FASTQs**
- SE/PE classification is derived from **FASTQ headers**, not filenames
- Supported header styles:
  - Illumina / CASAVA (` 1:` / ` 2:`)
  - ENA / older Illumina (`/1`, `/2`)
  - **SRA-export FASTQ** with missing mate markers
- SRA-style PE handling:
  - Files grouped by identical normalized first-read name
  - Exactly two files per key allowed
  - Mate assignment is deterministic by manifest order
  - Header keys and read counts must match
  - Ambiguous groupings are rejected
- fastp runs per unit (or unit-chunk), then pools outputs
- Pooled outputs are separated into:
  - **Unmerged PE** (R1 + R2)
  - **SE-like reads** (merged + unpaired + true SE)
- Optional pre-fastp chunking safety valve
- Optional **per-unit pilot limiting before fastp**
- Optional post-pooling mapping chunking (restartable)
- **Restartable mapping (per-chunk resume):** already-produced valid chunk BAMs are not remapped
- **Disk-space friendly mode:** optionally deletes chunk/mapchunk FASTQs as soon as they are safe to remove
- Automatic use of pigz if available
- Selectable duplicate handling: hybrid DeDup/samtools by default, or the v0.1 samtools-only policy
- Separate PE and SE-like duplicate counts and duplication rates
- Optional mapDamage for any library type, automatically parallelized by contig when `--threads` > 1
- Optional export of unmapped reads as FASTQ (`samtools fastq`)

---

## Installation

Clone the repository and make the script executable:

    git clone https://github.com/BiodiversityExtinction/PlainMap.git
    cd PlainMap
    chmod +x plainmap.sh

PlainMap is a single Bash script. No compilation required.

---

## Install with Conda

From the PlainMap directory:

```bash
conda env create -f environment.yml
conda activate plainmap
bash plainmap.sh --help
```

This environment includes the default mapping tools and pins the fastp, BWA,
samtools, and DeDup versions used in the evaluation. PlainMap itself remains a
Bash script; no separate installation is needed. Optional mapDamage and
AdaptClean are not included and must be installed separately if enabled.
The environment was tested on Linux x86-64. It uses Python 3.10 rather than
the Python 3.13.5 used in the original evaluation because the samtools 1.13
Conda package cannot be solved with Python 3.13.

## Dependencies

All tools must be available on PATH unless overridden via flags.

### Required

- bash (must be bash, not sh)
- fastp
- bwa
- samtools
- python3
- GNU split (must support `--filter`)
- gzip and zcat

### Optional (recommended)

- pigz (used automatically if present for parallel compression and integrity tests)

### Optional (only if enabled)

- [DeDup](https://github.com/apeltzer/DeDup) 0.12.9 or later (required for the default `--dedup-method hybrid`; not required for `samtools`)
- mapDamage (required only if `--run-mapdamage` is used)
- [AdaptClean](https://github.com/bilalbioinfo/AdaptClean) (required only with `--run-adaptclean`)

DeDup can be installed from Bioconda:

    conda install -c conda-forge -c bioconda dedup

AdaptClean currently builds from its source repository with `make`, producing
`bin/AdaptClean`. Provide that executable with `--adaptclean /path/to/bin/AdaptClean`
if it is not on `PATH`. PlainMap does not bundle AdaptClean.

---

## Input Format

### Manifest file

- One FASTQ.gz path per line
- Absolute or relative paths supported
- Relative paths resolved relative to the manifest directory
- Comments allowed with `#`
- Filenames do not need to follow a naming convention
- Manifest may contain mixed SE and PE
- Manifest may contain SRA-export FASTQs lacking mate markers

Example:

    Kumamoto_S1_L001_R1_001.fastq.gz
    Kumamoto_S1_L001_R2_001.fastq.gz
    Kumamoto_S1_L002_R1_001.fastq.gz
    Kumamoto_S1_L002_R2_001.fastq.gz
    single_end_lane.fastq.gz

---

## How PlainMap Determines SE vs PE

PlainMap inspects the **first FASTQ header line** of each file.

### Header classification

- Header contains ` 1:` → R1
- Header contains ` 2:` → R2
- Token ending with `/1` → R1
- Token ending with `/2` → R2

### SRA-export FASTQ (no mate markers)

Some SRA-export FASTQs contain **identical headers** for R1 and R2.

In this case:

- Files are grouped by identical normalized first-read name
- If exactly two files share the same key:
  - First file in the manifest is assigned R1
  - Second file is assigned R2
- Strict sanity checks are enforced:
  - Normalized first-read names must match
  - Read counts must match
- If more than two files share the same key, PlainMap refuses to guess

---

## Pipeline Overview

1. Preflight checks (tools, gzip support, split features)
2. Manifest parsing and path resolution
3. Header-based SE/PE classification and unit construction
4. Optional pre-fastp chunking (safety valve)
5. Optional per-unit pilot limiting (before fastp)
6. fastp trimming (JSON output only)
7. Pooling of trimmed outputs:
   - Unmerged PE kept separate
   - SE-like reads pooled together
8. Optional post-fastp mapping chunking
9. Mapping:
   - modern: bwa mem
   - ancient-fast: bwa aln + samse (SE-like only)
   - ancient-all: bwa aln + sampe / samse
10. Merge of mapping outputs (overwrite-safe)
11. Duplicate removal:
   - `hybrid` (default): DeDup `--merged` for SE-like reads and samtools `fixmate + markdup` for unmerged PE reads
   - `samtools`: v0.1-compatible samtools `markdup` processing of the combined BAM
12. Final sorting, indexing, and read group assignment
13. Optional mapDamage
14. Coverage and depth statistics
15. Fragment-aware summary statistics
16. Cleanup (optional retention of intermediates)

---

## Usage Examples

Common defaults (you can usually omit these):

- Resume is ON by default (`--no-resume` disables it)
- Delete-as-you-go is ON by default (`--no-delete-as-you-go` disables it)
- mapDamage is OFF by default (`--run-mapdamage` enables it)
- Unmapped FASTQ export is OFF by default (`--emit-unmapped-fastq` enables it)
- Hybrid duplicate handling is used by default (`--dedup-method samtools` selects the v0.1-compatible policy)
- Ancient modes use `bwa aln -o 2` by default (`--gap-opens 1` restores the earlier default)

PlainMap records the effective configuration of each run in
`SAMPLE.plainmap.run-config.tsv`. Reusing the same output directory and prefix requires
an exact configuration match. Changing any parameter, normalized input/output path,
manifest content, reference identity, PlainMap script, or tool command stops the run
before preprocessing or mapping and reports the changed fields. Argument order and an
explicit default versus an omitted default do not cause a mismatch. Use `--reset` to
intentionally discard the previous run and accept a new configuration.

If preprocessing is interrupted before its checkpoint completes, PlainMap reconstructs
the full set of input chunks from the original FASTQs and reruns preprocessing from the
beginning. Completed mapping chunks are reused separately after a mapping interruption.

Runs created by older PlainMap versions do not have a configuration file and cannot be
safely resumed. Use `--reset`, or select a new prefix, for the first run with configuration
tracking.

### Mapping modes

- `modern`: maps unmerged paired-end reads as paired-end with `bwa mem`, and maps SE-like reads as single-end.
- `ancient-fast`: maps only SE-like reads (merged reads, unpaired reads, and true single-end reads) with `bwa aln + samse`.
- `ancient-all`: maps both unmerged paired-end reads with `bwa aln + sampe` and SE-like reads with `bwa aln + samse`.

### Optional Residual Adapter Cleanup

`--run-adaptclean` applies the external AdaptClean executable to **true single-end FASTQ
units only**, after fastp and before reads are pooled. It does not process merged reads,
unpaired reads derived from paired-end units, or paired-end mates. PlainMap scans each
original SE FASTQ once for cycle inference (an extra input scan) and uses its
longest raw read length as the candidate cycle length for every chunk from that
file. Shorter raw reads are allowed. The candidate
must occur in at least two reads and at least 0.01% of the file; otherwise PlainMap
warns and skips AdaptClean for that file while continuing with fastp. PlainMap takes
a five-base adapter prefix from each chunk's fastp report. If fastp cannot identify
an adapter, specify `--adapter-r1 SEQ` or leave AdaptClean disabled. The post-cleanup
`--minlength` filter and read counts are recorded in `SAMPLE.adaptclean.tsv`; `NA`
for cycle length means uncertain inference, while `NA` for adapter prefix with a
numeric cycle means that chunk had no full-length reads remaining after fastp.

AdaptClean selects reads that still have the **full raw sequencing-cycle length after
fastp**. It can remove a genuine terminal base when a full-length read happens to end
with the first adapter base (and similarly for longer chance matches). It is therefore
off by default; compare results with and without it before using cleaned data for
downstream inference. If a FASTQ mixes multiple sequencing-cycle lengths, only reads
at the inferred maximum are eligible for AdaptClean; shorter cycles cannot be
distinguished reliably from pre-trimmed reads. It does not establish that a read is
endogenous.

## Parameters

Required parameters:

| Parameter | Description |
|---|---|
| `--manifest FILE` | Manifest with one FASTQ(.gz) path per line. |
| `--prefix STRING` | Sample/output prefix used in file names and read group tags. |
| `--ref FILE` | Reference FASTA. |
| `--outdir DIR` | Output directory. |

Optional parameters:

| Parameter | Default | Description |
|---|---:|---|
| `--library-type modern\|ancient-fast\|ancient-all` | `modern` | Mapping mode and aligner strategy. |
| `--threads INT` | `1` | Threads used by fastp, bwa, samtools, and pigz where applicable. |
| `--minlength INT` | `30` | Minimum read length after trimming (`fastp -l`). |
| `--mismatch FLOAT` | `0.01` | Mismatch rate for `bwa aln -n` (ancient-fast/ancient-all). |
| `--gap-opens INT` | `2` | Maximum gap opens for `bwa aln -o` (ancient-fast/ancient-all). Use `1` to reproduce the v0.1-v0.2.1 mapping default. |
| `--mapq INT` | `20` | MAPQ threshold applied during BAM generation. |
| `--dedup-method hybrid\|samtools` | `hybrid` | Duplicate-removal policy. `samtools` applies the combined-BAM v0.1 policy for comparable ongoing analyses. |
| `--emit-unmapped-fastq` | off | Export unmapped reads to `SAMPLE.unmapped.*.fastq.gz` without changing final BAM/stats filtering. |
| `--max-reads-per-chunk INT` | `0` | Pre/post-fastp chunking size cap (`0` disables chunk splitting). |
| `--pilot-fragments INT` | `0` | Per-unit read/fragment cap before fastp (`0` disables pilot mode). |
| `--adapter-r1 SEQ` | auto | Explicit adapter for read 1. |
| `--adapter-r2 SEQ` | auto | Explicit adapter for read 2. |
| `--run-adaptclean` | off | After fastp, run AdaptClean on true single-end reads only. Infers one cycle length per original SE file from a supported maximum raw length, then uses the first five bases of fastp's detected adapter. Skips files with an uncertain cycle length; an undetected adapter is an error. Reads shortened below `--minlength` are removed. |
| `--trim-only` | off | Run trimming/pooling and exit before mapping. |
| `--run-mapdamage` | off | Run mapDamage on final BAM. With `--threads` > 1 and more than one contig has mapped reads, the BAM is split by whole contig and scanned by mapDamage in parallel, then results are merged before final plotting; mapDamage itself is never modified, and the count tables are numerically equivalent to a single-process run. |
| `--mapdamage-serial` | off | Force the original single-process mapDamage behavior even when `--threads` > 1. |
| `--tmpdir DIR` | auto | Custom temp/scratch directory. |
| `--keep-intermediate` | off | Keep `<outdir>/<prefix>/work` after successful completion. |
| `--no-resume` | off | Disable checkpoint resume behavior. |
| `--no-delete-as-you-go` | off | Keep chunk/mapchunk FASTQs until cleanup. |
| `--dry-run` | off | Validate settings and print planned workflow only. |
| `--validate` | off | Validate tools and gzip integrity of manifest FASTQs, then exit. |
| `--reset` | off | Clear prior outputs/checkpoints and the recorded configuration, then restart using the current settings. |

Tool override parameters:

| Parameter | Description |
|---|---|
| `--fastp CMD` | fastp executable/command to use. |
| `--adaptclean CMD` | AdaptClean executable/command to use when `--run-adaptclean` is enabled. |
| `--bwa CMD` | bwa executable/command to use. |
| `--samtools CMD` | samtools executable/command to use. |
| `--dedup CMD` | DeDup executable/command to use. This is a tool-path override, not a method selector. |
| `--mapdamage CMD` | mapDamage executable/command to use. |
| `--python CMD` | Python executable/command to use. |

Modern DNA:

    bash plainmap.sh \
      --manifest manifest.txt \
      --prefix SAMPLE1 \
      --ref reference.fa \
      --outdir results \
      --threads 16

Ancient-fast DNA:

    bash plainmap.sh \
      --manifest manifest.txt \
      --prefix SAMPLE1 \
      --ref reference.fa \
      --outdir results \
      --library-type ancient-fast \
      --threads 16

Ancient-all DNA:

    bash plainmap.sh \
      --manifest manifest.txt \
      --prefix SAMPLE1 \
      --ref reference.fa \
      --outdir results \
      --library-type ancient-all \
      --threads 16

Samtools-only duplicate removal compatible with the v0.1 duplicate policy:

    bash plainmap.sh \
      --manifest manifest.txt \
      --prefix SAMPLE1 \
      --ref reference.fa \
      --outdir results \
      --dedup-method samtools \
      --threads 16

Earlier ancient mapping and duplicate-removal defaults:

    bash plainmap.sh \
      --manifest manifest.txt \
      --prefix SAMPLE1 \
      --ref reference.fa \
      --outdir results \
      --library-type ancient-all \
      --gap-opens 1 \
      --dedup-method samtools \
      --threads 16

Trim-only mode:

    bash plainmap.sh \
      --manifest manifest.txt \
      --prefix SAMPLE1 \
      --ref reference.fa \
      --outdir results \
      --trim-only

Validation mode:

    bash plainmap.sh \
      --manifest manifest.txt \
      --prefix SAMPLE1 \
      --ref reference.fa \
      --outdir results \
      --validate

Pilot run:

    bash plainmap.sh \
      --manifest manifest.txt \
      --prefix SAMPLE1 \
      --ref reference.fa \
      --outdir results \
      --pilot-fragments 10000000

Export unmapped reads as FASTQ:

    bash plainmap.sh \
      --manifest manifest.txt \
      --prefix SAMPLE1 \
      --ref reference.fa \
      --outdir results \
      --emit-unmapped-fastq

`--emit-unmapped-fastq` does not change final BAM filtering or summary statistics.
It must be included in the initial command because all parameters must remain unchanged
when resuming. To enable it for an existing prefix, restart with `--reset`.

---

## Output

Main outputs:

- SAMPLE.final.RG.bam
- SAMPLE.final.RG.bam.bai
- SAMPLE.plainmap.stats.tsv
- SAMPLE.plainmap.coverage.tsv
- SAMPLE.plainmap.run-config.tsv
- SAMPLE_plainmap.log
- SAMPLE/dedup_reports/ (hybrid mode only; DeDup log, histogram, JSON, and captured stdout when SE-like alignments are present)

Optional outputs:

- SAMPLE_mapdamage/ (if enabled)
- SAMPLE.unmapped.R1.fastq.gz (if `--emit-unmapped-fastq`)
- SAMPLE.unmapped.R2.fastq.gz (if `--emit-unmapped-fastq`)
- SAMPLE.unmapped.SE.fastq.gz (if `--emit-unmapped-fastq`)

---

## Summary Statistics

Fragment counts treat an intact paired-end pair as one fragment. Read-length,
mapped-base, and coverage statistics are read- or base-level measures.
Pilot mode affects all counts.

Reported fields include:

- raw reads and fragments
- trimmed, merged, and unpaired counts
- mapped reads and fragments (pre/post deduplication)
- endogenous fraction
- aggregate duplicate-fragment count and duplication rate
- separate PE and SE-like pre/post-deduplication fragment counts and rates
- duplicate-removal method recorded for each alignment class
- average read length
- mapped base pairs
- genome-wide average depth
- genome-wide percent covered

Selected fields are calculated as follows:

- `trimmed_fragments`: fragments retained after fastp and optional AdaptClean.
  An unmerged pair, a merged read, or a single read each counts as one; each
  surviving unpaired mate also counts as one.
- `endog`: unique mapped fragments after MAPQ filtering and duplicate removal,
  divided by raw input fragments. Multiply by 100 for a percentage; this is
  not an independent estimate of biological origin.
- `duprate` (also `pe_duprate` and `se_like_duprate`): duplicate fragments
  divided by mapped fragments before duplicate removal, for the relevant class.
- `avg_readlen`: mean length of the `SEQ` field for retained primary mapped
  reads in the final BAM.
- `mapped_bp`: sum of those `SEQ` lengths, not CIGAR-aligned reference bases.
- `avg_depth` and `pct_covered`: reference-length-weighted mean depth and
  covered-base percentage from `samtools coverage` at the configured MAPQ threshold.

---

## Evaluation

The [evaluation archive](manuscript_evaluation/README.md) contains the input
manifest, runnable protocols, tool-version and parameter information, and
reproduction notes. It separates the original v0.1.0 evaluation from later
tests of larger workloads and optional features:

| Evaluation | Data and question | Code and public materials |
|---|---|---|
| Original mapping modes and restart | Three million pilot fragments from a heterogeneous Japanese sea lion (Z3) manifest containing SE and PE FASTQs with SRA- and ENA-style headers; compare `modern`, `ancient-fast`, and `ancient-all`, then interrupt and resume a chunked `modern` run. | v0.1.0; [manifest](manuscript_evaluation/Z3_manifest_test.txt), [mapping-mode runner](manuscript_evaluation/run_mapping_mode.sh), [restart runner](manuscript_evaluation/run_restart_chunk_evaluation.sh), [reported results](manuscript_evaluation/evaluation_results.tsv), and [restart checks](manuscript_evaluation/restart_checks.tsv). |
| Larger Z3 workload | 27,943,970 input fragments across the same three manifest units; compare unchunked processing with one-million-fragment chunks under `modern` mapping. | v0.1.0 and, separately, a v0.3.0-equivalent pre-release snapshot using hybrid duplicate removal; [input preparation](manuscript_evaluation/prepare_z3_10m_inputs.sh) and [benchmark runner](manuscript_evaluation/run_large_z3_benchmark.sh). Results from the two duplicate policies are not directly interchangeable. |
| Chunking diagnostics | Compare five- and ten-million-fragment chunks with unchunked runs; repeat unchunked runs and inspect stage checksums to locate count differences. A controlled fastp audit compares automatic with fixed adapters. | Pre-release snapshot for the larger-run diagnostics; v0.1.0 for the [fastp audit](manuscript_evaluation/run_fastp_chunk_audit.sh). See the [chunk-size sweep](manuscript_evaluation/run_z3_chunk_size_sweep.sh), [mapping diagnostic](manuscript_evaluation/run_z3_mapping_diagnostic.sh), and [repeat runner](manuscript_evaluation/run_z3_unchunked_repeat.sh). |
| Short-read ancient sample | Fifty million pilot fragments from a single-end cave hyena library mapped in `ancient-fast` mode to striped hyena, with AdaptClean off versus on after fastp. | v0.3.0-equivalent pre-release snapshot; [runner](manuscript_evaluation/run_cavehyena_50m_adaptclean.sh) and [summary script](manuscript_evaluation/summarize_cavehyena_50m.sh). |

The original result tables are committed to Git. The large FASTQs, reference
genomes, BAMs, and output tables from the later runs are **not** in the
repository; the archive documents how to regenerate them and records the
evaluated script snapshot. The November 5 snapshot differs from the v0.3.0
release only in its version label and comments. Optional AdaptClean is not
part of the Conda environment above.

These tests check input handling, execution, chunking, and restart behavior.
Mapped-fragment counts do not establish alignment accuracy or biological
endogeneity. Chunking can change fastp's automatic adapter inference, so
small differences in downstream counts are possible; the controlled fixed-
adapter audit is described in the [evaluation notes](manuscript_evaluation/README.md).

---

## Duplicate Handling

PlainMap provides two duplicate-removal policies:

- `hybrid` (default): merged, unpaired, and true single-end reads are pooled as SE-like reads and processed by DeDup with `--merged`; unmerged paired-end reads are processed separately with samtools `fixmate` and `markdup`. The classes are recombined before downstream processing.
- `samtools`: the combined mapped BAM is processed using the same name-sort/fixmate/coordinate-sort/markdup sequence used by PlainMap v0.1.0. Use this mode for new samples that must remain methodologically comparable with v0.1 analyses.

Using DeDup `--merged` is well suited to merged reads whose aligned ends represent both molecular ends. For true SE reads shorter than the original molecule, the aligned 3' end may not be the molecular end; DeDup can therefore leave some duplicates unresolved. The selected method is recorded in the run configuration and summary statistics, and cannot be changed when resuming. The immutable earlier release tags remain available for exact software-version reproduction. Ancient analyses requiring the v0.1-v0.2.1 mapping default should specify `--gap-opens 1`; v0.2.2 and later default to `2`.

## Limitations (v0.3.0)

- Single read group per sample
- No QC plotting
- No contamination estimation
- No variant calling

---

## Development Status

PlainMap is **v0.3.0** and under active development.

Issues, feedback, and pull requests are welcome.

---

## Citation

If you use PlainMap, please cite the repository metadata in [CITATION.cff](CITATION.cff).

## Contributing

Bug reports, feature requests, and pull requests are welcome. Please see [CONTRIBUTING.md](CONTRIBUTING.md) for reporting guidelines and the suggested smoke test.
