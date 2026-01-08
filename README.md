# PlainMap

PlainMap is a transparent, failure-aware mapping pipeline for ancient and modern DNA.
It is designed to simplify read mapping for non-expert users while remaining robust,
inspectable, and restartable on HPC systems.

PlainMap deliberately avoids workflow engines and black-box abstractions.
Every command executed is visible, logged, and reproducible.

Current version: v0.1 (preliminary)

---

## Design Philosophy

PlainMap prioritises:

- transparency over abstraction
- robustness over silent failure
- restartability over monolithic execution
- simplicity over feature bloat

It is intended as a reference implementation of best-practice read mapping,
particularly suited to ancient, historical, or otherwise problematic datasets.

---

## Installation

### Basic installation

Clone the repository:

    git clone https://github.com/yourusername/plainmap.git
    cd plainmap

Make the script executable:

    chmod +x plainmap.sh

No further installation is required.

PlainMap is a single Bash script and does not require compilation, containers,
or workflow engines.

---

## Dependencies

PlainMap assumes required tools are callable from the system PATH, unless overridden.

### Required tools

- bash (must be run with bash, not sh)
- fastp
- bwa
- samtools
- python3
- gzip
- zcat
- GNU split (must support --filter)

### Optional (ancient DNA only)

- mapDamage (required only if -library-type ancient is used)

---

### Optional: Conda installation of dependencies

PlainMap does not require Conda, but Conda can be used to install dependencies:

    conda create -n plainmap \
      fastp bwa samtools coreutils python \
      -c bioconda -c conda-forge
    conda activate plainmap

For ancient DNA:

    conda install -c bioconda mapdamage

---

## Pipeline Overview

PlainMap performs the following steps in order:

1. Input FASTQ handling
   - Reads all FASTQ.gz paths listed in the manifest.
   - Resolves relative paths relative to the manifest directory.
   - Concatenates multiple FASTQs per sample into single R1/R2 files.
   - Detects R1/R2 orientation from FASTQ headers.

2. Optional chunking
   - Deterministic splitting of FASTQs into fixed-size chunks.
   - Each chunk processed independently and merged.
   - Enables safe restarts on HPC systems.

3. Adapter and quality trimming
   - fastp trimming with poly-G trimming enabled.
   - Minimum read length enforced.
   - Corrupt or truncated FASTQs cause immediate failure.

4. Reference preparation
   - BWA index created if missing.
   - Index creation protected by filesystem lock.

5. Read mapping
   - modern: bwa mem
   - ancient: bwa aln + bwa samse
   - Mapping quality filtering applied.

6. Merge and duplicate removal
   - Per-chunk BAMs merged.
   - Duplicates removed using samtools markdup -r.
   - Appropriate logic for SE, PE, and ancient data.

7. Final BAM processing
   - Coordinate sort and index.
   - Read groups added using samtools addreplacerg.

8. Optional mapDamage
   - Runs mapDamage if -library-type ancient.

9. Coverage and summary statistics
   - Fast per-contig coverage via samtools coverage.
   - Genome-wide average coverage computed as length-weighted mean depth.
   - Optional slow depth statistics via samtools depth.
   - Average read length and mapped base pairs computed.
   - One-line per-sample stats TSV written.

10. Cleanup
    - Intermediate files removed by default.
    - --keep-intermediate retains all work files.

---

## Input Format

### Manifest file

- One FASTQ.gz path per line
- Absolute or relative paths supported
- Comments allowed with #

Example:

    # Sample Kumamoto
    Kumamoto_S1_L001_R1_001.fastq.gz
    Kumamoto_S1_L001_R2_001.fastq.gz
    Kumamoto_S1_L002_R1_001.fastq.gz
    Kumamoto_S1_L002_R2_001.fastq.gz

File naming conventions do not matter.

---

## Validation and Dry-Run Modes

PlainMap provides two non-mapping execution modes for safety and inspection.
These modes serve different purposes and are not interchangeable.

### --dry-run

Dry-run mode prints the planned execution without running the pipeline.

Behaviour:
- No trimming or mapping is performed
- No FASTQs are read or tested
- No output files are created
- No gzip integrity checks are run

Dry-run prints:
- resolved absolute paths
- detected sequencing mode (SE/PE)
- planned pipeline steps

This mode is intended for:
- checking command-line arguments
- verifying path resolution
- inspecting the execution plan
- preparing HPC job scripts

Example:

    bash plainmap.sh \
      -manifest manifest.txt \
      -prefix SAMPLE1 \
      -ref reference.fa \
      -outdir results \
      --dry-run

---

### --validate

Validation mode actively checks inputs and tools, but does not perform trimming
or mapping.

Behaviour:
- No trimming or mapping is performed
- All required tools are checked for availability
- Every FASTQ listed in the manifest is tested with gzip -t
- Corrupt or truncated FASTQs cause immediate failure
- No mapping outputs are created

This mode is intended for:
- preflight validation before long runs
- detecting corrupted or incomplete FASTQs
- checking shared filesystems and staging areas
- automated checks in CI or HPC environments

Example:

    bash plainmap.sh \
      -manifest manifest.txt \
      -prefix SAMPLE1 \
      -ref reference.fa \
      -outdir results \
      --validate

In short:
- --dry-run answers “what would happen?”
- --validate answers “is my data and environment safe to run?”

## Usage

Modern DNA:

    bash plainmap.sh \
      -manifest manifest.txt \
      -prefix SAMPLE1 \
      -ref reference.fa \
      -outdir results \
      -t 16

Ancient DNA:

    bash plainmap.sh \
      -manifest manifest.txt \
      -prefix SAMPLE1 \
      -ref reference.fa \
      -outdir results \
      -library-type ancient \
      -t 16

Trim-only mode:

    bash plainmap.sh \
      -manifest manifest.txt \
      -prefix SAMPLE1 \
      -ref reference.fa \
      -outdir results \
      --trim-only

---

## Output

Main outputs:

- SAMPLE.final.RG.bam
- SAMPLE.final.RG.bam.bai
- SAMPLE.plainmap.stats.tsv
- SAMPLE.plainmap.coverage.tsv
- SAMPLE_plainmap.log

Optional:

- SAMPLE.plainmap.depth.tsv (with --depth-stats)
- SAMPLE_mapdamage/ (ancient DNA only)

---

## Summary Statistics

The stats TSV includes:

- raw and trimmed read counts
- merged and unpaired reads
- mapped reads (pre/post deduplication)
- endogenous fraction
- duplication rate
- average mapped read length
- mapped base pairs
- genome-wide average coverage

---

## Limitations (v0.1)

- Single read group per sample
- No QC plots
- No contamination estimation
- No variant calling
- MAPQ threshold fixed in script

---

## Development Status

PlainMap is v0.1 and under active development.

Issues, feedback, and pull requests are welcome.
