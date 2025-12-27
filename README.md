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

## Dependencies

PlainMap requires the following tools to be available in your PATH:

- bash (must be bash, not sh)
- gzip
- GNU split (with --filter support; part of coreutils)
- fastp
- bwa
- samtools

Notes:

- GNU split is required for deterministic FASTQ chunking.
  On macOS, install via Homebrew:

    brew install coreutils

- No Java, Picard, Docker, Singularity, Snakemake, or Nextflow is required.

---

## Installation

Clone the repository:

    git clone https://github.com/yourusername/plainmap.git
    cd plainmap

Make the script executable:

    chmod +x plainmap.sh

No further installation is required.

---

## Pipeline Overview

PlainMap performs the following steps in order:

1. **Input FASTQ handling**
   - All FASTQ files listed in the manifest are read.
   - If multiple FASTQs are provided per sample (e.g. across lanes or libraries), they are concatenated into a single R1 and R2 file.
   - Read orientation (R1/R2) is detected from FASTQ headers; file naming conventions are ignored.

2. **Adapter and quality trimming**
   - Reads are processed with `fastp`.
   - Adapter trimming is always enabled:
     - By default, adapters are auto-detected by `fastp`.
     - User-specified adapters can be provided if known.
   - Poly-G tail trimming (`fastp -g`) is enabled by default to support BGI/DNBseq data.
   - Reads shorter than the specified minimum length are removed.
   - Low-quality bases are trimmed according to `fastp` defaults.
   - Robust error checking ensures that truncated or corrupt FASTQs cause the pipeline to stop.

3. **Reference preparation**
   - The reference genome is indexed with BWA if required.
   - Index creation is protected by a filesystem lock to prevent corruption when multiple jobs use the same reference concurrently.

4. **Read mapping**
   - Trimmed reads are mapped to the reference genome using BWA:
     - `bwa mem` for modern DNA (paired-end mapping)
     - `bwa aln` for ancient DNA (merged reads, single-end logic)
   - Mapping quality filtering is applied.
   - BAM files are coordinate-sorted.

5. **Chunking and restartability (if enabled)**
   - For very large datasets, FASTQs can be split deterministically into fixed-size chunks.
   - Each chunk is processed independently.
   - Results are merged after mapping.
   - This enables safe restarts after wall-time limits or job interruption.

6. **Duplicate removal**
   - Duplicate reads are removed (not just marked) using `samtools markdup`.
   - Duplicate removal is performed globally after merging all chunks.

7. **Final BAM processing**
   - The final BAM file is sorted and indexed.
   - Read groups are added using `samtools addreplacerg`.

8. **Logging and cleanup**
   - A detailed log file is written, recording:
     - executed commands
     - timestamps
     - tool versions
     - fastp stderr output
     - total wall time
   - By default, all intermediate files are removed to minimise disk usage.
   - Optional flags allow intermediate files to be retained.

9. **Optional trim-only mode**
   - When `--trim-only` is specified, the pipeline stops after trimming and validation.
   - This mode is intended for checking adapter behaviour and read-length distributions before mapping.

## Input Format

### Manifest file

PlainMap uses a manifest file listing input FASTQ files:

- One FASTQ.gz path per line
- Absolute or relative paths supported
- Comments allowed with #

Example manifest.txt:

    # Sample Kumamoto
    Kumamoto_S1_L001_R1_001.fastq.gz
    Kumamoto_S1_L001_R2_001.fastq.gz
    Kumamoto_S1_L002_R1_001.fastq.gz
    Kumamoto_S1_L002_R2_001.fastq.gz

File naming conventions do not matter.
Read orientation (R1/R2) is detected from FASTQ headers.

---

## Usage

### Modern DNA (default)

    bash plainmap.sh \
      -manifest manifest.txt \
      -prefix SAMPLE1 \
      -ref reference.fa \
      -outdir results \
      -t 16

### Ancient DNA

    bash plainmap.sh \
      -manifest manifest.txt \
      -prefix SAMPLE1 \
      -ref reference.fa \
      -outdir results \
      -library-type ancient \
      -t 16

### Adapter trimming check (recommended for new datasets)

    bash plainmap.sh \
      -manifest manifest.txt \
      -prefix SAMPLE1 \
      -ref reference.fa \
      -outdir results \
      --trim-only

This runs adapter and poly-G trimming only, performs all validation checks,
and exits before mapping. It is recommended for new sequencing platforms
(e.g. BGI/DNBseq).

---

## Command-line Parameters

### Required parameters

-manifest FILE  
Text file listing input FASTQ.gz files (one per line).

-prefix STRING  
Sample name used for output files and read groups.

-ref FILE  
Reference genome FASTA file. A BWA index will be created automatically if missing.

-outdir DIR  
Output directory.

---

### Library type

-library-type modern|ancient

- modern (default)  
  Uses bwa mem for paired-end mapping.

- ancient  
  Uses bwa aln with merged reads and single-end duplicate handling.

Default: modern

---

### Threading

-t, --threads INT  
Total number of threads available to the pipeline.
Threads are internally balanced between BWA and samtools.

Default: 1

---

### Read filtering and mapping

-minlength INT  
Minimum read length after trimming.  
Default: 30

-mismatch FLOAT  
Mismatch parameter passed to bwa aln -n (ancient DNA only).  
Default: 0.01

MAPQ  
Minimum mapping quality applied during BAM filtering.  
Value (both modes): 20

---

### Adapter trimming

PlainMap relies on fastp for adapter trimming.

- Adapter trimming is always enabled
- By default, adapters are auto-detected by fastp
- Poly-G tail trimming (fastp -g) is enabled by default to support BGI/DNBseq data

Explicit adapter sequences may be provided when known:

--adapter-r1 SEQ  
--adapter-r2 SEQ  

These override auto-detection and are recommended for:

- non-Illumina platforms (e.g. BGI/DNBseq)
- short or degraded DNA fragments
- custom or non-standard library preparations

---

### Chunking large datasets

-max-reads-per-chunk INT  

Maximum number of reads per FASTQ chunk.

- Enables deterministic splitting of very large FASTQs
- Prevents memory and wall-time failures
- Allows partial restart from completed chunks
- Ensures no read appears in more than one chunk

Set to 0 to disable chunking.  
Default: 0

---

### Temporary files

--tmpdir DIR  

Custom temporary directory (e.g. local scratch).
A sample-specific subdirectory will be created automatically.

Default: <outdir>/tmp

---

### Cleanup behaviour

By default, PlainMap removes all intermediate files once the final BAM
has been successfully created.

--keep-intermediate  

Keeps all intermediate files, including:

- concatenated FASTQs
- chunk FASTQs
- per-chunk BAMs
- temporary files

Default: disabled

---

### Execution control

--resume  
Resume from completed checkpoints if present (default).

--no-resume  
Ignore checkpoints and rerun all steps from scratch.

--dry-run  
Print commands without executing them.

--validate  
Validate tools and input FASTQs, then exit without mapping.

--trim-only  
Run trimming and validation only, then exit before mapping.

---

### Tool overrides

--fastp CMD  
--bwa CMD  
--samtools CMD  

Example:

    --bwa /path/to/bwa

---

## Failure-safe Features

PlainMap is explicitly designed to fail early and loudly when something is wrong.

FASTQ integrity checking:
- All input FASTQs are tested with gzip -t
- Corrupt or truncated files cause immediate failure

Robust fastp error handling:
- fastp stderr is captured and inspected
- Known failure patterns (e.g. igzip errors, truncated input) cause the pipeline to stop
- Empty or missing fastp outputs are treated as fatal errors

Deterministic chunking:
- FASTQs are split by line count (4 lines per read)
- No random subsampling
- No read duplication across chunks
- Safe to restart after partial completion

Atomic BWA indexing:
- Reference indexing is protected by a filesystem lock
- Multiple concurrent jobs using the same reference will not corrupt the index

Restartability:
- Checkpoints are written after major steps
- Jobs can be safely resumed after wall-time limits or node failure

---

## Output

Main outputs:

    SAMPLE.final.RG.bam
    SAMPLE.final.RG.bam.bai

Log file:

    SAMPLE_plainmap.log

---

## Duplicate Handling

- Duplicates are removed, not just marked
- Duplicate removal is performed globally after merging all chunks
- Ancient DNA uses single-end duplicate logic
- Modern DNA uses name-sorted duplicate removal

---

## Read Groups

Read groups are added using:

    samtools addreplacerg

This avoids Java/Picard dependencies and ensures compatibility with downstream tools
such as FreeBayes, bcftools, and GATK.

---

## Limitations (v0.1)

PlainMap v0.1 is intentionally conservative.

- Single read group per sample
- No built-in QC plots
- No contamination estimation
- No variant calling

---

## When to Use PlainMap

PlainMap is ideal if you:

- want full visibility into mapping steps
- need robust restartability on HPC systems
- are working with ancient or historical DNA
- are teaching or learning mapping workflows
- want a simple, auditable reference pipeline

If you need a full end-to-end framework (damage profiling, contamination estimation,
variant calling), consider workflow-based pipelines such as nf-core/eager.

---

## Citation

If you use PlainMap, please cite:

Author(s). PlainMap: a transparent, failure-aware mapping pipeline for ancient and modern DNA. Journal, Year.

A DOI will be added upon publication.

---

## Development Status

PlainMap is currently v0.1 and under active development.

Feedback, issues, and pull requests are welcome.
