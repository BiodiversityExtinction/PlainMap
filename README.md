# PlainMap

PlainMap is a mapping pipeline for ancient and modern DNA designed to “just run” on HPC systems with minimal fuss.

Current version: **v0.1 (preliminary, under active development)**

---

## Installation

Clone the repository:

    git clone https://github.com/yourusername/plainmap.git
    cd plainmap

Make the script executable:

    chmod +x plainmap.sh

PlainMap is a single Bash script.

---

## Dependencies

PlainMap assumes required tools are callable from the system `PATH`, unless overridden with flags.

### Required tools

- bash (must be run with bash, not sh)
- fastp
- bwa
- samtools
- python3
- zcat
- GNU split (must support `--filter`)
- gzip

### Optional (recommended)

- pigz (used automatically if present for parallel compression/decompression and integrity tests)

### Optional (ancient DNA only)

- mapDamage (required only if `-library-type ancient` is used)

---

### Optional: Conda installation of dependencies

PlainMap does not require Conda, but Conda can be used to install dependencies:

    conda create -n plainmap \
      fastp bwa samtools coreutils python pigz \
      -c bioconda -c conda-forge
    conda activate plainmap

For ancient DNA:

    conda install -c bioconda mapdamage

---

## Input Format

### Manifest file

- One `FASTQ.gz` path per line
- Absolute or relative paths supported
- Relative paths are resolved relative to the manifest directory
- Comments allowed with `#`
- Filenames do **not** need to follow a naming convention
- The manifest may contain a **mix of SE and PE** FASTQs

Example:

    # Sample Kumamoto
    Kumamoto_S1_L001_R1_001.fastq.gz
    Kumamoto_S1_L001_R2_001.fastq.gz
    Kumamoto_S1_L002_R1_001.fastq.gz
    Kumamoto_S1_L002_R2_001.fastq.gz

---

## How PlainMap decides SE vs PE (including ENA headers)

PlainMap inspects the first header line in each FASTQ and classifies read direction as:

- Illumina/CASAVA style: header contains ` 1:` (R1) or ` 2:` (R2)
- ENA style: header token ends with `/1` (R1) or `/2` (R2)

PlainMap supports **mixed SE + PE** datasets in the same manifest.

---

## Pipeline Overview

PlainMap performs the following steps:

1. **Input handling**
   - Reads all FASTQ.gz paths from the manifest
   - Resolves paths to absolute paths
   - Classifies reads as R1/R2 using FASTQ header inspection (Illumina + ENA styles)
   - Groups files into “units”:
     - SE unit = one R1 FASTQ
     - PE unit = matching R1/R2 FASTQs (lane-splits must have matching counts)

2. **Adaptive pre-fastp chunking (safety valve)**
   - Controlled by `-max-reads-per-chunk INT`
   - For each unit:
     - If unit read count > INT: split **before fastp** into chunks of up to INT reads
     - If unit read count ≤ INT: do not chunk (processed as a single chunk)
   - If `-max-reads-per-chunk 0`, this is disabled

3. **Adapter/quality trimming (fastp, JSON only)**
   - fastp runs on each unit (or on each unit-chunk if pre-fastp chunking triggered)
   - Minimum read length enforced (`-minlength`)
   - Poly-G trimming enabled
   - fastp HTML output is **not generated** (JSON only)
   - For PE:
     - adapter detection is enabled automatically unless adapters are explicitly provided

4. **Pooling after trimming**
   - Trimmed outputs are appended into pooled trimmed FASTQs (concatenated gzip members)
   - This avoids mixing adapters *before trimming*, while still enabling fast mapping (map once)

5. **Reference preparation**
   - BWA index created if missing
   - Index creation protected by a filesystem lock

6. **Read mapping (pooled; map once)**
   - **modern**:
     - PE: `bwa mem` on pooled trimmed R1/R2
     - SE: `bwa mem` on pooled SE reads + pooled PE rescue reads (merged/unpaired)
   - **ancient**:
     - SE: `bwa aln` + `bwa samse`
     - PE: fastp merges reads first; pooled merged reads map as SE via `bwa aln` + `bwa samse`
   - MAPQ filter applied during BAM conversion (`-q`, default 20)

7. **Merge and duplicate removal**
   - All pooled mapping outputs merged
   - Duplicate removal using `samtools markdup -r`
   - For modern runs that include PE data, fixmate workflow is used before markdup

8. **Final BAM processing**
   - Coordinate sort and index
   - Read groups added using `samtools addreplacerg`

9. **Optional mapDamage**
   - Runs if `-library-type ancient`

10. **Coverage/depth + summary statistics**
   - Per-contig coverage and mean depth via `samtools coverage`
   - Genome-wide:
     - `avg_depth` = length-weighted mean of per-contig mean depth
     - `pct_covered` = length-weighted percent of bases covered (>=1x)
   - Fragment-aware mapping counts:
     - `mapped_fragments_*` counts templates/fragments (PE counted via mapped read1; SE counted via non-paired reads)

11. **Cleanup**
   - Intermediate files removed by default
   - `--keep-intermediate` retains the entire work directory

---

## Validation and Dry-Run Modes

### --dry-run

Dry-run prints the plan without running the pipeline.

Behaviour:
- No trimming or mapping is performed
- No FASTQs are tested
- No output files are created
- No gzip integrity checks are run

Example:

    bash plainmap.sh \
      -manifest manifest.txt \
      -prefix SAMPLE1 \
      -ref reference.fa \
      -outdir results \
      --dry-run

---

### --validate

Validation checks tools and validates that all FASTQs are readable and not truncated.

Behaviour:
- No trimming or mapping is performed
- All required tools are checked
- Every FASTQ is tested with `pigz -t` if pigz exists, otherwise `gzip -t`
- Corrupt FASTQs cause immediate failure

Example:

    bash plainmap.sh \
      -manifest manifest.txt \
      -prefix SAMPLE1 \
      -ref reference.fa \
      -outdir results \
      --validate

---

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

Enable pre-fastp chunking safety valve (recommended for extremely large single inputs):

    bash plainmap.sh \
      -manifest manifest.txt \
      -prefix SAMPLE1 \
      -ref reference.fa \
      -outdir results \
      -t 16 \
      -max-reads-per-chunk 50000000

---

## Output

Main outputs:

- `SAMPLE.final.RG.bam`
- `SAMPLE.final.RG.bam.bai`
- `SAMPLE.plainmap.stats.tsv`
- `SAMPLE.plainmap.coverage.tsv`
- `SAMPLE_plainmap.log`

Optional:

- `SAMPLE_mapdamage/` (ancient DNA only)

---

## Summary Statistics (TSV)

The stats file includes:

- `raw_R1_reads`, `raw_R2_reads`
- `raw_fragments` (fragment/template count; defined as total raw R1 reads across all units)
- `trimmed_fragments`
- `merged_reads`, `unpaired_reads`
- `mapped_reads_all`, `mapped_reads_unique`
- `mapped_fragments_all`, `mapped_fragments_unique` (**fragment-aware**)
- `endog` (mapped_fragments_unique / raw_fragments)
- `duprate` (fragment-space duplication rate)
- `avg_readlen`, `mapped_bp`
- `avg_depth` (genome-wide length-weighted mean depth)
- `pct_covered` (genome-wide percent bases covered >=1x)

---

## Limitations (v0.1)

- Single read group per sample
- No QC plotting
- No contamination estimation
- No variant calling
- MAPQ threshold is currently fixed in the script (default 20)

---

## Development Status

PlainMap is **v0.1** and under active development.

Issues, feedback, and pull requests are welcome.
