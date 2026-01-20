# PlainMap

PlainMap is a mapping pipeline for ancient and modern DNA designed to **just run** on HPC systems with minimal fuss.

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
     - PE unit = matching R1/R2 FASTQs  
       (lane-splits must be provided as matching R1/R2 lists; PlainMap refuses unsynchronized PE)

2. **Adaptive pre-fastp chunking (safety valve)**
   - Controlled by `-max-reads-per-chunk INT`
   - For each unit:
     - If unit read count > INT: split **before fastp** into chunks of up to INT reads
     - If unit read count ≤ INT: do not chunk (processed as a single chunk)
   - If `-max-reads-per-chunk 0`, this is disabled

3. **Optional pre-fastp pilot limiting**
   - Controlled by `--pilot-fragments INT`
   - If set (>0), PlainMap limits the number of raw fragments processed **before fastp**
   - Applied **per raw chunk**:
     - SE: first `INT` reads
     - PE: first `INT` pairs (R1 and R2 limited symmetrically)
   - This is intended for fast “pilot” comparisons without trimming/mapping the full dataset

4. **Adapter/quality trimming (fastp, JSON only)**
   - fastp runs per unit (or per unit-chunk if pre-fastp chunking is active)
   - Minimum read length enforced (`-minlength`)
   - Poly-G trimming enabled
   - fastp HTML output is **not generated** (JSON only)
   - For PE:
     - adapter detection is enabled automatically unless adapters are explicitly provided

5. **Pooling after trimming**
   - Trimmed outputs are appended into pooled trimmed FASTQs (concatenated gzip members)
   - This avoids mixing adapters *before trimming*, while still enabling fast mapping (map once)

6. **Reference preparation**
   - BWA index created if missing
   - Index creation protected by a filesystem lock

7. **Read mapping (pooled; map once)**
   - **modern**
     - PE: `bwa mem` on pooled trimmed R1/R2
     - SE: `bwa mem` on pooled SE reads + pooled PE rescue reads (merged/unpaired)
   - **ancient**
     - fastp merges reads first; pooled merged reads map as SE via `bwa aln` + `bwa samse`
     - pooled SE reads map as SE via `bwa aln` + `bwa samse`
   - MAPQ filter applied during BAM conversion (default MAPQ 20)

8. **Merge and duplicate removal**
   - All pooled mapping outputs merged
   - Duplicate removal using `samtools markdup -r`
   - For modern runs that include PE data, fixmate workflow is used before markdup

9. **Final BAM processing**
   - Coordinate sort and index
   - Read groups added using `samtools addreplacerg`

10. **Optional mapDamage**
    - Runs if `-library-type ancient`

11. **Coverage/depth + summary statistics**
    - Per-contig coverage and mean depth via `samtools coverage`
    - Genome-wide:
      - `avg_depth` = length-weighted mean of per-contig mean depth
      - `pct_covered` = length-weighted percent of bases covered (>=1x)
    - Fragment-aware mapping counts:
      - `mapped_fragments_*` counts templates/fragments  
        (PE counted via mapped read1; SE counted via non-paired reads)

12. **Cleanup**
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

Pilot run (limit processing before fastp; useful for comparing ancient vs modern settings quickly):

    bash plainmap.sh \
      -manifest manifest.txt \
      -prefix SAMPLE1 \
      -ref reference.fa \
      -outdir results \
      -t 16 \
      --pilot-fragments 10000000

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

PlainMap produces a single tab-separated summary file per sample.  
Each row corresponds to one PlainMap run (including pilot runs).

All statistics are **fragment-aware**, meaning paired-end reads are counted as a single sequencing fragment (template) where appropriate. When pilot mode is enabled, all values refer only to the subsampled input.

---

### Column definitions

- **sample**  
  Sample prefix provided with `-prefix`.

- **library_type**  
  Mapping mode used:
  - `modern`: `bwa mem`
  - `ancient`: `bwa aln` + `samse`

- **seq_mode**  
  Sequencing mode detected from the manifest:
  - `SE` (single-end only)
  - `PE` (paired-end present; may include mixed SE + PE inputs)

- **pilot_fragments**  
  Number of raw fragments processed per chunk in pilot mode.  
  A value of `0` indicates pilot mode is disabled.

- **max_reads_per_chunk**  
  Maximum reads per chunk used for adaptive pre-fastp chunking.  
  A value of `0` indicates chunking is disabled.

---

### Raw input statistics

- **raw_R1_reads**  
  Total number of raw R1 reads processed (pilot-aware).

- **raw_R2_reads**  
  Total number of raw R2 reads processed (pilot-aware).  
  This value is `0` for purely single-end datasets.

- **raw_fragments**  
  Total number of raw sequencing fragments (templates), defined as:
  - paired-end data: number of R1 reads
  - single-end data: number of reads

  This value is the denominator for fragment-level statistics and endogenous fraction.

---

### Trimming and merging statistics

- **trimmed_fragments**  
  Number of fragments remaining after fastp filtering and minimum length enforcement.

- **merged_reads**  
  Number of paired-end read pairs successfully merged by fastp.  
  Merged reads are treated as single-end fragments for mapping.

- **unpaired_reads**  
  Number of reads rescued as single-end after paired-end processing  
  (e.g. orphaned reads or failed merges).

- **pct_fragments_remaining**  
  Percentage of raw fragments retained after trimming, calculated as:  
  `trimmed_fragments / raw_fragments × 100`.

- **trimmed_over_raw**  
  Ratio of trimmed fragments to raw fragments.

---

### Mapping statistics (read-level)

- **mapped_reads_all**  
  Total number of mapped reads before duplicate removal.

- **mapped_reads_unique**  
  Number of mapped reads remaining after duplicate removal.

---

### Mapping statistics (fragment-level)

- **mapped_fragments_all**  
  Number of mapped fragments before duplicate removal, defined as:
  - paired-end data: mapped read1 counts as one fragment
  - single-end data: each mapped read counts as one fragment

- **mapped_fragments_unique**  
  Number of unique mapped fragments after duplicate removal.

- **endog**  
  Endogenous fraction, calculated as:  
  `mapped_fragments_unique / raw_fragments`.

- **duprate**  
  Fragment-level duplication rate, calculated from `samtools markdup` duplicate flags as:  
  `duplicate_mapped_fragments / mapped_fragments_all`,  
  where duplicate fragments are those flagged with SAM bit `0x400` (PCR/optical duplicate).


---

### Alignment and coverage statistics

- **avg_readlen**  
  Average length (bp) of mapped reads after duplicate removal.

- **mapped_bp**  
  Total number of mapped base pairs, estimated as:  
  `mapped_reads_all × avg_readlen`.

- **avg_depth**  
  Genome-wide average read depth, calculated as a length-weighted mean  
  of per-contig mean depth reported by `samtools coverage`.

- **pct_covered**  
  Genome-wide percentage of reference bases covered by at least one read (≥1×),  
  calculated as a length-weighted average across contigs.

---

### Notes

- All statistics are pilot-aware when pilot mode is enabled.
- Fragment-level statistics should be preferred over read-level statistics for paired-end and ancient DNA datasets.
- Coverage and depth statistics are computed on the final duplicate-removed BAM.

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
