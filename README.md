# PlainMap

PlainMap is a transparent, failure-aware mapping pipeline for **ancient, historical, and modern DNA** designed to **just run** on HPC systems with minimal fuss.

PlainMap is implemented as a **single Bash script** and is intentionally conservative: when inputs are ambiguous, it fails loudly rather than guessing.

Current version: **v0.1 (active development)**

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
- Automatic use of pigz if available
- Fragment-aware duplication rate derived from samtools markdup
- Optional mapDamage for any library type

---

## Installation

Clone the repository:

    git clone https://github.com/yourusername/plainmap.git
    cd plainmap

Make the script executable:

    chmod +x plainmap.sh

PlainMap is a single Bash script. No compilation required.

---

## Dependencies

All tools must be available on PATH unless overridden via flags.

### Required

- bash (must be bash, not sh)
- fastp
- bwa
- samtools
- python3
- GNU split (must support --filter)
- gzip and zcat

### Optional (recommended)

- pigz (used automatically if present for parallel compression and integrity tests)

### Optional (only if enabled)

- mapDamage (required only if --run-mapdamage is used)

---

## Input Format

### Manifest file

- One FASTQ.gz path per line
- Absolute or relative paths supported
- Relative paths resolved relative to the manifest directory
- Comments allowed with #
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
   - ancient: bwa aln + samse (SE-like only)
   - historical: bwa aln + sampe/samse
10. Merge of mapping outputs
11. Duplicate marking and filtering (samtools markdup)
12. Final sorting, indexing, and read group assignment
13. Optional mapDamage
14. Coverage and depth statistics
15. Fragment-aware summary statistics
16. Cleanup (optional retention of intermediates)

---

## Usage Examples

Modern DNA:

    bash plainmap.sh \
      --manifest manifest.txt \
      --prefix SAMPLE1 \
      --ref reference.fa \
      --outdir results \
      --threads 16

Ancient DNA:

    bash plainmap.sh \
      --manifest manifest.txt \
      --prefix SAMPLE1 \
      --ref reference.fa \
      --outdir results \
      --library-type ancient \
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

---

## Output

Main outputs:

- SAMPLE.final.RG.bam
- SAMPLE.final.RG.bam.bai
- SAMPLE.plainmap.stats.tsv
- SAMPLE.plainmap.coverage.tsv
- SAMPLE_plainmap.log

Optional outputs:

- SAMPLE_mapdamage/ (if enabled)

---

## Summary Statistics

All statistics are **fragment-aware**.

Paired-end reads are counted as one fragment where appropriate.
Pilot mode affects all counts.

Reported fields include:

- raw reads and fragments
- trimmed, merged, and unpaired counts
- mapped reads and fragments (pre/post deduplication)
- endogenous fraction
- duplication rate from duplicate flags
- average read length
- mapped base pairs
- genome-wide average depth
- genome-wide percent covered

---

## Limitations (v0.1)

- Single read group per sample
- No QC plotting
- No contamination estimation
- No variant calling

---

## Development Status

PlainMap is **v0.1** and under active development.

Issues, feedback, and pull requests are welcome.
