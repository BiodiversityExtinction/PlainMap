# PlainMap

PlainMap is a transparent, failure-aware mapping pipeline for ancient and modern DNA, designed to simplify high-throughput read mapping for non-expert users while remaining robust, inspectable, and restartable on HPC systems.

PlainMap deliberately avoids workflow engines and black-box abstractions. Every command executed is visible, logged, and reproducible.

## Installation

Clone the repository:

```
git clone https://github.com/yourusername/plainmap.git
cd plainmap
```

Make the script executable:

```
chmod +x plainmap.sh
```

No further installation is required.

## Input Format

### Manifest file

PlainMap uses a manifest file listing input FASTQ files:

- One FASTQ.gz path per line
- Absolute or relative paths supported
- Comments allowed with `#`

Example `manifest.txt`:

```
# Sample Kumamoto
Kumamoto_S1_L001_R1_001.fastq.gz
Kumamoto_S1_L001_R2_001.fastq.gz
Kumamoto_S1_L002_R1_001.fastq.gz
Kumamoto_S1_L002_R2_001.fastq.gz
```

File naming conventions do not matter — read orientation (R1/R2) is detected from FASTQ headers.

## Usage

### Ancient DNA

```
bash plainmap.sh \
  -manifest manifest.txt \
  -prefix SAMPLE1 \
  -ref reference.fa \
  -outdir results \
  -library-type ancient \
  -t 16
```

### Modern DNA

```
bash plainmap.sh \
  -manifest manifest.txt \
  -prefix SAMPLE1 \
  -ref reference.fa \
  -outdir results \
  -library-type modern \
  -t 24
```

## Command-line Parameters

### Required parameters

```
-manifest FILE
```
Text file listing input FASTQ.gz files (one per line).

```
-prefix STRING
```
Sample name / prefix used for output files and read groups.

```
-ref FILE
```
Reference genome FASTA file. BWA index will be created if missing.

```
-outdir DIR
```
Output directory.

---

### Library type

```
-library-type ancient|modern
```
Type of DNA being mapped.

- `ancient`  
  Uses `bwa aln`, merged reads, and single-end duplicate handling.

- `modern`  
  Uses `bwa mem`, paired-end and single-end mapping.

Default:
```
ancient
```

---

### Threading

```
-t, --threads INT
```
Total number of threads available to the pipeline.

Threads are internally balanced between BWA, samtools view, and samtools sort.

Default:
```
1
```

---

### Read filtering and mapping parameters

```
-minlength INT
```
Minimum read length after adapter trimming.

Default:
```
30
```

```
-mismatch FLOAT
```
Mismatch parameter passed to `bwa aln -n` (ancient DNA only).

Default:
```
0.01
```

---

### Chunking large datasets

```
-max-reads-per-chunk INT
```
Maximum number of reads per FASTQ chunk.

- Enables deterministic splitting of large FASTQs
- Prevents memory and wall-time failures
- Allows partial restart from completed chunks

Set to `0` to disable chunking.

Default:
```
0
```

---

### Cleanup behaviour

```
-clean normal|aggressive
```

- `normal`  
  Removes temporary files but keeps intermediate BAMs.

- `aggressive`  
  Removes chunk FASTQs, chunk BAMs, and concatenated FASTQs.

Default:
```
normal
```

---

### Execution control

```
--resume
```
Resume from completed checkpoints if present.

Default:
```
enabled
```

```
--no-resume
```
Ignore checkpoints and rerun all steps from scratch.

```
--dry-run
```
Print commands without executing them. Checkpoints are ignored.

```
--validate
```
Validate tools and input FASTQs, then exit without mapping.

---

### Tool overrides

PlainMap assumes tools are available in `PATH`, but paths can be overridden:

```
--fastp CMD
--bwa CMD
--samtools CMD
```

Example:

```
--bwa /path/to/bwa
```

---

### Help

```
-h, --help
```
Print usage information and exit.

---

## Chunking Large Datasets

For very large datasets, PlainMap can split FASTQs into fixed-size chunks by read count:

```
-max-reads-per-chunk 300000000
```

Chunking is deterministic:

- no random sampling
- no read duplication across chunks
- safe restart from partially completed runs
- global deduplication after merging

## Restarting After Failure

PlainMap automatically creates checkpoint files during execution.

If a job:

- runs out of wall time
- is cancelled
- crashes due to node failure

You can simply re-run the same command:

```
bash plainmap.sh [same arguments]
```

Completed steps will be skipped automatically.

To force a full re-run:

```
--no-resume
```

## Validation Mode

To check inputs and dependencies without running mapping:

```
bash plainmap.sh \
  -manifest manifest.txt \
  -prefix SAMPLE1 \
  -ref reference.fa \
  -outdir results \
  --validate
```

Validation checks:

- tool availability
- FASTQ presence
- gzip integrity (detects truncated downloads)

## Output

### Main outputs

Final BAM with read groups:

```
SAMPLE.final.RG.bam
SAMPLE.final.RG.bam.bai
```

Log file:

```
SAMPLE_plainmap.log
```

Log contents include:

- executed commands
- per-step wall times
- total wall time
- tool versions
- SLURM job ID (if available)

## Duplicate Handling

- Duplicates are removed, not just marked
- Global duplicate removal is performed after merging all chunks
- Ancient DNA uses single-end duplicate logic
- Modern DNA handles paired-end and single-end reads appropriately

## Read Groups

Read groups are added using:

```
samtools addreplacerg
```

This avoids Java and Picard dependencies and ensures compatibility with downstream tools such as FreeBayes, bcftools, and GATK.

## Cleanup Options

Default cleanup (`normal`) removes temporary files but keeps intermediate BAMs.

To aggressively reduce disk usage:

```
-clean aggressive
```

This removes:

- chunk FASTQs
- chunk BAMs
- concatenated FASTQs

## Limitations (v0.1)

PlainMap v0.1 is intentionally conservative.

Current limitations:

- Single read group per sample
- No built-in QC plots
- No contamination estimation
- No variant calling

These are design choices, not oversights.

## When to Use PlainMap

PlainMap is ideal if you:

- want full visibility into mapping steps
- need robust restartability on HPC systems
- are working with ancient or historical DNA
- are teaching or learning mapping workflows
- want a simple, auditable reference pipeline

If you need a full end-to-end framework (damage profiling, contamination, imputation), consider workflow-based pipelines such as nf-core/eager.

## Citation

If you use PlainMap, please cite:

Author(s). PlainMap: a transparent, failure-aware mapping pipeline for ancient and modern DNA. Journal, Year.

A DOI will be added upon publication.

## Development Status

PlainMap is currently v0.1 and under active development.

Feedback, issues, and pull requests are welcome.
