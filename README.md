# PlainMap

PlainMap is a transparent, failure-aware mapping pipeline for ancient and modern DNA, designed to simplify high-throughput read mapping for non-expert users while remaining robust, inspectable, and restartable on HPC systems.

PlainMap deliberately avoids workflow engines and black-box abstractions. Every command executed is visible, logged, and reproducible.

---

## Key Features

- Transparent execution  
  All mapping commands are explicit and logged — no hidden workflow logic.

- Failure-aware and restartable  
  Jobs can resume from intermediate steps after cancellation, wall-time limits, or node failures.

- Deterministic chunking  
  Very large FASTQ files can be split into fixed-size chunks without random subsampling or read duplication.

- Ancient and modern DNA support  
  - Ancient DNA: `bwa aln` with merged reads and single-end duplicate handling  
  - Modern DNA: `bwa mem` with paired-end and single-end handling

- Minimal dependencies  
  No Java, no containers, no workflow DSLs.

- Robust FASTQ validation  
  Detects truncated or corrupted gzip files and halts safely.

- HPC-friendly  
  Designed for SLURM-style environments and parallel job submission.

---

## Design Philosophy

PlainMap is intentionally plain.

Its goal is not to be a feature-complete end-to-end framework, but a reference implementation of best-practice read mapping that is:

- easy to understand  
- easy to debug  
- easy to audit  
- easy to teach  

This makes PlainMap particularly suitable for:
- MSc and PhD students
- wet-lab researchers
- method development and bias testing
- historical and ancient DNA datasets with non-standard properties

---

## Requirements

PlainMap requires the following tools to be available in your `PATH`:

- `bash` (not `sh`)
- `gzip`
- GNU `split` (with `--filter` support; part of coreutils)
- `fastp`
- `bwa`
- `samtools`

Note: PlainMap uses `split --filter`, which requires GNU coreutils.  
On macOS, install via Homebrew:

```bash
brew install coreutils
