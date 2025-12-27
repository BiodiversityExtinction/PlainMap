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

PlainMap requires the following tools to be available in your `PATH`:

- `bash` (must be bash, not sh)
- `gzip`
- GNU `split` (with `--filter` support; part of coreutils)
- `fastp`
- `bwa`
- `samtools`

Notes:

- GNU `split` is required for deterministic FASTQ chunking.
  On macOS, install via Homebrew:
