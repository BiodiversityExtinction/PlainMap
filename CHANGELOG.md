# Changelog

All notable changes to PlainMap will be documented here.

## v0.1.0 - 2026-07-07

Initial public release.

### Added

- Modern, ancient-fast, and ancient-all mapping modes.
- Header-based single-end and paired-end FASTQ classification.
- SRA-export paired-end handling for FASTQs with identical R1/R2 headers.
- Optional pilot subsampling before trimming.
- Optional pre-fastp and post-fastp chunking for large datasets.
- Restartable per-chunk mapping.
- Optional unmapped-read FASTQ export.
- Optional mapDamage integration.
- Fragment-aware summary statistics.
- Toy example dataset for smoke testing.
- Citation metadata and MIT license.

### Changed

- Renamed mapping modes to `ancient-fast` and `ancient-all`.
- Final BAM generation explicitly excludes unmapped reads while preserving optional unmapped FASTQ export behavior.

### Fixed

- Improved handling of SRA paired-end FASTQ files with matching read names across R1 and R2.
- Added retry logic around merged BAM validation.
- Added strict resume configuration tracking. PlainMap now refuses to reuse outputs or checkpoints if any effective parameter, input identity, script, or tool command has changed, and reports the differences.
- Added the exact duplicate-fragment count to summary statistics for reproducibility checks.
