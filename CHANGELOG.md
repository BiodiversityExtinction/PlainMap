# Changelog

All notable changes to PlainMap will be documented here.

## v0.2.1 - 2026-10-07

### Added

- Added `--dedup-method hybrid|samtools`. The default `hybrid` policy retains v0.2.0 behavior, while `samtools` reproduces the combined-BAM duplicate-removal sequence used by v0.1.0.
- Recorded the selected duplicate-removal policy in run configuration, logs, and summary statistics.
- Added regression coverage for both duplicate-removal policies in all mapping modes and for resume rejection when the policy changes.

### Changed

- DeDup is required only for `--dedup-method hybrid`; samtools-only runs no longer require DeDup to be installed.

## v0.2.0 - 2026-09-24

### Added

- mapDamage is now automatically parallelized by splitting the final BAM by whole contig when `--threads` > 1 and more than one contig has mapped reads, running independent mapDamage scans in parallel and merging the resulting tables before a final plotting pass. mapDamage's own code is never modified, and the merged count tables are numerically equivalent to a single-process run. Added `--mapdamage-serial` to force the original single-process behavior.
- Added DeDup as a required dependency and `--dedup CMD` as an executable-path override.
- Added separate PE and SE-like pre/post-deduplication fragment counts, duplicate counts, duplication rates, and method fields to summary statistics.
- Added persistent DeDup log, histogram, JSON, and stdout reports under `SAMPLE/dedup_reports/`.

### Changed

- Duplicate removal now uses DeDup with `--merged` for all SE-like alignments (merged, unpaired, and true SE) and samtools `fixmate + markdup` for unmerged PE alignments. These classes are deduplicated separately and recombined before downstream processing.
- Updated the integration tests to exercise a known SE duplicate with the real DeDup interface and verify class-specific statistics.

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
- Directed fastp HTML reports to the sample report directory instead of the launch directory.
- Added the reproducible manuscript evaluation manifest, commands, validation protocol, and reported results.
