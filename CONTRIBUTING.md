# Contributing to PlainMap

PlainMap is intended to stay lightweight, transparent in behavior, and easy to run on HPC systems. Contributions that improve robustness, documentation, portability, or reviewer/user testability are welcome.

## Reporting Bugs

Please open a GitHub Issue and include:

- the PlainMap command you ran
- the PlainMap log file, or the relevant error section
- the `plainmap.sh --help` version if available
- operating system or cluster environment
- versions of `fastp`, `bwa`, `samtools`, and `bash`
- whether the input is single-end, paired-end, mixed, or SRA-export FASTQ
- a minimal reproducible example if possible

Avoid uploading private sequencing data. If a bug depends on FASTQ headers, a few redacted header lines are usually enough.

## Feature Requests

Please open an Issue before large feature work. PlainMap deliberately has a narrow scope: preprocessing, read mapping, duplicate handling, final BAM generation, optional unmapped-read export, optional mapDamage, and summary statistics.

Features are more likely to fit if they:

- improve robustness for heterogeneous FASTQ inputs
- help compare mapping strategies reproducibly
- preserve restartable execution
- avoid adding heavy configuration systems
- keep downstream analysis choices outside PlainMap

## Pull Requests

Before opening a pull request:

```bash
bash -n plainmap.sh
gzip -t example/reads_R1.fastq.gz example/reads_R2.fastq.gz
cd example
bash ../plainmap.sh --manifest manifest.txt --prefix toy --ref reference.fa --outdir output --threads 1 --library-type modern
```

If your change affects outputs, please describe the expected behavior change and update the README or example documentation where relevant.

## Coding Style

- Keep the script Bash-only unless there is a strong reason otherwise.
- Prefer explicit validation and clear error messages.
- Preserve existing command-line behavior unless the change is intentionally documented.
- Keep optional dependencies optional.
- Do not assume FASTQ filename conventions; PlainMap should infer pairing from headers and validated manifest structure.
