# PlainMap Toy Example

This directory contains a tiny paired-end dataset for a quick PlainMap smoke test.

Run from this directory:

```bash
bash ../plainmap.sh \
  --manifest manifest.txt \
  --prefix toy \
  --ref reference.fa \
  --outdir output \
  --threads 1 \
  --library-type modern
```

Expected main outputs:

- `output/toy.final.RG.bam`
- `output/toy.final.RG.bam.bai`
- `output/toy.plainmap.stats.tsv`
- `output/toy.plainmap.coverage.tsv`
- `output/toy_plainmap.log`

The example is intentionally minimal and is meant only to verify installation and command-line behavior.
