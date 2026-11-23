**Table X. PlainMap evaluation using heterogeneous Japanese sea lion reads and short-read cave hyena data.**

| Experiment | Condition | Chunk size | Trimmed fragments | Merged fragments | Unique mapped fragments | Wall time (min) | Peak RSS (GiB) |
|---|---|---:|---:|---:|---:|---:|---:|
| Z3 library modes | modern | None | 27,928,803 | 15,036,758 | 70,348 | 9.23 | 5.94 |
| Z3 library modes | ancient-fast | None | 27,928,803 | 15,036,758 | 39,983 | 24.42 | 5.37 |
| Z3 library modes | ancient-all | None | 27,928,803 | 15,036,758 | 40,837 | 27.08 | 5.85 |
| Z3 chunking | modern | 1 million | 27,928,803 | 15,036,758 | 71,006 | 16.29 | 4.85 |
| Z3 chunking | modern | 5 million | 27,928,803 | 15,036,758 | 70,341 | 12.43 | 5.85 |
| Z3 chunking | modern | 10 million | 27,928,803 | 15,036,758 | 70,328 | 11.89 | 5.87 |
| Z3 chunking | ancient-fast | 1 million | 27,928,803 | 15,036,758 | 39,983 | 50.68 | 4.81 |
| Z3 chunking | ancient-fast | 10 million | 27,928,803 | 15,036,758 | 39,983 | 43.20 | 4.93 |
| Z3 chunking | ancient-all | 1 million | 27,928,803 | 15,036,758 | 40,863 | 55.74 | 4.97 |
| Z3 chunking | ancient-all | 10 million | 27,928,803 | 15,036,758 | 40,837 | 48.02 | 5.83 |
| Z3 restart | modern, resumed | 1 million | 27,928,803 | 15,036,758 | 71,006 | 8.32* | 4.85* |
| Cave hyena | AdaptClean off | None | 27,966,526 | 0 | 15,146,740 | 319.38 | 5.97 |
| Cave hyena | AdaptClean on | None | 27,966,526 | 0 | 15,156,531 | 442.18 | 5.93 |

**Legend.** Z3 comprised 27,943,970 raw fragments from one single-end and two
paired-end sequencing units, mapped to the California sea lion reference
(GCF_009762305.2). SRA-derived `SRR25607793` mates had identical headers and
were correctly inferred as paired-end. Cave hyena sample Ccsp043 comprised
50 million single-end pilot fragments mapped to the striped hyena reference
(GCF_003009895.1); its zero merged-fragment count reflects single-end input.
All runs used 16 threads, a 30-bp minimum retained length, MAPQ 20, and hybrid
duplicate removal. Unique mapped fragments were counted after MAPQ filtering
and duplicate removal. Peak RSS is maximum process memory, not total job memory.
The unchunked Z3 rows used the immutable `v0.3.0` tag; direct-chunk and restart
runs used PlainMap commit `f637f70` (the script still reports `0.3.0`). The
cave hyena runs used a frozen v0.3.0-equivalent snapshot. Unchunked Z3 runs
used `compute02`, modern direct-chunk runs used `compute04`, and ancient-mode
direct-chunk runs used `node12`. Timings are descriptive and should not be
compared directly across those groups. Mapped counts do not establish alignment
accuracy or biological endogeneity. *The restart reused three of 46 completed
mapping chunks and reproduced the uninterrupted one-million-fragment run's BAM
records and statistics exactly. Its time and RSS describe only the resumed
attempt, not the interrupted-plus-resumed total.*

Source data: [Z3 unchunked modes](z3_sra_v030_results.tsv),
[modern direct chunking](z3_modern_direct_chunk_results.tsv),
[ancient direct chunking](z3_ancient_direct_chunk_results.tsv),
[modern restart checks](z3_modern_direct_chunk_restart_checks.tsv), and
[cave hyena](cavehyena_50m_adaptclean_results.tsv).
