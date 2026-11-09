#!/usr/bin/env python3
"""Combine validated PlainMap mode runs into one manuscript-ready TSV."""

import csv
import pathlib
import sys


def read_stats(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected one statistics row in {path}")
    return rows[0]


def read_times(path):
    with path.open() as stream:
        return dict(line.rstrip("\n").split("\t", 1) for line in stream)


def read_provenance(path):
    with path.open() as stream:
        return dict(line.rstrip("\n").split("\t", 1) for line in stream)


root = pathlib.Path(sys.argv[1])
columns = (
    "library_type", "raw_fragments", "trimmed_fragments", "merged_reads",
    "unpaired_reads", "mapped_fragments_all", "mapped_fragments_unique",
    "duplicate_fragments", "pe_duplicate_fragments",
    "se_like_duplicate_fragments", "endog", "duprate",
    "wall_seconds", "user_cpu_seconds", "system_cpu_seconds",
    "max_process_rss_kb",
)
rows = []
shared_provenance = None
shared_counts = None
for mode in ("modern", "ancient-fast", "ancient-all"):
    prefix = f"Z3_28m_{mode.replace('-', '_')}"
    out = root / mode
    stats = read_stats(out / f"{prefix}.plainmap.stats.tsv")
    times = read_times(out / "time.tsv")
    provenance = read_provenance(out / "provenance.tsv")
    if stats["library_type"] != mode or times["exit_status"] != "0":
        raise ValueError(f"Invalid or mismatched run: {mode}")
    if provenance["mode"] != mode:
        raise ValueError(f"Provenance mode mismatch: {mode}")
    common = tuple(provenance[key] for key in (
        "plainmap_sha256", "manifest_sha256", "reference", "threads",
        "pilot_fragments_per_unit", "minlength", "mapq",
        "max_reads_per_chunk", "dedup_method", "fastp_version",
        "bwa_version", "samtools_version", "dedup_version",
    ))
    counts = tuple(stats[key] for key in (
        "raw_fragments", "trimmed_fragments", "merged_reads", "unpaired_reads",
    ))
    if shared_provenance is None:
        shared_provenance, shared_counts = common, counts
    elif common != shared_provenance or counts != shared_counts:
        raise ValueError(f"Input or configuration differs across modes: {mode}")
    rows.append({key: (stats if key in stats else times)[key] for key in columns})

with (root / "summary.tsv").open("w") as stream:
    writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
