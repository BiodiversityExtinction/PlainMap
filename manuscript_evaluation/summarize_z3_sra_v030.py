#!/usr/bin/env python3
"""Collect validated SRA-derived Z3 evaluations into one TSV."""

import csv
import pathlib
import sys


def read_row(path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 1:
        raise SystemExit(f"expected one statistics row: {path}")
    return rows[0]


def read_pairs(path):
    with path.open(newline="") as stream:
        return dict(csv.reader(stream, delimiter="\t"))


def collect(root):
    cases = (
        ("mode_modern", "modern", 0, "modes/modern", "Z3_28m_modern", "time.tsv", "full"),
        ("mode_ancient_fast", "ancient-fast", 0, "modes/ancient-fast", "Z3_28m_ancient_fast", "time.tsv", "full"),
        ("mode_ancient_all", "ancient-all", 0, "modes/ancient-all", "Z3_28m_ancient_all", "time.tsv", "full"),
        ("chunk_unchunked", "modern", 0, "chunking_baseline/unchunked", "Z3_large_modern", "time.tsv", "full"),
        ("chunk_1m", "modern", 1000000, "chunking_baseline/chunked", "Z3_large_modern", "time.tsv", "full"),
        ("chunk_5m", "modern", 5000000, "chunk_size_sweep/size_5000000", "Z3_large_modern", "time.tsv", "full"),
        ("chunk_10m", "modern", 10000000, "chunk_size_sweep/size_10000000", "Z3_large_modern", "time.tsv", "full"),
        ("restart_resume", "modern", 1000000, "restart/resumed", "Z3_large_modern", "../resume_time.tsv", "remaining_only"),
    )
    provenance = read_pairs(root / "provenance.tsv")
    checks = read_pairs(root / "restart/restart_checks.tsv")
    for key in ("resumed_vs_uninterrupted_alignments", "resumed_vs_uninterrupted_statistics", "final_bam_quickcheck"):
        if checks.get(key) != "PASS":
            raise SystemExit(f"restart check failed: {key}")

    fields = (
        "test", "library_type", "max_reads_per_chunk", "node", "raw_fragments",
        "trimmed_fragments", "merged_reads", "unpaired_reads", "mapped_fragments_all",
        "mapped_fragments_unique", "duplicate_fragments", "duprate", "wall_seconds",
        "user_cpu_seconds", "system_cpu_seconds", "max_process_rss_kb", "time_scope",
    )
    rows = []
    for name, mode, chunk_size, subdir, prefix, time_name, scope in cases:
        directory = root / subdir
        stats = read_row(directory / f"{prefix}.plainmap.stats.tsv")
        timing = read_pairs(directory / time_name)
        if stats["library_type"] != mode or int(stats["max_reads_per_chunk"]) != chunk_size:
            raise SystemExit(f"mode or chunk-size mismatch: {name}")
        if timing.get("exit_status") != "0":
            raise SystemExit(f"nonzero exit status: {name}")
        if int(stats["raw_fragments"]) != 27943970:
            raise SystemExit(f"wrong input fragment count: {name}")
        row = {
            "test": name,
            "library_type": mode,
            "max_reads_per_chunk": chunk_size,
            "node": provenance["node"],
            "time_scope": scope,
        }
        for key in fields:
            if key in stats:
                row[key] = stats[key]
            elif key in timing:
                row[key] = timing[key]
        row["max_process_rss_kb"] = timing.get("max_process_rss_kb", timing.get("max_rss_kb"))
        rows.append(row)

    if rows[0]["mapped_fragments_unique"] != rows[3]["mapped_fragments_unique"]:
        raise SystemExit("modern mode and unchunked baseline disagree")
    if rows[4]["mapped_fragments_unique"] != rows[7]["mapped_fragments_unique"]:
        raise SystemExit("resumed and uninterrupted chunked counts disagree")

    output = root / "results.tsv"
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(output)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: summarize_z3_sra_v030.py EVAL_ROOT")
    collect(pathlib.Path(sys.argv[1]))
