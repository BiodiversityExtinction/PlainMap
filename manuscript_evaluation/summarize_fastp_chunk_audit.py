#!/usr/bin/env python3
"""Compare fragment fates and sequences across a controlled fastp chunking audit."""

import collections
import gzip
import hashlib
import json
import pathlib
import sys


def reads(path):
    with gzip.open(path, "rt") as stream:
        while header := stream.readline():
            sequence = stream.readline().rstrip("\n\r")
            plus = stream.readline()
            quality = stream.readline().rstrip("\n\r")
            if not plus or len(quality) != len(sequence):
                raise ValueError(f"Invalid FASTQ record in {path}")
            name = header.split()[0].removeprefix("@").removesuffix("/1").removesuffix("/2")
            yield name, sequence


def fragment_fates(root, label):
    raw = root / label / "Z3_sra_chunk_audit" / "work" / "raw"
    prefix = "Z3_sra_chunk_audit.pool"
    r1 = raw / f"{prefix}.PE.R1.fastq.gz"
    r2 = raw / f"{prefix}.PE.R2.fastq.gz"
    se = raw / f"{prefix}.SE.all.fastq.gz"
    fates = {}
    mate_names = set()
    for name, sequence in reads(r1):
        if name in fates:
            raise ValueError(f"Duplicate R1 fragment ID in {label}: {name}")
        fates[name] = ("PE", len(sequence), hashlib.sha1(sequence.encode()).hexdigest())
    for name, _ in reads(r2):
        if name in mate_names:
            raise ValueError(f"Duplicate R2 fragment ID in {label}: {name}")
        mate_names.add(name)
    if mate_names != set(fates):
        raise ValueError(f"R1/R2 output mate IDs differ in {label}")
    for name, sequence in reads(se):
        if name in fates:
            raise ValueError(f"Fragment appears in both PE and SE-like output in {label}: {name}")
        fates[name] = ("SE-like", len(sequence), hashlib.sha1(sequence.encode()).hexdigest())
    return fates


def adapter_counts(root, label):
    reports = root / label / "Z3_sra_chunk_audit" / "fastp_reports"
    counts = collections.Counter()
    for path in reports.glob("*.json"):
        data = json.loads(path.read_text())
        counts[data["adapter_cutting"].get("read1_adapter_sequence", "<none>")] += 1
    return counts


def main(root):
    labels = ("auto_unchunked", "auto_chunked", "fixed_unchunked", "fixed_chunked")
    fates = {label: fragment_fates(root, label) for label in labels}
    with (root / "fragment_fates.tsv").open("w") as out:
        out.write("comparison\tfragment_id\tunchunked_class\tchunked_class\tunchunked_length\tchunked_length\tsequence_differs\n")
        for prefix in ("auto", "fixed"):
            unchunked = fates[f"{prefix}_unchunked"]
            chunked = fates[f"{prefix}_chunked"]
            summary = collections.Counter()
            for name in sorted(unchunked.keys() | chunked.keys()):
                before = unchunked.get(name)
                after = chunked.get(name)
                if before == after:
                    continue
                if before is None:
                    summary["only_chunked"] += 1
                elif after is None:
                    summary["only_unchunked"] += 1
                elif before[0] != after[0]:
                    summary["class_changed"] += 1
                else:
                    summary["sequence_or_length_changed"] += 1
                out.write(
                    f"{prefix}\t{name}\t{before[0] if before else 'absent'}\t"
                    f"{after[0] if after else 'absent'}\t"
                    f"{before[1] if before else 'NA'}\t{after[1] if after else 'NA'}\t"
                    f"{int(bool(before and after and before[2] != after[2]))}\n"
                )
            summary["unchunked_fragments"] = len(unchunked)
            summary["chunked_fragments"] = len(chunked)
            for key in ("unchunked_fragments", "chunked_fragments", "only_unchunked", "only_chunked", "class_changed", "sequence_or_length_changed"):
                print(f"{prefix}\t{key}\t{summary[key]}")
    with (root / "audit_summary.tsv").open("w") as out:
        out.write("case\tadapter_r1\treports\tfragments\n")
        for label in labels:
            for adapter, count in sorted(adapter_counts(root, label).items()):
                out.write(f"{label}\t{adapter}\t{count}\t{len(fates[label])}\n")


if __name__ == "__main__":
    main(pathlib.Path(sys.argv[1]))
