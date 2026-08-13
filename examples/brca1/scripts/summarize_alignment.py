#!/usr/bin/env python3
"""Summarize the fixed BRCA1 alignment and trimAl sensitivity profiles."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def fail(message: str) -> None:
    raise SystemExit(message)


def read_alignment(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []

    def store() -> None:
        if name is None:
            return
        sequence = "".join(chunks).upper()
        if not sequence:
            fail(f"empty aligned sequence in {path}: {name}")
        if name in records:
            fail(f"duplicate tip ID in {path}: {name}")
        records[name] = sequence

    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                store()
                name = line[1:].split()[0]
                chunks = []
            elif name is None:
                fail(f"sequence before first FASTA header: {path}")
            else:
                chunks.append(line)
    store()
    if len(records) < 4:
        fail(f"alignment has fewer than four sequences: {path}")
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) != 1:
        fail(f"aligned sequences have unequal lengths: {path}")
    return records


def parse_alignment(value: str) -> tuple[str, Path]:
    if "=" not in value:
        fail("--alignment values must use LABEL=PATH")
    label, raw_path = value.split("=", 1)
    if not label or not raw_path:
        fail("--alignment values must use non-empty LABEL=PATH")
    return label, Path(raw_path)


def column_stats(records: dict[str, str]) -> tuple[float, float, float, int, int]:
    sequences = list(records.values())
    n_sequences = len(sequences)
    n_columns = len(sequences[0])
    occupancies: list[float] = []
    invariant = 0
    informative = 0
    for column_index in range(n_columns):
        residues = [sequence[column_index] for sequence in sequences]
        nongaps = [residue for residue in residues if residue not in {"-", "."}]
        occupancies.append(len(nongaps) / n_sequences)
        counts: dict[str, int] = {}
        for residue in nongaps:
            counts[residue] = counts.get(residue, 0) + 1
        if len(counts) == 1 and nongaps:
            invariant += 1
        if sum(1 for count in counts.values() if count >= 2) >= 2:
            informative += 1
    return (
        sum(occupancies) / n_columns,
        min(occupancies),
        max(occupancies),
        invariant,
        informative,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", action="append", required=True, type=parse_alignment)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    if args.output_dir.exists():
        fail(f"refusing to overwrite output directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True)

    loaded: list[tuple[str, Path, dict[str, str]]] = []
    labels: set[str] = set()
    for label, path in args.alignment:
        if label in labels:
            fail(f"duplicate alignment label: {label}")
        labels.add(label)
        loaded.append((label, path, read_alignment(path)))

    raw_label, _raw_path, raw_records = loaded[0]
    if raw_label != "raw":
        fail("the first alignment must be labeled raw")
    raw_tips = set(raw_records)
    raw_columns = len(next(iter(raw_records.values())))

    alignment_rows: list[dict[str, str]] = []
    sequence_rows: list[dict[str, str]] = []
    for label, path, records in loaded:
        if set(records) != raw_tips:
            missing = sorted(raw_tips - set(records))
            extra = sorted(set(records) - raw_tips)
            fail(f"tip-set mismatch for {label}; missing={missing}; extra={extra}")
        columns = len(next(iter(records.values())))
        mean_occupancy, min_occupancy, max_occupancy, invariant, informative = column_stats(records)
        alignment_rows.append({
            "alignment_id": label,
            "path": path.as_posix(),
            "sequences": str(len(records)),
            "columns": str(columns),
            "retained_vs_raw": f"{columns / raw_columns:.6f}",
            "mean_occupancy": f"{mean_occupancy:.6f}",
            "min_occupancy": f"{min_occupancy:.6f}",
            "max_occupancy": f"{max_occupancy:.6f}",
            "invariant_columns": str(invariant),
            "parsimony_informative_columns": str(informative),
        })
        for tip_id, sequence in records.items():
            ungapped = sequence.replace("-", "").replace(".", "")
            leading = len(sequence) - len(sequence.lstrip("-."))
            trailing = len(sequence) - len(sequence.rstrip("-."))
            sequence_rows.append({
                "alignment_id": label,
                "tip_id": tip_id,
                "ungapped_length": str(len(ungapped)),
                "gap_fraction": f"{1 - len(ungapped) / columns:.6f}",
                "leading_gaps": str(leading),
                "trailing_gaps": str(trailing),
            })

    with (args.output_dir / "alignment_qc.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(alignment_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(alignment_rows)
    with (args.output_dir / "sequence_qc.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(sequence_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sequence_rows)


if __name__ == "__main__":
    main()
