#!/usr/bin/env python3
"""Audit whether trimAl profiles retain annotated BRCA1 domains.

The script maps each trimmed alignment back to the raw alignment by comparing
complete alignment-column vectors across the exact same tip set.  A mapping is
accepted only when it is the *unique* in-order column-subset mapping.  This is
deliberately stricter than matching each sequence independently: when repeated
column vectors make the provenance of a retained column ambiguous, the audit
fails instead of choosing an arbitrary raw column.

Domain coordinates come from ``candidate_qc.tsv`` and are 1-based inclusive
coordinates on each accession's ungapped protein sequence.  The output has one
row per accession, trim profile, domain type, and interval.  Validation errors
produce no report.  A completed report containing any interval below the
retention threshold is written with ``status=fail`` and the process exits 2.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


COORDINATE_PATTERN = re.compile(r"([1-9][0-9]*)-([1-9][0-9]*)")
REQUIRED_QC_COLUMNS = frozenset(
    {
        "accession",
        "subject_length",
        "ring_coordinates",
        "brct_coordinates",
        "domain_qc_status",
    }
)
OUTPUT_COLUMNS = (
    "accession",
    "profile",
    "domain",
    "interval_index",
    "residue_start",
    "residue_end",
    "raw_msa_start",
    "raw_msa_end",
    "total_residues",
    "retained_residues",
    "retained_fraction",
    "minimum_retained_fraction",
    "status",
)


class AuditError(RuntimeError):
    """Raised for fail-closed input, mapping, and domain validation errors."""


@dataclass(frozen=True)
class DomainInterval:
    domain: str
    index: int
    start: int
    end: int


@dataclass(frozen=True)
class CandidateDomains:
    accession: str
    subject_length: int
    intervals: tuple[DomainInterval, ...]


def fail(message: str) -> None:
    raise AuditError(message)


def parse_profile(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("profiles must use LABEL=PATH")
    label, raw_path = value.split("=", 1)
    if not label or not raw_path:
        raise argparse.ArgumentTypeError("profiles must use non-empty LABEL=PATH")
    if label == "raw":
        raise argparse.ArgumentTypeError("the reserved profile label 'raw' is not allowed")
    if any(character in label for character in "\t\r\n"):
        raise argparse.ArgumentTypeError("profile labels cannot contain tabs or newlines")
    return label, Path(raw_path)


def parse_fraction(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a number between 0 and 1") from exc
    if not math.isfinite(parsed) or not 0 <= parsed <= 1:
        raise argparse.ArgumentTypeError("must be a finite number between 0 and 1")
    return parsed


def read_alignment(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    chunks: list[str] = []

    def store() -> None:
        if identifier is None:
            return
        sequence = "".join(chunks).upper().replace(".", "-")
        if not sequence:
            fail(f"empty aligned sequence in {path}: {identifier}")
        invalid = sorted(set(sequence) - set("ABCDEFGHIJKLMNOPQRSTUVWXYZ-"))
        if invalid:
            fail(
                f"invalid aligned-sequence character(s) in {path} for "
                f"{identifier}: {''.join(invalid)}"
            )
        if identifier in records:
            fail(f"duplicate FASTA identifier in {path}: {identifier}")
        records[identifier] = sequence

    try:
        with path.open(encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    store()
                    header = line[1:].strip()
                    identifier = header.split()[0] if header else None
                    if identifier is None:
                        fail(f"empty FASTA identifier in {path}")
                    chunks = []
                elif identifier is None:
                    fail(f"sequence data precedes the first FASTA header in {path}")
                else:
                    chunks.append(line)
        store()
    except (OSError, UnicodeDecodeError) as exc:
        fail(f"cannot read alignment {path}: {exc}")

    if not records:
        fail(f"alignment contains no records: {path}")
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) != 1:
        fail(f"aligned sequences have unequal lengths in {path}")
    tip_order = sorted(records)
    for column_index in range(next(iter(lengths))):
        if all(records[tip][column_index] == "-" for tip in tip_order):
            fail(f"all-gap column {column_index + 1} in alignment {path}")
    return records


def parse_coordinate_field(
    value: str, *, accession: str, domain: str, subject_length: int
) -> tuple[DomainInterval, ...]:
    if not value:
        fail(f"{accession}: empty {domain}_coordinates")
    intervals: list[DomainInterval] = []
    previous_end = 0
    for index, token in enumerate(value.split(";"), start=1):
        stripped = token.strip()
        match = COORDINATE_PATTERN.fullmatch(stripped)
        if match is None:
            fail(
                f"{accession}: malformed {domain} interval {stripped!r}; "
                "expected 1-based START-END"
            )
        start, end = (int(part) for part in match.groups())
        if start > end:
            fail(f"{accession}: reversed {domain} interval {start}-{end}")
        if end > subject_length:
            fail(
                f"{accession}: {domain} interval {start}-{end} exceeds "
                f"subject_length {subject_length}"
            )
        if start <= previous_end:
            fail(
                f"{accession}: {domain} intervals must be strictly ordered and "
                f"non-overlapping; interval {start}-{end} follows residue {previous_end}"
            )
        previous_end = end
        intervals.append(DomainInterval(domain, index, start, end))
    return tuple(intervals)


def read_candidate_domains(path: Path) -> list[CandidateDomains]:
    candidates: list[CandidateDomains] = []
    seen: set[str] = set()
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                fail(f"candidate QC table has no header: {path}")
            missing = sorted(REQUIRED_QC_COLUMNS - set(reader.fieldnames))
            if missing:
                fail("candidate QC table is missing columns: " + ", ".join(missing))
            for line_number, raw_row in enumerate(reader, start=2):
                if None in raw_row:
                    fail(f"candidate QC row {line_number} has extra columns")
                accession = raw_row["accession"]
                if not accession or accession != accession.strip():
                    fail(f"candidate QC row {line_number} has an invalid accession")
                if accession in seen:
                    fail(f"duplicate candidate QC accession: {accession}")
                seen.add(accession)
                if raw_row["domain_qc_status"].strip().lower() != "pass":
                    fail(
                        f"{accession}: domain_qc_status is not pass: "
                        f"{raw_row['domain_qc_status']!r}"
                    )
                try:
                    subject_length = int(raw_row["subject_length"])
                except ValueError as exc:
                    raise AuditError(
                        f"{accession}: subject_length must be a positive integer"
                    ) from exc
                if subject_length <= 0:
                    fail(f"{accession}: subject_length must be a positive integer")
                ring = parse_coordinate_field(
                    raw_row["ring_coordinates"],
                    accession=accession,
                    domain="ring",
                    subject_length=subject_length,
                )
                brct = parse_coordinate_field(
                    raw_row["brct_coordinates"],
                    accession=accession,
                    domain="brct",
                    subject_length=subject_length,
                )
                candidates.append(CandidateDomains(accession, subject_length, ring + brct))
    except (OSError, UnicodeDecodeError) as exc:
        fail(f"cannot read candidate QC table {path}: {exc}")
    if not candidates:
        fail(f"candidate QC table contains no records: {path}")
    return candidates


def alignment_columns(records: dict[str, str], tip_order: Sequence[str]) -> list[tuple[str, ...]]:
    length = len(records[tip_order[0]])
    return [
        tuple(records[tip][column_index] for tip in tip_order)
        for column_index in range(length)
    ]


def unique_column_subset_map(
    raw_columns: Sequence[tuple[str, ...]],
    trimmed_columns: Sequence[tuple[str, ...]],
    *,
    profile: str,
) -> tuple[int, ...]:
    """Return zero-based raw indices for a unique in-order subsequence map."""

    earliest: list[int] = []
    trimmed_index = 0
    for raw_index, column in enumerate(raw_columns):
        if trimmed_index < len(trimmed_columns) and column == trimmed_columns[trimmed_index]:
            earliest.append(raw_index)
            trimmed_index += 1
    if trimmed_index != len(trimmed_columns):
        fail(
            f"profile {profile!r} is not an in-order column subset of the raw alignment; "
            f"first unmatched trimmed column is {trimmed_index + 1}"
        )

    latest_reversed: list[int] = []
    trimmed_index = len(trimmed_columns) - 1
    for raw_index in range(len(raw_columns) - 1, -1, -1):
        if trimmed_index >= 0 and raw_columns[raw_index] == trimmed_columns[trimmed_index]:
            latest_reversed.append(raw_index)
            trimmed_index -= 1
    latest = list(reversed(latest_reversed))
    if len(latest) != len(trimmed_columns):
        fail(f"internal error while deriving the latest map for profile {profile!r}")

    for trimmed_index, (first, last) in enumerate(zip(earliest, latest), start=1):
        if first != last:
            fail(
                f"profile {profile!r} has an ambiguous raw-column map at trimmed "
                f"column {trimmed_index}: raw columns {first + 1} and {last + 1} "
                "are both valid; refusing to guess"
            )
    return tuple(earliest)


def residue_to_raw_column(sequence: str, *, accession: str) -> dict[int, int]:
    mapping: dict[int, int] = {}
    residue_number = 0
    for column_index, character in enumerate(sequence, start=1):
        if character != "-":
            residue_number += 1
            mapping[residue_number] = column_index
    if not mapping:
        fail(f"{accession}: raw alignment sequence contains no residues")
    return mapping


def build_report(
    raw_records: dict[str, str],
    profiles: Sequence[tuple[str, dict[str, str]]],
    candidates: Sequence[CandidateDomains],
    minimum_fraction: float,
) -> list[dict[str, str]]:
    raw_tips = set(raw_records)
    qc_tips = {candidate.accession for candidate in candidates}
    if qc_tips != raw_tips:
        missing = sorted(raw_tips - qc_tips)
        extra = sorted(qc_tips - raw_tips)
        fail(
            "candidate QC/raw-alignment tip-set mismatch; "
            f"missing_from_qc={missing}; extra_in_qc={extra}"
        )

    tip_order = sorted(raw_records)
    raw_columns = alignment_columns(raw_records, tip_order)
    retained_raw_columns: dict[str, set[int]] = {}
    for label, records in profiles:
        profile_tips = set(records)
        if profile_tips != raw_tips:
            missing = sorted(raw_tips - profile_tips)
            extra = sorted(profile_tips - raw_tips)
            fail(
                f"profile {label!r}/raw-alignment tip-set mismatch; "
                f"missing={missing}; extra={extra}"
            )
        trimmed_columns = alignment_columns(records, tip_order)
        mapping = unique_column_subset_map(raw_columns, trimmed_columns, profile=label)
        retained_raw_columns[label] = {index + 1 for index in mapping}

    rows: list[dict[str, str]] = []
    for candidate in candidates:
        residue_map = residue_to_raw_column(
            raw_records[candidate.accession], accession=candidate.accession
        )
        if len(residue_map) != candidate.subject_length:
            fail(
                f"{candidate.accession}: ungapped raw-alignment length "
                f"{len(residue_map)} does not equal subject_length "
                f"{candidate.subject_length}"
            )
        for label, _records in profiles:
            retained_columns = retained_raw_columns[label]
            for interval in candidate.intervals:
                interval_columns = {
                    residue_map[residue]
                    for residue in range(interval.start, interval.end + 1)
                }
                total = interval.end - interval.start + 1
                retained = len(interval_columns & retained_columns)
                fraction = retained / total
                rows.append(
                    {
                        "accession": candidate.accession,
                        "profile": label,
                        "domain": interval.domain,
                        "interval_index": str(interval.index),
                        "residue_start": str(interval.start),
                        "residue_end": str(interval.end),
                        "raw_msa_start": str(residue_map[interval.start]),
                        "raw_msa_end": str(residue_map[interval.end]),
                        "total_residues": str(total),
                        "retained_residues": str(retained),
                        "retained_fraction": f"{fraction:.6f}",
                        "minimum_retained_fraction": f"{minimum_fraction:.6f}",
                        "status": "pass" if fraction >= minimum_fraction else "fail",
                    }
                )
    return rows


def write_report(path: Path, rows: Sequence[dict[str, str]]) -> None:
    if path.exists():
        fail(f"refusing to overwrite output report: {path}")
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
        )
        try:
            with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=OUTPUT_COLUMNS,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(rows)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary_name, path)
        except BaseException:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass
            raise
    except OSError as exc:
        fail(f"cannot write output report {path}: {exc}")


def run(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Fail-closed audit of annotated domain retention after alignment trimming."
    )
    parser.add_argument("--raw-alignment", required=True, type=Path)
    parser.add_argument(
        "--trimmed-alignment",
        action="append",
        required=True,
        type=parse_profile,
        metavar="LABEL=PATH",
    )
    parser.add_argument("--candidate-qc", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--minimum-retained-fraction",
        default=0.8,
        type=parse_fraction,
        metavar="FRACTION",
    )
    args = parser.parse_args(argv)

    labels = [label for label, _path in args.trimmed_alignment]
    duplicate_labels = sorted({label for label in labels if labels.count(label) > 1})
    if duplicate_labels:
        fail("duplicate trimmed-alignment label(s): " + ", ".join(duplicate_labels))

    raw_records = read_alignment(args.raw_alignment)
    profiles = [
        (label, read_alignment(path)) for label, path in args.trimmed_alignment
    ]
    candidates = read_candidate_domains(args.candidate_qc)
    rows = build_report(
        raw_records, profiles, candidates, args.minimum_retained_fraction
    )
    write_report(args.output, rows)

    failures = [row for row in rows if row["status"] == "fail"]
    if failures:
        examples = ", ".join(
            f"{row['accession']}:{row['profile']}:{row['domain']}"
            f"[{row['interval_index']}]={row['retained_fraction']}"
            for row in failures[:5]
        )
        if len(failures) > 5:
            examples += f", ... ({len(failures)} intervals total)"
        print(
            "ERROR: domain-retention threshold failed; " + examples,
            file=sys.stderr,
        )
        return 2

    print(
        f"Domain-retention audit passed: {len(rows)} intervals across "
        f"{len(profiles)} trim profile(s); report={args.output}"
    )
    return 0


def main() -> None:
    try:
        raise SystemExit(run())
    except AuditError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc


if __name__ == "__main__":
    main()
