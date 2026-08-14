#!/usr/bin/env python3
"""Fail-closed BRCA1 similarity and terminal-domain QC summarizer.

The companion ``run_brca1_qc.pbs`` script writes BLAST outfmt 6 with these
fourteen columns (and no header)::

    qseqid sseqid nident length mismatch gapopen qstart qend sstart send
    evalue bitscore qlen slen

InterProScan input is its standard TSV output.  Only Pfam and SMART signature
rows are accepted.  A candidate passes domain QC only when it has an
N-terminal RING-related signature and at least two *distinct* C-terminal
BRCT-related signature intervals.  Keyword recognition is supplemented by
well-known signature/integrated accessions, rather than depending on one
database identifier.

The script never edits its inputs.  It refuses to overwrite ``candidate_qc.tsv``
or ``candidates.tsv`` in the requested output directory.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


BLAST_COLUMNS = (
    "qseqid",
    "sseqid",
    "nident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
)

INTERPRO_COLUMNS = (
    "protein_accession",
    "sequence_md5",
    "sequence_length",
    "analysis",
    "signature_accession",
    "signature_description",
    "start",
    "end",
    "score",
    "status",
    "date",
    "interpro_accession",
    "interpro_description",
    "go_terms",
    "pathways",
)

REQUIRED_CANDIDATE_COLUMNS = {
    "accession",
    "relation",
    "is_canonical",
    "is_fragment",
    "query_coverage",
    "sequence_length",
    "bitscore",
    "evalue",
    "target_coverage",
    "percent_identity",
    "alignment_length",
    "domain_architecture",
    "notes",
}

QC_COLUMNS = (
    "accession",
    "query_accession",
    "hsp_count",
    "query_length",
    "subject_length",
    "query_covered_residues",
    "subject_covered_residues",
    "query_coverage",
    "target_coverage",
    "percent_identity",
    "alignment_length",
    "min_evalue",
    "sum_bitscore",
    "interpro_signature_hits",
    "n_terminal_ring_signature_hits",
    "ring_coordinates",
    "ring_evidence",
    "c_terminal_brct_signature_hits",
    "brct_repeat_count",
    "brct_coordinates",
    "brct_evidence",
    "domain_qc_status",
)

RING_ACCESSIONS = {
    "PF00097",  # Pfam zf-C3HC4/RING
    "SM00184",  # SMART RING
    "IPR001841",  # InterPro RING-type zinc finger
}
BRCT_ACCESSIONS = {
    "PF00533",  # Pfam BRCT
    "SM00292",  # SMART BRCT
    "IPR001357",  # InterPro BRCT domain
}

RING_TEXT = re.compile(
    r"(?:\bring(?:[- ]?(?:type|finger|domain))?\b|"
    r"\bzf[-_ ]?c3hc4\b|\bc3hc4(?:[- ]type)?\b|"
    r"zinc[- ]finger[^;]*\bring\b)",
    re.IGNORECASE,
)
BRCT_TEXT = re.compile(
    r"(?:\bbrct\b|\bbrca1[- ]?c[- ]?termin(?:al|us)\b)",
    re.IGNORECASE,
)


class QCError(RuntimeError):
    """Raised when an input or candidate fails closed validation."""


@dataclass(frozen=True)
class HSP:
    query: str
    subject: str
    nident: int
    alignment_length: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    evalue: float
    bitscore: float
    query_length: int
    subject_length: int


@dataclass(frozen=True)
class DomainHit:
    accession: str
    sequence_length: int
    analysis: str
    signature_accession: str
    signature_description: str
    start: int
    end: int
    interpro_accession: str
    interpro_description: str

    @property
    def midpoint(self) -> float:
        return (self.start + self.end) / 2

    @property
    def evidence(self) -> str:
        label = self.signature_description or self.interpro_description or "unnamed signature"
        return (
            f"{self.analysis}:{self.signature_accession}:{label}"
            f"@{self.start}-{self.end}"
        )

    def identifiers_and_text(self) -> tuple[set[str], str]:
        identifiers = {self.signature_accession.upper()}
        if self.interpro_accession:
            identifiers.add(self.interpro_accession.upper())
        text = " ".join(
            part
            for part in (
                self.signature_accession,
                self.signature_description,
                self.interpro_accession,
                self.interpro_description,
            )
            if part
        )
        return identifiers, text


def fail(message: str) -> None:
    raise QCError(message)


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier = ""
    chunks: list[str] = []

    def store() -> None:
        if not identifier:
            return
        sequence = "".join(chunks).upper()
        if not sequence or not re.fullmatch(r"[A-Z]+", sequence):
            fail(f"invalid amino-acid sequence for {identifier} in {path}")
        if identifier in records:
            fail(f"duplicate FASTA identifier in {path}: {identifier}")
        records[identifier] = sequence

    try:
        with path.open(encoding="utf-8") as handle:
            for raw in handle:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    store()
                    header = line[1:].strip()
                    identifier = header.split()[0] if header else ""
                    if not identifier:
                        fail(f"empty FASTA identifier in {path}")
                    chunks = []
                else:
                    if not identifier:
                        fail(f"sequence data precedes the first FASTA header in {path}")
                    chunks.append(line)
        store()
    except (OSError, UnicodeDecodeError) as exc:
        fail(f"cannot read FASTA {path}: {exc}")
    if not records:
        fail(f"FASTA contains no records: {path}")
    return records


def parse_candidates(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                fail(f"candidate table has no header: {path}")
            missing = sorted(REQUIRED_CANDIDATE_COLUMNS - set(reader.fieldnames))
            if missing:
                fail("candidate table is missing columns: " + ", ".join(missing))
            rows: list[dict[str, str]] = []
            seen: set[str] = set()
            for line_number, raw in enumerate(reader, start=2):
                if None in raw:
                    fail(f"candidate table row {line_number} has extra columns")
                row = {field: raw.get(field, "") for field in reader.fieldnames}
                accession = row["accession"]
                if not accession or accession != accession.strip():
                    fail(f"candidate table row {line_number} has an invalid accession")
                if accession in seen:
                    fail(f"duplicate candidate accession: {accession}")
                seen.add(accession)
                try:
                    sequence_length = int(row["sequence_length"])
                except ValueError as exc:
                    raise QCError(
                        f"{accession}: sequence_length must be an integer"
                    ) from exc
                if sequence_length <= 0:
                    fail(f"{accession}: sequence_length must be positive")
                rows.append(row)
    except (OSError, UnicodeDecodeError) as exc:
        fail(f"cannot read candidate table {path}: {exc}")
    if not rows:
        fail(f"candidate table contains no records: {path}")
    return list(reader.fieldnames), rows


def parse_nonnegative_int(value: str, context: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise QCError(f"{context} must be an integer: {value!r}") from exc
    if parsed < 0:
        fail(f"{context} must be non-negative: {value!r}")
    return parsed


def parse_nonnegative_float(value: str, context: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise QCError(f"{context} must be numeric: {value!r}") from exc
    if not math.isfinite(parsed) or parsed < 0:
        fail(f"{context} must be finite and non-negative: {value!r}")
    return parsed


def parse_blast(path: Path, accessions: set[str]) -> list[HSP]:
    hsps: list[HSP] = []
    try:
        with path.open(encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, start=1):
                line = raw.rstrip("\r\n")
                if not line or line.startswith("#"):
                    continue
                values = line.split("\t")
                if len(values) != len(BLAST_COLUMNS):
                    fail(
                        f"BLAST row {line_number} has {len(values)} columns; "
                        f"expected {len(BLAST_COLUMNS)}"
                    )
                row = dict(zip(BLAST_COLUMNS, values))
                if row["sseqid"] not in accessions:
                    fail(f"BLAST row {line_number} has an unknown subject: {row['sseqid']}")
                nident = parse_nonnegative_int(row["nident"], f"BLAST row {line_number} nident")
                alignment_length = parse_nonnegative_int(
                    row["length"], f"BLAST row {line_number} length"
                )
                query_start = parse_nonnegative_int(
                    row["qstart"], f"BLAST row {line_number} qstart"
                )
                query_end = parse_nonnegative_int(row["qend"], f"BLAST row {line_number} qend")
                subject_start = parse_nonnegative_int(
                    row["sstart"], f"BLAST row {line_number} sstart"
                )
                subject_end = parse_nonnegative_int(row["send"], f"BLAST row {line_number} send")
                query_length = parse_nonnegative_int(row["qlen"], f"BLAST row {line_number} qlen")
                subject_length = parse_nonnegative_int(row["slen"], f"BLAST row {line_number} slen")
                if alignment_length == 0 or query_length == 0 or subject_length == 0:
                    fail(f"BLAST row {line_number} contains a zero length")
                if nident > alignment_length:
                    fail(f"BLAST row {line_number} has nident greater than alignment length")
                if not (1 <= query_start <= query_end <= query_length):
                    fail(f"BLAST row {line_number} has invalid protein query coordinates")
                if not (1 <= subject_start <= subject_end <= subject_length):
                    fail(f"BLAST row {line_number} has invalid protein subject coordinates")
                hsps.append(
                    HSP(
                        query=row["qseqid"],
                        subject=row["sseqid"],
                        nident=nident,
                        alignment_length=alignment_length,
                        query_start=query_start,
                        query_end=query_end,
                        subject_start=subject_start,
                        subject_end=subject_end,
                        evalue=parse_nonnegative_float(
                            row["evalue"], f"BLAST row {line_number} evalue"
                        ),
                        bitscore=parse_nonnegative_float(
                            row["bitscore"], f"BLAST row {line_number} bitscore"
                        ),
                        query_length=query_length,
                        subject_length=subject_length,
                    )
                )
    except (OSError, UnicodeDecodeError) as exc:
        fail(f"cannot read BLAST table {path}: {exc}")
    if not hsps:
        fail(f"BLAST table contains no HSPs: {path}")
    return hsps


def parse_interpro(path: Path, sequences: dict[str, str]) -> list[DomainHit]:
    hits: list[DomainHit] = []
    try:
        with path.open(encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, start=1):
                line = raw.rstrip("\r\n")
                if not line or line.startswith("#"):
                    continue
                values = line.split("\t")
                if len(values) < 11 or len(values) > len(INTERPRO_COLUMNS):
                    fail(
                        f"InterProScan row {line_number} has {len(values)} columns; "
                        "expected 11 to 15 standard TSV columns"
                    )
                values.extend([""] * (len(INTERPRO_COLUMNS) - len(values)))
                row = dict(zip(INTERPRO_COLUMNS, values))
                accession = row["protein_accession"]
                if accession not in sequences:
                    fail(f"InterProScan row {line_number} has an unknown accession: {accession}")
                if row["analysis"] not in {"Pfam", "SMART"}:
                    fail(
                        f"InterProScan row {line_number} uses unexpected analysis "
                        f"{row['analysis']!r}; expected Pfam or SMART"
                    )
                if row["status"] != "T":
                    fail(f"InterProScan row {line_number} has non-true status {row['status']!r}")
                sequence_length = parse_nonnegative_int(
                    row["sequence_length"], f"InterProScan row {line_number} sequence length"
                )
                if sequence_length != len(sequences[accession]):
                    fail(
                        f"{accession}: InterProScan length {sequence_length} does not match "
                        f"FASTA length {len(sequences[accession])}"
                    )
                start = parse_nonnegative_int(row["start"], f"InterProScan row {line_number} start")
                end = parse_nonnegative_int(row["end"], f"InterProScan row {line_number} end")
                if not (1 <= start <= end <= sequence_length):
                    fail(f"InterProScan row {line_number} has invalid coordinates")
                hits.append(
                    DomainHit(
                        accession=accession,
                        sequence_length=sequence_length,
                        analysis=row["analysis"],
                        signature_accession=row["signature_accession"],
                        signature_description="" if row["signature_description"] == "-" else row["signature_description"],
                        start=start,
                        end=end,
                        interpro_accession="" if row["interpro_accession"] == "-" else row["interpro_accession"],
                        interpro_description="" if row["interpro_description"] == "-" else row["interpro_description"],
                    )
                )
    except (OSError, UnicodeDecodeError) as exc:
        fail(f"cannot read InterProScan table {path}: {exc}")
    if not hits:
        fail(f"InterProScan table contains no signature hits: {path}")
    return hits


def interval_union_length(intervals: Iterable[tuple[int, int]]) -> int:
    ordered = sorted(intervals)
    if not ordered:
        return 0
    total = 0
    current_start, current_end = ordered[0]
    for start, end in ordered[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start + 1
            current_start, current_end = start, end
    return total + current_end - current_start + 1


def related_to(hit: DomainHit, known: set[str], pattern: re.Pattern[str]) -> bool:
    identifiers, text = hit.identifiers_and_text()
    return bool(identifiers & known or pattern.search(text))


def overlap_fraction(left: DomainHit, right: DomainHit) -> float:
    overlap = max(0, min(left.end, right.end) - max(left.start, right.start) + 1)
    shortest = min(left.end - left.start + 1, right.end - right.start + 1)
    return overlap / shortest


def cluster_overlapping_hits(hits: Sequence[DomainHit]) -> list[list[DomainHit]]:
    """Group database signatures for the same repeat without merging tandem repeats."""

    clusters: list[list[DomainHit]] = []
    for hit in sorted(hits, key=lambda item: (item.start, item.end, item.analysis)):
        matching = [
            cluster
            for cluster in clusters
            if any(overlap_fraction(hit, existing) >= 0.35 for existing in cluster)
        ]
        if not matching:
            clusters.append([hit])
            continue
        primary = matching[0]
        primary.append(hit)
        for extra in matching[1:]:
            primary.extend(extra)
            clusters.remove(extra)
    return sorted(clusters, key=lambda cluster: min(item.start for item in cluster))


def cluster_coordinates(clusters: Sequence[Sequence[DomainHit]]) -> str:
    return ";".join(
        f"{min(hit.start for hit in cluster)}-{max(hit.end for hit in cluster)}"
        for cluster in clusters
    )


def evidence_text(hits: Sequence[DomainHit]) -> str:
    return ";".join(sorted({hit.evidence for hit in hits}))


def stable_decimal(value: float, places: int = 6) -> str:
    rendered = f"{value:.{places}f}".rstrip("0").rstrip(".")
    return rendered or "0"


def stable_evalue(value: float) -> str:
    return "0" if value == 0 else f"{value:.6g}"


def summarize_subject(
    accession: str,
    query_accession: str,
    sequence_length: int,
    hsps: Sequence[HSP],
    hits: Sequence[DomainHit],
    ring_max_midpoint_fraction: float,
    brct_min_midpoint_fraction: float,
) -> dict[str, str]:
    if not hsps:
        fail(f"{accession}: BLAST reported no HSP; similarity QC is incomplete")
    query_ids = {hsp.query for hsp in hsps}
    if query_ids != {query_accession}:
        fail(f"{accession}: BLAST query IDs do not exactly match {query_accession}")
    query_lengths = {hsp.query_length for hsp in hsps}
    subject_lengths = {hsp.subject_length for hsp in hsps}
    if len(query_lengths) != 1 or len(subject_lengths) != 1:
        fail(f"{accession}: BLAST reports inconsistent sequence lengths across HSPs")
    query_length = next(iter(query_lengths))
    subject_length = next(iter(subject_lengths))
    if subject_length != sequence_length:
        fail(
            f"{accession}: BLAST subject length {subject_length} does not match "
            f"FASTA length {sequence_length}"
        )
    aligned_total = sum(hsp.alignment_length for hsp in hsps)
    if aligned_total <= 0:
        fail(f"{accession}: BLAST summed alignment length is zero")
    query_covered = interval_union_length((hsp.query_start, hsp.query_end) for hsp in hsps)
    subject_covered = interval_union_length(
        (hsp.subject_start, hsp.subject_end) for hsp in hsps
    )

    ring_hits = [
        hit
        for hit in hits
        if related_to(hit, RING_ACCESSIONS, RING_TEXT)
        and hit.midpoint <= sequence_length * ring_max_midpoint_fraction
    ]
    brct_hits = [
        hit
        for hit in hits
        if related_to(hit, BRCT_ACCESSIONS, BRCT_TEXT)
        and hit.midpoint >= sequence_length * brct_min_midpoint_fraction
    ]
    ring_clusters = cluster_overlapping_hits(ring_hits)
    brct_clusters = cluster_overlapping_hits(brct_hits)
    problems: list[str] = []
    if not ring_clusters:
        problems.append("no N-terminal RING-related signature")
    if len(brct_clusters) < 2:
        problems.append(
            f"only {len(brct_clusters)} distinct C-terminal BRCT-related signature interval(s)"
        )
    if problems:
        fail(f"{accession}: terminal-domain QC failed: " + "; ".join(problems))

    return {
        "accession": accession,
        "query_accession": query_accession,
        "hsp_count": str(len(hsps)),
        "query_length": str(query_length),
        "subject_length": str(subject_length),
        "query_covered_residues": str(query_covered),
        "subject_covered_residues": str(subject_covered),
        "query_coverage": stable_decimal(query_covered / query_length),
        "target_coverage": stable_decimal(subject_covered / subject_length),
        "percent_identity": stable_decimal(
            100 * sum(hsp.nident for hsp in hsps) / aligned_total, places=3
        ),
        "alignment_length": str(aligned_total),
        "min_evalue": stable_evalue(min(hsp.evalue for hsp in hsps)),
        "sum_bitscore": stable_decimal(sum(hsp.bitscore for hsp in hsps), places=3),
        "interpro_signature_hits": str(len(hits)),
        "n_terminal_ring_signature_hits": str(len(ring_hits)),
        "ring_coordinates": cluster_coordinates(ring_clusters),
        "ring_evidence": evidence_text(ring_hits),
        "c_terminal_brct_signature_hits": str(len(brct_hits)),
        "brct_repeat_count": str(len(brct_clusters)),
        "brct_coordinates": cluster_coordinates(brct_clusters),
        "brct_evidence": evidence_text(brct_hits),
        "domain_qc_status": "pass",
    }


def update_candidate(row: dict[str, str], summary: dict[str, str], query: str) -> dict[str, str]:
    updated = dict(row)
    updated.update(
        {
            "query_coverage": summary["query_coverage"],
            "target_coverage": summary["target_coverage"],
            "percent_identity": summary["percent_identity"],
            "alignment_length": summary["alignment_length"],
            "bitscore": summary["sum_bitscore"],
            "evalue": summary["min_evalue"],
            "is_fragment": "false",
            "domain_architecture": (
                "InterProScan Pfam+SMART pass: N-terminal RING "
                f"({summary['ring_coordinates']}); tandem C-terminal BRCT "
                f"({summary['brct_coordinates']})"
            ),
        }
    )
    # RefSeq NP_ accessions are curated records, but that prefix alone is not
    # evidence that a non-query sequence is the organism's canonical isoform.
    if row["accession"] != query:
        updated["is_canonical"] = "false"
    note = (
        "QC: BLAST union coverage and Pfam+SMART terminal architecture passed "
        f"(RING {summary['ring_coordinates']}; BRCT {summary['brct_coordinates']})."
    )
    updated["notes"] = f"{row['notes'].rstrip()} {note}".strip()
    return updated


def bounded_fraction(value: str, name: str) -> float:
    parsed = parse_nonnegative_float(value, name)
    if not 0 < parsed < 1:
        fail(f"{name} must be greater than 0 and less than 1")
    return parsed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-fasta", type=Path, required=True)
    parser.add_argument("--candidate-table", type=Path, required=True)
    parser.add_argument("--blast-tsv", type=Path, required=True)
    parser.add_argument("--interpro-tsv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--ring-max-midpoint-fraction",
        default="0.25",
        help="largest sequence-relative midpoint accepted as N-terminal (default: 0.25)",
    )
    parser.add_argument(
        "--brct-min-midpoint-fraction",
        default="0.65",
        help="smallest sequence-relative midpoint accepted as C-terminal (default: 0.65)",
    )
    args = parser.parse_args()

    try:
        ring_fraction = bounded_fraction(
            args.ring_max_midpoint_fraction, "--ring-max-midpoint-fraction"
        )
        brct_fraction = bounded_fraction(
            args.brct_min_midpoint_fraction, "--brct-min-midpoint-fraction"
        )
        sequences = parse_fasta(args.candidate_fasta)
        fieldnames, candidates = parse_candidates(args.candidate_table)
        rows_by_accession = {row["accession"]: row for row in candidates}
        if set(rows_by_accession) != set(sequences):
            only_table = sorted(set(rows_by_accession) - set(sequences))
            only_fasta = sorted(set(sequences) - set(rows_by_accession))
            fail(
                "candidate table/FASTA accession mismatch; "
                f"table-only={only_table}, FASTA-only={only_fasta}"
            )
        for accession, row in rows_by_accession.items():
            if int(row["sequence_length"]) != len(sequences[accession]):
                fail(
                    f"{accession}: candidate-table length {row['sequence_length']} does not "
                    f"match FASTA length {len(sequences[accession])}"
                )
        self_rows = [row for row in candidates if row["relation"] == "self"]
        if len(self_rows) != 1:
            fail(f"expected exactly one relation=self query row; observed {len(self_rows)}")
        query_accession = self_rows[0]["accession"]

        hsps = parse_blast(args.blast_tsv, set(sequences))
        query_lengths = {hsp.query_length for hsp in hsps}
        if query_lengths != {len(sequences[query_accession])}:
            fail(
                "BLAST query length does not match the relation=self candidate sequence: "
                f"{query_lengths} versus {len(sequences[query_accession])}"
            )
        hits = parse_interpro(args.interpro_tsv, sequences)
        hsps_by_subject: dict[str, list[HSP]] = {accession: [] for accession in sequences}
        hits_by_accession: dict[str, list[DomainHit]] = {accession: [] for accession in sequences}
        for hsp in hsps:
            hsps_by_subject[hsp.subject].append(hsp)
        for hit in hits:
            hits_by_accession[hit.accession].append(hit)

        summaries: list[dict[str, str]] = []
        updated_candidates: list[dict[str, str]] = []
        for row in candidates:
            accession = row["accession"]
            summary = summarize_subject(
                accession,
                query_accession,
                len(sequences[accession]),
                hsps_by_subject[accession],
                hits_by_accession[accession],
                ring_fraction,
                brct_fraction,
            )
            summaries.append(summary)
            updated_candidates.append(update_candidate(row, summary, query_accession))

        args.out_dir.mkdir(parents=True, exist_ok=True)
        qc_path = args.out_dir / "candidate_qc.tsv"
        candidate_path = args.out_dir / "candidates.tsv"
        for output in (qc_path, candidate_path):
            if output.exists():
                fail(f"refusing to overwrite existing output: {output}")
        with qc_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=QC_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(summaries)
        with candidate_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(updated_candidates)
    except QCError as exc:
        parser.exit(2, f"BRCA1 QC failed: {exc}\n")

    print(f"Validated {len(candidates)} BRCA1 candidates.")
    print(f"QC audit: {qc_path}")
    print(f"Updated planner table: {candidate_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
