#!/usr/bin/env python3
"""Resolve scientific names against local NCBI Taxonomy dump files.

The resolver is deliberately local and strict. It performs case- and
punctuation-sensitive equality against ``names.dmp`` records whose name class
is exactly ``scientific name``. It never downloads a taxonomy snapshot, never
uses fuzzy matching, and never chooses the first result when a name is
ambiguous.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import tempfile
from datetime import datetime
from urllib.parse import urlsplit
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


VERSION = "0.1.0"
OFFICIAL_TAXONOMY_ROOT = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
OFFICIAL_CURRENT_ARCHIVE = (
    f"{OFFICIAL_TAXONOMY_ROOT}new_taxdump/new_taxdump.tar.gz"
)
ASCII_TAXON_ID = re.compile(r"[0-9]+")
NCBI_ARCHIVE_PATH = re.compile(
    r"/pub/taxonomy/(?:new_taxdump/new_taxdump\.(?:tar\.gz|zip)|"
    r"taxdump_archive/new_taxdump_[0-9]{4}-[0-9]{2}-[0-9]{2}\.zip)"
)
RESOLUTION_COLUMNS = (
    "record_id",
    "input_name",
    "requested_taxon_id",
    "matched_name",
    "name_class",
    "taxon_id",
    "parent_taxon_id",
    "rank",
    "unique_name",
    "status",
    "snapshot",
    "source_url",
    "retrieved_at",
    "names_sha256",
    "nodes_sha256",
    "resolver_version",
)


class TaxonomyError(Exception):
    """A stable, user-facing taxonomy validation failure."""

    def __init__(self, code: str, message: str, *, details: Sequence[Mapping[str, str]] = ()):
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = [dict(item) for item in details]


@dataclass(frozen=True)
class NameRequest:
    record_id: str
    input_name: str
    expected_taxon_id: str = ""


@dataclass(frozen=True)
class NameMatch:
    taxon_id: str
    matched_name: str
    unique_name: str
    name_class: str


@dataclass(frozen=True)
class Node:
    taxon_id: str
    parent_taxon_id: str
    rank: str


@dataclass(frozen=True)
class TaxonomyResolution:
    rows: Tuple[Mapping[str, str], ...]
    names_sha256: str
    nodes_sha256: str
    snapshot: str
    source_url: str
    retrieved_at: str
    resolver_version: str = VERSION


def _decode_line(raw_line: bytes, path: Path, line_number: int) -> str:
    try:
        return raw_line.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise TaxonomyError(
            "MALFORMED_DMP",
            f"{path.name}:{line_number}: record is not valid UTF-8.",
        ) from exc


def _parse_dmp_fields(raw_line: bytes, path: Path, line_number: int, minimum: int) -> List[str]:
    text = _decode_line(raw_line, path, line_number)
    if text.endswith("\t|\r\n"):
        payload = text[:-4]
    elif text.endswith("\t|\n"):
        payload = text[:-3]
    elif text.endswith("\t|"):
        payload = text[:-2]
    else:
        raise TaxonomyError(
            "MALFORMED_DMP",
            f"{path.name}:{line_number}: expected the official NCBI record terminator TAB|.",
        )
    fields = payload.split("\t|\t")
    if len(fields) < minimum:
        raise TaxonomyError(
            "MALFORMED_DMP",
            f"{path.name}:{line_number}: expected at least {minimum} fields, found {len(fields)}.",
        )
    return fields


def _validate_requests(requests: Sequence[NameRequest]) -> None:
    if not requests:
        raise TaxonomyError("EMPTY_REQUEST", "At least one organism name is required.")
    seen_ids: set[str] = set()
    for request in requests:
        if not request.record_id:
            raise TaxonomyError("INVALID_REQUEST", "Every record needs a non-empty record_id.")
        if request.record_id in seen_ids:
            raise TaxonomyError("DUPLICATE_RECORD_ID", f"Duplicate record_id: {request.record_id}")
        seen_ids.add(request.record_id)
        if not request.input_name:
            raise TaxonomyError("INVALID_REQUEST", f"{request.record_id}: organism name is empty.")
        if request.input_name != request.input_name.strip():
            raise TaxonomyError(
                "NONEXACT_INPUT_NAME",
                f"{request.record_id}: leading or trailing whitespace prevents exact matching; "
                f"retry explicitly with {request.input_name.strip()!r} if that is intended.",
            )
        if request.expected_taxon_id and not ASCII_TAXON_ID.fullmatch(
            request.expected_taxon_id
        ):
            raise TaxonomyError(
                "INVALID_TAXON_ID",
                f"{request.record_id}: expected TaxID must contain digits only.",
            )


def _scan_names(
    path: Path, wanted_names: Iterable[str]
) -> Tuple[Dict[str, List[NameMatch]], Dict[str, List[NameMatch]], str]:
    wanted = set(wanted_names)
    exact: Dict[str, List[NameMatch]] = {name: [] for name in wanted}
    secondary: Dict[str, List[NameMatch]] = {name: [] for name in wanted}
    digest = hashlib.sha256()
    try:
        handle = path.open("rb")
    except OSError as exc:
        raise TaxonomyError("DMP_IO_ERROR", f"Cannot open names.dmp: {exc}") from exc
    with handle:
        for line_number, raw_line in enumerate(handle, start=1):
            digest.update(raw_line)
            fields = _parse_dmp_fields(raw_line, path, line_number, 4)
            taxon_id, name_txt, unique_name, name_class = fields[:4]
            if not ASCII_TAXON_ID.fullmatch(taxon_id):
                raise TaxonomyError(
                    "MALFORMED_DMP", f"{path.name}:{line_number}: tax_id is not numeric."
                )
            if name_txt not in wanted:
                continue
            match = NameMatch(taxon_id, name_txt, unique_name, name_class)
            if name_class == "scientific name":
                exact[name_txt].append(match)
            else:
                secondary[name_txt].append(match)
    return exact, secondary, digest.hexdigest()


def _scan_nodes(path: Path, wanted_taxon_ids: Iterable[str]) -> Tuple[Dict[str, Node], str]:
    wanted = set(wanted_taxon_ids)
    nodes: Dict[str, Node] = {}
    digest = hashlib.sha256()
    try:
        handle = path.open("rb")
    except OSError as exc:
        raise TaxonomyError("DMP_IO_ERROR", f"Cannot open nodes.dmp: {exc}") from exc
    with handle:
        for line_number, raw_line in enumerate(handle, start=1):
            digest.update(raw_line)
            fields = _parse_dmp_fields(raw_line, path, line_number, 3)
            taxon_id, parent_taxon_id, rank = fields[:3]
            if not ASCII_TAXON_ID.fullmatch(taxon_id) or not ASCII_TAXON_ID.fullmatch(
                parent_taxon_id
            ):
                raise TaxonomyError(
                    "MALFORMED_DMP",
                    f"{path.name}:{line_number}: tax_id and parent tax_id must be numeric.",
                )
            if taxon_id not in wanted:
                continue
            node = Node(taxon_id, parent_taxon_id, rank)
            if taxon_id in nodes and nodes[taxon_id] != node:
                raise TaxonomyError(
                    "DUPLICATE_NODE", f"{path.name}: conflicting records for TaxID {taxon_id}."
                )
            nodes[taxon_id] = node
    return nodes, digest.hexdigest()


def resolve_exact_scientific_names(
    requests: Sequence[NameRequest],
    *,
    names_path: Path,
    nodes_path: Path,
    snapshot: str,
    source_url: str,
    retrieved_at: str,
) -> TaxonomyResolution:
    """Resolve all requests or fail without guessing or partial success."""

    _validate_requests(requests)
    if not snapshot or not source_url or not retrieved_at:
        raise TaxonomyError(
            "MISSING_PROVENANCE", "snapshot, source_url, and retrieved_at are required."
        )
    for field_name, value in (
        ("snapshot", snapshot),
        ("source_url", source_url),
        ("retrieved_at", retrieved_at),
    ):
        if value != value.strip():
            raise TaxonomyError(
                "INVALID_PROVENANCE",
                f"{field_name} has leading or trailing whitespace; provenance is never normalized silently.",
            )
    parsed_source = urlsplit(source_url)
    if (
        parsed_source.scheme != "https"
        or parsed_source.netloc != "ftp.ncbi.nlm.nih.gov"
        or not parsed_source.path.startswith("/pub/taxonomy/")
        or not NCBI_ARCHIVE_PATH.fullmatch(parsed_source.path)
        or parsed_source.query
        or parsed_source.fragment
    ):
        raise TaxonomyError(
            "UNOFFICIAL_SOURCE_URL",
            "source_url must be an exact HTTPS archive URL under "
            f"{OFFICIAL_TAXONOMY_ROOT}",
        )
    try:
        parsed_retrieved_at = datetime.fromisoformat(retrieved_at.replace("Z", "+00:00"))
    except ValueError as exc:
        raise TaxonomyError(
            "INVALID_RETRIEVED_AT",
            "retrieved_at must be an ISO 8601 UTC date or timestamp.",
        ) from exc
    date_only = bool(re.fullmatch(r"[0-9]{4}-[0-9]{2}-[0-9]{2}", retrieved_at))
    if not date_only and parsed_retrieved_at.tzinfo is None:
        raise TaxonomyError(
            "INVALID_RETRIEVED_AT",
            "retrieved_at timestamps must include Z or +00:00.",
        )
    if parsed_retrieved_at.tzinfo is not None and parsed_retrieved_at.utcoffset().total_seconds() != 0:
        raise TaxonomyError(
            "INVALID_RETRIEVED_AT",
            "retrieved_at must use UTC when a time-zone offset is present.",
        )
    exact, secondary, names_sha256 = _scan_names(
        names_path, (request.input_name for request in requests)
    )

    failures: List[Mapping[str, str]] = []
    chosen: Dict[str, NameMatch] = {}
    for name in sorted(exact):
        by_taxon = {match.taxon_id: match for match in exact[name]}
        if len(by_taxon) == 1:
            chosen[name] = next(iter(by_taxon.values()))
            continue
        if not by_taxon:
            alias_diagnostics = ";".join(
                sorted({f"{item.name_class}:{item.taxon_id}" for item in secondary[name]})
            )
            failures.append(
                {
                    "code": "UNRESOLVED_SCIENTIFIC_NAME",
                    "input_name": name,
                    "diagnostic_secondary_matches": alias_diagnostics,
                }
            )
            continue
        failures.append(
            {
                "code": "AMBIGUOUS_SCIENTIFIC_NAME",
                "input_name": name,
                "candidate_taxon_ids": ";".join(sorted(by_taxon, key=int)),
                "unique_names": ";".join(
                    sorted({match.unique_name for match in by_taxon.values() if match.unique_name})
                ),
            }
        )
    if failures:
        codes = sorted({failure["code"] for failure in failures})
        raise TaxonomyError(
            "+".join(codes),
            "Exact scientific-name resolution failed: "
            + "; ".join(
                f"{failure['input_name']!r}={failure['code']}" for failure in failures
            ),
            details=failures,
        )

    nodes, nodes_sha256 = _scan_nodes(nodes_path, (match.taxon_id for match in chosen.values()))
    missing_nodes = sorted(set(match.taxon_id for match in chosen.values()) - set(nodes), key=int)
    if missing_nodes:
        raise TaxonomyError(
            "NODE_MISSING",
            "Resolved TaxID missing from nodes.dmp in the same snapshot: " + ", ".join(missing_nodes),
        )

    rows: List[Mapping[str, str]] = []
    mismatches: List[Mapping[str, str]] = []
    for request in requests:
        match = chosen[request.input_name]
        node = nodes[match.taxon_id]
        if request.expected_taxon_id and request.expected_taxon_id != match.taxon_id:
            mismatches.append(
                {
                    "code": "TAXON_ID_MISMATCH",
                    "record_id": request.record_id,
                    "input_name": request.input_name,
                    "requested_taxon_id": request.expected_taxon_id,
                    "resolved_taxon_id": match.taxon_id,
                }
            )
        rows.append(
            {
                "record_id": request.record_id,
                "input_name": request.input_name,
                "requested_taxon_id": request.expected_taxon_id,
                "matched_name": match.matched_name,
                "name_class": match.name_class,
                "taxon_id": match.taxon_id,
                "parent_taxon_id": node.parent_taxon_id,
                "rank": node.rank,
                "unique_name": match.unique_name,
                "status": "resolved-exact-scientific-name",
                "snapshot": snapshot,
                "source_url": source_url,
                "retrieved_at": retrieved_at,
                "names_sha256": names_sha256,
                "nodes_sha256": nodes_sha256,
                "resolver_version": VERSION,
            }
        )
    if mismatches:
        raise TaxonomyError(
            "TAXON_ID_MISMATCH",
            "Resolved scientific name disagrees with a supplied TaxID: "
            + "; ".join(
                f"{item['record_id']}={item['requested_taxon_id']}->{item['resolved_taxon_id']}"
                for item in mismatches
            ),
            details=mismatches,
        )
    return TaxonomyResolution(
        rows=tuple(rows),
        names_sha256=names_sha256,
        nodes_sha256=nodes_sha256,
        snapshot=snapshot,
        source_url=source_url,
        retrieved_at=retrieved_at,
        resolver_version=VERSION,
    )


def write_resolution_tsv(path: Path, resolution: TaxonomyResolution) -> None:
    if path.exists():
        raise TaxonomyError("OUTPUT_EXISTS", f"Refusing to overwrite existing output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    handle_fd, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent))
    try:
        with os.fdopen(handle_fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=RESOLUTION_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(resolution.rows)
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except OSError:
            pass
        raise


def read_name_requests(path: Path, id_column: str, name_column: str, taxid_column: str) -> List[NameRequest]:
    try:
        handle = path.open("r", encoding="utf-8", newline="")
    except OSError as exc:
        raise TaxonomyError("INPUT_IO_ERROR", f"Cannot open input TSV: {exc}") from exc
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        required = [id_column, name_column]
        if taxid_column:
            required.append(taxid_column)
        missing = [column for column in required if column not in fields]
        if missing:
            raise TaxonomyError(
                "MISSING_COLUMNS", "Input TSV is missing columns: " + ", ".join(missing)
            )
        return [
            NameRequest(
                record_id=row[id_column],
                input_name=row[name_column],
                expected_taxon_id=row[taxid_column] if taxid_column else "",
            )
            for row in reader
        ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Resolve exact scientific names using local NCBI names.dmp and nodes.dmp files."
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("--names", required=True, help="Path to names.dmp from one NCBI snapshot.")
    parser.add_argument("--nodes", required=True, help="Path to nodes.dmp from the same snapshot.")
    parser.add_argument("--snapshot", required=True, help="Recorded snapshot date or archive label.")
    parser.add_argument(
        "--source-url", required=True, help="Exact NCBI archive URL used for this snapshot."
    )
    parser.add_argument("--retrieved-at", required=True, help="UTC retrieval date/time or date.")
    parser.add_argument("--out", required=True, help="New taxonomy_resolution.tsv path.")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--name", help="Resolve one exact scientific name.")
    source.add_argument("--input", help="Batch input TSV.")
    parser.add_argument("--expected-taxid", default="", help="Expected TaxID for --name.")
    parser.add_argument("--record-id", default="query", help="Record ID for --name.")
    parser.add_argument("--id-column", default="accession", help="Batch stable-ID column.")
    parser.add_argument("--name-column", default="species", help="Batch scientific-name column.")
    parser.add_argument("--taxid-column", default="taxon_id", help="Batch expected-TaxID column.")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.name is not None:
            requests = [NameRequest(args.record_id, args.name, args.expected_taxid)]
        else:
            requests = read_name_requests(
                Path(args.input), args.id_column, args.name_column, args.taxid_column
            )
        resolution = resolve_exact_scientific_names(
            requests,
            names_path=Path(args.names),
            nodes_path=Path(args.nodes),
            snapshot=args.snapshot,
            source_url=args.source_url,
            retrieved_at=args.retrieved_at,
        )
        output = Path(args.out)
        write_resolution_tsv(output, resolution)
        print(
            json.dumps(
                {
                    "status": "resolved",
                    "records": len(resolution.rows),
                    "match_mode": "exact-scientific-name",
                    "output": str(output),
                    "names_sha256": resolution.names_sha256,
                    "nodes_sha256": resolution.nodes_sha256,
                    "network_calls": 0,
                },
                sort_keys=True,
            )
        )
        return 0
    except TaxonomyError as exc:
        print(f"ERROR [{exc.code}] {exc.message}", file=sys.stderr)
        if exc.details:
            print(json.dumps({"details": exc.details}, ensure_ascii=False, sort_keys=True), file=sys.stderr)
        return 2
    except OSError as exc:
        print(f"ERROR [IO_ERROR] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
