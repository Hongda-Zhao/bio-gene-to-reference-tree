#!/usr/bin/env python3
"""Export complete NCBI lineages as post-plan, non-decision-bearing evidence.

The input ``taxonomy_resolution.tsv`` remains authoritative for the already
completed exact-name gate.  This script does not select, reject, reorder, or
modify sequences.  It verifies that ``names.dmp`` and ``nodes.dmp`` have the
exact SHA-256 digests recorded by every resolution row, validates each focal
node against that row, and follows parent links to TaxID 1 without guessing.

The full nodes table is indexed in a temporary SQLite database so the standard-
library-only implementation does not need to hold the current multi-million-
node NCBI taxonomy in Python objects.  The database is discarded after the
artifact is built and no temporary path is written to the result.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sqlite3
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


VERSION = "0.1.0"
ARTIFACT_ROLE = "post-plan-non-decision-bearing-taxonomy-evidence"
ASCII_TAXON_ID = re.compile(r"[0-9]+")
SHA256 = re.compile(r"[0-9a-f]{64}")
REQUIRED_RESOLUTION_COLUMNS = frozenset(
    {
        "record_id",
        "matched_name",
        "name_class",
        "taxon_id",
        "parent_taxon_id",
        "rank",
        "status",
        "snapshot",
        "source_url",
        "retrieved_at",
        "names_sha256",
        "nodes_sha256",
    }
)
OUTPUT_COLUMNS = (
    "accession",
    "taxon_id",
    "scientific_name",
    "rank",
    "parent_taxon_id",
    "lineage_node_count",
    "lineage_taxon_ids_json",
    "lineage_scientific_names_json",
    "lineage_ranks_json",
    "ranked_lineage_json",
    "artifact_role",
    "snapshot",
    "source_url",
    "retrieved_at",
    "names_sha256",
    "nodes_sha256",
    "lineage_generator_version",
)


class LineageError(RuntimeError):
    """A fail-closed taxonomy-lineage validation error."""


@dataclass(frozen=True)
class ResolutionRow:
    accession: str
    matched_name: str
    taxon_id: str
    parent_taxon_id: str
    rank: str
    snapshot: str
    source_url: str
    retrieved_at: str
    names_sha256: str
    nodes_sha256: str


@dataclass(frozen=True)
class LineageNode:
    taxon_id: str
    parent_taxon_id: str
    rank: str
    scientific_name: str


def fail(message: str) -> None:
    raise LineageError(message)


def parse_dmp_fields(
    raw_line: bytes, path: Path, line_number: int, minimum_fields: int
) -> list[str]:
    try:
        text = raw_line.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise LineageError(
            f"{path.name}:{line_number}: record is not valid UTF-8"
        ) from exc
    if text.endswith("\t|\r\n"):
        payload = text[:-4]
    elif text.endswith("\t|\n"):
        payload = text[:-3]
    elif text.endswith("\t|"):
        payload = text[:-2]
    else:
        fail(
            f"{path.name}:{line_number}: expected the official NCBI record "
            "terminator TAB|"
        )
    fields = payload.split("\t|\t")
    if len(fields) < minimum_fields:
        fail(
            f"{path.name}:{line_number}: expected at least {minimum_fields} "
            f"fields, found {len(fields)}"
        )
    return fields


def checked_text(value: str, context: str) -> str:
    if not value or value != value.strip() or any(char in value for char in "\t\r\n"):
        fail(f"{context} must be non-empty and free of surrounding whitespace")
    return value


def checked_taxon_id(value: str, context: str) -> str:
    if not ASCII_TAXON_ID.fullmatch(value) or int(value) < 1:
        fail(f"{context} must be a positive ASCII-decimal TaxID: {value!r}")
    return value


def read_resolution(path: Path) -> list[ResolutionRow]:
    rows: list[ResolutionRow] = []
    seen_accessions: set[str] = set()
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                fail(f"taxonomy resolution table has no header: {path}")
            missing = sorted(REQUIRED_RESOLUTION_COLUMNS - set(reader.fieldnames))
            if missing:
                fail("taxonomy resolution table is missing columns: " + ", ".join(missing))
            for line_number, raw_row in enumerate(reader, start=2):
                if None in raw_row:
                    fail(f"taxonomy resolution row {line_number} has extra columns")
                accession = checked_text(
                    raw_row["record_id"], f"taxonomy resolution row {line_number} record_id"
                )
                if accession in seen_accessions:
                    fail(f"duplicate taxonomy resolution record_id: {accession}")
                seen_accessions.add(accession)
                if raw_row["name_class"] != "scientific name":
                    fail(
                        f"{accession}: name_class must be exactly 'scientific name', "
                        f"found {raw_row['name_class']!r}"
                    )
                if raw_row["status"] != "resolved-exact-scientific-name":
                    fail(
                        f"{accession}: resolution status is not resolved-exact-scientific-name: "
                        f"{raw_row['status']!r}"
                    )
                names_sha256 = raw_row["names_sha256"]
                nodes_sha256 = raw_row["nodes_sha256"]
                if not SHA256.fullmatch(names_sha256):
                    fail(f"{accession}: names_sha256 is not a lowercase SHA-256 digest")
                if not SHA256.fullmatch(nodes_sha256):
                    fail(f"{accession}: nodes_sha256 is not a lowercase SHA-256 digest")
                rows.append(
                    ResolutionRow(
                        accession=accession,
                        matched_name=checked_text(
                            raw_row["matched_name"], f"{accession} matched_name"
                        ),
                        taxon_id=checked_taxon_id(raw_row["taxon_id"], f"{accession} taxon_id"),
                        parent_taxon_id=checked_taxon_id(
                            raw_row["parent_taxon_id"], f"{accession} parent_taxon_id"
                        ),
                        rank=checked_text(raw_row["rank"], f"{accession} rank"),
                        snapshot=checked_text(raw_row["snapshot"], f"{accession} snapshot"),
                        source_url=checked_text(raw_row["source_url"], f"{accession} source_url"),
                        retrieved_at=checked_text(
                            raw_row["retrieved_at"], f"{accession} retrieved_at"
                        ),
                        names_sha256=names_sha256,
                        nodes_sha256=nodes_sha256,
                    )
                )
    except (OSError, UnicodeDecodeError) as exc:
        fail(f"cannot read taxonomy resolution table {path}: {exc}")
    if not rows:
        fail(f"taxonomy resolution table contains no records: {path}")

    for attribute in (
        "snapshot",
        "source_url",
        "retrieved_at",
        "names_sha256",
        "nodes_sha256",
    ):
        values = {getattr(row, attribute) for row in rows}
        if len(values) != 1:
            fail(f"taxonomy resolution rows disagree on {attribute}")
    return rows


def load_nodes(path: Path, connection: sqlite3.Connection) -> str:
    digest = hashlib.sha256()
    connection.executescript(
        """
        PRAGMA journal_mode = OFF;
        PRAGMA synchronous = OFF;
        PRAGMA temp_store = MEMORY;
        CREATE TABLE nodes (
            taxon_id INTEGER PRIMARY KEY,
            parent_taxon_id INTEGER NOT NULL,
            rank TEXT NOT NULL
        ) WITHOUT ROWID;
        """
    )
    insert = "INSERT INTO nodes(taxon_id, parent_taxon_id, rank) VALUES (?, ?, ?)"
    count = 0
    try:
        with path.open("rb") as handle, connection:
            for line_number, raw_line in enumerate(handle, start=1):
                digest.update(raw_line)
                fields = parse_dmp_fields(raw_line, path, line_number, 3)
                taxon_id = checked_taxon_id(fields[0], f"{path.name}:{line_number} tax_id")
                parent_taxon_id = checked_taxon_id(
                    fields[1], f"{path.name}:{line_number} parent tax_id"
                )
                rank = checked_text(fields[2], f"{path.name}:{line_number} rank")
                try:
                    connection.execute(insert, (int(taxon_id), int(parent_taxon_id), rank))
                except sqlite3.IntegrityError as exc:
                    raise LineageError(
                        f"{path.name}:{line_number}: duplicate TaxID {taxon_id}"
                    ) from exc
                count += 1
    except OSError as exc:
        fail(f"cannot read nodes.dmp {path}: {exc}")
    if count == 0:
        fail(f"nodes.dmp contains no records: {path}")
    return digest.hexdigest()


def trace_lineages(
    rows: Sequence[ResolutionRow], connection: sqlite3.Connection
) -> dict[str, list[tuple[str, str, str]]]:
    """Return root-to-focal tuples of TaxID, parent TaxID, and rank."""

    lineages: dict[str, list[tuple[str, str, str]]] = {}
    for resolution in rows:
        reverse: list[tuple[str, str, str]] = []
        seen: set[str] = set()
        current = resolution.taxon_id
        while True:
            if current in seen:
                cycle = " -> ".join([node[0] for node in reverse] + [current])
                fail(f"{resolution.accession}: cycle in nodes.dmp parent chain: {cycle}")
            seen.add(current)
            result = connection.execute(
                "SELECT parent_taxon_id, rank FROM nodes WHERE taxon_id = ?",
                (int(current),),
            ).fetchone()
            if result is None:
                fail(
                    f"{resolution.accession}: TaxID {current} is missing from nodes.dmp "
                    "while traversing to root"
                )
            parent, rank = str(result[0]), str(result[1])
            reverse.append((current, parent, rank))
            if current == "1":
                if parent != "1":
                    fail(f"TaxID 1 must be self-parented, found parent TaxID {parent}")
                break
            if parent == current:
                fail(
                    f"{resolution.accession}: non-root TaxID {current} is self-parented"
                )
            current = parent

        lineage = list(reversed(reverse))
        focal_taxon, focal_parent, focal_rank = lineage[-1]
        if focal_taxon != resolution.taxon_id:
            fail(f"internal lineage endpoint error for {resolution.accession}")
        if focal_parent != resolution.parent_taxon_id:
            fail(
                f"{resolution.accession}: resolution parent_taxon_id "
                f"{resolution.parent_taxon_id} disagrees with nodes.dmp {focal_parent}"
            )
        if focal_rank != resolution.rank:
            fail(
                f"{resolution.accession}: resolution rank {resolution.rank!r} "
                f"disagrees with nodes.dmp {focal_rank!r}"
            )
        lineages[resolution.accession] = lineage
    return lineages


def load_scientific_names(path: Path, wanted_taxon_ids: Iterable[str]) -> tuple[dict[str, str], str]:
    wanted = set(wanted_taxon_ids)
    scientific_names: dict[str, str] = {}
    digest = hashlib.sha256()
    count = 0
    try:
        with path.open("rb") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                digest.update(raw_line)
                fields = parse_dmp_fields(raw_line, path, line_number, 4)
                taxon_id = checked_taxon_id(fields[0], f"{path.name}:{line_number} tax_id")
                name_text, name_class = fields[1], fields[3]
                if name_class != "scientific name" or taxon_id not in wanted:
                    continue
                checked_text(name_text, f"{path.name}:{line_number} scientific name")
                if taxon_id in scientific_names:
                    fail(
                        f"{path.name}:{line_number}: duplicate scientific-name record "
                        f"for TaxID {taxon_id}"
                    )
                scientific_names[taxon_id] = name_text
                count += 1
    except OSError as exc:
        fail(f"cannot read names.dmp {path}: {exc}")
    if count == 0:
        fail(f"names.dmp contains no requested scientific-name records: {path}")
    missing = sorted(wanted - set(scientific_names), key=int)
    if missing:
        fail(
            "lineage TaxIDs missing a unique scientific name in names.dmp: "
            + ", ".join(missing)
        )
    return scientific_names, digest.hexdigest()


def canonical_json(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def build_output_rows(
    resolutions: Sequence[ResolutionRow],
    lineages: dict[str, list[tuple[str, str, str]]],
    scientific_names: dict[str, str],
) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    for resolution in sorted(resolutions, key=lambda row: row.accession):
        lineage = lineages[resolution.accession]
        nodes = [
            LineageNode(
                taxon_id=taxon_id,
                parent_taxon_id=parent_taxon_id,
                rank=rank,
                scientific_name=scientific_names[taxon_id],
            )
            for taxon_id, parent_taxon_id, rank in lineage
        ]
        focal = nodes[-1]
        if focal.scientific_name != resolution.matched_name:
            fail(
                f"{resolution.accession}: matched_name {resolution.matched_name!r} "
                f"disagrees with names.dmp {focal.scientific_name!r}"
            )
        ranked_lineage = [
            {
                "taxon_id": node.taxon_id,
                "scientific_name": node.scientific_name,
                "rank": node.rank,
            }
            for node in nodes
        ]
        output.append(
            {
                "accession": resolution.accession,
                "taxon_id": focal.taxon_id,
                "scientific_name": focal.scientific_name,
                "rank": focal.rank,
                "parent_taxon_id": focal.parent_taxon_id,
                "lineage_node_count": str(len(nodes)),
                "lineage_taxon_ids_json": canonical_json(
                    [node.taxon_id for node in nodes]
                ),
                "lineage_scientific_names_json": canonical_json(
                    [node.scientific_name for node in nodes]
                ),
                "lineage_ranks_json": canonical_json([node.rank for node in nodes]),
                "ranked_lineage_json": canonical_json(ranked_lineage),
                "artifact_role": ARTIFACT_ROLE,
                "snapshot": resolution.snapshot,
                "source_url": resolution.source_url,
                "retrieved_at": resolution.retrieved_at,
                "names_sha256": resolution.names_sha256,
                "nodes_sha256": resolution.nodes_sha256,
                "lineage_generator_version": VERSION,
            }
        )
    return output


def write_output(path: Path, rows: Sequence[dict[str, str]]) -> None:
    if path.exists():
        fail(f"refusing to overwrite output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
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
        fail(f"cannot write output {path}: {exc}")


def run(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Export complete frozen-NCBI lineages as non-decision-bearing evidence."
    )
    parser.add_argument("--taxonomy-resolution", required=True, type=Path)
    parser.add_argument("--names-dmp", required=True, type=Path)
    parser.add_argument("--nodes-dmp", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)

    if args.output.exists():
        fail(f"refusing to overwrite output: {args.output}")
    resolutions = read_resolution(args.taxonomy_resolution)

    with tempfile.TemporaryDirectory(prefix="taxonomy-lineage-") as temporary:
        database = Path(temporary) / "nodes.sqlite3"
        with sqlite3.connect(database) as connection:
            nodes_sha256 = load_nodes(args.nodes_dmp, connection)
            expected_nodes_sha256 = resolutions[0].nodes_sha256
            if nodes_sha256 != expected_nodes_sha256:
                fail(
                    "nodes.dmp SHA-256 disagrees with taxonomy_resolution.tsv: "
                    f"expected {expected_nodes_sha256}, observed {nodes_sha256}"
                )
            lineages = trace_lineages(resolutions, connection)

        lineage_taxon_ids = {
            taxon_id
            for lineage in lineages.values()
            for taxon_id, _parent_taxon_id, _rank in lineage
        }
        scientific_names, names_sha256 = load_scientific_names(
            args.names_dmp, lineage_taxon_ids
        )
        expected_names_sha256 = resolutions[0].names_sha256
        if names_sha256 != expected_names_sha256:
            fail(
                "names.dmp SHA-256 disagrees with taxonomy_resolution.tsv: "
                f"expected {expected_names_sha256}, observed {names_sha256}"
            )

    rows = build_output_rows(resolutions, lineages, scientific_names)
    write_output(args.output, rows)
    print(
        f"Wrote {len(rows)} complete frozen-taxonomy lineages to {args.output}; "
        f"artifact_role={ARTIFACT_ROLE}"
    )
    return 0


def main() -> None:
    try:
        raise SystemExit(run())
    except (LineageError, sqlite3.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc


if __name__ == "__main__":
    main()
