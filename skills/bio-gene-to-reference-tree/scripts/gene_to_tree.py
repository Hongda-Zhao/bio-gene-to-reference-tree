#!/usr/bin/env python3
"""Compile a deterministic review bundle for an auditable protein gene tree.

The helper deliberately performs no network requests and launches no external
programs while planning.  A host agent may resolve a query and acquire candidate
records with authorized tools, then hand the materialized local bundle to this
script for deterministic selection, command planning, iTOL annotation, and
provenance capture.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple

from ncbi_taxonomy import (
    NameRequest,
    TaxonomyError,
    TaxonomyResolution,
    resolve_exact_scientific_names,
    write_resolution_tsv,
)


VERSION = "0.3.0"
OUTPUT_SCHEMA_VERSION = "0.3"
PROTEIN_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWYBXZJUO")
SAFE_ID = re.compile(r"^[A-Za-z0-9_.:-]+$")
SAFE_PROFILE_ID = re.compile(r"^[A-Za-z0-9_-]+$")
ASCII_TAXON_ID = re.compile(r"[0-9]+")
RELATION_PRIORITY = {
    "one2one_ortholog": 0,
    "ortholog": 1,
    "coortholog": 2,
    "homolog": 3,
    "paralog": 4,
    "self": 5,
}
REQUIRED_COLUMNS = (
    "accession",
    "taxon_id",
    "species",
    "role",
    "relation",
    "is_reviewed",
    "is_canonical",
    "is_fragment",
    "query_coverage",
    "sequence_length",
    "bitscore",
    "evalue",
    "source_db",
    "source_release",
    "retrieved_at",
    "clade",
)
OPTIONAL_COLUMNS = (
    "accession_version",
    "gene_name",
    "protein_name",
    "lineage",
    "analysis_group",
    "target_coverage",
    "percent_identity",
    "alignment_length",
    "orthology_source",
    "orthology_evidence",
    "domain_architecture",
    "cluster_id",
    "cluster_representative",
    "outgroup_rationale",
    "retrieval_query_id",
    "notes",
)
BASE_OUTPUT_COLUMNS = REQUIRED_COLUMNS + OPTIONAL_COLUMNS + ("sequence_sha256",)
ITOL_COLORS = {
    "study": "#E69F00",
    "expanded": "#009E73",
    "outgroup": "#999999",
}


class WorkflowError(Exception):
    """A user-facing validation or workflow error with a stable code."""

    def __init__(self, code: str, message: str, exit_code: int = 2) -> None:
        super().__init__(message)
        self.code = code
        self.message = message
        self.exit_code = exit_code


@dataclass(frozen=True)
class Candidate:
    accession: str
    taxon_id: str
    species: str
    role: str
    relation: str
    is_reviewed: bool
    is_canonical: bool
    is_fragment: bool
    query_coverage: float
    sequence_length: int
    bitscore: float
    evalue: float
    source_db: str
    source_release: str
    retrieved_at: str
    clade: str
    sequence: str
    accession_version: str = ""
    gene_name: str = ""
    protein_name: str = ""
    lineage: str = ""
    analysis_group: str = ""
    target_coverage: str = ""
    percent_identity: str = ""
    alignment_length: str = ""
    orthology_source: str = ""
    orthology_evidence: str = ""
    domain_architecture: str = ""
    cluster_id: str = ""
    cluster_representative: str = ""
    outgroup_rationale: str = ""
    retrieval_query_id: str = ""
    notes: str = ""

    @property
    def sequence_sha256(self) -> str:
        return sha256_text(self.sequence)

    @property
    def analysis_role(self) -> str:
        if self.analysis_group:
            return self.analysis_group
        if self.relation == "self":
            return "study"
        if self.role == "outgroup":
            return "outgroup"
        return "expanded"

    def rank_key(self) -> Tuple[Any, ...]:
        return (
            RELATION_PRIORITY[self.relation],
            -int(self.is_reviewed),
            -int(self.is_canonical),
            -self.query_coverage,
            -self.bitscore,
            self.accession,
        )

    def output_row(self) -> Dict[str, str]:
        return {
            "accession": self.accession,
            "taxon_id": self.taxon_id,
            "species": self.species,
            "role": self.role,
            "relation": self.relation,
            "is_reviewed": bool_text(self.is_reviewed),
            "is_canonical": bool_text(self.is_canonical),
            "is_fragment": bool_text(self.is_fragment),
            "query_coverage": stable_number(self.query_coverage),
            "sequence_length": str(self.sequence_length),
            "bitscore": stable_number(self.bitscore),
            "evalue": stable_number(self.evalue),
            "source_db": self.source_db,
            "source_release": self.source_release,
            "retrieved_at": self.retrieved_at,
            "clade": self.clade,
            "accession_version": self.accession_version,
            "gene_name": self.gene_name,
            "protein_name": self.protein_name,
            "lineage": self.lineage,
            "analysis_group": self.analysis_role,
            "target_coverage": self.target_coverage,
            "percent_identity": self.percent_identity,
            "alignment_length": self.alignment_length,
            "orthology_source": self.orthology_source,
            "orthology_evidence": self.orthology_evidence,
            "domain_architecture": self.domain_architecture,
            "cluster_id": self.cluster_id,
            "cluster_representative": self.cluster_representative,
            "outgroup_rationale": self.outgroup_rationale,
            "retrieval_query_id": self.retrieval_query_id,
            "notes": self.notes,
            "sequence_sha256": self.sequence_sha256,
        }

    def fingerprint_record(self) -> Dict[str, Any]:
        row: Dict[str, Any] = dict(self.output_row())
        return row


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def stable_number(value: float) -> str:
    return format(value, ".15g")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_text(text: str) -> str:
    return sha256_bytes(text.encode("utf-8"))


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise WorkflowError("INPUT_READ_ERROR", f"Cannot read {path.name}: {exc}") from exc
    return digest.hexdigest()


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def read_json(path: Path) -> Dict[str, Any]:
    try:
        raw = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise WorkflowError("REQUEST_READ_ERROR", f"Cannot read request {path}: {exc}") from exc
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise WorkflowError(
            "INVALID_REQUEST_JSON",
            f"Request is not valid JSON at line {exc.lineno}, column {exc.colno}: {exc.msg}",
        ) from exc
    if not isinstance(value, dict):
        raise WorkflowError("INVALID_REQUEST", "The request root must be a JSON object.")
    return value


def mapping_at(parent: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = parent.get(key)
    if not isinstance(value, dict):
        raise WorkflowError("INVALID_REQUEST", f"'{key}' must be an object.")
    return value


def string_at(parent: Mapping[str, Any], key: str, *, allow_empty: bool = False) -> str:
    value = parent.get(key)
    if not isinstance(value, str) or (not allow_empty and not value.strip()):
        raise WorkflowError("INVALID_REQUEST", f"'{key}' must be a non-empty string.")
    return value.strip()


def bool_at(parent: Mapping[str, Any], key: str) -> bool:
    value = parent.get(key)
    if not isinstance(value, bool):
        raise WorkflowError("INVALID_REQUEST", f"'{key}' must be true or false.")
    return value


def int_at(parent: Mapping[str, Any], key: str, minimum: int = 0) -> int:
    value = parent.get(key)
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise WorkflowError("INVALID_REQUEST", f"'{key}' must be an integer >= {minimum}.")
    return value


def float_at(parent: Mapping[str, Any], key: str, minimum: float = 0.0) -> float:
    value = parent.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise WorkflowError("INVALID_REQUEST", f"'{key}' must be a number.")
    result = float(value)
    if not math.isfinite(result) or result < minimum:
        raise WorkflowError("INVALID_REQUEST", f"'{key}' must be finite and >= {minimum}.")
    return result


def resolve_input_path(request_path: Path, configured: str, label: str) -> Path:
    candidate = Path(configured).expanduser()
    if not candidate.is_absolute():
        candidate = request_path.parent / candidate
    candidate = candidate.resolve()
    if not candidate.is_file():
        raise WorkflowError("MISSING_INPUT", f"{label} does not exist or is not a file: {configured}")
    return candidate


def parse_fasta(path: Path, label: str) -> List[Tuple[str, str, str]]:
    records: List[Tuple[str, str, str]] = []
    current_header: str | None = None
    chunks: List[str] = []

    def finish_record() -> None:
        nonlocal current_header, chunks
        if current_header is None:
            return
        identifier = current_header.split()[0]
        if not SAFE_ID.fullmatch(identifier):
            raise WorkflowError(
                "UNSAFE_SEQUENCE_ID",
                f"{label} FASTA identifier '{identifier}' contains unsupported characters.",
            )
        sequence = "".join(chunks).replace(" ", "").replace("\t", "").upper()
        if not sequence:
            raise WorkflowError("EMPTY_SEQUENCE", f"{label} FASTA record '{identifier}' has no sequence.")
        invalid = sorted(set(sequence) - PROTEIN_ALPHABET)
        if invalid:
            shown = "".join(invalid)
            raise WorkflowError(
                "INVALID_PROTEIN_SEQUENCE",
                f"{label} FASTA record '{identifier}' contains unsupported symbols: {shown}",
            )
        records.append((identifier, current_header, sequence))
        current_header = None
        chunks = []

    try:
        with path.open("r", encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    finish_record()
                    current_header = line[1:].strip()
                    if not current_header:
                        raise WorkflowError(
                            "INVALID_FASTA", f"{label} FASTA has an empty header at line {line_number}."
                        )
                else:
                    if current_header is None:
                        raise WorkflowError(
                            "INVALID_FASTA",
                            f"{label} FASTA has sequence data before its first header at line {line_number}.",
                        )
                    chunks.append(line)
        finish_record()
    except UnicodeDecodeError as exc:
        raise WorkflowError("INVALID_TEXT_ENCODING", f"{label} FASTA must be UTF-8 text.") from exc
    except OSError as exc:
        raise WorkflowError("INPUT_READ_ERROR", f"Cannot read {label} FASTA: {exc}") from exc

    if not records:
        raise WorkflowError("EMPTY_FASTA", f"{label} FASTA contains no records.")
    ids = [record[0] for record in records]
    duplicates = sorted({identifier for identifier in ids if ids.count(identifier) > 1})
    if duplicates:
        raise WorkflowError("DUPLICATE_SEQUENCE_ID", f"Duplicate {label} FASTA IDs: {', '.join(duplicates)}")
    return records


def parse_bool(value: str, accession: str, field: str) -> bool:
    normalized = value.strip().lower()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    raise WorkflowError("INVALID_CANDIDATE_TABLE", f"{accession}: {field} must be true or false.")


def parse_finite_float(value: str, accession: str, field: str, minimum: float = 0.0) -> float:
    try:
        result = float(value)
    except ValueError as exc:
        raise WorkflowError("INVALID_CANDIDATE_TABLE", f"{accession}: {field} must be numeric.") from exc
    if not math.isfinite(result) or result < minimum:
        raise WorkflowError(
            "INVALID_CANDIDATE_TABLE", f"{accession}: {field} must be finite and >= {minimum}."
        )
    return result


def normalize_optional_float(
    value: str,
    accession: str,
    field: str,
    *,
    minimum: float = 0.0,
    maximum: float | None = None,
) -> str:
    if not value:
        return ""
    result = parse_finite_float(value, accession, field, minimum)
    if maximum is not None and result > maximum:
        raise WorkflowError(
            "INVALID_CANDIDATE_TABLE",
            f"{accession}: {field} must be <= {stable_number(maximum)}.",
        )
    return stable_number(result)


def parse_candidate_table(path: Path, sequences: Mapping[str, str]) -> List[Candidate]:
    candidates: List[Candidate] = []
    seen: set[str] = set()
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise WorkflowError("INVALID_CANDIDATE_TABLE", "Candidate table has no header.")
            missing = [column for column in REQUIRED_COLUMNS if column not in reader.fieldnames]
            if missing:
                raise WorkflowError(
                    "INVALID_CANDIDATE_TABLE",
                    "Candidate table is missing columns: " + ", ".join(missing),
                )
            for row_number, raw in enumerate(reader, start=2):
                if None in raw:
                    raise WorkflowError(
                        "INVALID_CANDIDATE_TABLE", f"Candidate table row {row_number} has extra columns."
                    )
                for exact_field in ("species", "taxon_id"):
                    exact_value = raw.get(exact_field) or ""
                    if exact_value != exact_value.strip():
                        raise WorkflowError(
                            "NONEXACT_TAXONOMY_FIELD",
                            f"Candidate row {row_number} field '{exact_field}' has leading or trailing "
                            "whitespace; exact taxonomy matching never normalizes it silently.",
                        )
                row = {
                    key: (raw.get(key) or "").strip()
                    for key in REQUIRED_COLUMNS + OPTIONAL_COLUMNS
                }
                accession = row["accession"]
                if not SAFE_ID.fullmatch(accession):
                    raise WorkflowError(
                        "UNSAFE_SEQUENCE_ID", f"Candidate row {row_number} has unsafe accession '{accession}'."
                    )
                if accession in seen:
                    raise WorkflowError("DUPLICATE_ACCESSION", f"Duplicate candidate accession: {accession}")
                seen.add(accession)
                if accession not in sequences:
                    raise WorkflowError(
                        "CANDIDATE_SEQUENCE_MISMATCH",
                        f"Candidate table accession '{accession}' is absent from the candidate FASTA.",
                    )
                if not ASCII_TAXON_ID.fullmatch(row["taxon_id"]):
                    raise WorkflowError(
                        "INVALID_CANDIDATE_TABLE", f"{accession}: taxon_id must contain only digits."
                    )
                for field in ("species", "source_db"):
                    if not row[field]:
                        raise WorkflowError(
                            "INVALID_CANDIDATE_TABLE", f"{accession}: {field} must not be empty."
                        )
                if row["role"] not in {"ingroup", "outgroup"}:
                    raise WorkflowError(
                        "INVALID_CANDIDATE_TABLE", f"{accession}: role must be ingroup or outgroup."
                    )
                if row["relation"] not in RELATION_PRIORITY:
                    raise WorkflowError(
                        "INVALID_CANDIDATE_TABLE", f"{accession}: unsupported relation '{row['relation']}'."
                    )
                if row["analysis_group"] not in {"", "study", "expanded", "outgroup"}:
                    raise WorkflowError(
                        "INVALID_CANDIDATE_TABLE",
                        f"{accession}: analysis_group must be study, expanded, outgroup, or empty.",
                    )
                if row["role"] == "outgroup" and row["analysis_group"] not in {"", "outgroup"}:
                    raise WorkflowError(
                        "INVALID_CANDIDATE_TABLE",
                        f"{accession}: an outgroup cannot use analysis_group='{row['analysis_group']}'.",
                    )
                if row["role"] == "ingroup" and row["analysis_group"] == "outgroup":
                    raise WorkflowError(
                        "INVALID_CANDIDATE_TABLE",
                        f"{accession}: an ingroup cannot use analysis_group='outgroup'.",
                    )
                if row["relation"] == "self" and row["analysis_group"] not in {"", "study"}:
                    raise WorkflowError(
                        "INVALID_CANDIDATE_TABLE",
                        f"{accession}: a self/query record must use analysis_group='study'.",
                    )
                try:
                    sequence_length = int(row["sequence_length"])
                except ValueError as exc:
                    raise WorkflowError(
                        "INVALID_CANDIDATE_TABLE", f"{accession}: sequence_length must be an integer."
                    ) from exc
                if sequence_length <= 0 or sequence_length != len(sequences[accession]):
                    raise WorkflowError(
                        "SEQUENCE_LENGTH_MISMATCH",
                        f"{accession}: table length {sequence_length} does not match FASTA length "
                        f"{len(sequences[accession])}.",
                    )
                candidates.append(
                    Candidate(
                        accession=accession,
                        taxon_id=row["taxon_id"],
                        species=row["species"],
                        role=row["role"],
                        relation=row["relation"],
                        is_reviewed=parse_bool(row["is_reviewed"], accession, "is_reviewed"),
                        is_canonical=parse_bool(row["is_canonical"], accession, "is_canonical"),
                        is_fragment=parse_bool(row["is_fragment"], accession, "is_fragment"),
                        query_coverage=parse_finite_float(
                            row["query_coverage"], accession, "query_coverage"
                        ),
                        sequence_length=sequence_length,
                        bitscore=parse_finite_float(row["bitscore"], accession, "bitscore"),
                        evalue=parse_finite_float(row["evalue"], accession, "evalue"),
                        source_db=row["source_db"],
                        source_release=row["source_release"],
                        retrieved_at=row["retrieved_at"],
                        clade=row["clade"],
                        sequence=sequences[accession],
                        accession_version=row["accession_version"],
                        gene_name=row["gene_name"],
                        protein_name=row["protein_name"],
                        lineage=row["lineage"],
                        analysis_group=row["analysis_group"],
                        target_coverage=normalize_optional_float(
                            row["target_coverage"], accession, "target_coverage", maximum=1.0
                        ),
                        percent_identity=normalize_optional_float(
                            row["percent_identity"], accession, "percent_identity", maximum=100.0
                        ),
                        alignment_length=normalize_optional_float(
                            row["alignment_length"], accession, "alignment_length"
                        ),
                        orthology_source=row["orthology_source"],
                        orthology_evidence=row["orthology_evidence"],
                        domain_architecture=row["domain_architecture"],
                        cluster_id=row["cluster_id"],
                        cluster_representative=row["cluster_representative"],
                        outgroup_rationale=row["outgroup_rationale"],
                        retrieval_query_id=row["retrieval_query_id"],
                        notes=row["notes"],
                    )
                )
                if candidates[-1].query_coverage > 1:
                    raise WorkflowError(
                        "INVALID_CANDIDATE_TABLE",
                        f"{accession}: query_coverage must be between 0 and 1.",
                    )
    except UnicodeDecodeError as exc:
        raise WorkflowError("INVALID_TEXT_ENCODING", "Candidate table must be UTF-8 text.") from exc
    except OSError as exc:
        raise WorkflowError("INPUT_READ_ERROR", f"Cannot read candidate table: {exc}") from exc

    extra_fasta = sorted(set(sequences) - seen)
    if extra_fasta:
        raise WorkflowError(
            "CANDIDATE_SEQUENCE_MISMATCH",
            "Candidate FASTA IDs absent from the table: " + ", ".join(extra_fasta),
        )
    if not candidates:
        raise WorkflowError("EMPTY_CANDIDATE_SET", "Candidate bundle contains no candidates.")

    roles_by_taxon: Dict[str, set[str]] = {}
    for candidate in candidates:
        roles_by_taxon.setdefault(candidate.taxon_id, set()).add(candidate.role)
    conflicts = sorted(taxon for taxon, roles in roles_by_taxon.items() if len(roles) > 1)
    if conflicts:
        raise WorkflowError(
            "CONFLICTING_TAXON_ROLES",
            "A taxon cannot be both ingroup and outgroup: " + ", ".join(conflicts),
        )
    return candidates


def optional_mapping(parent: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = parent.get(key)
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise WorkflowError("INVALID_REQUEST", f"'{key}' must be an object.")
    return value


def optional_string(parent: Mapping[str, Any], key: str, default: str = "") -> str:
    value = parent.get(key, default)
    if not isinstance(value, str):
        raise WorkflowError("INVALID_REQUEST", f"'{key}' must be a string.")
    return value.strip()


def optional_bool(parent: Mapping[str, Any], key: str, default: bool) -> bool:
    value = parent.get(key, default)
    if not isinstance(value, bool):
        raise WorkflowError("INVALID_REQUEST", f"'{key}' must be true or false.")
    return value


def optional_int(parent: Mapping[str, Any], key: str, default: int, minimum: int = 0) -> int:
    value = parent.get(key, default)
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise WorkflowError("INVALID_REQUEST", f"'{key}' must be an integer >= {minimum}.")
    return value


def optional_float(
    parent: Mapping[str, Any], key: str, default: float, minimum: float = 0.0
) -> float:
    value = parent.get(key, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise WorkflowError("INVALID_REQUEST", f"'{key}' must be a number.")
    result = float(value)
    if not math.isfinite(result) or result < minimum:
        raise WorkflowError("INVALID_REQUEST", f"'{key}' must be finite and >= {minimum}.")
    return result


def string_list(parent: Mapping[str, Any], key: str, default: Sequence[str]) -> List[str]:
    value = parent.get(key, list(default))
    if not isinstance(value, list) or not all(isinstance(item, str) and item.strip() for item in value):
        raise WorkflowError("INVALID_REQUEST", f"'{key}' must be an array of non-empty strings.")
    return [item.strip() for item in value]


def normalize_request(request: Mapping[str, Any]) -> Dict[str, Any]:
    request_schema = string_at(request, "schema_version")
    if request_schema not in {"0.1", "0.2"}:
        raise WorkflowError(
            "UNSUPPORTED_SCHEMA_VERSION", "Supported request schema versions are '0.1' and '0.2'."
        )
    project_id = string_at(request, "project_id")
    objective = string_at(request, "objective")
    if objective not in {"ortholog-tree", "homolog-context", "within-species"}:
        raise WorkflowError(
            "INVALID_OBJECTIVE",
            "objective must be ortholog-tree, homolog-context, or within-species.",
        )
    sequence_context = optional_string(request, "sequence_context", "cellular")
    if sequence_context not in {"cellular", "viral"}:
        raise WorkflowError("INVALID_SEQUENCE_CONTEXT", "sequence_context must be cellular or viral.")

    query = mapping_at(request, "query")
    query_kind = string_at(query, "kind")
    allowed_query_kinds = {"protein-fasta", "accession", "protein-name", "gene-symbol"}
    if query_kind not in allowed_query_kinds:
        raise WorkflowError(
            "OFFLINE_QUERY_KIND_UNSUPPORTED",
            "query.kind must be protein-fasta, accession, protein-name, or gene-symbol; "
            "all routes must provide a resolved local protein FASTA before planning.",
        )
    if query_kind != "protein-fasta" and (not query.get("id") or not query.get("path")):
        raise WorkflowError(
            "QUERY_RESOLUTION_REQUIRED",
            "Accession and name routes must be resolved to a local protein FASTA with query.id and query.path.",
        )
    query_id = string_at(query, "id")
    if not SAFE_ID.fullmatch(query_id):
        raise WorkflowError("UNSAFE_SEQUENCE_ID", f"Query ID '{query_id}' contains unsupported characters.")
    raw_query_organism = query.get("organism")
    if isinstance(raw_query_organism, str) and raw_query_organism != raw_query_organism.strip():
        raise WorkflowError(
            "NONEXACT_TAXONOMY_FIELD",
            "query.organism has leading or trailing whitespace; exact taxonomy matching never "
            "normalizes it silently.",
        )
    query_organism = string_at(query, "organism")
    query_taxon_id = optional_string(query, "taxon_id")
    if query_taxon_id and not ASCII_TAXON_ID.fullmatch(query_taxon_id):
        raise WorkflowError("INVALID_REQUEST", "query.taxon_id must contain only digits when provided.")
    if query_kind in {"protein-name", "gene-symbol"} and not (query_organism or query_taxon_id):
        raise WorkflowError(
            "AMBIGUOUS_NAME_WITHOUT_TAXON", "A name or symbol requires a source organism or TaxID."
        )

    taxon_scope = mapping_at(request, "taxon_scope")
    privacy = mapping_at(request, "privacy")
    references = mapping_at(request, "references")
    selection = mapping_at(request, "selection")
    alignment = mapping_at(request, "alignment")
    tree = mapping_at(request, "tree")

    strategy = string_at(references, "strategy")
    if strategy not in {"ortholog-first", "homolog-first", "local-bundle"}:
        raise WorkflowError(
            "INVALID_REFERENCE_STRATEGY",
            "strategy must be ortholog-first, homolog-first, or local-bundle.",
        )
    if objective == "ortholog-tree" and strategy not in {"ortholog-first", "local-bundle"}:
        raise WorkflowError(
            "INCOMPATIBLE_REFERENCE_STRATEGY",
            "ortholog-tree requires references.strategy='ortholog-first' or a curated local-bundle.",
        )

    normalized_selection = {
        "min_query_coverage": float_at(selection, "min_query_coverage"),
        "min_target_coverage": optional_float(selection, "min_target_coverage", 0.0),
        "min_length_ratio": float_at(selection, "min_length_ratio"),
        "max_length_ratio": float_at(selection, "max_length_ratio"),
        "max_per_taxon": int_at(selection, "max_per_taxon", 1),
        "max_references": int_at(selection, "max_references", 1),
        "min_ingroup_taxa": int_at(selection, "min_ingroup_taxa", 1),
        "require_outgroup": bool_at(selection, "require_outgroup"),
        "outgroup_count": int_at(selection, "outgroup_count", 0),
        "allow_paralogs": bool_at(selection, "allow_paralogs"),
    }
    if (
        normalized_selection["min_query_coverage"] > 1
        or normalized_selection["min_target_coverage"] > 1
    ):
        raise WorkflowError("INVALID_REQUEST", "query and target coverage thresholds must be <= 1.")
    if normalized_selection["min_length_ratio"] > normalized_selection["max_length_ratio"]:
        raise WorkflowError("INVALID_REQUEST", "min_length_ratio must not exceed max_length_ratio.")
    if normalized_selection["require_outgroup"] and normalized_selection["outgroup_count"] < 1:
        raise WorkflowError("INVALID_REQUEST", "require_outgroup=true requires outgroup_count >= 1.")
    if not normalized_selection["require_outgroup"] and normalized_selection["outgroup_count"] != 0:
        raise WorkflowError("INVALID_REQUEST", "require_outgroup=false requires outgroup_count=0.")
    ingroup_scope = string_at(taxon_scope, "ingroup")
    outgroup_scope = string_at(taxon_scope, "outgroup", allow_empty=True)
    if normalized_selection["require_outgroup"] and not outgroup_scope:
        raise WorkflowError("INVALID_REQUEST", "A required outgroup needs taxon_scope.outgroup.")
    minimum_capacity = normalized_selection["min_ingroup_taxa"]
    if normalized_selection["require_outgroup"]:
        minimum_capacity += normalized_selection["outgroup_count"]
    if normalized_selection["max_references"] < minimum_capacity:
        raise WorkflowError(
            "INVALID_REQUEST",
            "max_references is too small for the ingroup minimum plus required outgroups.",
        )

    clustering_raw = optional_mapping(request, "clustering")
    clustering_mode = optional_string(clustering_raw, "mode", "off")
    if clustering_mode not in {"off", "auto", "precomputed"}:
        raise WorkflowError("INVALID_CLUSTERING_PLAN", "clustering.mode must be off, auto, or precomputed.")
    clustering_min_seq_id = optional_float(clustering_raw, "min_seq_id", 0.95)
    clustering_coverage = optional_float(clustering_raw, "coverage", 0.8)
    if clustering_min_seq_id > 1 or clustering_coverage > 1:
        raise WorkflowError(
            "INVALID_CLUSTERING_PLAN", "min_seq_id and coverage must be between 0 and 1."
        )
    clustering_cov_mode = optional_int(clustering_raw, "coverage_mode", 0)
    if clustering_cov_mode not in {0, 1, 2, 3, 4, 5}:
        raise WorkflowError("INVALID_CLUSTERING_PLAN", "coverage_mode must be an MMseqs2 mode 0..5.")
    normalized_clustering = {
        "mode": clustering_mode,
        "tool": optional_string(clustering_raw, "tool", "mmseqs2"),
        "trigger_min_sequences": optional_int(clustering_raw, "trigger_min_sequences", 200, 2),
        "algorithm": optional_string(clustering_raw, "algorithm", "easy-linclust"),
        "min_seq_id": clustering_min_seq_id,
        "coverage": clustering_coverage,
        "coverage_mode": clustering_cov_mode,
        "threads": optional_int(clustering_raw, "threads", 4, 1),
        "preserve_analysis_groups": string_list(
            clustering_raw, "preserve_analysis_groups", ["study", "outgroup"]
        ),
    }
    if normalized_clustering["tool"] != "mmseqs2" or normalized_clustering["algorithm"] not in {
        "easy-linclust",
        "easy-cluster",
    }:
        raise WorkflowError(
            "INVALID_CLUSTERING_PLAN",
            "v0.2 clustering supports mmseqs2 easy-linclust or easy-cluster.",
        )
    if set(normalized_clustering["preserve_analysis_groups"]) != {"study", "outgroup"}:
        raise WorkflowError(
            "INVALID_CLUSTERING_PLAN",
            "Study and outgroup sequences must both be preserved outside MMseqs2 clustering.",
        )

    alignment_tool = string_at(alignment, "tool")
    alignment_mode = string_at(alignment, "mode")
    if alignment_tool != "mafft" or alignment_mode not in {"auto", "linsi", "ginsi", "einsi"}:
        raise WorkflowError(
            "UNSUPPORTED_ALIGNMENT_PLAN",
            "v0.2 supports MAFFT modes auto, linsi, ginsi, and einsi.",
        )
    normalized_alignment = {
        "tool": alignment_tool,
        "mode": alignment_mode,
        "threads": optional_int(alignment, "threads", 4, 1),
    }

    if request_schema == "0.1":
        trimming_raw: Mapping[str, Any] = {
            "enabled": bool_at(alignment, "trim"),
            "tool": "trimal",
            "primary_profile": "balanced",
            "profiles": [{"id": "balanced", "gap_threshold": 0.9}],
            "min_retained_fraction": 0.3,
            "compare_topologies": False,
        }
    else:
        trimming_raw = optional_mapping(request, "trimming")
    trimming_enabled = optional_bool(trimming_raw, "enabled", False)
    raw_profiles = trimming_raw.get("profiles", [])
    if not isinstance(raw_profiles, list):
        raise WorkflowError("INVALID_TRIMMING_PLAN", "trimming.profiles must be an array.")
    trim_profiles: List[Dict[str, Any]] = []
    seen_profiles: set[str] = set()
    for index, raw_profile in enumerate(raw_profiles):
        if not isinstance(raw_profile, dict):
            raise WorkflowError("INVALID_TRIMMING_PLAN", f"trimming profile {index} must be an object.")
        profile_id = string_at(raw_profile, "id")
        if not SAFE_PROFILE_ID.fullmatch(profile_id) or profile_id in seen_profiles:
            raise WorkflowError(
                "INVALID_TRIMMING_PLAN", f"Invalid or duplicate trimming profile id '{profile_id}'."
            )
        seen_profiles.add(profile_id)
        threshold = float_at(raw_profile, "gap_threshold")
        if threshold > 1:
            raise WorkflowError("INVALID_TRIMMING_PLAN", "gap_threshold must be between 0 and 1.")
        trim_profiles.append({"id": profile_id, "gap_threshold": threshold})
    primary_profile = optional_string(trimming_raw, "primary_profile")
    if trimming_enabled and (not trim_profiles or primary_profile not in seen_profiles):
        raise WorkflowError(
            "INVALID_TRIMMING_PLAN",
            "Enabled trimming requires profiles and a primary_profile matching one profile id.",
        )
    min_retained_fraction = optional_float(trimming_raw, "min_retained_fraction", 0.3)
    if min_retained_fraction > 1:
        raise WorkflowError("INVALID_TRIMMING_PLAN", "min_retained_fraction must be <= 1.")
    normalized_trimming = {
        "enabled": trimming_enabled,
        "tool": optional_string(trimming_raw, "tool", "trimal"),
        "primary_profile": primary_profile,
        "profiles": trim_profiles,
        "min_retained_fraction": min_retained_fraction,
        "compare_topologies": optional_bool(trimming_raw, "compare_topologies", True),
    }
    if normalized_trimming["tool"] != "trimal":
        raise WorkflowError("INVALID_TRIMMING_PLAN", "v0.2 supports trimming.tool='trimal'.")

    if request_schema == "0.1":
        tree_mode = "accurate"
        tree_tool = string_at(tree, "tool")
        support_method = "ultrafast"
        support_replicates = int_at(tree, "ufboot", 1000)
        sh_alrt = int_at(tree, "sh_alrt", 1000)
        bnni = bool_at(tree, "bnni")
    else:
        tree_mode = string_at(tree, "mode")
        tree_tool = string_at(tree, "tool")
        support = mapping_at(tree, "support")
        support_method = string_at(support, "method")
        support_replicates = optional_int(support, "replicates", 1000, 0)
        sh_alrt = optional_int(support, "sh_alrt", 1000, 0)
        bnni = optional_bool(support, "bnni", True)
    if tree_mode not in {"fast", "accurate"}:
        raise WorkflowError("UNSUPPORTED_TREE_PLAN", "tree.mode must be fast or accurate.")
    if tree_mode == "fast":
        if tree_tool.lower() != "fasttree" or support_method != "sh-like-local":
            raise WorkflowError(
                "UNSUPPORTED_TREE_PLAN",
                "Fast mode requires tool=fasttree and support.method=sh-like-local.",
            )
    else:
        if tree_tool != "iqtree2" or support_method not in {"ultrafast", "standard-bootstrap"}:
            raise WorkflowError(
                "UNSUPPORTED_TREE_PLAN",
                "Accurate mode requires iqtree2 with ultrafast or standard-bootstrap support.",
            )
        if support_replicates < 1000 or sh_alrt < 1000:
            raise WorkflowError(
                "UNSUPPORTED_TREE_PLAN", "IQ-TREE2 support and SH-aLRT require at least 1000 replicates."
            )
        if support_method == "standard-bootstrap" and bnni:
            raise WorkflowError(
                "UNSUPPORTED_TREE_PLAN", "-bnni applies to UFBoot; set support.bnni=false for standard bootstrap."
            )
    tree_model = optional_string(tree, "model", "WAG" if tree_mode == "fast" else "MFP").upper()
    if tree_mode == "fast" and tree_model not in {"JTT", "WAG", "LG"}:
        raise WorkflowError(
            "UNSUPPORTED_TREE_PLAN", "FastTree protein model must be JTT, WAG, or LG."
        )
    if tree_mode == "accurate" and not re.fullmatch(r"[A-Za-z0-9+._,/-]+", tree_model):
        raise WorkflowError("UNSUPPORTED_TREE_PLAN", "tree.model contains unsupported characters.")
    rooting = string_at(tree, "rooting")
    if rooting not in {"unrooted", "outgroup"}:
        raise WorkflowError(
            "UNSUPPORTED_ROOTING", "rooting must be unrooted or outgroup; midpoint is not automatic."
        )
    if rooting == "outgroup" and not normalized_selection["require_outgroup"]:
        raise WorkflowError(
            "INVALID_REQUEST", "tree.rooting='outgroup' requires selection.require_outgroup=true."
        )
    normalized_tree = {
        "mode": tree_mode,
        "tool": tree_tool,
        "model": tree_model,
        "support_method": support_method,
        "support_replicates": support_replicates,
        "sh_alrt": sh_alrt,
        "bnni": bnni,
        "threads": optional_int(tree, "threads", 4, 1),
        "seed": optional_int(tree, "seed", 12345, 1),
        "rooting": rooting,
    }

    taxonomy_raw = optional_mapping(request, "taxonomy")
    if "taxonomy" in request:
        allowed_taxonomy_keys = {
            "enabled",
            "source",
            "match_mode",
            "names_dmp",
            "nodes_dmp",
            "snapshot",
            "source_url",
            "retrieved_at",
        }
        unknown_taxonomy_keys = sorted(set(taxonomy_raw) - allowed_taxonomy_keys)
        if unknown_taxonomy_keys:
            raise WorkflowError(
                "INVALID_TAXONOMY_PLAN",
                "Unknown taxonomy fields: " + ", ".join(unknown_taxonomy_keys),
            )
        if "enabled" not in taxonomy_raw:
            raise WorkflowError(
                "INVALID_TAXONOMY_PLAN",
                "A taxonomy object must explicitly set enabled to true or false.",
            )
    taxonomy_enabled = optional_bool(taxonomy_raw, "enabled", False)
    normalized_taxonomy = {
        "enabled": taxonomy_enabled,
        "source": optional_string(taxonomy_raw, "source", "ncbi-taxdump"),
        "match_mode": optional_string(
            taxonomy_raw, "match_mode", "exact-scientific-name"
        ),
        "names_dmp": optional_string(taxonomy_raw, "names_dmp"),
        "nodes_dmp": optional_string(taxonomy_raw, "nodes_dmp"),
        "snapshot": optional_string(taxonomy_raw, "snapshot"),
        "source_url": optional_string(taxonomy_raw, "source_url"),
        "retrieved_at": optional_string(taxonomy_raw, "retrieved_at"),
    }
    if taxonomy_enabled:
        if normalized_taxonomy["source"] != "ncbi-taxdump":
            raise WorkflowError(
                "INVALID_TAXONOMY_PLAN", "taxonomy.source must be 'ncbi-taxdump'."
            )
        if normalized_taxonomy["match_mode"] != "exact-scientific-name":
            raise WorkflowError(
                "INVALID_TAXONOMY_PLAN",
                "taxonomy.match_mode must be 'exact-scientific-name'; fuzzy and alias matching "
                "cannot assign TaxIDs automatically.",
            )
        missing_taxonomy = [
            key
            for key in ("names_dmp", "nodes_dmp", "snapshot", "source_url", "retrieved_at")
            if not normalized_taxonomy[key]
        ]
        if missing_taxonomy:
            raise WorkflowError(
                "INVALID_TAXONOMY_PLAN",
                "Enabled taxonomy validation requires: " + ", ".join(missing_taxonomy),
            )
        if not normalized_taxonomy["source_url"].startswith(
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
        ):
            raise WorkflowError(
                "INVALID_TAXONOMY_PLAN",
                "taxonomy.source_url must identify an official NCBI Taxonomy FTP archive.",
            )

    itol_raw = optional_mapping(request, "itol")
    colors_raw = optional_mapping(itol_raw, "colors")
    colors = {
        role: optional_string(colors_raw, role, default)
        for role, default in ITOL_COLORS.items()
    }
    if any(not re.fullmatch(r"#[0-9A-Fa-f]{6}", color) for color in colors.values()):
        raise WorkflowError("INVALID_ITOL_PLAN", "Every iTOL role color must be a six-digit hex value.")
    if len({color.upper() for color in colors.values()}) != 3:
        raise WorkflowError("INVALID_ITOL_PLAN", "Study, expanded, and outgroup colors must be distinct.")
    normalized_itol = {
        "enabled": optional_bool(itol_raw, "enabled", True),
        "dataset_label": optional_string(itol_raw, "dataset_label", "Sequence roles"),
        "colors": colors,
        "generate_ranges_after_tree_qc": optional_bool(
            itol_raw, "generate_ranges_after_tree_qc", False
        ),
    }
    if any(character in normalized_itol["dataset_label"] for character in "\t\r\n"):
        raise WorkflowError("INVALID_ITOL_PLAN", "iTOL dataset_label must be a single TSV-safe line.")

    literature_raw = optional_mapping(request, "literature")
    normalized_literature = {
        "enabled": optional_bool(literature_raw, "enabled", True),
        "years_back": optional_int(literature_raw, "years_back", 10, 1),
        "include_foundational": optional_bool(literature_raw, "include_foundational", True),
        "taxon_fallback_ranks": string_list(
            literature_raw, "taxon_fallback_ranks", ["species", "genus", "family", "order"]
        ),
        "sources": string_list(
            literature_raw,
            "sources",
            ["primary literature", "Open Tree of Life", "authoritative taxonomy"],
        ),
    }

    return {
        "schema_version": OUTPUT_SCHEMA_VERSION,
        "request_schema_version": request_schema,
        "project_id": project_id,
        "objective": objective,
        "sequence_context": sequence_context,
        "query": {
            "kind": query_kind,
            "path": string_at(query, "path"),
            "id": query_id,
            "original_value": optional_string(query, "original_value", query_id),
            "organism": query_organism,
            "taxon_id": query_taxon_id,
            "gene_name": optional_string(query, "gene_name"),
            "protein_name": optional_string(query, "protein_name"),
            "source_db": optional_string(query, "source_db"),
            "source_release": optional_string(query, "source_release"),
            "retrieved_at": optional_string(query, "retrieved_at"),
        },
        "taxon_scope": {"ingroup": ingroup_scope, "outgroup": outgroup_scope},
        "privacy": {
            "remote_search_allowed": bool_at(privacy, "remote_search_allowed"),
            "unpublished_sequence": optional_bool(privacy, "unpublished_sequence", False),
        },
        "references": {
            "strategy": strategy,
            "candidate_table": string_at(references, "candidate_table"),
            "candidate_fasta": string_at(references, "candidate_fasta"),
            "discovery_tiers": string_list(
                references,
                "discovery_tiers",
                ["curated orthologs", "RefSeq protein", "Swiss-Prot", "UniProtKB/nr", "profile/domain"],
            ),
        },
        "selection": normalized_selection,
        "clustering": normalized_clustering,
        "alignment": normalized_alignment,
        "trimming": normalized_trimming,
        "tree": normalized_tree,
        "taxonomy": normalized_taxonomy,
        "itol": normalized_itol,
        "literature": normalized_literature,
    }


def allowed_relations(config: Mapping[str, Any]) -> set[str]:
    allowed = {"one2one_ortholog", "ortholog", "coortholog"}
    if config["objective"] in {"homolog-context", "within-species"}:
        allowed.add("homolog")
    if config["selection"]["allow_paralogs"]:
        allowed.add("paralog")
    return allowed


def clustering_evaluation(
    candidates: Sequence[Candidate], config: Mapping[str, Any]
) -> Tuple[Dict[str, Any], List[str]]:
    policy = config["clustering"]
    expanded = [candidate for candidate in candidates if candidate.analysis_role == "expanded"]
    triggered = policy["mode"] == "precomputed" or (
        policy["mode"] == "auto" and len(expanded) >= policy["trigger_min_sequences"]
    )
    complete = bool(expanded) and all(candidate.cluster_id for candidate in expanded)
    blockers: List[str] = []
    if triggered and not complete:
        blockers.append("MMSEQS_CLUSTERING_REQUIRED")
    if policy["mode"] == "off":
        status = "disabled"
    elif not triggered:
        status = "not-needed"
    elif complete:
        status = "precomputed"
    else:
        status = "pending-upstream"
    command = {
        "id": "cluster-expanded-candidates",
        "stage": "cluster-candidates",
        "tool": "mmseqs2",
        "argv": [
            "mmseqs",
            policy["algorithm"],
            "expanded_candidates.faa",
            "clusters",
            "mmseqs_tmp",
            "--min-seq-id",
            stable_number(policy["min_seq_id"]),
            "-c",
            stable_number(policy["coverage"]),
            "--cov-mode",
            str(policy["coverage_mode"]),
            "--threads",
            str(policy["threads"]),
        ],
        "inputs": ["expanded_candidates.faa"],
        "outputs": ["clusters_cluster.tsv", "clusters_rep_seq.fasta"],
        "status": "planned" if triggered and not complete else "skipped",
        "executed": False,
    }
    return (
        {
            **policy,
            "expanded_candidate_count": len(expanded),
            "triggered": triggered,
            "cluster_ids_complete": complete,
            "status": status,
            "protected_analysis_groups": ["study", "outgroup"],
            "command": command,
            "replan_required_after_execution": triggered and not complete,
        },
        blockers,
    )


def select_candidates(
    candidates: Sequence[Candidate], query_id: str, query_sequence: str, config: Mapping[str, Any]
) -> Tuple[List[Candidate], List[Tuple[Candidate, List[str]]], List[str]]:
    policy = config["selection"]
    permitted = allowed_relations(config)
    reasons: Dict[str, List[str]] = {}
    eligible: List[Candidate] = []

    query_candidates = [
        candidate
        for candidate in candidates
        if candidate.accession == query_id or candidate.relation == "self"
    ]
    if len(query_candidates) != 1:
        raise WorkflowError(
            "QUERY_CANDIDATE_COUNT",
            "Candidate bundle must contain exactly one self record matching the query ID.",
        )
    query_candidate = query_candidates[0]
    if query_candidate.accession != query_id or query_candidate.relation != "self":
        raise WorkflowError(
            "INVALID_QUERY_CANDIDATE",
            "The query candidate must use the query accession and relation='self'.",
        )
    if query_candidate.role != "ingroup":
        raise WorkflowError("INVALID_QUERY_CANDIDATE", "The query self record must use role='ingroup'.")
    if query_candidate.sequence != query_sequence:
        raise WorkflowError(
            "QUERY_SEQUENCE_MISMATCH",
            "The query self record in the candidate FASTA differs from the query FASTA.",
        )
    query_length = len(query_sequence)

    for candidate in candidates:
        if candidate is query_candidate:
            continue
        current: List[str] = []
        if candidate.is_fragment:
            current.append("FRAGMENT_FLAG")
        if candidate.query_coverage < policy["min_query_coverage"]:
            current.append("LOW_QUERY_COVERAGE")
        if (
            candidate.target_coverage
            and float(candidate.target_coverage) < policy["min_target_coverage"]
        ):
            current.append("LOW_TARGET_COVERAGE")
        ratio = candidate.sequence_length / query_length
        if ratio < policy["min_length_ratio"] or ratio > policy["max_length_ratio"]:
            current.append("LENGTH_RATIO_OUT_OF_RANGE")
        if candidate.relation not in permitted and candidate.relation != "self":
            current.append("RELATION_NOT_ALLOWED")
        if candidate.analysis_role == "study":
            eligible.append(candidate)
            continue
        if current:
            reasons[candidate.accession] = current
        else:
            eligible.append(candidate)

    clustering_plan, _cluster_blockers = clustering_evaluation(candidates, config)
    if clustering_plan["triggered"] and clustering_plan["cluster_ids_complete"]:
        by_cluster: Dict[str, List[Candidate]] = {}
        protected: List[Candidate] = []
        for candidate in eligible:
            if candidate.analysis_role in {"study", "outgroup"}:
                protected.append(candidate)
            else:
                by_cluster.setdefault(candidate.cluster_id, []).append(candidate)
        clustered: List[Candidate] = list(protected)
        for cluster_id in sorted(by_cluster):
            ranked_cluster = sorted(by_cluster[cluster_id], key=Candidate.rank_key)
            clustered.append(ranked_cluster[0])
            for candidate in ranked_cluster[1:]:
                reasons[candidate.accession] = ["MMSEQS_CLUSTER_REDUNDANT"]
        eligible = clustered

    protected_studies = sorted(
        (
            candidate
            for candidate in eligible
            if candidate.analysis_role == "study" and candidate is not query_candidate
        ),
        key=lambda candidate: candidate.accession,
    )
    quota_eligible = [candidate for candidate in eligible if candidate.analysis_role != "study"]
    representatives: List[Candidate] = []
    by_taxon: Dict[str, List[Candidate]] = {}
    for candidate in quota_eligible:
        by_taxon.setdefault(candidate.taxon_id, []).append(candidate)
    for taxon_id in sorted(by_taxon):
        ranked = sorted(by_taxon[taxon_id], key=Candidate.rank_key)
        keep = ranked[: policy["max_per_taxon"]]
        representatives.extend(keep)
        for candidate in ranked[policy["max_per_taxon"] :]:
            reasons[candidate.accession] = ["PER_TAXON_LIMIT"]

    outgroups = sorted(
        (candidate for candidate in representatives if candidate.role == "outgroup"),
        key=Candidate.rank_key,
    )
    ingroups = [
        candidate
        for candidate in representatives
        if candidate.role == "ingroup" and candidate.analysis_role == "expanded"
    ]
    wanted_outgroups = policy["outgroup_count"] if policy["require_outgroup"] else 0
    chosen_outgroups = outgroups[:wanted_outgroups]
    for candidate in outgroups[wanted_outgroups:]:
        reasons[candidate.accession] = ["OUTGROUP_LIMIT"]

    capacity = max(0, policy["max_references"] - len(chosen_outgroups))
    by_clade: Dict[str, List[Candidate]] = {}
    for candidate in ingroups:
        by_clade.setdefault(candidate.clade, []).append(candidate)
    for clade in by_clade:
        by_clade[clade].sort(key=Candidate.rank_key)

    balanced: List[Candidate] = []
    clades = sorted(by_clade)
    while len(balanced) < capacity:
        added = False
        for clade in clades:
            if by_clade[clade] and len(balanced) < capacity:
                balanced.append(by_clade[clade].pop(0))
                added = True
        if not added:
            break
    for remaining in by_clade.values():
        for candidate in remaining:
            reasons[candidate.accession] = ["MAX_REFERENCE_LIMIT"]

    selected = [query_candidate]
    selected.extend(protected_studies)
    selected.extend(sorted(balanced, key=lambda c: (c.clade, c.taxon_id, c.accession)))
    selected.extend(sorted(chosen_outgroups, key=lambda c: (c.clade, c.taxon_id, c.accession)))
    selected_ids = {candidate.accession for candidate in selected}
    rejected = [
        (candidate, reasons[candidate.accession])
        for candidate in sorted(candidates, key=lambda c: c.accession)
        if candidate.accession not in selected_ids
    ]

    blockers: List[str] = []
    ingroup_taxa = {
        candidate.taxon_id
        for candidate in selected
        if candidate.role == "ingroup" and candidate.relation != "self"
    }
    selected_outgroups = [candidate for candidate in selected if candidate.role == "outgroup"]
    if config["objective"] == "within-species":
        ingroup_references = [
            candidate
            for candidate in selected
            if candidate.role == "ingroup" and candidate.relation != "self"
        ]
        if len(ingroup_references) < policy["min_ingroup_taxa"]:
            blockers.append("INSUFFICIENT_INGROUP_SEQUENCES")
    elif len(ingroup_taxa) < policy["min_ingroup_taxa"]:
        blockers.append("INSUFFICIENT_INGROUP_TAXA")
    if policy["require_outgroup"] and len(selected_outgroups) < policy["outgroup_count"]:
        blockers.append("MISSING_REQUIRED_OUTGROUP")
    if config["objective"] == "within-species":
        if len(selected) < 4:
            blockers.append("INSUFFICIENT_TOTAL_SEQUENCES")
    elif len({candidate.taxon_id for candidate in selected}) < 4:
        blockers.append("INSUFFICIENT_TOTAL_TAXA")
    if (
        config["request_schema_version"] == "0.2"
        and policy["require_outgroup"]
        and any(not candidate.outgroup_rationale for candidate in selected_outgroups)
    ):
        blockers.append("OUTGROUP_RATIONALE_REQUIRED")
    return selected, rejected, blockers


def taxonomy_audit_policy(
    config: Mapping[str, Any], resolution: TaxonomyResolution | None
) -> Dict[str, Any]:
    taxonomy = config["taxonomy"]
    if not taxonomy["enabled"]:
        return {
            "enabled": False,
            "status": "not-requested",
            "source": "ncbi-taxdump",
            "match_mode": "exact-scientific-name",
        }
    if resolution is None:
        raise WorkflowError("TAXONOMY_VALIDATION_MISSING", "Enabled taxonomy validation has no evidence.")
    return {
        "enabled": True,
        "status": "validated",
        "source": taxonomy["source"],
        "match_mode": taxonomy["match_mode"],
        "snapshot": resolution.snapshot,
        "source_url": resolution.source_url,
        "retrieved_at": resolution.retrieved_at,
        "names_sha256": resolution.names_sha256,
        "nodes_sha256": resolution.nodes_sha256,
        "records_validated": len(resolution.rows),
        "resolver_version": resolution.resolver_version,
        "resolution_artifact": "taxonomy_resolution.tsv",
    }


def build_fingerprint(
    config: Mapping[str, Any],
    query_sequence: str,
    candidates: Sequence[Candidate],
    taxonomy_resolution: TaxonomyResolution | None,
) -> Dict[str, Any]:
    return {
        "schema_version": config["schema_version"],
        "request_schema_version": config["request_schema_version"],
        "objective": config["objective"],
        "sequence_context": config["sequence_context"],
        "query": {
            "kind": config["query"]["kind"],
            "id": config["query"]["id"],
            "organism": config["query"]["organism"],
            "taxon_id": config["query"]["taxon_id"],
            "source_db": config["query"]["source_db"],
            "sequence_sha256": sha256_text(query_sequence),
        },
        "taxon_scope": config["taxon_scope"],
        "reference_strategy": config["references"]["strategy"],
        "selection": config["selection"],
        "clustering": config["clustering"],
        "alignment": config["alignment"],
        "trimming": config["trimming"],
        "tree": config["tree"],
        "taxonomy": taxonomy_audit_policy(config, taxonomy_resolution),
        "itol": config["itol"],
        "literature": config["literature"],
        "candidates": [
            candidate.fingerprint_record()
            for candidate in sorted(candidates, key=lambda item: item.accession)
        ],
    }


def build_commands(
    config: Mapping[str, Any], clustering_plan: Mapping[str, Any], *, downstream_blocked: bool
) -> List[Dict[str, Any]]:
    alignment = config["alignment"]
    mafft_modes = {
        "auto": ["--auto"],
        "linsi": ["--localpair", "--maxiterate", "1000"],
        "ginsi": ["--globalpair", "--maxiterate", "1000"],
        "einsi": ["--genafpair", "--maxiterate", "1000"],
    }
    mafft_args = [*mafft_modes[alignment["mode"]], "--thread", str(alignment["threads"])]
    mafft_args.append("reference_set.faa")
    commands: List[Dict[str, Any]] = []
    downstream_command_status = "blocked" if downstream_blocked else "planned"
    if clustering_plan["command"]["status"] == "planned":
        commands.append(dict(clustering_plan["command"]))
    commands.append(
        {
            "id": "align-proteins",
            "stage": "align-and-qc",
            "tool": "mafft",
            "argv": ["mafft", *mafft_args],
            "inputs": ["reference_set.faa"],
            "outputs": ["alignment.raw.faa"],
            "stdout": "alignment.raw.faa",
            "status": downstream_command_status,
            "executed": False,
        }
    )
    tree_input = "alignment.raw.faa"
    trimming = config["trimming"]
    if trimming["enabled"]:
        for profile in trimming["profiles"]:
            output_name = f"alignment.trimmed.{profile['id']}.faa"
            commands.append(
                {
                    "id": f"trim-{profile['id']}",
                    "stage": "trim-and-qc",
                    "tool": "trimal",
                    "argv": [
                        "trimal",
                        "-in",
                        "alignment.raw.faa",
                        "-out",
                        output_name,
                        "-gt",
                        stable_number(profile["gap_threshold"]),
                    ],
                    "inputs": ["alignment.raw.faa"],
                    "outputs": [output_name],
                    "stdout": None,
                    "status": downstream_command_status,
                    "executed": False,
                    "threshold_semantics": (
                        "Retain columns whose fraction of non-gap residues meets the configured threshold."
                    ),
                }
            )
        tree_input = f"alignment.trimmed.{trimming['primary_profile']}.faa"

    tree = config["tree"]
    if tree["mode"] == "fast":
        model_args = {"JTT": [], "WAG": ["-wag"], "LG": ["-lg"]}[tree["model"]]
        commands.append(
            {
                "id": "infer-fast-tree",
                "stage": "infer-tree",
                "tool": "fasttree",
                "argv": ["FastTree", *model_args, "-gamma", tree_input],
                "inputs": [tree_input],
                "outputs": ["gene-tree.fast.unrooted.nwk"],
                "stdout": "gene-tree.fast.unrooted.nwk",
                "status": downstream_command_status,
                "executed": False,
                "support_semantics": "SH-like local support; not a global bootstrap.",
            }
        )
    else:
        iqtree_args = ["-s", tree_input, "-m", tree["model"]]
        if tree["support_method"] == "ultrafast":
            iqtree_args.extend(["-B", str(tree["support_replicates"])])
            if tree["bnni"]:
                iqtree_args.append("-bnni")
            support_semantics = "UFBoot2; interpret >=95 separately from standard bootstrap."
        else:
            iqtree_args.extend(["-b", str(tree["support_replicates"])])
            support_semantics = "Standard nonparametric bootstrap; not UFBoot2."
        iqtree_args.extend(
            [
                "-alrt",
                str(tree["sh_alrt"]),
                "-T",
                str(tree["threads"]),
                "-seed",
                str(tree["seed"]),
                "--prefix",
                "gene-tree",
            ]
        )
        commands.append(
            {
                "id": "infer-accurate-tree",
                "stage": "infer-tree",
                "tool": "iqtree2",
                "argv": ["iqtree2", *iqtree_args],
                "inputs": [tree_input],
                "outputs": ["gene-tree.treefile", "gene-tree.contree", "gene-tree.iqtree"],
                "stdout": None,
                "status": downstream_command_status,
                "executed": False,
                "support_semantics": support_semantics,
            }
        )
    return commands


def build_plan(
    run_id: str,
    config: Mapping[str, Any],
    query_sequence: str,
    selected: Sequence[Candidate],
    rejected: Sequence[Tuple[Candidate, Sequence[str]]],
    blockers: Sequence[str],
    input_hashes: Mapping[str, Mapping[str, str]],
    taxonomy_resolution: TaxonomyResolution | None,
) -> Dict[str, Any]:
    all_candidates = list(selected) + [candidate for candidate, _reason_codes in rejected]
    clustering_plan, clustering_blockers = clustering_evaluation(all_candidates, config)
    all_blockers = list(dict.fromkeys([*blockers, *clustering_blockers]))
    downstream_status = "blocked" if all_blockers else "planned"
    outgroup_ids = [candidate.accession for candidate in selected if candidate.role == "outgroup"]
    retained_ingroup_clades = {
        candidate.clade
        for candidate in selected
        if candidate.role == "ingroup" and candidate.relation != "self"
    }
    candidate_ingroup_clades = {
        candidate.clade
        for candidate in selected
        if candidate.role == "ingroup" and candidate.relation != "self"
    }
    candidate_ingroup_clades.update(
        candidate.clade
        for candidate, _reason_codes in rejected
        if candidate.role == "ingroup" and candidate.relation != "self"
    )
    unsampled_clades = sorted(candidate_ingroup_clades - retained_ingroup_clades)
    commands = build_commands(config, clustering_plan, downstream_blocked=bool(all_blockers))
    warnings = [
        "No network request or external executable was run.",
        "Similarity alone does not establish orthology.",
        "A supported gene-tree topology is not automatically a species tree.",
        "Choose a nearby homologous sister-lineage outgroup; never choose the most distant hit automatically.",
    ]
    if config["request_schema_version"] == "0.1":
        warnings.append("Request schema 0.1 was migrated in memory; use schema 0.2 for new projects.")
    if config["sequence_context"] == "viral":
        warnings.append(
            "Viral gene trees require explicit checks for recombination, reassortment, segmentation, and mosaic ancestry."
        )
    if any(not candidate.clade for candidate in all_candidates):
        warnings.append("One or more candidates lack clade metadata; taxonomic balancing may be incomplete.")
    if any(not candidate.lineage for candidate in all_candidates):
        warnings.append("One or more candidates lack a recorded taxonomic lineage.")
    if config["selection"]["min_target_coverage"] > 0 and any(
        not candidate.target_coverage for candidate in all_candidates
    ):
        warnings.append(
            "Target coverage is missing for one or more candidates; the target-coverage threshold could not be applied to them."
        )
    study_qc_accessions = sorted(
        candidate.accession
        for candidate in selected
        if candidate.analysis_role == "study"
        and candidate.relation != "self"
        and (
            candidate.is_fragment
            or candidate.query_coverage < config["selection"]["min_query_coverage"]
            or candidate.sequence_length / len(query_sequence) < config["selection"]["min_length_ratio"]
            or candidate.sequence_length / len(query_sequence) > config["selection"]["max_length_ratio"]
            or candidate.relation not in allowed_relations(config)
        )
    )
    if study_qc_accessions:
        warnings.append(
            "Focal study sequences were preserved despite QC or relationship warnings: "
            + ", ".join(study_qc_accessions)
            + "."
        )
    extreme_trim_profiles = [
        profile["id"]
        for profile in config["trimming"]["profiles"]
        if profile["gap_threshold"] <= 0.1 or profile["gap_threshold"] >= 0.98
    ]
    if extreme_trim_profiles:
        warnings.append(
            "Extreme trimAl profiles are sensitivity analyses, not automatic defaults: "
            + ", ".join(extreme_trim_profiles)
            + "."
        )
    if unsampled_clades:
        warnings.append(
            "Candidate ingroup clades with no retained reference: " + ", ".join(unsampled_clades) + "."
        )
    plan: Dict[str, Any] = {
        "schema_version": OUTPUT_SCHEMA_VERSION,
        "request_schema_version": config["request_schema_version"],
        "project_id": config["project_id"],
        "run_id": run_id,
        "state": (
            "pending-clustering"
            if "MMSEQS_CLUSTERING_REQUIRED" in all_blockers
            else ("blocked" if all_blockers else "pending-reference-approval")
        ),
        "created_at": None,
        "offline": True,
        "mode": "offline-dry-run",
        "objective": config["objective"],
        "sequence_context": config["sequence_context"],
        "query": {
            "id": config["query"]["id"],
            "original_kind": config["query"]["kind"],
            "original_value": config["query"]["original_value"],
            "organism": config["query"]["organism"],
            "taxon_id": config["query"]["taxon_id"],
            "length": len(query_sequence),
            "sequence_sha256": sha256_text(query_sequence),
        },
        "query_resolution": {
            "status": "materialized-local-record",
            "source_db": config["query"]["source_db"],
            "source_release": config["query"]["source_release"],
            "retrieved_at": config["query"]["retrieved_at"],
            "gene_name": config["query"]["gene_name"],
            "protein_name": config["query"]["protein_name"],
            "live_lookup_performed_by_planner": False,
        },
        "taxon_scope": config["taxon_scope"],
        "taxonomy_plan": taxonomy_audit_policy(config, taxonomy_resolution),
        "privacy": config["privacy"],
        "selection_parameters": config["selection"],
        "reference_discovery": {
            "strategy": config["references"]["strategy"],
            "fallback_tiers": config["references"]["discovery_tiers"],
            "status": "provided-local-candidate-bundle",
        },
        "clustering_plan": clustering_plan,
        "alignment_plan": {
            **config["alignment"],
            "raw_alignment": "alignment.raw.faa",
            "qc_required_before_trimming": True,
        },
        "trimming_plan": {
            **config["trimming"],
            "preserve_raw_alignment": True,
            "approval_required_before_tree": True,
        },
        "tree_plan": {
            **config["tree"],
            "tree_claim": "gene-tree",
            "support_is_repeatability_not_correctness": True,
            "preserve_unrooted": True,
        },
        "annotation_plan": {
            "itol_colorstrip": "itol_roles.txt" if config["itol"]["enabled"] else None,
            "local_ggtree_renderer": {
                "script": "scripts/render_tree_ggtree.R",
                "engine": "ggtree+ggplot2",
                "status": "post-tree-local",
                "joins": "exact tip_id equality",
                "default_outputs": ["SVG", "PDF", "settings TSV"],
                "reroots_tree": False,
            },
            "colors": config["itol"]["colors"],
            "range_dataset": (
                "post-tree-only after contiguity/monophyly review"
                if config["itol"]["generate_ranges_after_tree_qc"]
                else "disabled"
            ),
            "sequence_metadata": "sequence_metadata.tsv",
        },
        "literature_plan": {
            **config["literature"],
            "status": "pending-live-evidence-search" if config["literature"]["enabled"] else "disabled",
            "exact_taxon_then_escalate": True,
            "gene_tree_species_tree_comparison": "qualitative; discordance is not automatically an error",
        },
        "input_hashes": input_hashes,
        "scientific_scope": {
            "objective": config["objective"],
            "claim": "gene-tree",
            "not_a_species_tree": True,
        },
        "selection_summary": {
            "query_id": config["query"]["id"],
            "query_length": len(query_sequence),
            "selected_accessions": [candidate.accession for candidate in selected],
            "selected_ingroup_taxa": len(
                {
                    candidate.taxon_id
                    for candidate in selected
                    if candidate.role == "ingroup" and candidate.relation != "self"
                }
            ),
            "selected_outgroups": outgroup_ids,
            "selected_study_sequences": [
                candidate.accession for candidate in selected if candidate.analysis_role == "study"
            ],
            "unsampled_candidate_ingroup_clades": unsampled_clades,
            "rejected_count": len(rejected),
            "blockers": list(all_blockers),
        },
        "selected_accessions": [candidate.accession for candidate in selected],
        "rejected_accessions_and_reasons": [
            {"accession": candidate.accession, "reason_codes": list(reason_codes)}
            for candidate, reason_codes in rejected
        ],
        "outgroup_accessions": outgroup_ids,
        "candidate_and_selection_counts": {
            "candidate_records_including_query": len(selected) + len(rejected),
            "selected_records_including_query": len(selected),
            "selected_references_excluding_query": max(0, len(selected) - 1),
            "selected_references_excluding_study": len(
                [candidate for candidate in selected if candidate.analysis_role != "study"]
            ),
            "rejected_records": len(rejected),
        },
        "hard_stops": list(all_blockers),
        "approval": None,
        "decision_gates": [
            {
                "id": "reference-and-outgroup",
                "status": "blocked" if all_blockers else "pending",
                "requires": [
                    "selected references",
                    "paralog policy",
                    "outgroup rationale",
                    "taxonomy evidence when enabled",
                ],
            },
            {
                "id": "alignment-and-trimming-qc",
                "status": "not-reached",
                "requires": ["raw MSA QC", "trim profile comparison", "primary alignment choice"],
            },
        ],
        "stages": [
            {"id": "intake", "status": "completed"},
            {"id": "resolve-query", "status": "completed", "method": "local protein FASTA"},
            {
                "id": "discover-candidates",
                "status": "completed",
                "method": "provided local candidate bundle; no database query",
            },
            {
                "id": "cluster-candidates",
                "status": clustering_plan["status"],
                "protected": ["study", "outgroup"],
            },
            {"id": "select-references", "status": "blocked" if all_blockers else "completed"},
            {"id": "align-and-qc", "status": downstream_status},
            {"id": "trim-and-qc", "status": downstream_status},
            {"id": "infer-tree", "status": downstream_status, "tree_state": "unrooted"},
            {"id": "root-and-annotate", "status": downstream_status},
            {
                "id": "literature-context",
                "status": "pending" if config["literature"]["enabled"] else "skipped",
            },
        ],
        "planned_commands": commands,
        "rooting_plan": {
            "method": config["tree"]["rooting"],
            "outgroup_accessions": outgroup_ids,
            "preserve_unrooted": True,
            "executed": False,
        },
        "warnings": warnings,
        "candidate_semantic_hash": sha256_text(
            canonical_json(
                [
                    candidate.fingerprint_record()
                    for candidate in sorted(all_candidates, key=lambda item: item.accession)
                ]
            )
        ),
    }
    plan_hash_source = dict(plan)
    plan_hash_source.pop("created_at")
    plan_hash_source.pop("approval")
    plan_hash_source.pop("input_hashes")
    plan["plan_hash"] = sha256_text(canonical_json(plan_hash_source))
    return plan


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")


def write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def reference_fasta_text(
    selected: Sequence[Candidate],
) -> str:
    parts: List[str] = []
    for candidate in selected:
        parts.append(f">{candidate.accession}\n{wrap_sequence(candidate.sequence)}")
    return "\n".join(parts) + "\n"


def sequence_metadata_rows(
    selected: Sequence[Candidate], rejected: Sequence[Tuple[Candidate, Sequence[str]]]
) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for index, candidate in enumerate(selected, start=1):
        row = candidate.output_row()
        row.update(
            {
                "tip_id": candidate.accession,
                "analysis_role": candidate.analysis_role,
                "inclusion_status": "selected",
                "selection_order": str(index),
                "reason_codes": "",
            }
        )
        rows.append(row)
    for candidate, reason_codes in sorted(rejected, key=lambda item: item[0].accession):
        row = candidate.output_row()
        row.update(
            {
                "tip_id": candidate.accession,
                "analysis_role": candidate.analysis_role,
                "inclusion_status": "rejected",
                "selection_order": "",
                "reason_codes": ";".join(reason_codes),
            }
        )
        rows.append(row)
    return rows


def itol_colorstrip_text(selected: Sequence[Candidate], config: Mapping[str, Any]) -> str:
    colors = config["itol"]["colors"]
    label = config["itol"]["dataset_label"]
    lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR TAB",
        f"DATASET_LABEL\t{label}",
        "COLOR\t#000000",
        "STRIP_WIDTH\t25",
        "MARGIN\t5",
        "BORDER_WIDTH\t0",
        "SHOW_INTERNAL\t0",
        "LEGEND_TITLE\tSequence role",
        "LEGEND_SHAPES\t1\t1\t1",
        f"LEGEND_COLORS\t{colors['study']}\t{colors['expanded']}\t{colors['outgroup']}",
        "LEGEND_LABELS\tStudy\tExpanded\tOutgroup",
        "DATA",
    ]
    display = {"study": "Study", "expanded": "Expanded", "outgroup": "Outgroup"}
    for candidate in selected:
        role = candidate.analysis_role
        lines.append(f"{candidate.accession}\t{colors[role]}\t{display[role]}")
    return "\n".join(lines) + "\n"


def logical_input(role: str, filename: str, digest: str) -> Dict[str, str]:
    return {"role": role, "logical_path": f"inputs/{filename}", "sha256": digest}


def write_bundle(
    output_path: Path,
    request_path: Path,
    query_path: Path,
    candidate_fasta_path: Path,
    candidate_table_path: Path,
    config: Mapping[str, Any],
    query_sequence: str,
    selected: Sequence[Candidate],
    rejected: Sequence[Tuple[Candidate, Sequence[str]]],
    blockers: Sequence[str],
    run_id: str,
    taxonomy_resolution: TaxonomyResolution | None = None,
    taxonomy_names_path: Path | None = None,
    taxonomy_nodes_path: Path | None = None,
) -> None:
    if output_path.exists():
        raise WorkflowError("OUTPUT_EXISTS", f"Refusing to overwrite existing output path: {output_path}")
    parent = output_path.parent
    try:
        parent.mkdir(parents=True, exist_ok=True)
        temp_path = Path(tempfile.mkdtemp(prefix=f".{output_path.name}.tmp-", dir=str(parent)))
    except OSError as exc:
        raise WorkflowError("OUTPUT_CREATE_ERROR", f"Cannot prepare output directory: {exc}") from exc

    try:
        selected_rows: List[MutableMapping[str, str]] = []
        for index, candidate in enumerate(selected, start=1):
            row = candidate.output_row()
            row["selection_order"] = str(index)
            if candidate.relation == "self":
                row["decision_reason"] = "QUERY"
            elif candidate.analysis_role == "study":
                row["decision_reason"] = "SELECTED_STUDY"
            elif candidate.role == "outgroup":
                row["decision_reason"] = "SELECTED_OUTGROUP"
            else:
                row["decision_reason"] = "SELECTED_INGROUP"
            selected_rows.append(row)
        rejected_rows: List[MutableMapping[str, str]] = []
        for candidate, reason_codes in rejected:
            row = candidate.output_row()
            row["reason_codes"] = ";".join(reason_codes)
            rejected_rows.append(row)

        write_tsv(
            temp_path / "selected_references.tsv",
            BASE_OUTPUT_COLUMNS + ("selection_order", "decision_reason"),
            selected_rows,
        )
        write_tsv(
            temp_path / "rejected_references.tsv",
            BASE_OUTPUT_COLUMNS + ("reason_codes",),
            rejected_rows,
        )
        (temp_path / "reference_set.faa").write_text(
            reference_fasta_text(selected),
            encoding="utf-8",
        )
        metadata_columns = (
            "tip_id",
            "analysis_role",
            "inclusion_status",
            "selection_order",
            "reason_codes",
        ) + BASE_OUTPUT_COLUMNS
        write_tsv(
            temp_path / "sequence_metadata.tsv",
            metadata_columns,
            sequence_metadata_rows(selected, rejected),
        )
        if config["itol"]["enabled"]:
            (temp_path / "itol_roles.txt").write_text(
                itol_colorstrip_text(selected, config), encoding="utf-8"
            )
        if taxonomy_resolution is not None:
            write_resolution_tsv(temp_path / "taxonomy_resolution.tsv", taxonomy_resolution)
        input_hashes = {
            "request": {
                "logical_path": f"inputs/{request_path.name}",
                "sha256": file_sha256(request_path),
            },
            "query_fasta": {
                "logical_path": f"inputs/{query_path.name}",
                "sha256": file_sha256(query_path),
            },
            "candidate_fasta": {
                "logical_path": f"inputs/{candidate_fasta_path.name}",
                "sha256": file_sha256(candidate_fasta_path),
            },
            "candidate_table": {
                "logical_path": f"inputs/{candidate_table_path.name}",
                "sha256": file_sha256(candidate_table_path),
            },
        }
        if taxonomy_resolution is not None:
            if taxonomy_names_path is None or taxonomy_nodes_path is None:
                raise WorkflowError(
                    "TAXONOMY_VALIDATION_MISSING",
                    "Taxonomy evidence is missing its names.dmp or nodes.dmp input path.",
                )
            input_hashes["taxonomy_names"] = {
                "logical_path": "inputs/names.dmp",
                "sha256": taxonomy_resolution.names_sha256,
            }
            input_hashes["taxonomy_nodes"] = {
                "logical_path": "inputs/nodes.dmp",
                "sha256": taxonomy_resolution.nodes_sha256,
            }
        plan = build_plan(
            run_id,
            config,
            query_sequence,
            selected,
            rejected,
            blockers,
            input_hashes,
            taxonomy_resolution,
        )
        if plan["clustering_plan"]["command"]["status"] == "planned":
            expanded = [
                candidate
                for candidate in list(selected) + [item[0] for item in rejected]
                if candidate.analysis_role == "expanded"
            ]
            (temp_path / "expanded_candidates.faa").write_text(
                reference_fasta_text(sorted(expanded, key=lambda item: item.accession)),
                encoding="utf-8",
            )
        write_json(temp_path / "plan.json", plan)

        output_names: Tuple[str, ...] = (
            "selected_references.tsv",
            "rejected_references.tsv",
            "reference_set.faa",
            "sequence_metadata.tsv",
            "plan.json",
        )
        if config["itol"]["enabled"]:
            output_names += ("itol_roles.txt",)
        if taxonomy_resolution is not None:
            output_names += ("taxonomy_resolution.tsv",)
        if (temp_path / "expanded_candidates.faa").is_file():
            output_names += ("expanded_candidates.faa",)
        media_types = {
            "selected_references.tsv": "text/tab-separated-values",
            "rejected_references.tsv": "text/tab-separated-values",
            "reference_set.faa": "text/x-fasta",
            "sequence_metadata.tsv": "text/tab-separated-values",
            "itol_roles.txt": "text/plain",
            "taxonomy_resolution.tsv": "text/tab-separated-values",
            "expanded_candidates.faa": "text/x-fasta",
            "plan.json": "application/json",
        }
        request_digest = input_hashes["request"]["sha256"]
        database_provenance = [
            {"source_db": source, "source_release": release, "retrieved_at": retrieved}
            for source, release, retrieved in sorted(
                {
                    (candidate.source_db, candidate.source_release, candidate.retrieved_at)
                    for candidate in [item for item in selected]
                    + [item[0] for item in rejected]
                }
            )
        ]
        if taxonomy_resolution is not None:
            database_provenance.append(
                {
                    "source_db": "NCBI Taxonomy new_taxdump",
                    "source_release": taxonomy_resolution.snapshot,
                    "retrieved_at": taxonomy_resolution.retrieved_at,
                    "source_url": taxonomy_resolution.source_url,
                    "names_sha256": taxonomy_resolution.names_sha256,
                    "nodes_sha256": taxonomy_resolution.nodes_sha256,
                    "resolver_version": taxonomy_resolution.resolver_version,
                }
            )
        manifest_inputs = [
            logical_input("request", request_path.name, request_digest),
            logical_input("query_fasta", query_path.name, input_hashes["query_fasta"]["sha256"]),
            logical_input(
                "candidate_fasta",
                candidate_fasta_path.name,
                input_hashes["candidate_fasta"]["sha256"],
            ),
            logical_input(
                "candidate_table",
                candidate_table_path.name,
                input_hashes["candidate_table"]["sha256"],
            ),
        ]
        if taxonomy_resolution is not None:
            manifest_inputs.extend(
                [
                    logical_input("taxonomy_names", "names.dmp", taxonomy_resolution.names_sha256),
                    logical_input("taxonomy_nodes", "nodes.dmp", taxonomy_resolution.nodes_sha256),
                ]
            )
        manifest: Dict[str, Any] = {
            "schema_version": OUTPUT_SCHEMA_VERSION,
            "workflow_version": VERSION,
            "run_id": run_id,
            "workflow_state": plan["state"],
            "created_at": None,
            "offline": True,
            "mode": "offline-dry-run",
            "project_id": config["project_id"],
            "request_path_and_hash": {
                "logical_path": f"inputs/{request_path.name}",
                "sha256": request_digest,
            },
            "input_artifacts": manifest_inputs,
            "query": {
                "id": config["query"]["id"],
                "organism": config["query"]["organism"],
                "taxon_id": config["query"]["taxon_id"],
                "length": len(query_sequence),
                "sequence_sha256": sha256_text(query_sequence),
            },
            "database_provenance": database_provenance,
            "policy": {
                "objective": config["objective"],
                "sequence_context": config["sequence_context"],
                "reference_strategy": config["references"]["strategy"],
                "taxon_scope": config["taxon_scope"],
                "selection": config["selection"],
                "clustering": config["clustering"],
                "alignment": config["alignment"],
                "trimming": config["trimming"],
                "tree": config["tree"],
                "taxonomy": taxonomy_audit_policy(config, taxonomy_resolution),
                "itol": config["itol"],
                "literature": config["literature"],
            },
            "decisions": {
                "selected_accessions": [candidate.accession for candidate in selected],
                "rejected": [
                    {"accession": candidate.accession, "reason_codes": list(reason_codes)}
                    for candidate, reason_codes in rejected
                ],
                "blockers": list(plan["hard_stops"]),
            },
            "plan_hash": plan["plan_hash"],
            "approved_plan_hash": None,
            "tool_versions": {
                "mmseqs2": {"status": "not-inspected", "version": None},
                "mafft": {"status": "not-inspected", "version": None},
                "trimal": {"status": "not-inspected", "version": None},
                "fasttree": {"status": "not-inspected", "version": None},
                "iqtree2": {"status": "not-inspected", "version": None},
                "rscript": {"status": "not-inspected", "version": None},
            },
            "commands": plan["planned_commands"],
            "execution": {"network_calls": 0, "external_processes": 0},
            "warnings": list(plan["warnings"]),
            "errors": list(plan["hard_stops"]),
            "output_artifacts": [
                {
                    "logical_path": name,
                    "media_type": media_types[name],
                    "bytes": (temp_path / name).stat().st_size,
                    "sha256": file_sha256(temp_path / name),
                }
                for name in output_names
            ],
        }
        write_json(temp_path / "manifest.json", manifest)

        if output_path.exists():
            raise WorkflowError("OUTPUT_EXISTS", f"Output path appeared during the run: {output_path}")
        os.replace(temp_path, output_path)
    except Exception:
        shutil.rmtree(temp_path, ignore_errors=True)
        raise


def run_plan(args: argparse.Namespace) -> int:
    if not args.offline or not args.dry_run:
        raise WorkflowError(
            "LIVE_MODE_NOT_IMPLEMENTED",
            "Planning requires both --offline and --dry-run; it never contacts databases or runs tools.",
        )
    request_path = Path(args.request).expanduser().resolve()
    output_path = Path(args.out).expanduser().resolve()
    if output_path.exists():
        raise WorkflowError("OUTPUT_EXISTS", f"Refusing to overwrite existing output path: {output_path}")
    if not request_path.is_file():
        raise WorkflowError("MISSING_REQUEST", f"Request file does not exist: {args.request}")

    request = read_json(request_path)
    config = normalize_request(request)
    query_path = resolve_input_path(request_path, config["query"]["path"], "Query FASTA")
    candidate_fasta_path = resolve_input_path(
        request_path, config["references"]["candidate_fasta"], "Candidate FASTA"
    )
    candidate_table_path = resolve_input_path(
        request_path, config["references"]["candidate_table"], "Candidate table"
    )

    query_records = parse_fasta(query_path, "query")
    if len(query_records) != 1:
        raise WorkflowError("QUERY_RECORD_COUNT", "Query FASTA must contain exactly one protein record.")
    query_id, _query_header, query_sequence = query_records[0]
    if query_id != config["query"]["id"]:
        raise WorkflowError(
            "QUERY_ID_MISMATCH",
            f"Configured query ID '{config['query']['id']}' does not match FASTA ID '{query_id}'.",
        )

    candidate_records = parse_fasta(candidate_fasta_path, "candidate")
    candidate_sequences = {identifier: sequence for identifier, _header, sequence in candidate_records}
    candidates = parse_candidate_table(candidate_table_path, candidate_sequences)
    taxonomy_resolution: TaxonomyResolution | None = None
    taxonomy_names_path: Path | None = None
    taxonomy_nodes_path: Path | None = None
    if config["taxonomy"]["enabled"]:
        taxonomy_names_path = resolve_input_path(
            request_path, config["taxonomy"]["names_dmp"], "NCBI names.dmp"
        )
        taxonomy_nodes_path = resolve_input_path(
            request_path, config["taxonomy"]["nodes_dmp"], "NCBI nodes.dmp"
        )
        query_self = [
            candidate
            for candidate in candidates
            if candidate.accession == config["query"]["id"] or candidate.relation == "self"
        ]
        if len(query_self) != 1:
            raise WorkflowError(
                "QUERY_CANDIDATE_COUNT",
                "Candidate bundle must contain exactly one self record matching the query ID.",
            )
        if (
            query_self[0].accession != config["query"]["id"]
            or query_self[0].relation != "self"
        ):
            raise WorkflowError(
                "INVALID_QUERY_CANDIDATE",
                "The query candidate must use the query accession and relation='self'.",
            )
        if query_self[0].species != config["query"]["organism"]:
            raise WorkflowError(
                "QUERY_TAXONOMY_MISMATCH",
                "query.organism must exactly equal the query self-record species before "
                "NCBI Taxonomy validation.",
            )
        if config["query"]["taxon_id"] and query_self[0].taxon_id != config["query"]["taxon_id"]:
            raise WorkflowError(
                "QUERY_TAXONOMY_MISMATCH",
                "query.taxon_id must equal the query self-record TaxID before NCBI Taxonomy "
                "validation.",
            )
        try:
            taxonomy_resolution = resolve_exact_scientific_names(
                [
                    NameRequest(candidate.accession, candidate.species, candidate.taxon_id)
                    for candidate in candidates
                ],
                names_path=taxonomy_names_path,
                nodes_path=taxonomy_nodes_path,
                snapshot=config["taxonomy"]["snapshot"],
                source_url=config["taxonomy"]["source_url"],
                retrieved_at=config["taxonomy"]["retrieved_at"],
            )
        except TaxonomyError as exc:
            detail = (
                " Details: " + json.dumps(exc.details, ensure_ascii=False, sort_keys=True)
                if exc.details
                else ""
            )
            raise WorkflowError(exc.code, exc.message + detail) from exc
        query_resolution_row = next(
            row
            for row in taxonomy_resolution.rows
            if row["record_id"] == config["query"]["id"]
        )
        config["query"]["taxon_id"] = query_resolution_row["taxon_id"]
    selected, rejected, blockers = select_candidates(candidates, query_id, query_sequence, config)
    _clustering_plan, clustering_blockers = clustering_evaluation(candidates, config)
    all_blockers = list(dict.fromkeys([*blockers, *clustering_blockers]))
    fingerprint = build_fingerprint(config, query_sequence, candidates, taxonomy_resolution)
    run_id = "gtr-" + sha256_text(canonical_json(fingerprint))[:16]
    write_bundle(
        output_path,
        request_path,
        query_path,
        candidate_fasta_path,
        candidate_table_path,
        config,
        query_sequence,
        selected,
        rejected,
        all_blockers,
        run_id,
        taxonomy_resolution,
        taxonomy_names_path,
        taxonomy_nodes_path,
    )

    status = (
        "pending-clustering"
        if "MMSEQS_CLUSTERING_REQUIRED" in all_blockers
        else ("blocked" if all_blockers else "pending-reference-approval")
    )
    summary = {
        "run_id": run_id,
        "status": status,
        "output": str(output_path),
        "selected": len(selected),
        "rejected": len(rejected),
        "blockers": all_blockers,
        "network_calls": 0,
        "external_processes": 0,
    }
    print(json.dumps(summary, sort_keys=True))
    return 3 if all_blockers else 0


def run_doctor(args: argparse.Namespace) -> int:
    probes = {
        "mmseqs2": ["mmseqs", "version"],
        "mafft": ["mafft", "--version"],
        "trimal": ["trimal", "--version"],
        "fasttree": ["FastTree", "-help"],
        "iqtree2": ["iqtree2", "--version"],
        "rscript": ["Rscript", "--version"],
    }
    results: Dict[str, Dict[str, Any]] = {}
    for name, argv in probes.items():
        executable = shutil.which(argv[0])
        if executable is None:
            results[name] = {"status": "missing", "executable": None, "version": None}
            continue
        version: str | None = None
        status = "available"
        try:
            completed = subprocess.run(
                [executable, *argv[1:]],
                capture_output=True,
                text=True,
                timeout=5,
                check=False,
            )
            rendered = (completed.stdout or completed.stderr).strip().splitlines()
            version = rendered[0][:300] if rendered else None
            if completed.returncode not in {0, 1}:
                status = "probe-failed"
        except (OSError, subprocess.SubprocessError) as exc:
            status = "probe-failed"
            version = str(exc)
        results[name] = {
            "status": status,
            "executable": Path(executable).name,
            "version": version,
        }
    payload = {"workflow_version": VERSION, "tools": results}
    print(json.dumps(payload, indent=None if args.json else 2, sort_keys=True))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gene_to_tree.py",
        description="Compile an auditable gene-to-reference-tree review bundle or inspect local tools.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    subparsers = parser.add_subparsers(dest="command", required=True)
    plan = subparsers.add_parser("plan", help="Validate a local bundle and create a review plan.")
    plan.add_argument("--request", required=True, help="Path to request JSON.")
    plan.add_argument("--out", required=True, help="New output directory; existing paths are refused.")
    plan.add_argument("--offline", action="store_true", help="Assert that all inputs are local.")
    plan.add_argument("--dry-run", action="store_true", help="Plan only; launch no external tools.")
    plan.set_defaults(handler=run_plan)
    doctor = subparsers.add_parser(
        "doctor", help="Inspect whether optional local bioinformatics executables are available."
    )
    doctor.add_argument("--json", action="store_true", help="Emit compact JSON.")
    doctor.set_defaults(handler=run_doctor)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.handler(args))
    except WorkflowError as exc:
        print(f"ERROR [{exc.code}] {exc.message}", file=sys.stderr)
        return exc.exit_code
    except OSError as exc:
        print(f"ERROR [IO_ERROR] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
