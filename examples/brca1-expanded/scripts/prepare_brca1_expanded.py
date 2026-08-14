#!/usr/bin/env python3
"""Materialize a manifest-selected BRCA1 reference-review bundle.

This helper does not infer orthology. It validates a declared, taxonomically
balanced selection against one frozen NCBI Ortholog package, records the full
provider species/protein inventory, and writes the local inputs required for
BLAST/domain QC and the deterministic gene-tree planner.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


RETRIEVED_AT = "2026-08-13"
SOURCE_DB = "NCBI Gene / RefSeq"
ORTHOLOGY_SOURCE = "NCBI Orthologs"
ORTHOLOGY_EVIDENCE = (
    "NCBI Ortholog group membership for human GeneID 672; the provider method "
    "combines protein similarity, nucleotide alignment, and microsynteny; "
    "pair-specific evidence components were not exported"
)
DOMAIN_ARCHITECTURE = "pending Pfam+SMART terminal-domain QC"
SHA256 = re.compile(r"[0-9a-f]{64}")

MANIFEST_FIELDS = (
    "accession", "taxon_id", "species", "gene_id", "analysis_group",
    "clade", "is_reviewed", "selection_rationale",
)
CANDIDATE_FIELDS = (
    "accession", "taxon_id", "species", "role", "relation", "is_reviewed",
    "is_canonical", "is_fragment", "query_coverage", "sequence_length",
    "bitscore", "evalue", "source_db", "source_release", "retrieved_at",
    "clade", "accession_version", "gene_name", "protein_name", "lineage",
    "analysis_group", "target_coverage", "percent_identity",
    "alignment_length", "orthology_source", "orthology_evidence",
    "domain_architecture", "cluster_id", "cluster_representative",
    "outgroup_rationale", "retrieval_query_id", "notes",
)
INVENTORY_FIELDS = (
    "taxon_id", "species", "gene_id", "provider_protein_count",
    "minimum_sequence_length", "maximum_sequence_length",
    "within_declared_length_gate_count", "selected_accession",
    "rejected_accession", "selection_status", "reason_code", "evidence_note",
)
REJECTION_FIELDS = (
    "accession", "taxon_id", "species", "gene_id", "stage", "decision",
    "reason", "replacement_accession", "evidence",
)


class PreparationError(RuntimeError):
    pass


@dataclass(frozen=True)
class ProteinRecord:
    accession: str
    header: str
    sequence: str
    species: str
    gene_id: str
    isoform: str


def fail(message: str) -> None:
    raise PreparationError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != MANIFEST_FIELDS:
            fail("selection manifest header does not match the exact contract")
        rows = list(reader)
    if not rows:
        fail("selection manifest is empty")
    for line_number, row in enumerate(rows, start=2):
        if None in row or any(value != value.strip() or not value for value in row.values()):
            fail(f"selection manifest row {line_number} has empty, extra, or padded fields")
        if row["analysis_group"] not in {"study", "expanded", "outgroup"}:
            fail(f"row {line_number}: unsupported analysis_group")
        if row["is_reviewed"] not in {"true", "false"}:
            fail(f"row {line_number}: is_reviewed must be true or false")
        if not re.fullmatch(r"[A-Z]{2}_[0-9]+\.[0-9]+", row["accession"]):
            fail(f"row {line_number}: accession must be versioned NP_/XP_")
        expected_reviewed = str(row["accession"].startswith("NP_")).lower()
        if row["is_reviewed"] != expected_reviewed:
            fail(
                f"row {line_number}: is_reviewed must follow the conservative "
                "RefSeq NP_=true / XP_=false accession-class rule"
            )
        if not re.fullmatch(r"[0-9]+", row["taxon_id"]):
            fail(f"row {line_number}: invalid TaxID")
    for field in ("accession", "taxon_id", "species", "gene_id"):
        values = [row[field] for row in rows]
        if len(values) != len(set(values)):
            fail(f"selection manifest contains duplicate {field}")
    study = [row for row in rows if row["analysis_group"] == "study"]
    if len(study) != 1 or study[0]["accession"] != "NP_009225.1":
        fail("manifest must contain exactly one NP_009225.1 study record")
    if sum(row["analysis_group"] == "outgroup" for row in rows) < 2:
        fail("manifest must retain at least two outgroup candidates")
    return rows


def read_rejections(path: Path | None) -> list[dict[str, str]]:
    if path is None:
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != REJECTION_FIELDS:
            fail("rejected-candidate ledger header does not match the exact contract")
        rows = list(reader)
    for line_number, row in enumerate(rows, start=2):
        if None in row or any(value != value.strip() or not value for value in row.values()):
            fail(f"rejected-candidate row {line_number} has empty, extra, or padded fields")
        if row["decision"] != "rejected":
            fail(f"rejected-candidate row {line_number} must have decision=rejected")
        if not re.fullmatch(r"[A-Z]{2}_[0-9]+\.[0-9]+", row["accession"]):
            fail(f"rejected-candidate row {line_number} has an invalid accession")
    for field in ("accession", "taxon_id", "gene_id"):
        values = [row[field] for row in rows]
        if len(values) != len(set(values)):
            fail(f"rejected-candidate ledger contains duplicate {field}")
    return rows


def parse_fasta(path: Path) -> dict[str, ProteinRecord]:
    records: dict[str, ProteinRecord] = {}
    header = ""
    sequence: list[str] = []

    def store() -> None:
        if not header:
            return
        accession = header.split()[0]
        organism = re.search(r"\[organism=([^]]+)\]", header)
        gene = re.search(r"\[GeneID=([^]]+)\]", header)
        isoform = re.search(r"\[isoform=([^]]+)\]", header)
        joined = "".join(sequence).upper()
        if accession in records:
            fail(f"duplicate provider FASTA accession: {accession}")
        if not organism or not gene or not joined or not re.fullmatch(r"[A-Z]+", joined):
            fail(f"invalid provider FASTA record: {accession}")
        records[accession] = ProteinRecord(
            accession=accession,
            header=header,
            sequence=joined,
            species=organism.group(1),
            gene_id=gene.group(1),
            isoform=isoform.group(1) if isoform else "",
        )

    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                store()
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
    store()
    if not records:
        fail("provider protein FASTA is empty")
    return records


def load_reports(path: Path) -> dict[str, dict]:
    reports: dict[str, dict] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            record = json.loads(raw)
            gene_id = str(record.get("geneId") or "")
            if not gene_id or gene_id in reports:
                fail(f"missing or duplicate GeneID in data-report row {line_number}")
            reports[gene_id] = record
    if not reports:
        fail("provider data report is empty")
    return reports


def annotation_context(report: dict) -> str:
    annotations = report.get("annotations") or []
    if not annotations:
        return "NCBI gene record; no assembly annotation context exported"
    annotation = annotations[0]
    label = str(annotation.get("annotationName") or annotation.get("assemblyAccession") or "")
    return f"NCBI gene-record annotation context: {label}"


def write_fasta(path: Path, records: list[ProteinRecord]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(f">{record.accession} BRCA1 [organism={record.species}]\n")
            for index in range(0, len(record.sequence), 70):
                handle.write(record.sequence[index:index + 70] + "\n")


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein-fasta", type=Path, required=True)
    parser.add_argument("--data-report", type=Path, required=True)
    parser.add_argument("--selection-manifest", type=Path, required=True)
    parser.add_argument("--rejected-candidates", type=Path)
    parser.add_argument(
        "--externally-verified-provider-archive-sha256",
        dest="provider_archive_sha256",
        required=True,
        help=(
            "lowercase SHA-256 independently verified by the host against the "
            "frozen provider archive; this helper does not receive the archive itself"
        ),
    )
    parser.add_argument("--min-length-ratio", type=float, default=0.75)
    parser.add_argument("--max-length-ratio", type=float, default=1.10)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    if args.out.exists():
        fail(f"refusing to overwrite output directory: {args.out}")
    if not SHA256.fullmatch(args.provider_archive_sha256):
        fail("provider archive SHA-256 must be a lowercase digest")
    if not 0 < args.min_length_ratio <= args.max_length_ratio:
        fail("invalid length-ratio interval")

    manifest = read_manifest(args.selection_manifest)
    rejections = read_rejections(args.rejected_candidates)
    proteins = parse_fasta(args.protein_fasta)
    reports = load_reports(args.data_report)
    query_length = len(proteins["NP_009225.1"].sequence)
    selected_records: list[ProteinRecord] = []
    candidate_rows: list[dict[str, str]] = []
    provenance_rows: list[dict[str, str]] = []

    selected_accessions = {row["accession"] for row in manifest}
    rejected_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for rejected in rejections:
        if rejected["replacement_accession"] not in selected_accessions:
            fail(
                "rejected candidate names a replacement that is not selected: "
                f"{rejected['replacement_accession']}"
            )
        record = proteins.get(rejected["accession"])
        report = reports.get(rejected["gene_id"])
        if record is None or report is None:
            fail(f"rejected record is absent from provider package: {rejected['accession']}")
        if rejected["accession"] in selected_accessions:
            fail(f"accession cannot be both selected and rejected: {rejected['accession']}")
        if record.species != rejected["species"] or record.gene_id != rejected["gene_id"]:
            fail(f"rejected provider FASTA metadata mismatch: {rejected['accession']}")
        if (
            str(report.get("taxId")) != rejected["taxon_id"]
            or report.get("taxname") != rejected["species"]
        ):
            fail(f"rejected provider gene metadata mismatch: {rejected['accession']}")
        rejected_by_gene[rejected["gene_id"]].append(rejected)

    for selected in manifest:
        accession = selected["accession"]
        record = proteins.get(accession)
        report = reports.get(selected["gene_id"])
        if record is None or report is None:
            fail(f"selected record is absent from provider package: {accession}")
        if record.species != selected["species"] or record.gene_id != selected["gene_id"]:
            fail(f"provider FASTA metadata mismatch: {accession}")
        if str(report.get("taxId")) != selected["taxon_id"] or report.get("taxname") != selected["species"]:
            fail(f"provider gene metadata mismatch: {accession}")
        ratio = len(record.sequence) / query_length
        if ratio < args.min_length_ratio or ratio > args.max_length_ratio:
            fail(f"selected accession {accession} has out-of-policy length ratio {ratio:.6f}")
        if re.search(r"\b(fragment|partial)\b", record.header, flags=re.IGNORECASE):
            fail(f"selected accession is explicitly labelled fragment/partial: {accession}")
        analysis_group = selected["analysis_group"]
        relation = "self" if analysis_group == "study" else "ortholog"
        role = "outgroup" if analysis_group == "outgroup" else "ingroup"
        outgroup_rationale = (
            "Amphibian NCBI Ortholog group member outside Amniota; retain as a "
            "candidate only until combined outgroup, branch-length, and stability review"
            if analysis_group == "outgroup" else "not applicable"
        )
        candidate_rows.append({
            "accession": accession,
            "taxon_id": selected["taxon_id"],
            "species": selected["species"],
            "role": role,
            "relation": relation,
            "is_reviewed": selected["is_reviewed"],
            "is_canonical": str(accession == "NP_009225.1").lower(),
            "is_fragment": "false",
            "query_coverage": "",
            "sequence_length": str(len(record.sequence)),
            "bitscore": "",
            "evalue": "",
            "source_db": SOURCE_DB,
            "source_release": annotation_context(report),
            "retrieved_at": RETRIEVED_AT,
            "clade": selected["clade"],
            "accession_version": accession,
            "gene_name": str(report.get("symbol") or "BRCA1"),
            "protein_name": "BRCA1 DNA repair associated",
            "lineage": "",
            "analysis_group": analysis_group,
            "target_coverage": "",
            "percent_identity": "",
            "alignment_length": "",
            "orthology_source": "query" if relation == "self" else ORTHOLOGY_SOURCE,
            "orthology_evidence": (
                "NCBI GeneID 672 and RefSeq NP_009225.1" if relation == "self"
                else ORTHOLOGY_EVIDENCE
            ),
            "domain_architecture": DOMAIN_ARCHITECTURE,
            "cluster_id": "",
            "cluster_representative": "",
            "outgroup_rationale": outgroup_rationale,
            "retrieval_query_id": "datasets download gene gene-id 672 --ortholog all",
            "notes": selected["selection_rationale"],
        })
        provenance_rows.append({
            "accession": accession,
            "sequence_sha256": hashlib.sha256(record.sequence.encode("ascii")).hexdigest(),
            "taxon_id": selected["taxon_id"],
            "species": selected["species"],
            "gene_id": selected["gene_id"],
            "sequence_length": str(len(record.sequence)),
            "isoform_label": record.isoform,
            "gene_record_annotation_context": annotation_context(report),
            "selection_rationale": selected["selection_rationale"],
        })
        selected_records.append(record)

    selected_by_gene = {row["gene_id"]: row["accession"] for row in manifest}
    proteins_by_gene: dict[str, list[ProteinRecord]] = defaultdict(list)
    for record in proteins.values():
        if record.gene_id not in reports:
            fail(f"provider FASTA GeneID is absent from data report: {record.gene_id}")
        proteins_by_gene[record.gene_id].append(record)
    inventory_rows: list[dict[str, str]] = []
    for gene_id, report in reports.items():
        records = proteins_by_gene.get(gene_id, [])
        lengths = [len(record.sequence) for record in records]
        eligible = sum(
            args.min_length_ratio <= length / query_length <= args.max_length_ratio
            for length in lengths
        )
        selected_accession = selected_by_gene.get(gene_id, "")
        rejected_rows = rejected_by_gene.get(gene_id, [])
        rejected_accession = ";".join(row["accession"] for row in rejected_rows)
        if selected_accession:
            selection_status = "selected_for_reference_review"
            reason_code = "TAXONOMICALLY_BALANCED_REVIEW_SET"
            evidence_note = (
                "Selected for executable TaxID, similarity, and terminal-domain review; "
                "not yet approved for inference."
            )
        elif rejected_rows:
            selection_status = "rejected_during_reference_review"
            reason_code = ";".join(
                f"{row['stage'].upper()}_FAILED" for row in rejected_rows
            )
            evidence_note = " | ".join(
                f"{row['accession']}: {row['reason']} Replacement: "
                f"{row['replacement_accession']}. Evidence: {row['evidence']}"
                for row in rejected_rows
            )
        elif not records:
            selection_status = "not_selected"
            reason_code = "NO_PROVIDER_PROTEIN_RECORD"
            evidence_note = (
                "The provider gene record contained no protein FASTA record; no downstream "
                "QC failure is asserted."
            )
        elif eligible == 0:
            selection_status = "not_selected"
            reason_code = "NO_SEQUENCE_WITHIN_LENGTH_GATE"
            evidence_note = (
                "No provider protein for this species passed the declared 0.75-1.10 "
                "query-length gate; no terminal-domain failure is asserted."
            )
        else:
            selection_status = "not_selected"
            reason_code = "NOT_SELECTED_TAXONOMIC_BALANCE"
            evidence_note = (
                "Provider species inventoried but not selected for the balanced 50-tip "
                "review set; no QC failure is asserted."
            )
        inventory_rows.append({
            "taxon_id": str(report.get("taxId") or ""),
            "species": str(report.get("taxname") or ""),
            "gene_id": gene_id,
            "provider_protein_count": str(len(records)),
            "minimum_sequence_length": str(min(lengths)) if lengths else "",
            "maximum_sequence_length": str(max(lengths)) if lengths else "",
            "within_declared_length_gate_count": str(eligible),
            "selected_accession": selected_accession,
            "rejected_accession": rejected_accession,
            "selection_status": selection_status,
            "reason_code": reason_code,
            "evidence_note": evidence_note,
        })
    inventory_rows.sort(key=lambda row: (row["species"], row["gene_id"]))

    args.out.mkdir(parents=True)
    write_fasta(args.out / "candidates.faa", selected_records)
    write_fasta(args.out / "query.faa", [proteins["NP_009225.1"]])
    write_tsv(args.out / "candidates.pre-qc.tsv", CANDIDATE_FIELDS, candidate_rows)
    write_tsv(
        args.out / "candidate_provenance.tsv",
        tuple(provenance_rows[0]),
        provenance_rows,
    )
    write_tsv(args.out / "provider_species_inventory.tsv", INVENTORY_FIELDS, inventory_rows)
    summary = {
        "schema_version": "0.2",
        "status": "reference-review-candidate",
        "externally_verified_provider_archive_sha256": args.provider_archive_sha256,
        "provider_protein_fasta_sha256": sha256_file(args.protein_fasta),
        "provider_data_report_sha256": sha256_file(args.data_report),
        "provider_gene_record_count": len(reports),
        "provider_protein_record_count": len(proteins),
        "selected_tip_count": len(manifest),
        "study_count": sum(row["analysis_group"] == "study" for row in manifest),
        "expanded_count": sum(row["analysis_group"] == "expanded" for row in manifest),
        "outgroup_candidate_count": sum(row["analysis_group"] == "outgroup" for row in manifest),
        "unique_selected_taxon_count": len({row["taxon_id"] for row in manifest}),
        "length_ratio_gate": [args.min_length_ratio, args.max_length_ratio],
        "selection_manifest_sha256": sha256_file(args.selection_manifest),
        "rejected_candidate_count": len(rejections),
        "rejected_candidates_sha256": (
            sha256_file(args.rejected_candidates) if args.rejected_candidates else None
        ),
    }
    (args.out / "provider_scope_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        f"Materialized {len(manifest)} review candidates from "
        f"{len(reports)} provider genes / {len(proteins)} proteins."
    )
    return 0


def main() -> None:
    try:
        raise SystemExit(run())
    except (PreparationError, OSError, UnicodeError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=__import__("sys").stderr)
        raise SystemExit(2) from exc


if __name__ == "__main__":
    main()
