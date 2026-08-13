#!/usr/bin/env python3
"""Inventory provider proteins considered around the fixed BRCA1 example set.

This is an acquisition-audit helper, not a candidate selector.  It records every
protein in the frozen NCBI package for the 18 sampled species and two explicitly
screened extra taxa.  Unselected provider records are not represented as having
failed BLAST, domain, or tree QC unless an executable test actually established
that fact.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path


EXTRA_SCREENED_TAXA = (
    "Ornithorhynchus anatinus",
    "Sphenodon punctatus",
)

FIELDS = (
    "accession_version",
    "taxon_id",
    "species",
    "gene_id",
    "sequence_length",
    "length_ratio_to_query",
    "isoform_label",
    "archive_scope",
    "source_record_status",
    "selection_status",
    "reason_code",
    "evidence_note",
)


@dataclass(frozen=True)
class ProteinRecord:
    accession: str
    organism: str
    gene_id: str
    isoform: str
    sequence_length: int


def fail(message: str) -> None:
    raise SystemExit(message)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        fail(f"empty TSV: {path}")
    return rows


def parse_fasta(path: Path) -> list[ProteinRecord]:
    records: list[ProteinRecord] = []
    header = ""
    sequence_length = 0
    seen: set[str] = set()

    def store() -> None:
        nonlocal header, sequence_length
        if not header:
            return
        accession = header.split()[0]
        organism_match = re.search(r"\[organism=([^]]+)\]", header)
        gene_match = re.search(r"\[GeneID=([^]]+)\]", header)
        isoform_match = re.search(r"\[isoform=([^]]+)\]", header)
        if not organism_match or not gene_match:
            fail(f"missing organism or GeneID in FASTA header: {accession}")
        if accession in seen:
            fail(f"duplicate FASTA accession: {accession}")
        if sequence_length < 1:
            fail(f"empty FASTA sequence: {accession}")
        seen.add(accession)
        records.append(
            ProteinRecord(
                accession=accession,
                organism=organism_match.group(1),
                gene_id=gene_match.group(1),
                isoform=isoform_match.group(1) if isoform_match else "",
                sequence_length=sequence_length,
            )
        )

    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                store()
                header = line[1:]
                sequence_length = 0
            else:
                if not re.fullmatch(r"[A-Za-z*.-]+", line):
                    fail(f"invalid FASTA sequence line after: {header.split()[0]}")
                sequence_length += len(line.replace("-", "").replace(".", ""))
    store()
    return records


def load_gene_reports(path: Path) -> dict[str, dict[str, str]]:
    reports: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            record = json.loads(raw)
            gene_id = str(record.get("geneId") or "")
            taxon_id = str(record.get("taxId") or "")
            species = str(record.get("taxname") or "")
            if not gene_id or not taxon_id or not species:
                fail(f"incomplete data-report row {line_number}")
            if gene_id in reports:
                fail(f"duplicate GeneID in data report: {gene_id}")
            reports[gene_id] = {
                "taxon_id": taxon_id,
                "species": species,
            }
    return reports


def exact_scientific_taxon_id(path: Path, species: str) -> str:
    matches: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            fields = [field.strip() for field in raw.split("|")]
            if len(fields) >= 4 and fields[1] == species and fields[3] == "scientific name":
                matches.add(fields[0])
    if len(matches) != 1:
        fail(
            f"expected one exact NCBI scientific-name match for {species!r}, "
            f"found {sorted(matches)}"
        )
    return next(iter(matches))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein-fasta", type=Path, required=True)
    parser.add_argument("--data-report", type=Path, required=True)
    parser.add_argument("--names-dmp", type=Path, required=True)
    parser.add_argument("--selected-candidates", type=Path, required=True)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.output.exists():
        fail(f"refusing to overwrite output: {args.output}")

    selected_rows = read_tsv(args.selected_candidates)
    selected_by_accession = {row["accession"]: row for row in selected_rows}
    if len(selected_by_accession) != len(selected_rows):
        fail("selected candidate table contains duplicate accessions")
    selected_species = {row["species"] for row in selected_rows}
    if len(selected_species) != len(selected_rows):
        fail("fixed BRCA1 set must contain no more than one protein per species")

    request = json.loads(args.request.read_text(encoding="utf-8"))
    min_length_ratio = float(request["selection"]["min_length_ratio"])
    query_accession = str(request["query"]["id"])
    if query_accession not in selected_by_accession:
        fail(f"query accession absent from selected candidates: {query_accession}")
    query_length = int(selected_by_accession[query_accession]["sequence_length"])

    target_species = selected_species | set(EXTRA_SCREENED_TAXA)
    taxonomy_ids = {
        species: exact_scientific_taxon_id(args.names_dmp, species)
        for species in target_species
    }
    reports = load_gene_reports(args.data_report)
    proteins = parse_fasta(args.protein_fasta)
    source_records = [record for record in proteins if record.organism in target_species]
    source_accessions = {record.accession for record in source_records}
    missing_selected = sorted(set(selected_by_accession) - source_accessions)
    if missing_selected:
        fail(f"selected accessions absent from provider FASTA: {missing_selected}")

    output_rows: list[dict[str, str]] = []
    observed_species: set[str] = set()
    for record in source_records:
        report = reports.get(record.gene_id)
        if report is None:
            fail(f"FASTA GeneID absent from data report: {record.gene_id}")
        if report["species"] != record.organism:
            fail(f"organism mismatch for {record.accession}")
        if report["taxon_id"] != taxonomy_ids[record.organism]:
            fail(f"TaxID mismatch for {record.accession}")

        observed_species.add(record.organism)
        selected = record.accession in selected_by_accession
        length_ratio = record.sequence_length / query_length
        if selected:
            chosen = selected_by_accession[record.accession]
            if (
                chosen["species"] != record.organism
                or chosen["taxon_id"] != report["taxon_id"]
                or int(chosen["sequence_length"]) != record.sequence_length
            ):
                fail(f"selected metadata mismatch for {record.accession}")
            reason_code = "FIXED_ANALYSIS_SET"
            evidence_note = (
                "Member of the committed 18-tip set; downstream sequence and domain QC "
                "is reported in the separate QC artifacts."
            )
        elif record.organism in selected_species:
            reason_code = "NOT_PROMOTED_ONE_PER_SPECIES_SET"
            evidence_note = (
                "Provider record inventoried but not promoted into the fixed one-per-species "
                "set; it was not individually subjected to downstream BLAST, domain, or tree QC."
            )
        elif length_ratio < min_length_ratio:
            reason_code = "BELOW_CONFIGURED_LENGTH_RATIO"
            evidence_note = (
                f"Length ratio {length_ratio:.6f} is below request min_length_ratio "
                f"{min_length_ratio:.6f}; no domain-failure claim was made."
            )
        else:
            reason_code = "NOT_PROMOTED_EXTRA_TAXON"
            evidence_note = (
                "Extra taxon was screened at acquisition but not promoted into the fixed set; "
                "no downstream QC failure is asserted."
            )

        output_rows.append({
            "accession_version": record.accession,
            "taxon_id": report["taxon_id"],
            "species": record.organism,
            "gene_id": record.gene_id,
            "sequence_length": str(record.sequence_length),
            "length_ratio_to_query": f"{length_ratio:.6f}",
            "isoform_label": record.isoform,
            "archive_scope": (
                "selected_species" if record.organism in selected_species
                else "screened_extra_taxon"
            ),
            "source_record_status": "present",
            "selection_status": "selected" if selected else "not_selected",
            "reason_code": reason_code,
            "evidence_note": evidence_note,
        })

    for species in sorted(target_species - observed_species):
        output_rows.append({
            "accession_version": "",
            "taxon_id": taxonomy_ids[species],
            "species": species,
            "gene_id": "",
            "sequence_length": "",
            "length_ratio_to_query": "",
            "isoform_label": "",
            "archive_scope": (
                "selected_species" if species in selected_species else "screened_extra_taxon"
            ),
            "source_record_status": "absent",
            "selection_status": "not_selected",
            "reason_code": "NO_RECORD_IN_FROZEN_NCBI_PACKAGE",
            "evidence_note": (
                "Exact scientific-name TaxID was present in the frozen names.dmp, but no "
                "matching gene-report or protein-FASTA record was present in the frozen package."
            ),
        })

    output_rows.sort(
        key=lambda row: (
            row["species"],
            0 if row["selection_status"] == "selected" else 1,
            row["accession_version"],
        )
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output_rows)


if __name__ == "__main__":
    main()
