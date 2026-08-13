#!/usr/bin/env python3
"""Prepare the fixed, auditable BRCA1 README example from an NCBI package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


RETRIEVED_AT = "2026-08-13"
SOURCE_DB = "NCBI Gene / RefSeq"
ORTHOLOGY_SOURCE = "NCBI Orthologs"
ORTHOLOGY_EVIDENCE = "NCBI gene group 672; protein similarity, nucleotide alignment, and microsynteny"
DOMAIN_ARCHITECTURE = "pending Pfam+SMART terminal-domain QC"


SELECTED = (
    # accession, taxid, species, gene id, role, clade, reviewed, isoform note
    ("NP_009225.1", "9606", "Homo sapiens", "672", "study", "Primates", True, "RefSeq isoform 1; focal sequence"),
    ("NP_001038958.1", "9598", "Pan troglodytes", "449497", "expanded", "Primates", True, "RefSeq NP record"),
    ("NP_033894.3", "10090", "Mus musculus", "12189", "expanded", "Glires", True, "RefSeq NP record"),
    ("XP_017204702.2", "9986", "Oryctolagus cuniculus", "100347269", "expanded", "Glires", False, "full-length predicted X1 representative"),
    ("NP_001013434.1", "9615", "Canis lupus familiaris", "403437", "expanded", "Laurasiatheria", True, "RefSeq NP record"),
    ("NP_848668.1", "9913", "Bos taurus", "353120", "expanded", "Laurasiatheria", True, "RefSeq NP record"),
    ("XP_014595447.2", "9796", "Equus caballus", "100051990", "expanded", "Laurasiatheria", False, "full-length predicted X1 representative"),
    ("XP_003414318.3", "9785", "Loxodonta africana", "100653763", "expanded", "Afrotheria", False, "full-length predicted representative"),
    ("XP_058140293.1", "9361", "Dasypus novemcinctus", "101429125", "expanded", "Xenarthra", False, "full-length predicted X1 representative"),
    ("NP_001029141.1", "13616", "Monodelphis domestica", "554178", "expanded", "Marsupialia", True, "RefSeq NP record"),
    ("NP_989500.1", "9031", "Gallus gallus", "373983", "expanded", "Aves", True, "RefSeq NP record"),
    ("XP_072775070.1", "59729", "Taeniopygia guttata", "100224649", "expanded", "Aves", False, "predicted X1; X2 logged as an alternative"),
    ("XP_019406054.1", "8502", "Crocodylus porosus", "109320418", "expanded", "Crocodylia", False, "only protein in current NCBI package"),
    ("XP_023967135.2", "8478", "Chrysemys picta bellii", "101942702", "expanded", "Testudines", False, "full-length predicted X1 representative"),
    ("XP_008111382.1", "28377", "Anolis carolinensis", "100553919", "expanded", "Lepidosauria", False, "predicted X1 representative"),
    ("XP_026576759.1", "8673", "Pseudonaja textilis", "113449821", "expanded", "Lepidosauria", False, "shorter model retained only after terminal-domain QC"),
    ("NP_001107963.1", "8364", "Xenopus tropicalis", "733513", "outgroup", "Anura", True, "RefSeq NP amphibian outgroup"),
    ("XP_029429046.1", "194408", "Rhinatrema bivittatum", "115074072", "outgroup", "Gymnophiona", False, "full-length predicted X1 amphibian outgroup"),
)


CANDIDATE_FIELDS = (
    "accession", "taxon_id", "species", "role", "relation", "is_reviewed",
    "is_canonical", "is_fragment", "query_coverage", "sequence_length", "bitscore",
    "evalue", "source_db", "source_release", "retrieved_at", "clade",
    "accession_version", "gene_name", "protein_name", "lineage", "analysis_group",
    "target_coverage", "percent_identity", "alignment_length", "orthology_source",
    "orthology_evidence", "domain_architecture", "cluster_id", "cluster_representative",
    "outgroup_rationale", "retrieval_query_id", "notes",
)


def fail(message: str) -> None:
    raise SystemExit(message)


def parse_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header = ""
    sequence: list[str] = []

    def store() -> None:
        if not header:
            return
        accession = header.split()[0]
        if accession in records:
            fail(f"duplicate FASTA accession: {accession}")
        joined = "".join(sequence).upper()
        if not joined or not re.fullmatch(r"[A-Z*]+", joined) or "*" in joined:
            fail(f"invalid protein sequence for {accession}")
        records[accession] = (header, joined)

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
    return records


def load_reports(path: Path) -> dict[str, dict]:
    reports: dict[str, dict] = {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            record = json.loads(raw)
            reports[str(record["taxId"])] = record
    return reports


def source_release(report: dict) -> str:
    annotations = report.get("annotations") or []
    if not annotations:
        return "NCBI Gene record without assembly annotation label"
    annotation = annotations[0]
    return str(annotation.get("annotationName") or annotation.get("assemblyAccession") or "")


def write_fasta(path: Path, selected_records: list[tuple[tuple, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for definition, sequence in selected_records:
            accession, _taxid, species, _gene_id, _role, _clade, _reviewed, _note = definition
            handle.write(f">{accession} BRCA1 [organism={species}]\n")
            for index in range(0, len(sequence), 70):
                handle.write(sequence[index:index + 70] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein-fasta", type=Path, required=True)
    parser.add_argument("--data-report", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        fail(f"refusing to overwrite output directory: {args.out}")
    args.out.mkdir(parents=True)

    fasta = parse_fasta(args.protein_fasta)
    reports = load_reports(args.data_report)
    chosen: list[tuple[tuple, str]] = []
    rows: list[dict[str, str]] = []
    provenance: list[dict[str, str]] = []

    for definition in SELECTED:
        accession, taxid, species, gene_id, analysis_group, clade, reviewed, note = definition
        if accession not in fasta:
            fail(f"selected accession absent from package: {accession}")
        header, sequence = fasta[accession]
        report = reports.get(taxid)
        if report is None:
            fail(f"selected TaxID absent from data report: {taxid}")
        if report.get("taxname") != species or str(report.get("geneId")) != gene_id:
            fail(f"metadata mismatch for {accession}")
        if f"[organism={species}]" not in header or f"[GeneID={gene_id}]" not in header:
            fail(f"FASTA metadata mismatch for {accession}")
        role = "outgroup" if analysis_group == "outgroup" else "ingroup"
        relation = "self" if analysis_group == "study" else "one2one_ortholog"
        outgroup_rationale = (
            "Amphibian NCBI ortholog outside Amniota; root only if both amphibians form the approved root split"
            if role == "outgroup" else "not applicable"
        )
        rows.append({
            "accession": accession,
            "taxon_id": taxid,
            "species": species,
            "role": role,
            "relation": relation,
            "is_reviewed": str(reviewed).lower(),
            # RefSeq NP_ denotes a curated record, not necessarily the
            # organism's canonical isoform.  Canonical status is asserted only
            # for the human focal isoform and is otherwise left false.
            "is_canonical": str(accession == "NP_009225.1").lower(),
            "is_fragment": "false",
            "query_coverage": "",
            "sequence_length": str(len(sequence)),
            "bitscore": "",
            "evalue": "",
            "source_db": SOURCE_DB,
            "source_release": source_release(report),
            "retrieved_at": RETRIEVED_AT,
            "clade": clade,
            "accession_version": accession,
            "gene_name": str(report.get("symbol") or "BRCA1"),
            "protein_name": "BRCA1 DNA repair associated",
            "lineage": "",
            "analysis_group": analysis_group,
            "target_coverage": "",
            "percent_identity": "",
            "alignment_length": "",
            "orthology_source": ORTHOLOGY_SOURCE if relation != "self" else "query",
            "orthology_evidence": ORTHOLOGY_EVIDENCE if relation != "self" else "NCBI Gene 672 and RefSeq NP_009225.1",
            "domain_architecture": DOMAIN_ARCHITECTURE,
            "cluster_id": "",
            "cluster_representative": "",
            "outgroup_rationale": outgroup_rationale,
            "retrieval_query_id": "datasets download gene gene-id 672 --ortholog all",
            "notes": note,
        })
        provenance.append({
            "accession": accession,
            "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            "taxon_id": taxid,
            "species": species,
            "gene_id": gene_id,
            "sequence_length": str(len(sequence)),
            "annotation_release": source_release(report),
            "annotation_release_date": str((report.get("annotations") or [{}])[0].get("annotationReleaseDate") or ""),
            "decision": note,
        })
        chosen.append((definition, sequence))

    write_fasta(args.out / "candidates.faa", chosen)
    query = [item for item in chosen if item[0][0] == "NP_009225.1"]
    write_fasta(args.out / "query.faa", query)
    with (args.out / "candidates.pre-qc.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CANDIDATE_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with (args.out / "candidate_provenance.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = tuple(provenance[0])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(provenance)


if __name__ == "__main__":
    main()
