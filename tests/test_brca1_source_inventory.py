"""Audit the record-level source inventory for the BRCA1 worked example."""

from __future__ import annotations

import csv
import re
import unittest
from collections import Counter
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = REPOSITORY_ROOT / "examples" / "brca1"
INVENTORY = EXAMPLE / "inputs" / "source_protein_inventory.tsv"

EXPECTED_PRESENT_COUNTS = {
    "Anolis carolinensis": 2,
    "Bos taurus": 11,
    "Canis lupus familiaris": 12,
    "Chrysemys picta bellii": 5,
    "Crocodylus porosus": 1,
    "Dasypus novemcinctus": 6,
    "Equus caballus": 31,
    "Gallus gallus": 5,
    "Homo sapiens": 368,
    "Loxodonta africana": 2,
    "Monodelphis domestica": 4,
    "Mus musculus": 5,
    "Ornithorhynchus anatinus": 2,
    "Oryctolagus cuniculus": 4,
    "Pan troglodytes": 51,
    "Pseudonaja textilis": 1,
    "Rhinatrema bivittatum": 4,
    "Sphenodon punctatus": 0,
    "Taeniopygia guttata": 2,
    "Xenopus tropicalis": 2,
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class Brca1SourceInventoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.rows = read_tsv(INVENTORY)
        cls.candidates = read_tsv(EXAMPLE / "inputs" / "candidates.tsv")

    def test_inventory_is_lossless_for_declared_species_scope(self) -> None:
        present = [row for row in self.rows if row["source_record_status"] == "present"]
        absent = [row for row in self.rows if row["source_record_status"] == "absent"]
        self.assertEqual(len(present), 518)
        self.assertEqual(len(absent), 1)
        self.assertEqual(len(self.rows), 519)

        accessions = [row["accession_version"] for row in present]
        self.assertEqual(len(accessions), len(set(accessions)))
        self.assertTrue(all(re.fullmatch(r"[A-Z]{2}_[0-9]+\.[0-9]+", item) for item in accessions))

        counts = Counter(row["species"] for row in present)
        self.assertEqual(
            {species: counts.get(species, 0) for species in EXPECTED_PRESENT_COUNTS},
            EXPECTED_PRESENT_COUNTS,
        )
        self.assertEqual(set(counts) | {row["species"] for row in absent}, set(EXPECTED_PRESENT_COUNTS))

    def test_selected_rows_match_the_committed_fixed_set(self) -> None:
        selected = {
            row["accession_version"]: row
            for row in self.rows
            if row["selection_status"] == "selected"
        }
        candidates = {row["accession"]: row for row in self.candidates}
        self.assertEqual(set(selected), set(candidates))
        self.assertEqual(len({row["species"] for row in selected.values()}), 18)
        for accession, source in selected.items():
            candidate = candidates[accession]
            self.assertEqual(source["taxon_id"], candidate["taxon_id"])
            self.assertEqual(source["species"], candidate["species"])
            self.assertEqual(source["sequence_length"], candidate["sequence_length"])
            self.assertEqual(source["reason_code"], "FIXED_ANALYSIS_SET")

    def test_unpromoted_provider_records_do_not_claim_unrun_qc(self) -> None:
        alternatives = [
            row for row in self.rows
            if row["archive_scope"] == "selected_species"
            and row["selection_status"] == "not_selected"
        ]
        self.assertEqual(len(alternatives), 498)
        self.assertEqual(
            {row["reason_code"] for row in alternatives},
            {"NOT_PROMOTED_ONE_PER_SPECIES_SET"},
        )
        for row in alternatives:
            note = row["evidence_note"].lower()
            self.assertIn("not individually subjected", note)
            self.assertNotIn("failed", note)

        zebra_x2 = next(
            row for row in alternatives
            if row["accession_version"] == "XP_030112190.4"
        )
        self.assertEqual(zebra_x2["species"], "Taeniopygia guttata")
        self.assertEqual(zebra_x2["sequence_length"], "1804")
        self.assertEqual(zebra_x2["isoform_label"], "X2")

    def test_extra_taxa_have_only_evidence_supported_reasons(self) -> None:
        platypus = [
            row for row in self.rows
            if row["species"] == "Ornithorhynchus anatinus"
        ]
        self.assertEqual(
            {row["accession_version"] for row in platypus},
            {"XP_028930515.1", "XP_028930516.1"},
        )
        self.assertTrue(all(row["sequence_length"] == "1346" for row in platypus))
        self.assertTrue(all(row["length_ratio_to_query"] == "0.722491" for row in platypus))
        self.assertTrue(
            all(row["reason_code"] == "BELOW_CONFIGURED_LENGTH_RATIO" for row in platypus)
        )
        self.assertTrue(all("no domain-failure claim" in row["evidence_note"] for row in platypus))

        tuatara = [row for row in self.rows if row["species"] == "Sphenodon punctatus"]
        self.assertEqual(len(tuatara), 1)
        self.assertEqual(tuatara[0]["taxon_id"], "8508")
        self.assertEqual(tuatara[0]["accession_version"], "")
        self.assertEqual(tuatara[0]["source_record_status"], "absent")
        self.assertEqual(tuatara[0]["reason_code"], "NO_RECORD_IN_FROZEN_NCBI_PACKAGE")


if __name__ == "__main__":
    unittest.main()
