"""Tests for deterministic, non-decision-bearing taxonomy lineage export."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    REPOSITORY_ROOT
    / "examples"
    / "brca1"
    / "scripts"
    / "export_taxonomy_lineage.py"
)


class Brca1TaxonomyLineageTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory(
            prefix="brca1-taxonomy-lineage-test-"
        )
        self.root = Path(self.temporary_directory.name)
        self.names = self.root / "names.dmp"
        self.nodes = self.root / "nodes.dmp"
        self.names.write_text(
            "1\t|\troot\t|\t\t|\tscientific name\t|\n"
            "10\t|\tGenus alpha\t|\t\t|\tscientific name\t|\n"
            "11\t|\tSpecies alpha\t|\t\t|\tscientific name\t|\n"
            "11\t|\talias alpha\t|\t\t|\tsynonym\t|\n"
            "20\t|\tSpecies beta\t|\t\t|\tscientific name\t|\n",
            encoding="utf-8",
        )
        self.nodes.write_text(
            "1\t|\t1\t|\tno rank\t|\n"
            "10\t|\t1\t|\tgenus\t|\n"
            "11\t|\t10\t|\tspecies\t|\n"
            "20\t|\t1\t|\tspecies\t|\n",
            encoding="utf-8",
        )
        self.resolution = self.root / "taxonomy_resolution.tsv"
        self.write_resolution()

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    @staticmethod
    def sha256(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def write_resolution(
        self,
        *,
        names_sha256: str | None = None,
        nodes_sha256: str | None = None,
        species_alpha_parent: str = "10",
        species_alpha_rank: str = "species",
        species_alpha_name: str = "Species alpha",
    ) -> None:
        fields = (
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
        )
        with self.resolution.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            # Deliberately reverse lexical accession order.  Output is sorted
            # independently of the planner artifact's row order.
            writer.writerows(
                (
                    {
                        "record_id": "Z_PROTEIN",
                        "matched_name": species_alpha_name,
                        "name_class": "scientific name",
                        "taxon_id": "11",
                        "parent_taxon_id": species_alpha_parent,
                        "rank": species_alpha_rank,
                        "status": "resolved-exact-scientific-name",
                        "snapshot": "synthetic-2026-08-01",
                        "source_url": "https://ftp.ncbi.nlm.nih.gov/example.zip",
                        "retrieved_at": "2026-08-13T00:00:00Z",
                        "names_sha256": names_sha256 or self.sha256(self.names),
                        "nodes_sha256": nodes_sha256 or self.sha256(self.nodes),
                    },
                    {
                        "record_id": "A_PROTEIN",
                        "matched_name": "Species beta",
                        "name_class": "scientific name",
                        "taxon_id": "20",
                        "parent_taxon_id": "1",
                        "rank": "species",
                        "status": "resolved-exact-scientific-name",
                        "snapshot": "synthetic-2026-08-01",
                        "source_url": "https://ftp.ncbi.nlm.nih.gov/example.zip",
                        "retrieved_at": "2026-08-13T00:00:00Z",
                        "names_sha256": names_sha256 or self.sha256(self.names),
                        "nodes_sha256": nodes_sha256 or self.sha256(self.nodes),
                    },
                )
            )

    def run_script(
        self, *, output_name: str = "taxonomy_lineage.tsv"
    ) -> tuple[subprocess.CompletedProcess[str], Path]:
        output = self.root / output_name
        completed = subprocess.run(
            (
                sys.executable,
                str(SCRIPT),
                "--taxonomy-resolution",
                str(self.resolution),
                "--names-dmp",
                str(self.names),
                "--nodes-dmp",
                str(self.nodes),
                "--output",
                str(output),
            ),
            capture_output=True,
            text=True,
            check=False,
        )
        return completed, output

    def test_full_root_to_focal_lineage_and_deterministic_order(self) -> None:
        completed, output = self.run_script()
        self.assertEqual(completed.returncode, 0, completed.stderr)
        with output.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["accession"] for row in rows], ["A_PROTEIN", "Z_PROTEIN"])
        alpha = rows[1]
        self.assertEqual(alpha["lineage_node_count"], "3")
        self.assertEqual(
            json.loads(alpha["lineage_taxon_ids_json"]), ["1", "10", "11"]
        )
        self.assertEqual(
            json.loads(alpha["lineage_scientific_names_json"]),
            ["root", "Genus alpha", "Species alpha"],
        )
        self.assertEqual(
            json.loads(alpha["lineage_ranks_json"]),
            ["no rank", "genus", "species"],
        )
        ranked = json.loads(alpha["ranked_lineage_json"])
        self.assertEqual(ranked[1], {
            "taxon_id": "10",
            "scientific_name": "Genus alpha",
            "rank": "genus",
        })
        self.assertEqual(
            alpha["artifact_role"],
            "post-plan-non-decision-bearing-taxonomy-evidence",
        )
        self.assertNotIn(str(self.root), output.read_text(encoding="utf-8"))

    def test_snapshot_hash_mismatch_fails_without_output(self) -> None:
        self.write_resolution(nodes_sha256="0" * 64)
        completed, output = self.run_script()
        self.assertEqual(completed.returncode, 2)
        self.assertIn("nodes.dmp SHA-256 disagrees", completed.stderr)
        self.assertFalse(output.exists())

    def test_cycle_fails_without_output(self) -> None:
        self.nodes.write_text(
            "1\t|\t1\t|\tno rank\t|\n"
            "10\t|\t11\t|\tgenus\t|\n"
            "11\t|\t10\t|\tspecies\t|\n"
            "20\t|\t1\t|\tspecies\t|\n",
            encoding="utf-8",
        )
        self.write_resolution(species_alpha_parent="10")
        completed, output = self.run_script()
        self.assertEqual(completed.returncode, 2)
        self.assertIn("cycle in nodes.dmp parent chain", completed.stderr)
        self.assertFalse(output.exists())

    def test_resolution_parent_and_rank_are_revalidated(self) -> None:
        self.write_resolution(species_alpha_parent="1")
        completed, output = self.run_script(output_name="parent.tsv")
        self.assertEqual(completed.returncode, 2)
        self.assertIn("parent_taxon_id", completed.stderr)
        self.assertFalse(output.exists())

        self.write_resolution(species_alpha_rank="subspecies")
        completed, output = self.run_script(output_name="rank.tsv")
        self.assertEqual(completed.returncode, 2)
        self.assertIn("resolution rank", completed.stderr)
        self.assertFalse(output.exists())

    def test_every_lineage_node_requires_one_scientific_name(self) -> None:
        self.names.write_text(
            "1\t|\troot\t|\t\t|\tscientific name\t|\n"
            "11\t|\tSpecies alpha\t|\t\t|\tscientific name\t|\n"
            "20\t|\tSpecies beta\t|\t\t|\tscientific name\t|\n",
            encoding="utf-8",
        )
        self.write_resolution()
        completed, output = self.run_script()
        self.assertEqual(completed.returncode, 2)
        self.assertIn("missing a unique scientific name", completed.stderr)
        self.assertIn("10", completed.stderr)
        self.assertFalse(output.exists())

    def test_focal_matched_name_is_revalidated_and_output_is_not_overwritten(self) -> None:
        self.write_resolution(species_alpha_name="Wrong species")
        completed, output = self.run_script(output_name="wrong-name.tsv")
        self.assertEqual(completed.returncode, 2)
        self.assertIn("matched_name", completed.stderr)
        self.assertFalse(output.exists())

        self.write_resolution()
        completed, output = self.run_script(output_name="existing.tsv")
        self.assertEqual(completed.returncode, 0, completed.stderr)
        original = output.read_bytes()
        repeated, _ = self.run_script(output_name="existing.tsv")
        self.assertEqual(repeated.returncode, 2)
        self.assertIn("refusing to overwrite", repeated.stderr)
        self.assertEqual(output.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
