"""Contract tests for strict local NCBI Taxonomy resolution."""

from __future__ import annotations

import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SKILL_ROOT = REPOSITORY_ROOT / "skills" / "bio-gene-to-reference-tree"
SCRIPT = SKILL_ROOT / "scripts" / "ncbi_taxonomy.py"
FIXTURE = REPOSITORY_ROOT / "tests" / "fixtures" / "taxonomy"

SPEC = importlib.util.spec_from_file_location("ncbi_taxonomy", SCRIPT)
assert SPEC and SPEC.loader
TAXONOMY = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = TAXONOMY
SPEC.loader.exec_module(TAXONOMY)


class StrictTaxonomyTests(unittest.TestCase):
    def resolve(self, *requests):
        return TAXONOMY.resolve_exact_scientific_names(
            list(requests),
            names_path=FIXTURE / "names.dmp",
            nodes_path=FIXTURE / "nodes.dmp",
            snapshot="synthetic-2026-08-01",
            source_url=TAXONOMY.OFFICIAL_CURRENT_ARCHIVE,
            retrieved_at="2026-08-13T00:00:00Z",
        )

    def test_exact_scientific_name_and_taxid_resolve(self) -> None:
        resolution = self.resolve(TAXONOMY.NameRequest("QUERY_001", "Homo sapiens", "9606"))
        self.assertEqual(len(resolution.rows), 1)
        row = resolution.rows[0]
        self.assertEqual(row["taxon_id"], "9606")
        self.assertEqual(row["parent_taxon_id"], "9605")
        self.assertEqual(row["rank"], "species")
        self.assertEqual(row["name_class"], "scientific name")
        self.assertEqual(row["status"], "resolved-exact-scientific-name")
        self.assertEqual(row["resolver_version"], TAXONOMY.VERSION)
        self.assertRegex(resolution.names_sha256, r"^[0-9a-f]{64}$")
        self.assertRegex(resolution.nodes_sha256, r"^[0-9a-f]{64}$")

    def test_matching_is_case_sensitive_and_does_not_accept_aliases(self) -> None:
        for name in ("homo sapiens", "human"):
            with self.subTest(name=name), self.assertRaises(TAXONOMY.TaxonomyError) as caught:
                self.resolve(TAXONOMY.NameRequest("query", name))
            self.assertIn("UNRESOLVED_SCIENTIFIC_NAME", caught.exception.code)

    def test_ambiguous_exact_scientific_name_is_not_auto_selected(self) -> None:
        with self.assertRaises(TAXONOMY.TaxonomyError) as caught:
            self.resolve(TAXONOMY.NameRequest("query", "Duplicata exacta"))
        self.assertIn("AMBIGUOUS_SCIENTIFIC_NAME", caught.exception.code)
        self.assertEqual(caught.exception.details[0]["candidate_taxon_ids"], "10;11")

    def test_taxid_mismatch_fails_closed(self) -> None:
        with self.assertRaises(TAXONOMY.TaxonomyError) as caught:
            self.resolve(TAXONOMY.NameRequest("query", "Homo sapiens", "10090"))
        self.assertEqual(caught.exception.code, "TAXON_ID_MISMATCH")

    def test_leading_or_trailing_whitespace_is_not_silently_normalized(self) -> None:
        with self.assertRaises(TAXONOMY.TaxonomyError) as caught:
            self.resolve(TAXONOMY.NameRequest("query", " Homo sapiens"))
        self.assertEqual(caught.exception.code, "NONEXACT_INPUT_NAME")

    def test_batch_cli_writes_auditable_tsv_and_refuses_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "taxonomy_resolution.tsv"
            command = [
                sys.executable,
                str(SCRIPT),
                "--names",
                str(FIXTURE / "names.dmp"),
                "--nodes",
                str(FIXTURE / "nodes.dmp"),
                "--snapshot",
                "synthetic-2026-08-01",
                "--source-url",
                TAXONOMY.OFFICIAL_CURRENT_ARCHIVE,
                "--retrieved-at",
                "2026-08-13T00:00:00Z",
                "--input",
                str(FIXTURE / "organisms.tsv"),
                "--out",
                str(output),
            ]
            completed = subprocess.run(command, text=True, capture_output=True, check=False)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with output.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["taxon_id"] for row in rows], ["9606", "10090"])
            self.assertTrue(all(row["source_url"].startswith("https://ftp.ncbi.nlm.nih.gov/") for row in rows))

            repeated = subprocess.run(command, text=True, capture_output=True, check=False)
            self.assertEqual(repeated.returncode, 2)
            self.assertIn("OUTPUT_EXISTS", repeated.stderr)

    def test_malformed_dmp_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            malformed = Path(temporary) / "names.dmp"
            malformed.write_text("9606\tHomo sapiens\n", encoding="utf-8")
            with self.assertRaises(TAXONOMY.TaxonomyError) as caught:
                TAXONOMY.resolve_exact_scientific_names(
                    [TAXONOMY.NameRequest("query", "Homo sapiens")],
                    names_path=malformed,
                    nodes_path=FIXTURE / "nodes.dmp",
                    snapshot="synthetic",
                    source_url=TAXONOMY.OFFICIAL_CURRENT_ARCHIVE,
                    retrieved_at="2026-08-13",
                )
            self.assertEqual(caught.exception.code, "MALFORMED_DMP")

    def test_non_ncbi_or_directory_source_url_is_rejected(self) -> None:
        for source_url in (
            "https://example.org/new_taxdump.tar.gz",
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/",
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt",
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/new_taxdump.tar.gz",
            "http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz",
        ):
            with self.subTest(source_url=source_url), self.assertRaises(
                TAXONOMY.TaxonomyError
            ) as caught:
                TAXONOMY.resolve_exact_scientific_names(
                    [TAXONOMY.NameRequest("query", "Homo sapiens")],
                    names_path=FIXTURE / "names.dmp",
                    nodes_path=FIXTURE / "nodes.dmp",
                    snapshot="synthetic",
                    source_url=source_url,
                    retrieved_at="2026-08-13",
                )
            self.assertEqual(caught.exception.code, "UNOFFICIAL_SOURCE_URL")

    def test_retrieval_time_requires_iso8601_and_utc_offsets(self) -> None:
        for retrieved_at in (
            "yesterday",
            "2026-08-13T09:00:00+09:00",
            "2026-08-13T12:00:00",
            "2026-08-13 12:00:00",
        ):
            with self.subTest(retrieved_at=retrieved_at), self.assertRaises(
                TAXONOMY.TaxonomyError
            ) as caught:
                TAXONOMY.resolve_exact_scientific_names(
                    [TAXONOMY.NameRequest("query", "Homo sapiens")],
                    names_path=FIXTURE / "names.dmp",
                    nodes_path=FIXTURE / "nodes.dmp",
                    snapshot="synthetic",
                    source_url=TAXONOMY.OFFICIAL_CURRENT_ARCHIVE,
                    retrieved_at=retrieved_at,
                )
            self.assertEqual(caught.exception.code, "INVALID_RETRIEVED_AT")


if __name__ == "__main__":
    unittest.main()
