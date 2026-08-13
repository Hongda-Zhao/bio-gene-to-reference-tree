"""Contract tests for the executable BRCA1 candidate-QC example."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPOSITORY_ROOT / "examples" / "brca1" / "scripts" / "summarize_brca1_qc.py"


class Brca1QcTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory(prefix="brca1-qc-test-")
        self.root = Path(self.temporary_directory.name)
        sequence = "A" * 100
        (self.root / "candidates.faa").write_text(
            f">QUERY\n{sequence}\n>SUBJECT\n{sequence}\n", encoding="utf-8"
        )
        fields = (
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
        )
        with (self.root / "candidates.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(
                (
                    {
                        "accession": "QUERY",
                        "relation": "self",
                        "is_canonical": "true",
                        "is_fragment": "false",
                        "query_coverage": "1",
                        "sequence_length": "100",
                        "bitscore": "0",
                        "evalue": "0",
                        "target_coverage": "1",
                        "percent_identity": "",
                        "alignment_length": "",
                        "domain_architecture": "unverified",
                        "notes": "focal",
                    },
                    {
                        "accession": "SUBJECT",
                        "relation": "one2one_ortholog",
                        "is_canonical": "true",
                        "is_fragment": "false",
                        "query_coverage": "1",
                        "sequence_length": "100",
                        "bitscore": "0",
                        "evalue": "0",
                        "target_coverage": "1",
                        "percent_identity": "",
                        "alignment_length": "",
                        "domain_architecture": "unverified",
                        "notes": "candidate",
                    },
                )
            )
        # Two overlapping HSPs exercise union coverage: 1-60 plus 51-100 is
        # 100 covered residues, while alignment_length remains the auditable
        # sum of 110 aligned columns.
        (self.root / "blast.tsv").write_text(
            "QUERY\tQUERY\t60\t60\t0\t0\t1\t60\t1\t60\t0\t120\t100\t100\n"
            "QUERY\tQUERY\t50\t50\t0\t0\t51\t100\t51\t100\t1e-30\t100\t100\t100\n"
            "QUERY\tSUBJECT\t48\t60\t12\t0\t1\t60\t1\t60\t1e-20\t90\t100\t100\n"
            "QUERY\tSUBJECT\t35\t50\t15\t0\t51\t100\t51\t100\t2e-15\t70\t100\t100\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    @staticmethod
    def interpro_row(
        accession: str,
        analysis: str,
        signature: str,
        description: str,
        start: int,
        end: int,
        interpro: str = "-",
        interpro_description: str = "-",
    ) -> str:
        return "\t".join(
            (
                accession,
                "synthetic-md5",
                "100",
                analysis,
                signature,
                description,
                str(start),
                str(end),
                "1e-10",
                "T",
                "13-08-2026",
                interpro,
                interpro_description,
            )
        )

    def write_interpro(self, *, include_second_repeat: bool = True) -> Path:
        rows: list[str] = []
        for accession in ("QUERY", "SUBJECT"):
            rows.extend(
                (
                    self.interpro_row(
                        accession, "Pfam", "CUSTOM_RING", "RING finger domain", 4, 15
                    ),
                    self.interpro_row(
                        accession, "SMART", "SM00184", "RING", 5, 14
                    ),
                    self.interpro_row(
                        accession, "Pfam", "CUSTOM_BRCT_1", "BRCT domain", 70, 79
                    ),
                    self.interpro_row(
                        accession, "SMART", "SM00292", "BRCT", 71, 80
                    ),
                )
            )
            if include_second_repeat:
                rows.extend(
                    (
                        self.interpro_row(
                            accession,
                            "Pfam",
                            "CUSTOM_BRCT_2",
                            "BRCA1 C-terminal domain",
                            86,
                            96,
                        ),
                        self.interpro_row(
                            accession, "SMART", "SM00292", "BRCT", 85, 95
                        ),
                    )
                )
        path = self.root / "interpro.tsv"
        path.write_text("\n".join(rows) + "\n", encoding="utf-8")
        return path

    def run_script(self, interpro: Path) -> tuple[subprocess.CompletedProcess[str], Path]:
        output = self.root / "output"
        completed = subprocess.run(
            (
                sys.executable,
                str(SCRIPT),
                "--candidate-fasta",
                str(self.root / "candidates.faa"),
                "--candidate-table",
                str(self.root / "candidates.tsv"),
                "--blast-tsv",
                str(self.root / "blast.tsv"),
                "--interpro-tsv",
                str(interpro),
                "--out-dir",
                str(output),
            ),
            capture_output=True,
            text=True,
            check=False,
        )
        return completed, output

    def test_union_coverage_distinct_brct_repeats_and_honest_canonical_flag(self) -> None:
        completed, output = self.run_script(self.write_interpro())
        self.assertEqual(completed.returncode, 0, msg=completed.stderr)
        with (output / "candidate_qc.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            qc = {row["accession"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(qc["SUBJECT"]["query_covered_residues"], "100")
        self.assertEqual(qc["SUBJECT"]["subject_covered_residues"], "100")
        self.assertEqual(qc["SUBJECT"]["alignment_length"], "110")
        self.assertEqual(qc["SUBJECT"]["brct_repeat_count"], "2")
        self.assertEqual(qc["SUBJECT"]["domain_qc_status"], "pass")
        with (output / "candidates.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            candidates = {
                row["accession"]: row for row in csv.DictReader(handle, delimiter="\t")
            }
        self.assertEqual(candidates["QUERY"]["is_canonical"], "true")
        self.assertEqual(candidates["SUBJECT"]["is_canonical"], "false")
        self.assertEqual(candidates["SUBJECT"]["query_coverage"], "1")
        self.assertIn("tandem C-terminal BRCT", candidates["SUBJECT"]["domain_architecture"])

    def test_missing_second_brct_repeat_fails_without_outputs(self) -> None:
        completed, output = self.run_script(
            self.write_interpro(include_second_repeat=False)
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("only 1 distinct C-terminal BRCT", completed.stderr)
        self.assertFalse((output / "candidate_qc.tsv").exists())
        self.assertFalse((output / "candidates.tsv").exists())


if __name__ == "__main__":
    unittest.main()
