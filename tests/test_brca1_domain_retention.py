"""Tests for the fail-closed BRCA1 domain-retention audit."""

from __future__ import annotations

import csv
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
    / "audit_domain_retention.py"
)


class Brca1DomainRetentionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory(
            prefix="brca1-domain-retention-test-"
        )
        self.root = Path(self.temporary_directory.name)
        self.raw = self.root / "raw.faa"
        self.raw.write_text(
            ">A\nA-BCDEFG-H\n>B\nAQBCDEFRSH\n",
            encoding="utf-8",
        )
        self.balanced = self.root / "balanced.faa"
        self.balanced.write_text(
            ">B\nAQBCEFRS\n>A\nA-BCEFG-\n",
            encoding="utf-8",
        )
        self.strict = self.root / "strict.faa"
        self.strict.write_text(
            ">A\nABCEFH\n>B\nABCEFS\n",
            encoding="utf-8",
        )
        self.qc = self.root / "candidate_qc.tsv"
        with self.qc.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=(
                    "accession",
                    "subject_length",
                    "ring_coordinates",
                    "brct_coordinates",
                    "domain_qc_status",
                ),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(
                (
                    {
                        "accession": "A",
                        "subject_length": "8",
                        "ring_coordinates": "1-2",
                        "brct_coordinates": "6-8",
                        "domain_qc_status": "pass",
                    },
                    {
                        "accession": "B",
                        "subject_length": "10",
                        "ring_coordinates": "1-2",
                        "brct_coordinates": "8-10",
                        "domain_qc_status": "pass",
                    },
                )
            )

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def run_script(
        self,
        *profiles: tuple[str, Path],
        raw: Path | None = None,
        qc: Path | None = None,
        threshold: str = "0.8",
        output_name: str = "report.tsv",
    ) -> tuple[subprocess.CompletedProcess[str], Path]:
        output = self.root / output_name
        command = [
            sys.executable,
            str(SCRIPT),
            "--raw-alignment",
            str(raw or self.raw),
        ]
        for label, path in profiles:
            command.extend(("--trimmed-alignment", f"{label}={path}"))
        command.extend(
            (
                "--candidate-qc",
                str(qc or self.qc),
                "--minimum-retained-fraction",
                threshold,
                "--output",
                str(output),
            )
        )
        completed = subprocess.run(
            command, capture_output=True, text=True, check=False
        )
        return completed, output

    def test_unique_column_map_and_per_interval_counts(self) -> None:
        completed, output = self.run_script(
            ("balanced", self.balanced), threshold="0.6"
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        with output.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 4)
        indexed = {
            (row["accession"], row["domain"], row["interval_index"]): row
            for row in rows
        }
        self.assertEqual(indexed[("A", "ring", "1")]["raw_msa_start"], "1")
        self.assertEqual(indexed[("A", "ring", "1")]["raw_msa_end"], "3")
        self.assertEqual(indexed[("A", "brct", "1")]["total_residues"], "3")
        self.assertEqual(indexed[("A", "brct", "1")]["retained_residues"], "2")
        self.assertEqual(indexed[("A", "brct", "1")]["retained_fraction"], "0.666667")

    def test_threshold_failure_writes_explicit_failed_report_and_exits_two(self) -> None:
        completed, output = self.run_script(
            ("balanced", self.balanced), threshold="0.8"
        )
        # A's terminal H and B's terminal H/S are deleted in the balanced
        # profile, so both BRCT intervals fall below 0.8.
        self.assertEqual(completed.returncode, 2)
        self.assertIn("domain-retention threshold failed", completed.stderr)
        self.assertTrue(output.is_file())
        with output.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        failed = [row for row in rows if row["status"] == "fail"]
        self.assertEqual(
            {(row["accession"], row["domain"]) for row in failed},
            {("A", "brct"), ("B", "brct")},
        )

    def test_profile_tip_set_must_exactly_match(self) -> None:
        mismatch = self.root / "mismatch.faa"
        mismatch.write_text(">A\nA-BCEFG-\n>C\nAQBCEFRS\n", encoding="utf-8")
        completed, output = self.run_script(("mismatch", mismatch))
        self.assertEqual(completed.returncode, 2)
        self.assertIn("tip-set mismatch", completed.stderr)
        self.assertFalse(output.exists())

    def test_non_subset_column_fails_without_report(self) -> None:
        changed = self.root / "changed.faa"
        changed.write_text(">A\nA-ZCEFG-\n>B\nAQBCEFRS\n", encoding="utf-8")
        completed, output = self.run_script(("changed", changed))
        self.assertEqual(completed.returncode, 2)
        self.assertIn("not an in-order column subset", completed.stderr)
        self.assertFalse(output.exists())

    def test_ambiguous_column_map_fails_instead_of_guessing(self) -> None:
        ambiguous_raw = self.root / "ambiguous-raw.faa"
        ambiguous_raw.write_text(">A\nABAC\n>B\nABAC\n", encoding="utf-8")
        ambiguous_trimmed = self.root / "ambiguous-trimmed.faa"
        ambiguous_trimmed.write_text(">A\nAC\n>B\nAC\n", encoding="utf-8")
        ambiguous_qc = self.root / "ambiguous-qc.tsv"
        ambiguous_qc.write_text(
            "accession\tsubject_length\tring_coordinates\tbrct_coordinates\t"
            "domain_qc_status\n"
            "A\t4\t1-1\t3-4\tpass\n"
            "B\t4\t1-1\t3-4\tpass\n",
            encoding="utf-8",
        )
        completed, output = self.run_script(
            ("ambiguous", ambiguous_trimmed),
            raw=ambiguous_raw,
            qc=ambiguous_qc,
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("ambiguous raw-column map", completed.stderr)
        self.assertIn("refusing to guess", completed.stderr)
        self.assertFalse(output.exists())

    def test_domain_coordinates_and_ungapped_lengths_fail_closed(self) -> None:
        invalid_qc = self.root / "invalid-qc.tsv"
        invalid_qc.write_text(
            "accession\tsubject_length\tring_coordinates\tbrct_coordinates\t"
            "domain_qc_status\n"
            "A\t8\t1-2\t6-9\tpass\n"
            "B\t10\t1-2\t8-10\tpass\n",
            encoding="utf-8",
        )
        completed, output = self.run_script(
            ("balanced", self.balanced), qc=invalid_qc
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("exceeds subject_length", completed.stderr)
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
