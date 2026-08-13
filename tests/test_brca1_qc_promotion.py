"""Regression tests for the BRCA1 post-QC artifact promotion receipt."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = REPOSITORY_ROOT / "examples" / "brca1"
VERIFIER = EXAMPLE / "scripts" / "verify_brca1_qc_promotion.py"
COMMITTED_RECEIPT = EXAMPLE / "report" / "qc_promotion_receipt.json"
VERIFIED_AT = "2026-08-13T12:45:08Z"


def verifier_argv(output: Path, promoted_qc: str = "qc/candidate_qc.tsv") -> list[str]:
    return [
        sys.executable,
        "scripts/verify_brca1_qc_promotion.py",
        "--summarizer",
        "scripts/summarize_brca1_qc.py",
        "--candidate-fasta",
        "inputs/candidates.faa",
        "--candidate-table",
        "inputs/candidates.pre-qc.tsv",
        "--blast-tsv",
        "qc/blastp-human-vs-candidates.tsv",
        "--interpro-tsv",
        "qc/interproscan-pfam-smart.tsv",
        "--promoted-qc",
        promoted_qc,
        "--promoted-candidates",
        "inputs/candidates.tsv",
        "--verified-at-utc",
        VERIFIED_AT,
        "--output",
        str(output),
    ]


def ledger_records() -> list[dict[str, object]]:
    return [
        json.loads(line)
        for line in (EXAMPLE / "report" / "commands.jsonl")
        .read_text(encoding="utf-8")
        .splitlines()
        if line
    ]


class Brca1QcPromotionTests(unittest.TestCase):
    def test_committed_receipt_reproduces_exactly(self) -> None:
        with tempfile.TemporaryDirectory(prefix="brca1-qc-promotion-test-") as temporary:
            output = Path(temporary) / "receipt.json"
            completed = subprocess.run(
                verifier_argv(output),
                cwd=EXAMPLE,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(
                json.loads(output.read_text(encoding="utf-8")),
                json.loads(COMMITTED_RECEIPT.read_text(encoding="utf-8")),
            )

    def test_mismatched_promoted_output_fails_closed(self) -> None:
        # Keep the test file repository-relative because public receipt paths
        # deliberately reject absolute and parent-traversal paths.
        with tempfile.TemporaryDirectory(
            dir=EXAMPLE, prefix=".brca1-qc-promotion-test-"
        ) as temporary:
            temporary_path = Path(temporary)
            altered = temporary_path / "candidate_qc.tsv"
            altered.write_bytes((EXAMPLE / "qc" / "candidate_qc.tsv").read_bytes() + b"\n")
            output = temporary_path / "receipt.json"
            relative_altered = altered.relative_to(EXAMPLE).as_posix()
            completed = subprocess.run(
                verifier_argv(output, promoted_qc=relative_altered),
                cwd=EXAMPLE,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn("not byte-identical", completed.stderr)
            self.assertFalse(output.exists())

    def test_receipt_is_explicitly_post_hoc_not_invented_replay(self) -> None:
        receipt = json.loads(COMMITTED_RECEIPT.read_text(encoding="utf-8"))
        self.assertEqual(receipt["status"], "verified")
        self.assertEqual(
            receipt["record_type"],
            "post-hoc-qc-artifact-promotion-reconciliation",
        )
        self.assertIs(receipt["prospective_replay_claimed"], False)
        self.assertEqual(receipt["historical_copy_utility"], "not-recorded")
        self.assertEqual(receipt["historical_copy_time_utc"], "not-recorded")
        self.assertEqual(len(receipt["promotions"]), 2)
        self.assertTrue(
            all(
                item["verification"] == "byte-identical-deterministic-regeneration"
                for item in receipt["promotions"]
            )
        )

    def test_ledger_orders_dependencies_and_marks_immutable_plan_as_provenance(self) -> None:
        records = ledger_records()
        by_id = {record["id"]: record for record in records}
        order = {record["id"]: index for index, record in enumerate(records)}
        self.assertLess(order["summarize-candidate-qc"], order["verify-qc-output-promotion"])
        self.assertLess(order["verify-qc-output-promotion"], order["inventory-source-proteins"])
        self.assertLess(order["inventory-source-proteins"], order["compile-reference-plan"])
        self.assertLess(order["compile-reference-plan"], order["export-taxonomy-lineages"])

        inventory_argv = by_id["inventory-source-proteins"]["argv"]
        self.assertIn("inputs/candidates.tsv", inventory_argv)
        planner = by_id["compile-reference-plan"]
        self.assertEqual(planner["command_record_type"], "normalized-provenance-argv")
        self.assertIn("<immutable-review-output>", planner["argv"])
        self.assertIn("not replayable", planner["path_policy"])

    def test_verifier_uses_only_standard_library_modules(self) -> None:
        source = VERIFIER.read_text(encoding="utf-8")
        for forbidden in ("requests", "pandas", "Bio.", "urllib.request"):
            self.assertNotIn(forbidden, source)


if __name__ == "__main__":
    unittest.main()
