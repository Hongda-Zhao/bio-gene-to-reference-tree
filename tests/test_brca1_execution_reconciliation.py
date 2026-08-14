"""Integrity contract for the BRCA1 post-hoc execution reconciliation."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import unittest
from datetime import datetime
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = REPOSITORY_ROOT / "examples" / "brca1"
RECEIPT_PATH = EXAMPLE / "report" / "execution_reconciliation.json"
SCHEMA_PATH = EXAMPLE / "report" / "execution_reconciliation.schema.json"
PLAN_PATH = EXAMPLE / "review" / "plan.json"
LEDGER_PATH = EXAMPLE / "report" / "commands.jsonl"

COMMAND_PAIRS = (
    ("align-proteins", "align-einsi"),
    ("trim-permissive", "trim-permissive"),
    ("trim-balanced", "trim-balanced"),
    ("trim-strict", "trim-strict"),
    ("infer-accurate-tree", "infer-primary-tree"),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_argv_sha256(argv: list[str]) -> str:
    payload = json.dumps(
        argv,
        ensure_ascii=False,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def load_receipt() -> dict[str, object]:
    return json.loads(RECEIPT_PATH.read_text(encoding="utf-8"))


def load_ledger() -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    for raw in LEDGER_PATH.read_text(encoding="utf-8").splitlines():
        record = json.loads(raw)
        command_id = record["id"]
        if command_id in records:
            raise AssertionError(f"duplicate command id: {command_id}")
        records[command_id] = record
    return records


def load_key_value_receipt(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if any(len(row) != 2 for row in rows):
        raise AssertionError(f"invalid key/value receipt: {path}")
    return dict(rows)


class Brca1ExecutionReconciliationTests(unittest.TestCase):
    def test_schema_declares_nonretroactive_pending_and_final_states(self) -> None:
        schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
        self.assertEqual(schema["$schema"], "https://json-schema.org/draft/2020-12/schema")
        required = set(schema["required"])
        self.assertTrue(
            {
                "status",
                "prospective_command_approval",
                "scientific_equivalence_claimed",
                "plan_binding",
                "ledger_binding",
                "non_promoted_attempts",
                "command_reconciliations",
            }.issubset(required)
        )
        properties = schema["properties"]
        self.assertEqual(properties["prospective_command_approval"]["const"], False)
        self.assertEqual(properties["scientific_equivalence_claimed"]["const"], False)
        self.assertEqual(
            set(properties["status"]["enum"]),
            {"pending-final-iqtree-artifacts", "reconciled-post-hoc"},
        )

        receipt = load_receipt()
        self.assertEqual(receipt["schema_version"], "1.0")
        self.assertEqual(
            receipt["schema_path"], "report/execution_reconciliation.schema.json"
        )
        self.assertFalse(receipt["prospective_command_approval"])
        self.assertFalse(receipt["scientific_equivalence_claimed"])

    def test_plan_ledger_and_materialized_artifacts_are_hash_bound(self) -> None:
        receipt = load_receipt()
        plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
        plan_binding = receipt["plan_binding"]
        ledger_binding = receipt["ledger_binding"]

        self.assertEqual(plan_binding["path"], "review/plan.json")
        self.assertEqual(plan_binding["sha256"], sha256(PLAN_PATH))
        self.assertEqual(plan_binding["run_id"], plan["run_id"])
        self.assertEqual(plan_binding["plan_hash"], plan["plan_hash"])
        self.assertEqual(ledger_binding["path"], "report/commands.jsonl")
        self.assertEqual(ledger_binding["sha256"], sha256(LEDGER_PATH))

        approval = load_key_value_receipt(
            EXAMPLE / str(plan_binding["approval_receipt_path"])
        )
        self.assertEqual(approval["run_id"], plan_binding["run_id"])
        self.assertEqual(approval["plan_hash"], plan_binding["plan_hash"])
        self.assertEqual(approval["approved_at_utc"], plan_binding["approved_at_utc"])

        for name, binding in receipt["artifact_bindings"].items():
            with self.subTest(artifact=name):
                path = EXAMPLE / binding["path"]
                self.assertTrue(path.is_file(), path)
                self.assertEqual(binding["sha256"], sha256(path))

    def test_receipt_copies_exact_plan_and_ledger_argv_with_recomputed_hashes(self) -> None:
        receipt = load_receipt()
        plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
        planned = {record["id"]: record for record in plan["planned_commands"]}
        executed = load_ledger()
        reconciliations = {
            (record["planned_command_id"], record["executed_command_id"]): record
            for record in receipt["command_reconciliations"]
        }
        self.assertEqual(set(reconciliations), set(COMMAND_PAIRS))

        for planned_id, executed_id in COMMAND_PAIRS:
            with self.subTest(planned=planned_id, executed=executed_id):
                record = reconciliations[(planned_id, executed_id)]
                self.assertEqual(record["planned_argv"], planned[planned_id]["argv"])
                self.assertEqual(record["executed_argv"], executed[executed_id]["argv"])
                self.assertEqual(
                    record["planned_argv_sha256"],
                    canonical_argv_sha256(record["planned_argv"]),
                )
                self.assertEqual(
                    record["executed_argv_sha256"],
                    canonical_argv_sha256(record["executed_argv"]),
                )
                self.assertNotEqual(record["planned_argv"], record["executed_argv"])
                self.assertFalse(record["exact_argv_match"])
                self.assertFalse(record["prospectively_approved"])

                for binding_name in ("input_binding", "output_binding"):
                    binding = record[binding_name]
                    if binding["sha256"] is None:
                        continue
                    path = EXAMPLE / binding["path"]
                    self.assertTrue(path.is_file(), path)
                    self.assertEqual(binding["sha256"], sha256(path))

    def test_every_exact_difference_has_honest_scientific_classification(self) -> None:
        receipt = load_receipt()
        records = {record["stage"]: record for record in receipt["command_reconciliations"]}
        self.assertEqual(
            set(records),
            {"alignment", "trim-permissive", "trim-balanced", "trim-strict", "final-tree"},
        )

        mafft = records["alignment"]
        self.assertFalse(mafft["potentially_result_affecting"])
        self.assertEqual(
            [(item["field"], item["classification"]) for item in mafft["differences"]],
            [("input_path", "path-only-hash-bound")],
        )

        for stage, threshold in (
            ("trim-permissive", 0.1),
            ("trim-balanced", 0.5),
            ("trim-strict", 0.9),
        ):
            with self.subTest(stage=stage):
                record = records[stage]
                self.assertFalse(record["potentially_result_affecting"])
                self.assertEqual(record["scientific_parameters_unchanged"], {"gap_threshold": threshold})
                differences = {item["field"]: item for item in record["differences"]}
                self.assertEqual(set(differences), {"input_path", "output_path", "output_format"})
                self.assertEqual(differences["input_path"]["classification"], "path-only-hash-bound")
                self.assertEqual(differences["output_path"]["classification"], "output-naming-only")
                self.assertEqual(differences["output_format"]["planned"], None)
                self.assertEqual(differences["output_format"]["executed"], "-fasta")
                self.assertEqual(differences["output_format"]["classification"], "output-format-only")
                self.assertTrue(
                    all(not item["potentially_result_affecting"] for item in differences.values())
                )

        iqtree = records["final-tree"]
        self.assertFalse(iqtree["potentially_result_affecting"])
        self.assertEqual(
            iqtree["scientific_parameters_unchanged"],
            {
                "model_selection": "MFP",
                "ufboot2_replicates": 1000,
                "ufboot_bnni": True,
                "sh_alrt_replicates": 1000,
                "threads": 8,
                "seed": 20260813,
            },
        )
        differences = {item["field"]: item for item in iqtree["differences"]}
        self.assertEqual(set(differences), {"input_path", "output_prefix"})
        self.assertEqual(differences["input_path"]["classification"], "path-only-hash-bound")
        self.assertEqual(differences["output_prefix"]["classification"], "output-naming-only")
        self.assertTrue(
            all(not item["potentially_result_affecting"] for item in differences.values())
        )

    def test_pending_is_local_only_and_final_receipt_is_closed(self) -> None:
        receipt = load_receipt()
        status = receipt["status"]
        self.assertIn(status, {"pending-final-iqtree-artifacts", "reconciled-post-hoc"})
        if status == "pending-final-iqtree-artifacts":
            self.assertEqual(receipt["publication_disposition"], "pending")
            self.assertIsNone(receipt["reconciled_at_utc"])
            self.assertIsNone(receipt["reviewer_acknowledgement"])
            self.assertIsNone(receipt["actual_iqtree"]["reported_version"])
            for binding in receipt["actual_iqtree"].values():
                if isinstance(binding, dict):
                    self.assertIsNone(binding["sha256"])
            final_tree = next(
                item
                for item in receipt["command_reconciliations"]
                if item["stage"] == "final-tree"
            )
            self.assertIsNone(final_tree["output_binding"]["sha256"])
            if os.environ.get("CI"):
                self.fail(
                    "execution_reconciliation.json cannot remain pending in a public CI run"
                )
            return

        self.assertEqual(
            receipt["publication_disposition"],
            "accepted-with-post-hoc-protocol-deviation",
        )
        self.assertFalse(receipt["prospective_command_approval"])
        self.assertFalse(receipt["scientific_equivalence_claimed"])
        acknowledgement = receipt["reviewer_acknowledgement"]
        self.assertIsInstance(acknowledgement, str)
        self.assertIn("does not retroactively authorize", acknowledgement.lower())

        started = datetime.fromisoformat(
            receipt["executed_inference_started_at_utc"].replace("Z", "+00:00")
        )
        reconciled = datetime.fromisoformat(receipt["reconciled_at_utc"].replace("Z", "+00:00"))
        self.assertGreater(reconciled, started)

        actual = receipt["actual_iqtree"]
        self.assertIn("2.4.0", actual["reported_version"])
        self.assertEqual(actual["started_at_utc"], "2026-08-13T13:10:04Z")
        self.assertEqual(actual["exit_code"], 0)
        self.assertEqual(actual["wall_seconds"], 99)
        self.assertEqual(
            actual["input_alignment_sha256"],
            receipt["artifact_bindings"]["balanced_alignment"]["sha256"],
        )
        for name in (
            "version_receipt",
            "public_log",
            "public_treefile",
            "redactions_receipt",
            "raw_output_receipt",
        ):
            with self.subTest(final_binding=name):
                binding = actual[name]
                self.assertRegex(binding["sha256"], r"^[0-9a-f]{64}$")
                path = EXAMPLE / binding["path"]
                self.assertTrue(path.is_file(), path)
                self.assertEqual(binding["sha256"], sha256(path))

        version_text = (EXAMPLE / actual["version_receipt"]["path"]).read_text(
            encoding="utf-8"
        )
        self.assertIn("2.4.0", version_text)
        final_tree = next(
            item for item in receipt["command_reconciliations"] if item["stage"] == "final-tree"
        )
        self.assertEqual(
            final_tree["output_binding"]["sha256"],
            actual["public_treefile"]["sha256"],
        )
        ledger = load_ledger()
        self.assertIs(ledger["infer-primary-tree"]["executed"], True)
        self.assertEqual(ledger["infer-primary-tree"]["status"], "completed")
        self.assertEqual(ledger["infer-primary-tree"]["tool"], "iqtree2")
        self.assertEqual(ledger["infer-primary-tree"]["exit_code"], 0)

    def test_canceled_iqtree3_attempt_is_audited_but_never_promoted(self) -> None:
        receipt = load_receipt()
        attempts = receipt["non_promoted_attempts"]
        self.assertEqual(len(attempts), 1)
        attempt = attempts[0]
        self.assertEqual(attempt["tool"], "iqtree3")
        self.assertEqual(attempt["reported_module_version"], "3.1.3")
        self.assertEqual(attempt["started_at_utc"], "2026-08-13T09:12:23Z")
        self.assertEqual(attempt["exit_code"], 271)
        self.assertEqual(attempt["outcome"], "canceled")
        self.assertIs(attempt["promoted"], False)
        self.assertEqual(attempt["promoted_output_count"], 0)

        ledger = load_ledger()
        ledger_attempt = ledger[attempt["command_id"]]
        self.assertEqual(ledger_attempt["tool"], "iqtree3")
        self.assertEqual(ledger_attempt["exit_code"], 271)
        self.assertIs(ledger_attempt["promoted"], False)
        self.assertIn("non-promoted", ledger_attempt["execution_outcome"])

    def test_public_narrative_discloses_deviation_without_equivalence_claim(self) -> None:
        documents = (
            REPOSITORY_ROOT / "README.md",
            EXAMPLE / "README.md",
            EXAMPLE / "PROVENANCE.md",
            EXAMPLE / "report" / "report.md",
        )
        for path in documents:
            with self.subTest(document=path.relative_to(REPOSITORY_ROOT)):
                text = path.read_text(encoding="utf-8").lower()
                self.assertIn("post-hoc", text)
                self.assertIn("iqtree2", text)
                self.assertIn("iqtree3", text)
                self.assertIn("scientific equivalence", text)
                self.assertTrue(
                    "not retroactive" in text or "does not retroactively" in text,
                    path,
                )


if __name__ == "__main__":
    unittest.main()
