"""Black-box contract tests for the offline planning MVP.

The fixtures are deliberately synthetic.  These tests must never require network
access or the MAFFT/IQ-TREE executables; the planner is expected to describe
those future commands without executing them.
"""

from __future__ import annotations

import csv
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Iterable


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    REPOSITORY_ROOT
    / "skills"
    / "bio-gene-to-reference-tree"
    / "scripts"
    / "gene_to_tree.py"
)
FIXTURES = Path(__file__).resolve().parent / "fixtures" / "review"
ASSETS = REPOSITORY_ROOT / "skills" / "bio-gene-to-reference-tree" / "assets"
SECRET_SENTINEL = "offline-test-secret-must-not-appear"


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _values_for_key(value: Any, wanted: str) -> Iterable[Any]:
    if isinstance(value, dict):
        for key, child in value.items():
            if key == wanted:
                yield child
            yield from _values_for_key(child, wanted)
    elif isinstance(value, list):
        for child in value:
            yield from _values_for_key(child, wanted)


def _all_strings(value: Any) -> Iterable[str]:
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for key, child in value.items():
            yield str(key)
            yield from _all_strings(child)
    elif isinstance(value, list):
        for child in value:
            yield from _all_strings(child)


def _planned_argv(plan: dict[str, Any]) -> list[list[str]]:
    """Return command argv arrays while rejecting shell-string commands."""

    commands = plan.get("planned_commands", plan.get("commands"))
    if commands is None:
        # A stage-oriented plan is also acceptable, but its command is still an
        # argv array.  This keeps the contract independent of presentation.
        commands = [
            stage
            for stage in plan.get("stages", [])
            if isinstance(stage, dict)
            and ("argv" in stage or "command" in stage)
        ]
    if not isinstance(commands, list):
        raise AssertionError("plan commands must be a JSON array")

    result: list[list[str]] = []
    for entry in commands:
        if isinstance(entry, list):
            argv = entry
        elif isinstance(entry, dict):
            argv = entry.get("argv", entry.get("command"))
        else:
            raise AssertionError(f"invalid planned command entry: {entry!r}")
        if not isinstance(argv, list) or not argv:
            raise AssertionError(f"planned command is not a non-empty argv array: {entry!r}")
        if not all(isinstance(part, str) and part for part in argv):
            raise AssertionError(f"planned argv contains a non-string or empty part: {argv!r}")
        result.append(argv)
    return result


class OfflineWorkflowContractTests(unittest.TestCase):
    maxDiff = None

    def setUp(self) -> None:
        self._temporary_directory = tempfile.TemporaryDirectory(
            prefix="gene tree offline review "
        )
        self.temp_root = Path(self._temporary_directory.name)

    def tearDown(self) -> None:
        self._temporary_directory.cleanup()

    def _run_plan(
        self,
        request: str | Path = "request.json",
        *,
        out: Path | None = None,
        include_offline: bool = True,
        include_dry_run: bool = True,
    ) -> tuple[subprocess.CompletedProcess[str], Path]:
        output = out or self.temp_root / "new output"
        request_path = Path(request)
        if not request_path.is_absolute():
            request_path = FIXTURES / request_path
        command = [
            sys.executable,
            str(SCRIPT),
            "plan",
            "--request",
            str(request_path),
        ]
        if include_offline:
            command.append("--offline")
        if include_dry_run:
            command.append("--dry-run")
        command.extend(["--out", str(output)])

        environment = os.environ.copy()
        environment.update(
            {
                "GENE_TO_TREE_TEST_SECRET": SECRET_SENTINEL,
                "NCBI_API_KEY": SECRET_SENTINEL,
                "ENTREZ_EMAIL": "private-reviewer@example.invalid",
                # A network regression should fail promptly instead of reaching
                # a live database.  Correct offline code never consults these.
                "HTTP_PROXY": "http://127.0.0.1:9",
                "HTTPS_PROXY": "http://127.0.0.1:9",
                "ALL_PROXY": "http://127.0.0.1:9",
                "NO_PROXY": "",
                # The planner must not require or launch external executables in
                # dry-run mode.  Python itself is already addressed absolutely.
                "PATH": "",
            }
        )
        completed = subprocess.run(
            command,
            cwd=self.temp_root,
            env=environment,
            capture_output=True,
            text=True,
            timeout=15,
            check=False,
        )
        return completed, output

    def _successful_plan(
        self, request: str = "request.json", *, out: Path | None = None
    ) -> Path:
        completed, output = self._run_plan(request, out=out)
        self.assertEqual(
            completed.returncode,
            0,
            msg=f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}",
        )
        self.assertTrue(output.is_dir())
        return output

    def _write_v02_request(self, name: str, mutate: Any | None = None) -> Path:
        request = json.loads((ASSETS / "request.example.json").read_text(encoding="utf-8"))
        request["query"]["path"] = str(ASSETS / "query.example.faa")
        request["references"]["candidate_table"] = str(ASSETS / "candidates.example.tsv")
        request["references"]["candidate_fasta"] = str(ASSETS / "candidates.example.faa")
        if mutate is not None:
            mutate(request)
        request_path = self.temp_root / name
        request_path.write_text(json.dumps(request), encoding="utf-8")
        return request_path

    def test_review_fixture_selects_and_rejects_auditable_candidates(self) -> None:
        output = self._successful_plan()
        expected = json.loads((FIXTURES / "expected.json").read_text(encoding="utf-8"))

        selected = _read_tsv(output / "selected_references.tsv")
        rejected = _read_tsv(output / "rejected_references.tsv")
        self.assertTrue(selected, "selected reference audit must not be empty")
        self.assertTrue(rejected, "rejected reference audit must not be empty")
        self.assertIn("accession", selected[0])
        self.assertIn("accession", rejected[0])
        self.assertIn("reason_codes", rejected[0])

        self.assertEqual(
            {row["accession"] for row in selected},
            set(expected["selected_accessions"]),
        )
        rejected_by_accession = {row["accession"]: row for row in rejected}
        for accession, reason_code in expected["rejected_reason_codes"].items():
            self.assertIn(accession, rejected_by_accession)
            actual_codes = {
                code.strip()
                for code in re.split(",|;|\\|", rejected_by_accession[accession]["reason_codes"])
                if code.strip()
            }
            self.assertIn(reason_code, actual_codes)

        reference_headers = [
            line[1:].split()[0]
            for line in (output / "reference_set.faa").read_text(encoding="utf-8").splitlines()
            if line.startswith(">")
        ]
        self.assertEqual(len(reference_headers), len(set(reference_headers)))
        self.assertEqual(
            set(reference_headers),
            {"QUERY_001", *expected["selected_accessions"]},
        )

    def test_dry_run_writes_a_plan_but_no_alignment_or_tree(self) -> None:
        output = self._successful_plan()
        required = {
            "selected_references.tsv",
            "rejected_references.tsv",
            "reference_set.faa",
            "sequence_metadata.tsv",
            "itol_roles.txt",
            "plan.json",
            "manifest.json",
        }
        self.assertEqual(required, {path.name for path in output.iterdir()})

        forbidden_suffixes = {
            ".aln",
            ".aln-fasta",
            ".phy",
            ".nwk",
            ".tree",
            ".treefile",
            ".contree",
            ".iqtree",
            ".log",
        }
        unexpected = [
            path.relative_to(output).as_posix()
            for path in output.rglob("*")
            if path.is_file()
            and (
                path.suffix.lower() in forbidden_suffixes
                or "alignment" in path.name.lower()
                or "treefile" in path.name.lower()
            )
        ]
        self.assertEqual(unexpected, [])

    def test_plan_commands_are_argument_arrays(self) -> None:
        output = self._successful_plan()
        plan = json.loads((output / "plan.json").read_text(encoding="utf-8"))
        commands = _planned_argv(plan)
        self.assertGreaterEqual(len(commands), 2)
        executed_values = list(_values_for_key(plan, "executed"))
        self.assertTrue(executed_values and all(value is False for value in executed_values))

        executables = {Path(argv[0]).name.lower() for argv in commands}
        self.assertTrue(
            {"mafft", "iqtree2"}.issubset(executables),
            msg=f"planned executables were {sorted(executables)!r}",
        )
        for argv in commands:
            # Redirection belongs in structured stdout/stderr fields, never in
            # an argv token interpreted by a shell.
            self.assertFalse(any(token in {">", ">>", "|", ";", "&&"} for token in argv))

    def test_manifest_is_portable_and_does_not_capture_secrets(self) -> None:
        output = self._successful_plan()
        manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
        rendered = json.dumps(manifest, sort_keys=True)

        self.assertNotIn(SECRET_SENTINEL, rendered)
        self.assertNotIn("private-reviewer@example.invalid", rendered)
        self.assertNotIn(str(REPOSITORY_ROOT.resolve()), rendered)
        self.assertNotIn(str(self.temp_root.resolve()), rendered)

        absolute_paths = [
            item
            for item in _all_strings(manifest)
            if Path(item).is_absolute() or re.match(r"^[A-Za-z]:[\\\\/]", item)
        ]
        self.assertEqual(absolute_paths, [])

        offline_values = list(_values_for_key(manifest, "offline"))
        self.assertTrue(offline_values and all(value is True for value in offline_values))
        executed_values = list(_values_for_key(manifest, "executed"))
        self.assertTrue(executed_values and all(value is False for value in executed_values))
        self.assertIsInstance(manifest.get("run_id"), str)
        self.assertTrue(manifest["run_id"])

    def test_live_or_execute_mode_is_not_available_in_the_mvp(self) -> None:
        for offline, dry_run in ((False, True), (True, False), (False, False)):
            with self.subTest(offline=offline, dry_run=dry_run):
                output = self.temp_root / f"live-{offline}-{dry_run}"
                completed, _ = self._run_plan(
                    out=output,
                    include_offline=offline,
                    include_dry_run=dry_run,
                )
                self.assertNotEqual(completed.returncode, 0)
                self.assertFalse(output.exists())

    def test_existing_output_directory_is_refused_without_modification(self) -> None:
        output = self.temp_root / "already exists"
        output.mkdir()
        sentinel = output / "reviewer-note.txt"
        sentinel.write_text("keep me\n", encoding="utf-8")

        completed, _ = self._run_plan(out=output)
        self.assertNotEqual(completed.returncode, 0)
        self.assertEqual(sentinel.read_text(encoding="utf-8"), "keep me\n")
        self.assertEqual({path.name for path in output.iterdir()}, {sentinel.name})

    def test_gene_symbol_cannot_bypass_the_offline_input_contract(self) -> None:
        request = json.loads((FIXTURES / "request.json").read_text(encoding="utf-8"))
        request["query"] = {"kind": "gene-symbol", "symbol": "BRCA1"}
        request_path = self.temp_root / "symbol request.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")

        completed, output = self._run_plan(request_path)

        self.assertEqual(completed.returncode, 2)
        self.assertIn("QUERY_RESOLUTION_REQUIRED", completed.stderr)
        self.assertFalse(output.exists())

    def test_insufficient_taxa_produce_a_blocked_audit_bundle(self) -> None:
        request = json.loads((FIXTURES / "request.json").read_text(encoding="utf-8"))
        request["query"]["path"] = str(FIXTURES / "query.faa")
        request["references"]["candidate_table"] = str(FIXTURES / "candidates.tsv")
        request["references"]["candidate_fasta"] = str(FIXTURES / "candidates.faa")
        request["selection"]["min_ingroup_taxa"] = 3
        request["selection"]["max_references"] = 4
        request_path = self.temp_root / "insufficient request.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")

        completed, output = self._run_plan(request_path, out=self.temp_root / "blocked")

        self.assertEqual(completed.returncode, 3)
        self.assertTrue(output.is_dir())
        plan = json.loads((output / "plan.json").read_text(encoding="utf-8"))
        self.assertEqual(plan["state"], "blocked")
        self.assertIn("INSUFFICIENT_INGROUP_TAXA", plan["hard_stops"])
        self.assertTrue(all(command["executed"] is False for command in plan["planned_commands"]))
        manifest_text = (output / "manifest.json").read_text(encoding="utf-8")
        self.assertNotIn(str(FIXTURES), manifest_text)

    def test_semantically_identical_candidate_order_has_stable_run_id(self) -> None:
        first = self._successful_plan("request.json", out=self.temp_root / "ordered")
        second = self._successful_plan(
            "request-shuffled.json", out=self.temp_root / "shuffled"
        )
        first_manifest = json.loads((first / "manifest.json").read_text(encoding="utf-8"))
        second_manifest = json.loads((second / "manifest.json").read_text(encoding="utf-8"))

        self.assertEqual(first_manifest.get("run_id"), second_manifest.get("run_id"))
        first_plan = json.loads((first / "plan.json").read_text(encoding="utf-8"))
        second_plan = json.loads((second / "plan.json").read_text(encoding="utf-8"))
        self.assertEqual(first_plan["plan_hash"], second_plan["plan_hash"])
        first_selected = {row["accession"] for row in _read_tsv(first / "selected_references.tsv")}
        second_selected = {row["accession"] for row in _read_tsv(second / "selected_references.tsv")}
        self.assertEqual(first_selected, second_selected)

    def test_v02_example_generates_itol_and_complete_metadata(self) -> None:
        request_path = self._write_v02_request("v02 request.json")
        output = self._successful_plan(request_path, out=self.temp_root / "v02")

        plan = json.loads((output / "plan.json").read_text(encoding="utf-8"))
        self.assertEqual(plan["schema_version"], "0.2")
        self.assertEqual(plan["state"], "pending-reference-approval")
        self.assertEqual(
            [gate["id"] for gate in plan["decision_gates"]],
            ["reference-and-outgroup", "alignment-and-trimming-qc"],
        )

        metadata = _read_tsv(output / "sequence_metadata.tsv")
        candidates = _read_tsv(ASSETS / "candidates.example.tsv")
        self.assertEqual(len(metadata), len(candidates))
        self.assertEqual(len({row["accession"] for row in metadata}), len(metadata))
        by_accession = {row["accession"]: row for row in metadata}
        self.assertEqual(by_accession["QUERY_001"]["analysis_role"], "study")
        self.assertEqual(by_accession["MOUSE_CAN"]["analysis_role"], "expanded")
        self.assertEqual(by_accession["CIONA_OUT"]["analysis_role"], "outgroup")
        self.assertEqual(by_accession["FISH_PARALOG"]["inclusion_status"], "rejected")

        itol = (output / "itol_roles.txt").read_text(encoding="utf-8")
        self.assertTrue(itol.startswith("DATASET_COLORSTRIP\nSEPARATOR TAB\n"))
        self.assertIn("LEGEND_COLORS\t#E69F00\t#009E73\t#999999", itol)
        selected_ids = {row["accession"] for row in metadata if row["inclusion_status"] == "selected"}
        data_lines = itol.split("\nDATA\n", 1)[1].strip().splitlines()
        self.assertEqual({line.split("\t", 1)[0] for line in data_lines}, selected_ids)
        self.assertNotIn("FISH_PARALOG\t", itol)

        manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
        recorded_outputs = {item["logical_path"] for item in manifest["output_artifacts"]}
        self.assertEqual(recorded_outputs, {path.name for path in output.iterdir()} - {"manifest.json"})

    def test_v02_commands_use_trimal_and_support_aware_iqtree(self) -> None:
        request_path = self._write_v02_request("commands request.json")
        output = self._successful_plan(request_path, out=self.temp_root / "commands")
        plan = json.loads((output / "plan.json").read_text(encoding="utf-8"))
        commands = {command["id"]: command for command in plan["planned_commands"]}

        self.assertIn("--thread", commands["align-proteins"]["argv"])
        self.assertEqual(commands["trim-balanced"]["tool"], "trimal")
        self.assertEqual(commands["trim-balanced"]["argv"][-2:], ["-gt", "0.9"])
        iqtree = commands["infer-accurate-tree"]["argv"]
        self.assertIn("-B", iqtree)
        self.assertIn("-bnni", iqtree)
        self.assertIn("-alrt", iqtree)
        self.assertNotIn("-b", iqtree)
        self.assertNotIn("clipkit", json.dumps(plan).lower())

    def test_mmseqs_auto_trigger_blocks_and_preserves_study_outgroup(self) -> None:
        def trigger(request: dict[str, Any]) -> None:
            request["clustering"]["trigger_min_sequences"] = 2

        request_path = self._write_v02_request("cluster request.json", trigger)
        completed, output = self._run_plan(request_path, out=self.temp_root / "cluster")
        self.assertEqual(completed.returncode, 3)
        self.assertTrue((output / "expanded_candidates.faa").is_file())
        expanded_text = (output / "expanded_candidates.faa").read_text(encoding="utf-8")
        self.assertNotIn(">QUERY_001", expanded_text)
        self.assertNotIn(">CIONA_OUT", expanded_text)

        plan = json.loads((output / "plan.json").read_text(encoding="utf-8"))
        self.assertEqual(plan["state"], "pending-clustering")
        self.assertIn("MMSEQS_CLUSTERING_REQUIRED", plan["hard_stops"])
        command = plan["clustering_plan"]["command"]["argv"]
        self.assertIn("--min-seq-id", command)
        self.assertIn("-c", command)
        self.assertIn("--cov-mode", command)
        self.assertEqual(plan["clustering_plan"]["protected_analysis_groups"], ["study", "outgroup"])

    def test_fasttree_support_is_not_called_bootstrap(self) -> None:
        def fast(request: dict[str, Any]) -> None:
            request["trimming"]["enabled"] = False
            request["tree"] = {
                "mode": "fast",
                "tool": "fasttree",
                "model": "WAG",
                "support": {
                    "method": "sh-like-local",
                    "replicates": 0,
                    "sh_alrt": 0,
                    "bnni": False
                },
                "threads": 4,
                "seed": 12345,
                "rooting": "outgroup"
            }

        request_path = self._write_v02_request("fast request.json", fast)
        output = self._successful_plan(request_path, out=self.temp_root / "fast")
        plan = json.loads((output / "plan.json").read_text(encoding="utf-8"))
        command = next(item for item in plan["planned_commands"] if item["tool"] == "fasttree")
        self.assertIn("SH-like local support", command["support_semantics"])
        self.assertIn("not a global bootstrap", command["support_semantics"].lower())
        self.assertFalse(any(token in {"-B", "-b"} for token in command["argv"]))

    def test_standard_bootstrap_uses_lowercase_b_without_bnni(self) -> None:
        def standard(request: dict[str, Any]) -> None:
            request["tree"]["support"] = {
                "method": "standard-bootstrap",
                "replicates": 1000,
                "sh_alrt": 1000,
                "bnni": False,
            }

        request_path = self._write_v02_request("standard request.json", standard)
        output = self._successful_plan(request_path, out=self.temp_root / "standard")
        plan = json.loads((output / "plan.json").read_text(encoding="utf-8"))
        command = next(item for item in plan["planned_commands"] if item["tool"] == "iqtree2")
        self.assertIn("-b", command["argv"])
        self.assertNotIn("-B", command["argv"])
        self.assertNotIn("-bnni", command["argv"])

    def test_resolved_accession_route_preserves_original_input(self) -> None:
        def accession(request: dict[str, Any]) -> None:
            request["query"]["kind"] = "accession"
            request["query"]["original_value"] = "QUERY_001.1"

        request_path = self._write_v02_request("accession request.json", accession)
        output = self._successful_plan(request_path, out=self.temp_root / "accession")
        plan = json.loads((output / "plan.json").read_text(encoding="utf-8"))
        self.assertEqual(plan["query"]["original_kind"], "accession")
        self.assertEqual(plan["query"]["original_value"], "QUERY_001.1")
        self.assertEqual(plan["query_resolution"]["status"], "materialized-local-record")
        self.assertFalse(plan["query_resolution"]["live_lookup_performed_by_planner"])

    def test_precomputed_mmseqs_clusters_are_applied_deterministically(self) -> None:
        rows = _read_tsv(ASSETS / "candidates.example.tsv")
        cluster_ids = {
            "MOUSE_CAN": "cluster-1",
            "MOUSE_ALT": "cluster-1",
            "CHICKEN_OK": "cluster-2",
            "FROG_FRAGMENT": "cluster-3",
            "FISH_PARALOG": "cluster-4",
        }
        table_path = self.temp_root / "clustered candidates.tsv"
        fieldnames = list(rows[0])
        if "cluster_id" not in fieldnames:
            fieldnames.append("cluster_id")
        with table_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                row["cluster_id"] = cluster_ids.get(row["accession"], "")
                writer.writerow(row)

        def precomputed(request: dict[str, Any]) -> None:
            request["references"]["candidate_table"] = str(table_path)
            request["clustering"]["mode"] = "precomputed"

        request_path = self._write_v02_request("precomputed request.json", precomputed)
        output = self._successful_plan(request_path, out=self.temp_root / "precomputed")
        rejected = {row["accession"]: row["reason_codes"] for row in _read_tsv(output / "rejected_references.tsv")}
        self.assertEqual(rejected["MOUSE_ALT"], "MMSEQS_CLUSTER_REDUNDANT")
        plan = json.loads((output / "plan.json").read_text(encoding="utf-8"))
        self.assertEqual(plan["clustering_plan"]["status"], "precomputed")
        self.assertFalse(plan["clustering_plan"]["replan_required_after_execution"])

    def test_decision_bearing_itol_color_changes_plan_hash(self) -> None:
        first_request = self._write_v02_request("first color.json")
        first = self._successful_plan(first_request, out=self.temp_root / "first-color")

        def recolor(request: dict[str, Any]) -> None:
            request["itol"]["colors"]["study"] = "#F0A000"

        second_request = self._write_v02_request("second color.json", recolor)
        second = self._successful_plan(second_request, out=self.temp_root / "second-color")
        first_plan = json.loads((first / "plan.json").read_text(encoding="utf-8"))
        second_plan = json.loads((second / "plan.json").read_text(encoding="utf-8"))
        self.assertNotEqual(first_plan["plan_hash"], second_plan["plan_hash"])

    def test_additional_study_sequence_is_preserved_outside_reference_quotas(self) -> None:
        rows = _read_tsv(ASSETS / "candidates.example.tsv")
        table_path = self.temp_root / "multiple study.tsv"
        with table_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for row in rows:
                if row["accession"] == "MOUSE_ALT":
                    row["analysis_group"] = "study"
                writer.writerow(row)

        def multiple_study(request: dict[str, Any]) -> None:
            request["references"]["candidate_table"] = str(table_path)

        request_path = self._write_v02_request("multiple study request.json", multiple_study)
        output = self._successful_plan(request_path, out=self.temp_root / "multiple-study")
        selected = {row["accession"]: row for row in _read_tsv(output / "selected_references.tsv")}
        self.assertIn("MOUSE_ALT", selected)
        self.assertIn("MOUSE_CAN", selected)
        self.assertEqual(selected["MOUSE_ALT"]["decision_reason"], "SELECTED_STUDY")
        itol = (output / "itol_roles.txt").read_text(encoding="utf-8")
        self.assertIn("MOUSE_ALT\t#E69F00\tStudy", itol)
        plan = json.loads((output / "plan.json").read_text(encoding="utf-8"))
        self.assertEqual(plan["candidate_and_selection_counts"]["selected_references_excluding_study"], 3)

    def test_schema_files_are_valid_json(self) -> None:
        references = REPOSITORY_ROOT / "skills" / "bio-gene-to-reference-tree" / "references"
        request_schema = json.loads((references / "request-0.2.schema.json").read_text(encoding="utf-8"))
        plan_schema = json.loads((references / "plan-0.2.schema.json").read_text(encoding="utf-8"))
        self.assertEqual(request_schema["properties"]["schema_version"]["const"], "0.2")
        self.assertEqual(plan_schema["properties"]["schema_version"]["const"], "0.2")


if __name__ == "__main__":
    unittest.main()
