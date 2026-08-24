"""Black-box tests for capability snapshots and deterministic environment routing."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path
from typing import Any


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SKILL_ROOT = REPOSITORY_ROOT / "skills" / "bio-gene-to-reference-tree"
SCRIPT = SKILL_ROOT / "scripts" / "gene_to_tree.py"
TOOL_NAMES = (
    "ncbi-datasets",
    "blastp",
    "blastn",
    "blastx",
    "tblastn",
    "tblastx",
    "interproscan",
    "jackhmmer",
    "hhsearch",
    "foldseek",
    "mmseqs2",
    "mafft",
    "mafft-qinsi",
    "trimal",
    "pal2nal",
    "macse",
    "fasttree",
    "iqtree2",
    "iqtree3",
    "iqtree",
    "rscript",
)
R_PACKAGES = ("ape", "ggplot2", "ggtree", "openssl", "svglite")


def _probe(status: str = "missing", version: str | None = None) -> dict[str, Any]:
    return {
        "status": status,
        "executable": "tool" if status == "available" else None,
        "version": version,
    }


def _package(status: str = "missing") -> dict[str, Any]:
    return {"status": status, "version": "1.0" if status == "available" else None}


def _snapshot(
    *available_tools: str,
    kind: str = "local",
    available_packages: tuple[str, ...] = (),
    available_executors: tuple[str, ...] = (),
) -> dict[str, Any]:
    return {
        "schema_version": "0.1",
        "workflow_version": "0.4.0",
        "environment_id": f"test-{kind}",
        "environment_kind": kind,
        "probe_scope": "current-process-environment",
        "probe_mode": "version-command",
        "network": "not-probed",
        "environment": {
            "operating_system": "TestOS",
            "architecture": "test-arch",
            "python_version": "3.12.0",
            "cpu_count": 8,
            "memory_gb": 32.0,
            "scratch_free_gb": 100.0,
        },
        "tools": {
            name: _probe("available", "1.0")
            if name in available_tools
            else _probe()
            for name in TOOL_NAMES
        },
        "r_packages": {
            name: _package("available") if name in available_packages else _package()
            for name in R_PACKAGES
        },
        "executors": {
            name: _probe("available", "1.0")
            if name in available_executors
            else _probe()
            for name in ("ssh", "slurm", "pbs", "lsf")
        },
    }


def _profile() -> dict[str, Any]:
    return {
        "schema_version": "0.1",
        "profile_id": "test-profile",
        "intent": "accurate",
        "inputs": {
            "query_kind": "accession-or-name",
            "query_resolved": True,
            "query_location": "compute-target",
            "candidates_location": "compute-target",
            "cluster_mapping_location": "absent",
            "alignment_location": "absent",
            "tree_location": "absent",
            "tree_root_state": "none",
            "metadata_location": "absent",
            "candidate_count": 20,
            "taxdump_location": "compute-target",
            "sequence_database_location": "absent",
            "sequence_database_format": "none",
            "cached_literature": False,
        },
        "capabilities": {
            "compute_shell": True,
            "host_shell": True,
            "database_lookup": False,
            "literature_search": False,
            "file_transfer": False,
            "compute_tree_io": False,
            "host_tree_io": False,
            "web_upload": False,
        },
        "permissions": {
            "public_database_network": False,
            "remote_sequence_submission": False,
            "literature_network": False,
            "itol_upload": False,
        },
        "compute": {
            "target": "local",
            "scheduler": "none",
            "threads": 4,
            "memory_gb": 16,
        },
        "requirements": {
            "trimming": True,
            "exact_taxonomy": True,
            "rooted_tree": False,
            "visualization": "none",
            "current_literature": False,
            "report": False,
            "clustering_trigger": 200,
        },
    }


def _profile_v02(
    molecule: str = "protein",
    *,
    analysis_kind: str | None = None,
    coding_status: str | None = None,
    genetic_code: int | None = None,
) -> dict[str, Any]:
    """Return a current profile while preserving the legacy fixture separately."""
    profile = _profile()
    profile["schema_version"] = "0.2"
    if molecule == "protein":
        resolved_analysis = analysis_kind or "protein"
        resolved_coding_status = coding_status or "not-applicable"
        alignment_strategy = "mafft-protein"
        trimming_strategy = "trimal-columns"
    elif molecule in {"noncoding-dna", "noncoding-rna"}:
        resolved_analysis = analysis_kind or "nucleotide"
        resolved_coding_status = coding_status or "not-applicable"
        alignment_strategy = "mafft-nucleotide"
        trimming_strategy = "trimal-columns"
    elif molecule in {"coding-dna", "coding-rna"}:
        resolved_analysis = analysis_kind or "codon"
        resolved_coding_status = coding_status or "clean"
        genetic_code = 1 if genetic_code is None else genetic_code
        alignment_strategy = "translate-mafft-backtranslate"
        trimming_strategy = "protein-mask-to-codons"
    else:
        resolved_analysis = analysis_kind or "auto"
        resolved_coding_status = coding_status or "unknown"
        alignment_strategy = "auto"
        trimming_strategy = "none"
        profile["requirements"]["trimming"] = False

    profile["inputs"].update(
        {
            "query_molecule": molecule,
            "candidates_molecule": molecule,
            "candidate_bundle_kind": (
                "protein-fasta-metadata"
                if molecule == "protein"
                else (
                    "noncoding-nucleotide-fasta-metadata"
                    if molecule in {"noncoding-dna", "noncoding-rna"}
                    else (
                        "clean-cds-translations-metadata"
                        if molecule in {"coding-dna", "coding-rna"}
                        and resolved_coding_status == "clean"
                        else (
                            "disrupted-cds-metadata"
                            if molecule in {"coding-dna", "coding-rna"}
                            and resolved_coding_status == "disrupted"
                            else "unknown"
                        )
                    )
                )
            ),
            "alignment_kind": "absent",
            "tree_data_kind": "absent",
            "sequence_database_molecule": "none",
        }
    )
    profile["requirements"].update(
        {
            "analysis_kind": resolved_analysis,
            "genetic_code": genetic_code,
            "coding_status": resolved_coding_status,
            "alignment_strategy": alignment_strategy,
            "trimming_strategy": trimming_strategy,
            "rna_structure_aware": False,
        }
    )
    return profile


class EnvironmentRoutingTests(unittest.TestCase):
    maxDiff = None

    def setUp(self) -> None:
        self._temporary_directory = tempfile.TemporaryDirectory(prefix="environment route ")
        self.temp_root = Path(self._temporary_directory.name)

    def tearDown(self) -> None:
        self._temporary_directory.cleanup()

    def _route(
        self,
        profile: dict[str, Any],
        snapshot: dict[str, Any],
        *,
        host_snapshot: dict[str, Any] | None = None,
        output: Path | None = None,
    ) -> tuple[subprocess.CompletedProcess[str], dict[str, Any] | None]:
        profile_path = self.temp_root / "profile.json"
        snapshot_path = self.temp_root / "snapshot.json"
        profile_path.write_text(json.dumps(profile), encoding="utf-8")
        snapshot_path.write_text(json.dumps(snapshot), encoding="utf-8")
        command = [
            sys.executable,
            str(SCRIPT),
            "route",
            "--profile",
            str(profile_path),
            "--snapshot",
            str(snapshot_path),
            "--json",
        ]
        if host_snapshot is not None:
            host_snapshot_path = self.temp_root / "host-snapshot.json"
            host_snapshot_path.write_text(json.dumps(host_snapshot), encoding="utf-8")
            command.extend(("--host-snapshot", str(host_snapshot_path)))
        if output is not None:
            command.extend(("--out", str(output)))
        environment = os.environ.copy()
        environment.update(
            {
                "PATH": "",
                "HTTP_PROXY": "http://127.0.0.1:9",
                "HTTPS_PROXY": "http://127.0.0.1:9",
                "ALL_PROXY": "http://127.0.0.1:9",
                "NO_PROXY": "",
            }
        )
        completed = subprocess.run(
            command,
            cwd=self.temp_root,
            env=environment,
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        payload: dict[str, Any] | None = None
        if output is not None and output.is_file():
            payload = json.loads(output.read_text(encoding="utf-8"))
        elif completed.stdout.strip():
            payload = json.loads(completed.stdout)
        return completed, payload

    @staticmethod
    def _stage(decision: dict[str, Any], step: int) -> dict[str, Any]:
        return next(stage for stage in decision["stages"] if stage["step"] == step)

    def test_accurate_intent_never_silently_falls_back_to_fasttree(self) -> None:
        completed, decision = self._route(
            _profile(), _snapshot("mafft", "trimal", "fasttree")
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIsNotNone(decision)
        assert decision is not None
        self.assertEqual(decision["status"], "blocked")
        self.assertEqual(decision["selected_tree_mode"], "accurate")
        self.assertEqual(self._stage(decision, 8)["reason_code"], "IQTREE2_REQUIRED")
        self.assertIn(
            {"route": "quick", "decision": "rejected-explicit-accurate-no-silent-downgrade"},
            decision["considered_routes"],
        )

    def test_iqtree3_is_detected_but_not_substituted_for_iqtree2(self) -> None:
        completed, decision = self._route(
            _profile(), _snapshot("mafft", "trimal", "iqtree3")
        )
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(
            self._stage(decision, 8)["reason_code"],
            "IQTREE_ALTERNATE_VERSION_REQUIRES_APPROVAL",
        )
        self.assertIn(
            "IQTREE_ALTERNATE_VERSION_NOT_AUTO_SELECTED", decision["warnings"]
        )

    def test_auto_selects_quick_only_when_accurate_is_unavailable(self) -> None:
        profile = _profile()
        profile["intent"] = "auto"
        completed, decision = self._route(
            profile, _snapshot("mafft", "trimal", "fasttree")
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["status"], "ready")
        self.assertEqual(decision["selected_tree_mode"], "quick")
        self.assertEqual(self._stage(decision, 8)["route"], "fasttree")
        self.assertIn("AUTO_FALLBACK_TO_QUICK_TREE", decision["warnings"])

    def test_auto_falls_back_to_planning_when_compute_stack_is_absent(self) -> None:
        profile = _profile()
        profile["intent"] = "auto"
        completed, decision = self._route(profile, _snapshot())
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["status"], "planning-only")
        self.assertIsNone(decision["selected_tree_mode"])
        self.assertEqual(decision["required_steps"], [1, 2, 3, 4])
        self.assertEqual(self._stage(decision, 6)["status"], "skipped")

    def test_explicit_untrimmed_route_does_not_require_trimal(self) -> None:
        profile = _profile()
        profile["requirements"]["trimming"] = False
        completed, decision = self._route(profile, _snapshot("mafft", "iqtree2"))
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["status"], "ready")
        self.assertEqual(
            self._stage(decision, 7)["reason_code"], "TRIMMING_DISABLED_EXPLICITLY"
        )

    def test_untrimmed_route_does_not_bypass_a_missing_alignment(self) -> None:
        profile = _profile()
        profile["requirements"]["trimming"] = False
        completed, decision = self._route(profile, _snapshot("iqtree2"))
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(self._stage(decision, 6)["reason_code"], "MAFFT_REQUIRED")
        self.assertEqual(
            self._stage(decision, 7)["reason_code"],
            "UPSTREAM_ALIGNMENT_CAPABILITY_MISSING",
        )
        self.assertEqual(
            self._stage(decision, 8)["reason_code"],
            "UPSTREAM_ALIGNMENT_CAPABILITY_MISSING",
        )

    def test_existing_alignment_starts_at_alignment_review_without_query_inputs(self) -> None:
        profile = _profile()
        profile["inputs"].update(
            {
                "query_resolved": False,
                "query_location": "absent",
                "candidates_location": "absent",
                "alignment_location": "compute-target",
                "taxdump_location": "absent",
            }
        )
        completed, decision = self._route(
            profile, _snapshot("mafft", "trimal", "iqtree2")
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["required_steps"], [6, 7, 8])
        self.assertEqual(self._stage(decision, 1)["status"], "skipped")
        self.assertEqual(
            self._stage(decision, 6)["reason_code"],
            "ALIGNMENT_ALREADY_MATERIALIZED",
        )
        self.assertEqual(decision["selected_alignment_route"], "existing-alignment")
        self.assertNotIn("mafft", decision["used_tools"])

    def test_existing_alignment_bypasses_obsolete_clustering_trigger(self) -> None:
        profile = _profile()
        profile["intent"] = "auto"
        profile["inputs"]["alignment_location"] = "compute-target"
        profile["inputs"]["candidate_count"] = 500
        completed, decision = self._route(
            profile, _snapshot("trimal", "iqtree2")
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["selected_tree_mode"], "accurate")
        self.assertNotIn(5, decision["required_steps"])

    def test_auto_regenerates_when_existing_tree_cannot_reach_a_rooting_route(self) -> None:
        profile = _profile()
        profile["intent"] = "auto"
        profile["compute"].update({"target": "hpc", "scheduler": "pbs"})
        profile["inputs"].update(
            {
                "alignment_location": "compute-target",
                "tree_location": "host",
                "tree_root_state": "unrooted",
            }
        )
        profile["requirements"]["rooted_tree"] = True
        profile["capabilities"].update(
            {
                "compute_tree_io": True,
                "host_tree_io": False,
                "file_transfer": False,
            }
        )
        completed, decision = self._route(
            profile,
            _snapshot(
                "trimal",
                "iqtree2",
                kind="hpc",
                available_executors=("pbs",),
            ),
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["selected_tree_mode"], "accurate")
        self.assertEqual(self._stage(decision, 8)["status"], "ready")
        self.assertIn("then-root-on-compute-target", self._stage(decision, 8)["route"])

    def test_mmseqs2_is_required_only_after_the_declared_trigger(self) -> None:
        below = _profile()
        completed, decision = self._route(
            below, _snapshot("mafft", "trimal", "iqtree2")
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(self._stage(decision, 5)["status"], "skipped")

        triggered = _profile()
        triggered["inputs"]["candidate_count"] = 200
        completed, decision = self._route(
            triggered, _snapshot("mafft", "trimal", "iqtree2")
        )
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(self._stage(decision, 5)["reason_code"], "MMSEQS2_REQUIRED")

        precomputed = _profile()
        precomputed["inputs"]["candidate_count"] = 200
        precomputed["inputs"]["cluster_mapping_location"] = "compute-target"
        completed, decision = self._route(
            precomputed, _snapshot("mafft", "trimal", "iqtree2")
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(
            self._stage(decision, 5)["reason_code"],
            "CLUSTER_MAPPING_ALREADY_MATERIALIZED",
        )

    def test_auto_falls_back_to_planning_when_triggered_clustering_is_unavailable(self) -> None:
        profile = _profile()
        profile["intent"] = "auto"
        profile["inputs"]["candidate_count"] = 200
        completed, decision = self._route(
            profile, _snapshot("mafft", "trimal", "iqtree2", "fasttree")
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["status"], "planning-only")
        self.assertIsNone(decision["selected_tree_mode"])
        self.assertNotIn(5, decision["required_steps"])

    def test_local_visualization_requires_tree_metadata_and_every_r_package(self) -> None:
        profile = _profile()
        profile["intent"] = "visualization"
        profile["inputs"]["tree_location"] = "host"
        profile["inputs"]["tree_root_state"] = "unrooted"
        profile["inputs"]["metadata_location"] = "host"
        profile["requirements"]["visualization"] = "local"
        incomplete_packages = tuple(name for name in R_PACKAGES if name != "svglite")
        completed, decision = self._route(
            profile,
            _snapshot("rscript", available_packages=incomplete_packages),
        )
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(
            self._stage(decision, 9)["reason_code"], "VISUALIZATION_ROUTE_UNAVAILABLE"
        )

        profile["inputs"]["metadata_location"] = "absent"
        completed, decision = self._route(
            profile,
            _snapshot("rscript", available_packages=R_PACKAGES),
        )
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(self._stage(decision, 9)["reason_code"], "METADATA_INPUT_MISSING")

    def test_itol_permission_is_independent_from_public_database_lookup(self) -> None:
        profile = _profile()
        profile["intent"] = "visualization"
        profile["inputs"]["tree_location"] = "host"
        profile["inputs"]["tree_root_state"] = "unrooted"
        profile["inputs"]["metadata_location"] = "host"
        profile["requirements"]["visualization"] = "itol"
        profile["capabilities"]["web_upload"] = True
        profile["permissions"]["itol_upload"] = True
        self.assertFalse(profile["permissions"]["public_database_network"])
        completed, decision = self._route(profile, _snapshot())
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(self._stage(decision, 9)["route"], "authorized-itol-upload")
        self.assertFalse(decision["execution_authorized"])

    def test_exact_taxonomy_blocks_without_taxdump_or_acquisition_route(self) -> None:
        profile = _profile()
        profile["intent"] = "planning"
        profile["inputs"]["taxdump_location"] = "absent"
        completed, decision = self._route(profile, _snapshot())
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(self._stage(decision, 4)["reason_code"], "LOCAL_TAXDUMP_MISSING")

    def test_current_literature_is_routed_even_when_full_report_is_false(self) -> None:
        profile = _profile()
        profile["requirements"]["current_literature"] = True
        self.assertFalse(profile["requirements"]["report"])
        completed, decision = self._route(
            profile, _snapshot("mafft", "trimal", "iqtree2")
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertIn(10, decision["required_steps"])
        self.assertEqual(self._stage(decision, 10)["status"], "conditional")
        self.assertEqual(decision["status"], "ready-with-limitations")

    def test_current_literature_is_not_dropped_for_planning_or_visualization(self) -> None:
        planning = _profile()
        planning["intent"] = "planning"
        planning["requirements"]["current_literature"] = True
        completed, decision = self._route(planning, _snapshot())
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertIn(10, decision["required_steps"])

        visualization = _profile()
        visualization["intent"] = "visualization"
        visualization["inputs"]["tree_location"] = "host"
        visualization["inputs"]["tree_root_state"] = "unrooted"
        visualization["inputs"]["metadata_location"] = "host"
        visualization["requirements"]["visualization"] = "itol"
        visualization["requirements"]["current_literature"] = True
        visualization["capabilities"]["web_upload"] = True
        visualization["permissions"]["itol_upload"] = True
        completed, decision = self._route(visualization, _snapshot())
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(set(decision["required_steps"]), {9, 10})

    def test_network_denial_and_unpublished_submission_are_fail_closed(self) -> None:
        denied = _profile()
        denied["intent"] = "planning"
        denied["inputs"]["query_resolved"] = False
        denied["inputs"]["candidates_location"] = "absent"
        completed, decision = self._route(denied, _snapshot())
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(
            self._stage(decision, 1)["reason_code"], "QUERY_RESOLUTION_CAPABILITY_MISSING"
        )

        unpublished = deepcopy(denied)
        unpublished["inputs"]["query_kind"] = "unpublished-sequence"
        unpublished["capabilities"]["database_lookup"] = True
        unpublished["permissions"]["public_database_network"] = True
        completed, decision = self._route(unpublished, _snapshot())
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(
            self._stage(decision, 1)["reason_code"],
            "REMOTE_SEQUENCE_SUBMISSION_NOT_ALLOWED",
        )

        public_raw = deepcopy(denied)
        public_raw["inputs"]["query_kind"] = "public-sequence"
        public_raw["capabilities"]["database_lookup"] = True
        public_raw["permissions"]["public_database_network"] = True
        completed, decision = self._route(public_raw, _snapshot())
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(
            self._stage(decision, 1)["reason_code"],
            "REMOTE_SEQUENCE_SUBMISSION_NOT_ALLOWED",
        )

    def test_unknown_raw_sequence_may_search_locally_but_cannot_invent_source_identity(self) -> None:
        profile = _profile()
        profile["intent"] = "planning"
        profile["inputs"].update(
            {
                "query_kind": "unpublished-sequence",
                "query_resolved": False,
                "candidates_location": "absent",
                "sequence_database_location": "compute-target",
                "sequence_database_format": "blast",
            }
        )
        completed, decision = self._route(profile, _snapshot("blastp"))
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(self._stage(decision, 1)["status"], "conditional")
        self.assertEqual(
            self._stage(decision, 1)["reason_code"], "SOURCE_IDENTITY_REVIEW_REQUIRED"
        )
        self.assertEqual(self._stage(decision, 3)["status"], "ready")
        self.assertEqual(self._stage(decision, 4)["reason_code"], "UPSTREAM_INPUT_MISSING")

    def test_raw_sequence_search_requires_the_sequence_at_a_reachable_location(self) -> None:
        profile = _profile()
        profile["intent"] = "planning"
        profile["inputs"].update(
            {
                "query_kind": "unpublished-sequence",
                "query_resolved": False,
                "query_location": "absent",
                "candidates_location": "absent",
                "sequence_database_location": "compute-target",
                "sequence_database_format": "blast",
            }
        )
        profile["capabilities"]["database_lookup"] = True
        profile["permissions"].update(
            {
                "public_database_network": True,
                "remote_sequence_submission": True,
            }
        )
        completed, decision = self._route(profile, _snapshot("blastp"))
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(
            self._stage(decision, 1)["reason_code"], "QUERY_SEQUENCE_INPUT_MISSING"
        )

    def test_local_search_requires_a_database_compatible_with_the_tool(self) -> None:
        profile = _profile()
        profile["intent"] = "planning"
        profile["inputs"].update(
            {
                "query_kind": "unpublished-sequence",
                "query_resolved": False,
                "candidates_location": "absent",
                "sequence_database_location": "compute-target",
                "sequence_database_format": "mmseqs2",
            }
        )
        completed, decision = self._route(profile, _snapshot("blastp"))
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(
            self._stage(decision, 1)["reason_code"],
            "QUERY_RESOLUTION_CAPABILITY_MISSING",
        )

    def test_verified_pbs_target_yields_handoff_not_execution_claim(self) -> None:
        profile = _profile()
        profile["compute"].update({"target": "hpc", "scheduler": "pbs"})
        profile["capabilities"]["file_transfer"] = True
        completed, decision = self._route(
            profile,
            _snapshot(
                "mafft",
                "trimal",
                "iqtree2",
                kind="hpc",
                available_executors=("pbs",),
            ),
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["status"], "handoff-required")
        self.assertFalse(decision["execution_authorized"])
        self.assertEqual(decision["selected_route"], "hpc-offline")

    def test_hybrid_hpc_can_render_on_a_separate_host_snapshot(self) -> None:
        profile = _profile()
        profile["compute"].update({"target": "hpc", "scheduler": "pbs"})
        profile["capabilities"]["file_transfer"] = True
        profile["requirements"]["visualization"] = "local"
        compute_snapshot = _snapshot(
            "mafft",
            "trimal",
            "iqtree2",
            kind="hpc",
            available_executors=("pbs",),
        )
        host_snapshot = _snapshot(
            "rscript", available_packages=R_PACKAGES
        )
        completed, decision = self._route(
            profile, compute_snapshot, host_snapshot=host_snapshot
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(self._stage(decision, 9)["route"], "host-ggtree")
        self.assertEqual(
            decision["host_environment_id"], host_snapshot["environment_id"]
        )

    def test_host_itol_requires_tree_and_metadata_to_reach_host(self) -> None:
        profile = _profile()
        profile["intent"] = "visualization"
        profile["compute"].update({"target": "hpc", "scheduler": "pbs"})
        profile["inputs"]["tree_location"] = "compute-target"
        profile["inputs"]["tree_root_state"] = "unrooted"
        profile["inputs"]["metadata_location"] = "compute-target"
        profile["requirements"]["visualization"] = "itol"
        profile["capabilities"].update({"web_upload": True, "file_transfer": False})
        profile["permissions"]["itol_upload"] = True
        completed, decision = self._route(
            profile,
            _snapshot(kind="hpc", available_executors=("pbs",)),
        )
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(
            self._stage(decision, 9)["reason_code"],
            "VISUALIZATION_ROUTE_UNAVAILABLE",
        )

    def test_remote_target_is_not_a_handoff_when_all_required_work_is_on_host(self) -> None:
        profile = _profile()
        profile["intent"] = "visualization"
        profile["compute"].update({"target": "hpc", "scheduler": "pbs"})
        profile["inputs"]["tree_location"] = "host"
        profile["inputs"]["tree_root_state"] = "unrooted"
        profile["inputs"]["metadata_location"] = "host"
        profile["requirements"]["visualization"] = "local"
        completed, decision = self._route(
            profile,
            _snapshot(kind="hpc", available_executors=("pbs",)),
            host_snapshot=_snapshot(
                "rscript", available_packages=R_PACKAGES
            ),
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["status"], "ready")
        self.assertEqual(decision["selected_route"], "workstation-offline")
        self.assertIn("DECLARED_REMOTE_COMPUTE_TARGET_NOT_USED", decision["warnings"])

    def test_artifact_available_on_both_sites_does_not_force_remote_review(self) -> None:
        profile = _profile()
        profile["compute"].update({"target": "hpc", "scheduler": "pbs"})
        profile["inputs"]["tree_location"] = "both"
        profile["inputs"]["tree_root_state"] = "rooted"
        completed, decision = self._route(
            profile,
            _snapshot(kind="hpc", available_executors=("pbs",)),
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["selected_tree_mode"], "materialized")
        self.assertEqual(self._stage(decision, 8)["environment"], "host")
        self.assertEqual(decision["status"], "ready")

    def test_remote_planning_target_is_reported_as_a_handoff(self) -> None:
        profile = _profile()
        profile["intent"] = "planning"
        profile["compute"].update({"target": "hpc", "scheduler": "pbs"})
        completed, decision = self._route(
            profile,
            _snapshot(kind="hpc", available_executors=("pbs",)),
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["status"], "handoff-required")
        self.assertEqual(decision["selected_route"], "hpc-offline")

    def test_materialized_unrooted_tree_still_requires_a_rooting_route(self) -> None:
        profile = _profile()
        profile["intent"] = "visualization"
        profile["inputs"]["tree_location"] = "host"
        profile["inputs"]["tree_root_state"] = "unrooted"
        profile["inputs"]["metadata_location"] = "host"
        profile["requirements"]["rooted_tree"] = True
        profile["requirements"]["visualization"] = "itol"
        profile["capabilities"]["web_upload"] = True
        profile["permissions"]["itol_upload"] = True
        completed, decision = self._route(profile, _snapshot())
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertIn(8, decision["required_steps"])
        self.assertEqual(self._stage(decision, 8)["reason_code"], "ROOTING_TOOL_REQUIRED")

        completed, decision = self._route(
            profile,
            _snapshot("rscript", available_packages=("ape",)),
        )
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(self._stage(decision, 8)["reason_code"], "ROOTING_TOOL_REQUIRED")

        profile["capabilities"]["host_tree_io"] = True
        completed, decision = self._route(profile, _snapshot())
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(
            self._stage(decision, 8)["reason_code"],
            "EXISTING_TREE_ROOTING_AVAILABLE",
        )

    def test_auto_hpc_route_requires_the_selected_scheduler_launcher(self) -> None:
        profile = _profile()
        profile["intent"] = "auto"
        profile["compute"].update({"target": "hpc", "scheduler": "pbs"})
        profile["inputs"]["alignment_location"] = "compute-target"
        profile["requirements"]["trimming"] = False
        completed, decision = self._route(
            profile,
            _snapshot("iqtree2", "fasttree", kind="hpc"),
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["selected_tree_mode"], None)
        self.assertEqual(decision["status"], "handoff-required")

    def test_visualization_only_route_does_not_apply_tree_job_resource_limits(self) -> None:
        profile = _profile()
        profile["intent"] = "visualization"
        profile["compute"].update({"threads": 64, "memory_gb": 512})
        profile["inputs"]["tree_location"] = "host"
        profile["inputs"]["tree_root_state"] = "unrooted"
        profile["inputs"]["metadata_location"] = "host"
        profile["requirements"]["visualization"] = "itol"
        profile["capabilities"]["web_upload"] = True
        profile["permissions"]["itol_upload"] = True
        completed, decision = self._route(profile, _snapshot())
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["status"], "ready")
        self.assertEqual(
            decision["resource_assessment"]["status"], "not-applicable"
        )
        self.assertEqual(decision["resource_assessment"]["limitations"], [])

    def test_route_label_requires_permission_not_capability_alone(self) -> None:
        profile = _profile()
        profile["capabilities"]["database_lookup"] = True
        profile["capabilities"]["literature_search"] = True
        self.assertFalse(profile["permissions"]["public_database_network"])
        self.assertFalse(profile["permissions"]["literature_network"])
        completed, decision = self._route(
            profile, _snapshot("mafft", "trimal", "iqtree2")
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["selected_route"], "workstation-offline")

    def test_remote_materialized_inputs_must_be_location_scoped(self) -> None:
        profile = _profile()
        profile["compute"].update({"target": "hpc", "scheduler": "pbs"})
        profile["inputs"]["query_location"] = "host"
        profile["inputs"]["candidates_location"] = "host"
        profile["inputs"]["taxdump_location"] = "host"
        profile["capabilities"]["file_transfer"] = False
        completed, decision = self._route(
            profile,
            _snapshot(
                "mafft",
                "trimal",
                "iqtree2",
                kind="hpc",
                available_executors=("pbs",),
            ),
        )
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(
            self._stage(decision, 4)["reason_code"], "PLANNER_INPUT_LOCATION_UNAVAILABLE"
        )

    def test_selected_hpc_scheduler_must_be_observed(self) -> None:
        profile = _profile()
        profile["compute"].update({"target": "hpc", "scheduler": "pbs"})
        completed, decision = self._route(
            profile,
            _snapshot("mafft", "trimal", "iqtree2", kind="hpc"),
        )
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(
            self._stage(decision, 6)["reason_code"], "SCHEDULER_LAUNCHER_MISSING"
        )

    def test_resource_overcommit_is_visible_without_guessing_job_failure(self) -> None:
        profile = _profile()
        profile["compute"].update({"threads": 16, "memory_gb": 64})
        completed, decision = self._route(
            profile, _snapshot("mafft", "trimal", "iqtree2")
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["status"], "ready-with-limitations")
        self.assertEqual(decision["resource_assessment"]["status"], "review-required")
        self.assertEqual(
            set(decision["resource_assessment"]["limitations"]),
            {
                "REQUESTED_THREADS_EXCEED_OBSERVED_CPUS",
                "REQUESTED_MEMORY_EXCEEDS_OBSERVED_MEMORY",
            },
        )

    def test_unknown_resource_observations_require_review(self) -> None:
        profile = _profile()
        snapshot = _snapshot("mafft", "trimal", "iqtree2")
        snapshot["environment"]["cpu_count"] = None
        snapshot["environment"]["memory_gb"] = None
        completed, decision = self._route(profile, snapshot)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["resource_assessment"]["status"], "review-required")
        self.assertEqual(
            set(decision["resource_assessment"]["limitations"]),
            {"OBSERVED_CPU_COUNT_UNKNOWN", "OBSERVED_MEMORY_UNKNOWN"},
        )

    def test_scheduler_is_rejected_for_non_hpc_profile(self) -> None:
        profile = _profile()
        profile["compute"]["scheduler"] = "pbs"
        completed, decision = self._route(
            profile, _snapshot("mafft", "trimal", "iqtree2")
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIsNone(decision)
        self.assertIn("INVALID_ENVIRONMENT_PROFILE", completed.stderr)

    def test_required_nullable_profile_fields_cannot_be_omitted(self) -> None:
        for schema_version in ("0.1", "0.2"):
            with self.subTest(schema=schema_version, field="candidate_count"):
                profile = _profile() if schema_version == "0.1" else _profile_v02()
                del profile["inputs"]["candidate_count"]
                completed, decision = self._route(profile, _snapshot())
                self.assertEqual(completed.returncode, 2)
                self.assertIsNone(decision)
                self.assertIn("candidate_count is required", completed.stderr)
            with self.subTest(schema=schema_version, field="memory_gb"):
                profile = _profile() if schema_version == "0.1" else _profile_v02()
                del profile["compute"]["memory_gb"]
                completed, decision = self._route(profile, _snapshot())
                self.assertEqual(completed.returncode, 2)
                self.assertIsNone(decision)
                self.assertIn("memory_gb is required", completed.stderr)

    def test_current_profile_enum_values_are_not_silently_trimmed(self) -> None:
        for name, mutate in (
            ("profile-id", lambda value: value.update({"profile_id": " test-profile "})),
            ("intent", lambda value: value.update({"intent": " accurate "})),
            (
                "query-kind",
                lambda value: value["inputs"].update(
                    {"query_kind": " public-sequence "}
                ),
            ),
        ):
            with self.subTest(field=name):
                profile = _profile_v02()
                mutate(profile)
                completed, decision = self._route(profile, _snapshot())
                self.assertEqual(completed.returncode, 2)
                self.assertIsNone(decision)
                self.assertIn("leading or trailing whitespace", completed.stderr)

    def test_imported_snapshot_rejects_extra_fields_and_private_paths(self) -> None:
        profile = _profile()
        snapshot = _snapshot("mafft", "trimal", "iqtree2")
        snapshot["hostname"] = "private-host"
        completed, decision = self._route(profile, snapshot)
        self.assertEqual(completed.returncode, 2)
        self.assertIsNone(decision)
        self.assertIn("unsupported fields", completed.stderr)

        snapshot = _snapshot("mafft", "trimal", "iqtree2")
        snapshot["tools"]["iqtree2"]["version"] = "built at /private/home/user/iqtree2"
        completed, decision = self._route(profile, snapshot)
        self.assertEqual(completed.returncode, 2)
        self.assertIsNone(decision)
        self.assertIn("contains a host path", completed.stderr)

    def test_remote_kind_mismatch_and_missing_snapshot_fail_before_routing(self) -> None:
        profile = _profile()
        profile["compute"].update({"target": "hpc", "scheduler": "pbs"})
        completed, decision = self._route(profile, _snapshot(kind="local"))
        self.assertEqual(completed.returncode, 2)
        self.assertIsNone(decision)
        self.assertIn("ENVIRONMENT_KIND_MISMATCH", completed.stderr)

        profile_path = self.temp_root / "remote-profile.json"
        profile_path.write_text(json.dumps(profile), encoding="utf-8")
        completed = subprocess.run(
            [sys.executable, str(SCRIPT), "route", "--profile", str(profile_path)],
            cwd=self.temp_root,
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("REMOTE_SNAPSHOT_REQUIRED", completed.stderr)

    def test_route_hash_is_stable_and_separate_from_scientific_plan_hash(self) -> None:
        profile = _profile()
        snapshot = _snapshot("mafft", "trimal", "iqtree2")
        first_completed, first = self._route(profile, snapshot)
        second_completed, second = self._route(profile, snapshot)
        self.assertEqual(first_completed.returncode, 0)
        self.assertEqual(second_completed.returncode, 0)
        assert first is not None and second is not None
        self.assertEqual(first["route_hash"], second["route_hash"])
        self.assertNotIn("plan_hash", first)

        changed = deepcopy(snapshot)
        changed["tools"]["fasttree"] = _probe("available", "9.9")
        completed, third = self._route(profile, changed)
        self.assertEqual(completed.returncode, 0)
        assert third is not None
        self.assertNotEqual(first["route_hash"], third["route_hash"])

    def test_legacy_profile_remains_an_explicit_protein_route(self) -> None:
        completed, decision = self._route(
            _profile(), _snapshot("mafft", "trimal", "iqtree2")
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["schema_version"], "0.2")
        self.assertEqual(decision["source_profile_schema_version"], "0.1")
        self.assertEqual(decision["resolved_query_molecule"], "protein")
        self.assertEqual(decision["selected_analysis_kind"], "protein")
        self.assertEqual(decision["selected_alignment_route"], "mafft-protein")

    def test_noncoding_dna_uses_blastn_and_nucleotide_alignment(self) -> None:
        profile = _profile_v02("noncoding-dna")
        profile["inputs"].update(
            {
                "candidates_location": "absent",
                "candidates_molecule": "absent",
                "candidate_bundle_kind": "absent",
                "sequence_database_location": "compute-target",
                "sequence_database_format": "blast",
                "sequence_database_molecule": "nucleotide",
            }
        )
        completed, decision = self._route(
            profile, _snapshot("blastn", "mafft", "trimal", "iqtree2")
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["status"], "ready")
        self.assertEqual(decision["selected_search_mode"], "nucleotide-nucleotide")
        self.assertEqual(decision["selected_analysis_kind"], "nucleotide")
        self.assertEqual(decision["selected_alignment_route"], "mafft-nucleotide")
        self.assertEqual(self._stage(decision, 3)["route"], "local-sequence-database")
        self.assertEqual(self._stage(decision, 6)["reason_code"], "NUCLEOTIDE_MSA_AVAILABLE")
        self.assertIn("blastn", decision["used_tools"])
        self.assertNotIn("blastp", decision["used_tools"])

    def test_nucleotide_query_rejects_a_protein_database_route(self) -> None:
        profile = _profile_v02("noncoding-dna")
        profile["inputs"].update(
            {
                "candidates_location": "absent",
                "candidates_molecule": "absent",
                "candidate_bundle_kind": "absent",
                "sequence_database_location": "compute-target",
                "sequence_database_format": "blast",
                "sequence_database_molecule": "protein",
            }
        )
        completed, decision = self._route(
            profile, _snapshot("blastp", "mafft", "trimal", "iqtree2")
        )
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(decision["status"], "blocked")
        self.assertEqual(decision["selected_search_mode"], "none")
        self.assertIn("SEQUENCE_DATABASE_MOLECULE_MISMATCH", decision["warnings"])
        self.assertEqual(
            self._stage(decision, 3)["reason_code"],
            "CANDIDATE_DISCOVERY_CAPABILITY_MISSING",
        )
        self.assertNotIn("blastp", decision["used_tools"])

    def test_clean_cds_selects_protein_alignment_and_codon_backtranslation(self) -> None:
        profile = _profile_v02("coding-dna", genetic_code=11)
        completed, decision = self._route(
            profile, _snapshot("mafft", "trimal", "iqtree2")
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["status"], "ready")
        self.assertEqual(decision["selected_analysis_kind"], "codon")
        self.assertEqual(decision["selected_tree_data_kind"], "codon")
        self.assertEqual(decision["genetic_code"], 11)
        self.assertEqual(
            decision["selected_alignment_route"],
            "translate-mafft-trimal-backtranslate",
        )
        self.assertEqual(self._stage(decision, 6)["reason_code"], "CODON_ALIGNMENT_AVAILABLE")
        self.assertEqual(
            self._stage(decision, 7)["reason_code"],
            "CODON_SAFE_TRIMMING_AVAILABLE",
        )
        operations = {item["operation"] for item in decision["transformations"]}
        self.assertEqual(
            operations,
            {"validate-and-use-cds-translations", "backtranslate-protein-alignment"},
        )

    def test_disrupted_cds_requires_and_selects_macse_without_silent_mafft(self) -> None:
        profile = _profile_v02(
            "coding-dna", coding_status="disrupted", genetic_code=11
        )
        profile["requirements"].update(
            {
                "trimming": False,
                "alignment_strategy": "macse",
                "trimming_strategy": "none",
            }
        )
        completed, decision = self._route(
            profile, _snapshot("mafft", "macse", "iqtree2")
        )
        self.assertEqual(completed.returncode, 2, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["status"], "blocked")
        self.assertEqual(decision["selected_alignment_route"], "macse-codon-aware")
        self.assertEqual(self._stage(decision, 6)["status"], "conditional")
        self.assertEqual(self._stage(decision, 6)["reason_code"], "MACSE_REVIEW_REQUIRED")
        self.assertEqual(
            self._stage(decision, 7)["reason_code"],
            "MACSE_EXPORT_POLICY_REQUIRED",
        )
        self.assertEqual(self._stage(decision, 8)["status"], "blocked")
        self.assertEqual(
            self._stage(decision, 8)["reason_code"],
            "UPSTREAM_ALIGNMENT_CAPABILITY_MISSING",
        )
        self.assertIn("macse", decision["used_tools"])
        self.assertNotIn("mafft", decision["used_tools"])
        self.assertEqual(
            decision["transformations"],
            [
                {
                    "operation": "export-macse-frameshift-symbols",
                    "status": "review-required",
                    "reason_code": "MACSE_EXPORT_POLICY_REQUIRED",
                }
            ],
        )

    def test_structure_aware_rna_requires_the_dedicated_qinsi_executable(self) -> None:
        profile = _profile_v02("noncoding-rna")
        profile["requirements"]["rna_structure_aware"] = True

        missing, decision = self._route(
            profile, _snapshot("mafft", "trimal", "iqtree2")
        )
        self.assertEqual(missing.returncode, 2)
        assert decision is not None
        self.assertIsNone(decision["selected_alignment_route"])
        self.assertEqual(
            self._stage(decision, 6)["reason_code"], "MAFFT_QINSI_REQUIRED"
        )
        self.assertNotIn("mafft", decision["used_tools"])

        ready, decision = self._route(
            profile, _snapshot("mafft-qinsi", "trimal", "iqtree2")
        )
        self.assertEqual(ready.returncode, 0, ready.stderr)
        assert decision is not None
        self.assertEqual(decision["selected_alignment_route"], "mafft-qinsi")
        self.assertEqual(
            self._stage(decision, 6)["reason_code"],
            "RNA_STRUCTURE_MSA_AVAILABLE",
        )
        self.assertIn("mafft-qinsi", decision["used_tools"])
        self.assertNotIn("mafft", decision["used_tools"])

    def test_precomputed_alignment_request_blocks_when_artifact_is_absent(self) -> None:
        profile = _profile_v02("noncoding-dna")
        profile["requirements"]["alignment_strategy"] = "precomputed"
        completed, decision = self._route(
            profile, _snapshot("mafft", "trimal", "iqtree2")
        )
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertIsNone(decision["selected_alignment_route"])
        self.assertEqual(
            self._stage(decision, 6)["reason_code"],
            "PRECOMPUTED_ALIGNMENT_REQUIRED",
        )

    def test_coding_nucleotide_analysis_keeps_codon_alignment_but_uses_dna_tree(self) -> None:
        profile = _profile_v02(
            "coding-dna", analysis_kind="nucleotide", genetic_code=11
        )
        completed, decision = self._route(
            profile, _snapshot("mafft", "trimal", "iqtree2")
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        assert decision is not None
        self.assertEqual(decision["selected_alignment_kind"], "codon")
        self.assertEqual(decision["selected_tree_data_kind"], "nucleotide")
        self.assertEqual(decision["selected_tree_model_family"], "nucleotide")
        self.assertEqual(
            decision["selected_alignment_route"],
            "translate-mafft-trimal-backtranslate",
        )

    def test_candidate_bundle_kind_must_match_the_declared_molecule(self) -> None:
        profile = _profile_v02("coding-dna", genetic_code=11)
        profile["inputs"]["candidate_bundle_kind"] = (
            "noncoding-nucleotide-fasta-metadata"
        )
        completed, decision = self._route(
            profile, _snapshot("mafft", "trimal", "iqtree2")
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIsNone(decision)
        self.assertIn("INVALID_ENVIRONMENT_PROFILE", completed.stderr)

    def test_fasttree_codon_route_fails_closed(self) -> None:
        profile = _profile_v02("coding-dna", genetic_code=11)
        profile["intent"] = "quick"
        completed, decision = self._route(
            profile, _snapshot("mafft", "trimal", "fasttree")
        )
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(decision["selected_tree_mode"], "quick")
        self.assertEqual(
            self._stage(decision, 8)["reason_code"],
            "FASTTREE_CODON_MODEL_UNAVAILABLE",
        )
        self.assertNotIn("fasttree", decision["used_tools"])

    def test_unknown_molecule_blocks_before_alignment_selection(self) -> None:
        profile = _profile_v02("unknown")
        completed, decision = self._route(
            profile, _snapshot("blastp", "blastn", "mafft", "trimal", "iqtree2")
        )
        self.assertEqual(completed.returncode, 2)
        assert decision is not None
        self.assertEqual(decision["status"], "blocked")
        self.assertEqual(decision["resolved_query_molecule"], "unknown")
        self.assertEqual(decision["selected_analysis_kind"], "unknown")
        self.assertEqual(decision["selected_alignment_kind"], "unknown")
        self.assertIsNone(decision["selected_alignment_route"])
        self.assertEqual(
            self._stage(decision, 1)["reason_code"],
            "MOLECULE_CLASSIFICATION_REQUIRED",
        )
        self.assertEqual(self._stage(decision, 6)["status"], "blocked")

    def test_route_hash_changes_with_molecule_and_genetic_code(self) -> None:
        tools = _snapshot("mafft", "trimal", "iqtree2")
        dna_completed, dna = self._route(_profile_v02("noncoding-dna"), tools)
        rna_completed, rna = self._route(_profile_v02("noncoding-rna"), tools)
        code1_completed, code1 = self._route(
            _profile_v02("coding-dna", genetic_code=1), tools
        )
        code11_completed, code11 = self._route(
            _profile_v02("coding-dna", genetic_code=11), tools
        )
        self.assertEqual(
            (dna_completed.returncode, rna_completed.returncode), (0, 0)
        )
        self.assertEqual(
            (code1_completed.returncode, code11_completed.returncode), (0, 0)
        )
        assert dna is not None and rna is not None
        assert code1 is not None and code11 is not None
        self.assertNotEqual(dna["route_hash"], rna["route_hash"])
        self.assertNotEqual(code1["route_hash"], code11["route_hash"])

    def test_blocked_route_is_written_and_existing_output_is_refused(self) -> None:
        output = self.temp_root / "route.json"
        completed, decision = self._route(
            _profile(), _snapshot("mafft", "trimal", "fasttree"), output=output
        )
        self.assertEqual(completed.returncode, 2)
        self.assertTrue(output.is_file())
        assert decision is not None
        self.assertEqual(decision["status"], "blocked")

        completed, _decision = self._route(
            _profile(), _snapshot("mafft", "trimal", "fasttree"), output=output
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("OUTPUT_EXISTS", completed.stderr)

    def test_doctor_default_is_path_only_and_reports_packages_and_executors(self) -> None:
        environment = os.environ.copy()
        environment.update(
            {
                "PATH": "",
                "HTTP_PROXY": "http://127.0.0.1:9",
                "HTTPS_PROXY": "http://127.0.0.1:9",
                "ALL_PROXY": "http://127.0.0.1:9",
                "NO_PROXY": "",
            }
        )
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "doctor",
                "--environment-id",
                "empty-path",
                "--kind",
                "local",
                "--json",
            ],
            cwd=self.temp_root,
            env=environment,
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        snapshot = json.loads(completed.stdout)
        self.assertEqual(snapshot["network"], "not-probed")
        self.assertEqual(snapshot["probe_mode"], "path-only")
        self.assertEqual(snapshot["environment"]["python_version"], sys.version.split()[0])
        self.assertEqual(set(snapshot["r_packages"]), set(R_PACKAGES))
        self.assertEqual(set(snapshot["executors"]), {"ssh", "slurm", "pbs", "lsf"})
        self.assertEqual(set(snapshot["tools"]), set(TOOL_NAMES))
        self.assertTrue(
            {"blastn", "blastx", "tblastn", "tblastx", "pal2nal", "macse"}
            <= set(snapshot["tools"])
        )
        self.assertTrue(all(item["status"] == "missing" for item in snapshot["tools"].values()))
        self.assertNotIn(str(Path.home()), completed.stdout)

    @unittest.skipIf(os.name == "nt", "POSIX wrapper fixture; routing is tested on Windows")
    def test_doctor_path_only_does_not_launch_visible_commands(self) -> None:
        marker = self.temp_root / "launched"
        fake = self.temp_root / "mafft"
        fake.write_text(
            f"#!/bin/sh\n: > '{marker}'\necho 'MAFFT v7.525'\n",
            encoding="utf-8",
        )
        fake.chmod(0o755)
        environment = os.environ.copy()
        environment["PATH"] = str(self.temp_root)
        completed = subprocess.run(
            [sys.executable, str(SCRIPT), "doctor", "--json"],
            cwd=self.temp_root,
            env=environment,
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        snapshot = json.loads(completed.stdout)
        self.assertEqual(snapshot["tools"]["mafft"]["status"], "available")
        self.assertIsNone(snapshot["tools"]["mafft"]["version"])
        self.assertFalse(marker.exists())

    @unittest.skipIf(os.name == "nt", "POSIX wrapper fixture; routing is tested on Windows")
    def test_active_doctor_keeps_only_version_tokens_and_rejects_error_exit(self) -> None:
        mafft = self.temp_root / "mafft"
        mafft.write_text(
            "#!/bin/sh\necho 'MAFFT v7.525 built at /private/secret on private-host'\n",
            encoding="utf-8",
        )
        mafft.chmod(0o755)
        blastp = self.temp_root / "blastp"
        blastp.write_text("#!/bin/sh\necho 'blastp 2.15.0 error'\nexit 1\n", encoding="utf-8")
        blastp.chmod(0o755)
        environment = os.environ.copy()
        environment["PATH"] = str(self.temp_root)
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "doctor",
                "--run-version-probes",
                "--json",
            ],
            cwd=self.temp_root,
            env=environment,
            capture_output=True,
            text=True,
            timeout=20,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        snapshot = json.loads(completed.stdout)
        self.assertEqual(snapshot["probe_mode"], "version-command")
        self.assertEqual(snapshot["tools"]["mafft"]["version"], "7.525")
        self.assertEqual(snapshot["tools"]["blastp"]["status"], "probe-failed")
        self.assertNotIn("private", completed.stdout)

    @unittest.skipIf(os.name == "nt", "POSIX wrapper fixture; routing is tested on Windows")
    def test_active_doctor_accepts_official_mafft_version_only_output(self) -> None:
        for executable_name in ("mafft", "mafft-qinsi"):
            executable = self.temp_root / executable_name
            executable.write_text(
                "#!/bin/sh\necho 'v7.526 (2024/Apr/22)'\n", encoding="utf-8"
            )
            executable.chmod(0o755)
        environment = os.environ.copy()
        environment["PATH"] = str(self.temp_root)
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "doctor",
                "--run-version-probes",
                "--json",
            ],
            cwd=self.temp_root,
            env=environment,
            capture_output=True,
            text=True,
            timeout=20,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        snapshot = json.loads(completed.stdout)
        for name in ("mafft", "mafft-qinsi"):
            with self.subTest(tool=name):
                self.assertEqual(snapshot["tools"][name]["status"], "available")
                self.assertEqual(snapshot["tools"][name]["version"], "7.526")

    def test_environment_schemas_and_example_are_valid_json(self) -> None:
        profile = json.loads(
            (SKILL_ROOT / "assets" / "environment-profile.example.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(profile["schema_version"], "0.2")
        for name in (
            "environment-profile-0.1.schema.json",
            "environment-profile-0.2.schema.json",
            "environment-snapshot-0.1.schema.json",
            "route-decision-0.1.schema.json",
            "route-decision-0.2.schema.json",
        ):
            schema = json.loads((SKILL_ROOT / "references" / name).read_text(encoding="utf-8"))
            self.assertEqual(schema["$schema"], "https://json-schema.org/draft/2020-12/schema")


if __name__ == "__main__":
    unittest.main()
