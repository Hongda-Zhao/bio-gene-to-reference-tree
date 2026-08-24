"""Draft 2020-12 validation for molecule-aware requests, plans, and routes."""

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

try:
    from jsonschema import Draft202012Validator, FormatChecker
except ImportError:  # pragma: no cover - CI installs the pinned validator.
    Draft202012Validator = None  # type: ignore[assignment]
    FormatChecker = None  # type: ignore[assignment]

from tests import test_nucleotide_workflow as nucleotide_workflow
from tests.test_environment_routing import _profile, _profile_v02, _snapshot


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SKILL_ROOT = REPOSITORY_ROOT / "skills" / "bio-gene-to-reference-tree"
REFERENCES = SKILL_ROOT / "references"
ASSETS = SKILL_ROOT / "assets"
SCRIPT = SKILL_ROOT / "scripts" / "gene_to_tree.py"
SCHEMA_FILES = (
    "request-0.3.schema.json",
    "plan-0.4.schema.json",
    "environment-profile-0.2.schema.json",
    "route-decision-0.2.schema.json",
)


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON object key: {key}")
        result[key] = value
    return result


def _load_json(path: Path, *, reject_duplicates: bool = False) -> dict[str, Any]:
    kwargs = {"object_pairs_hook": _reject_duplicate_keys} if reject_duplicates else {}
    value = json.loads(path.read_text(encoding="utf-8"), **kwargs)
    if not isinstance(value, dict):
        raise TypeError(f"Expected a JSON object in {path}")
    return value


@unittest.skipIf(Draft202012Validator is None, "jsonschema is not installed")
class MoleculeSchemaTests(unittest.TestCase):
    maxDiff = None

    def setUp(self) -> None:
        self._temporary_directory = tempfile.TemporaryDirectory(
            prefix="molecule schema contracts "
        )
        self.temp_root = Path(self._temporary_directory.name)
        assert Draft202012Validator is not None
        assert FormatChecker is not None
        self.validators = {
            name: Draft202012Validator(
                _load_json(REFERENCES / name), format_checker=FormatChecker()
            )
            for name in SCHEMA_FILES
        }

    def tearDown(self) -> None:
        self._temporary_directory.cleanup()

    def _validate(self, schema_name: str, instance: dict[str, Any]) -> None:
        errors = sorted(
            self.validators[schema_name].iter_errors(instance),
            key=lambda error: list(error.absolute_path),
        )
        if errors:
            details = "\n".join(
                f"{'.'.join(map(str, error.absolute_path)) or '<root>'}: {error.message}"
                for error in errors
            )
            self.fail(f"{schema_name} rejected the instance:\n{details}")

    def _run_route(
        self, profile: dict[str, Any], snapshot: dict[str, Any], name: str
    ) -> tuple[subprocess.CompletedProcess[str], dict[str, Any]]:
        profile_path = self.temp_root / f"{name}.profile.json"
        snapshot_path = self.temp_root / f"{name}.snapshot.json"
        profile_path.write_text(json.dumps(profile), encoding="utf-8")
        snapshot_path.write_text(json.dumps(snapshot), encoding="utf-8")
        environment = os.environ.copy()
        environment["PATH"] = ""
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "route",
                "--profile",
                str(profile_path),
                "--snapshot",
                str(snapshot_path),
                "--json",
            ],
            cwd=self.temp_root,
            env=environment,
            capture_output=True,
            text=True,
            timeout=15,
            check=False,
        )
        self.assertTrue(completed.stdout.strip(), completed.stderr)
        decision = json.loads(completed.stdout)
        self.assertIsInstance(decision, dict)
        return completed, decision

    def test_schema_documents_are_duplicate_free_and_valid_draft_2020_12(self) -> None:
        assert Draft202012Validator is not None
        for name in SCHEMA_FILES:
            with self.subTest(schema=name):
                schema = _load_json(REFERENCES / name, reject_duplicates=True)
                Draft202012Validator.check_schema(schema)

    def test_checked_in_request_and_environment_profile_validate(self) -> None:
        self._validate(
            "request-0.3.schema.json",
            _load_json(ASSETS / "request.nucleotide.example.json"),
        )
        self._validate(
            "environment-profile-0.2.schema.json",
            _load_json(ASSETS / "environment-profile.example.json"),
        )

    def test_request_schema_covers_protein_dna_rna_and_clean_cds_routes(self) -> None:
        base = _load_json(ASSETS / "request.nucleotide.example.json")
        variants: dict[str, dict[str, Any]] = {"noncoding-dna": base}

        for source_encoding in ("rna-u", "dna-t"):
            request = deepcopy(base)
            request["molecule"].update(
                {"type": "noncoding-rna", "source_encoding": source_encoding}
            )
            variants[f"noncoding-rna-{source_encoding}"] = request

        protein = deepcopy(base)
        protein["query"].update(
            {"kind": "protein-fasta", "strand": "not-applicable"}
        )
        protein["molecule"].update(
            {
                "type": "protein",
                "analysis_kind": "protein",
                "coding_status": "not-applicable",
                "genetic_code": None,
                "alignment_strategy": "mafft-protein",
                "trimming_strategy": "trimal-columns",
                "rna_structure_aware": False,
                "source_encoding": "not-applicable",
            }
        )
        variants["protein"] = protein

        for molecule_type, source_encoding in (
            ("coding-dna", "not-applicable"),
            ("coding-rna", "rna-u"),
        ):
            for analysis_kind in ("codon", "nucleotide"):
                request = deepcopy(base)
                request["molecule"].update(
                    {
                        "type": molecule_type,
                        "analysis_kind": analysis_kind,
                        "coding_status": "clean",
                        "genetic_code": 11,
                        "alignment_strategy": "translate-mafft-backtranslate",
                        "trimming_strategy": "protein-mask-to-codons",
                        "rna_structure_aware": False,
                        "source_encoding": source_encoding,
                    }
                )
                request["references"]["translation_fasta"] = "translations.faa"
                variants[f"{molecule_type}-{analysis_kind}"] = request

        for name, request in variants.items():
            with self.subTest(route=name):
                self._validate("request-0.3.schema.json", request)

        standard_bootstrap = deepcopy(base)
        standard_bootstrap["tree"]["support"] = {
            "method": "standard-bootstrap",
            "replicates": 1000,
            "sh_alrt": 1000,
        }
        self._validate("request-0.3.schema.json", standard_bootstrap)

        invalid_standard_bnni = deepcopy(standard_bootstrap)
        invalid_standard_bnni["tree"]["support"]["bnni"] = True
        self.assertTrue(
            list(
                self.validators["request-0.3.schema.json"].iter_errors(
                    invalid_standard_bnni
                )
            )
        )

        invalid_accurate_model = deepcopy(base)
        invalid_accurate_model["tree"]["model"] = "LG+G4"
        self.assertTrue(
            list(
                self.validators["request-0.3.schema.json"].iter_errors(
                    invalid_accurate_model
                )
            )
        )

        invalid = deepcopy(base)
        invalid["unexpected_top_level"] = True
        self.assertTrue(
            list(self.validators["request-0.3.schema.json"].iter_errors(invalid))
        )
        invalid = deepcopy(base)
        invalid["query"]["retrieved_at"] = "not-a-date"
        self.assertTrue(
            list(self.validators["request-0.3.schema.json"].iter_errors(invalid))
        )
        invalid["query"]["retrieved_at"] = "2026-99-99"
        self.assertTrue(
            list(self.validators["request-0.3.schema.json"].iter_errors(invalid))
        )

    def test_generated_nucleotide_rna_and_codon_plans_validate(self) -> None:
        fixture = nucleotide_workflow.NucleotideWorkflowTests(
            methodName="test_noncoding_dna_example_uses_explicit_nucleotide_commands"
        )
        fixture.setUp()
        try:
            completed, output = fixture._run(ASSETS / "request.nucleotide.example.json")
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self._validate("plan-0.4.schema.json", _load_json(output / "plan.json"))

            untrimmed_request = _load_json(ASSETS / "request.nucleotide.example.json")
            untrimmed_request["query"]["path"] = str(
                ASSETS / "query.nucleotide.example.fna"
            )
            untrimmed_request["references"].update(
                {
                    "candidate_fasta": str(
                        ASSETS / "candidates.nucleotide.example.fna"
                    ),
                    "candidate_table": str(
                        ASSETS / "candidates.nucleotide.example.tsv"
                    ),
                }
            )
            untrimmed_request["molecule"]["trimming_strategy"] = "none"
            untrimmed_request["trimming"]["enabled"] = False
            untrimmed_path = fixture.temp_root / "schema-untrimmed-request.json"
            untrimmed_path.write_text(json.dumps(untrimmed_request), encoding="utf-8")
            completed, output = fixture._run(
                untrimmed_path, output_name="schema-untrimmed"
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            untrimmed_plan = _load_json(output / "plan.json")
            self.assertEqual(untrimmed_plan["trimming_plan"]["primary_profile"], "")
            self.assertEqual(untrimmed_plan["trimming_plan"]["profiles"], [])
            self._validate("plan-0.4.schema.json", untrimmed_plan)

            rna_sequences = {
                "QUERY_NT1": "ACGTACGTAA" * 6,
                "REF_NT1": "ACGTACGTAG" * 6,
                "REF_NT2": "ACGTACGCAA" * 6,
                "OUT_NT1": "ACGTTCGTAA" * 6,
            }
            query, candidates, table, _ = fixture._write_nucleotide_bundle(
                molecule_type="noncoding-rna", sequences=rna_sequences
            )
            request = fixture._request_for_bundle(
                "noncoding-rna", query, candidates, table
            )
            request["molecule"]["source_encoding"] = "dna-t"
            request_path = fixture.temp_root / "schema-rna-request.json"
            request_path.write_text(json.dumps(request), encoding="utf-8")
            completed, output = fixture._run(request_path, output_name="schema-rna")
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self._validate("plan-0.4.schema.json", _load_json(output / "plan.json"))

            cds_sequences = {
                "QUERY_NT1": "ATG" + "GCT" * 8 + "TAA",
                "REF_NT1": "ATG" + "GCC" * 8 + "TAA",
                "REF_NT2": "ATG" + "GCA" * 8 + "TAA",
                "OUT_NT1": "ATG" + "GTT" * 8 + "TAA",
            }
            translations = {
                "QUERY_NT1": "M" + "A" * 8,
                "REF_NT1": "M" + "A" * 8,
                "REF_NT2": "M" + "A" * 8,
                "OUT_NT1": "M" + "V" * 8,
            }
            query, candidates, table, translation = fixture._write_nucleotide_bundle(
                molecule_type="coding-dna",
                sequences=cds_sequences,
                translations=translations,
            )
            assert translation is not None
            for analysis_kind in ("codon", "nucleotide"):
                request = fixture._request_for_bundle(
                    "coding-dna", query, candidates, table, translation
                )
                request["molecule"].update(
                    {
                        "analysis_kind": analysis_kind,
                        "coding_status": "clean",
                        "genetic_code": 11,
                        "alignment_strategy": "translate-mafft-backtranslate",
                        "trimming_strategy": "protein-mask-to-codons",
                    }
                )
                request_path = fixture.temp_root / f"schema-{analysis_kind}.json"
                request_path.write_text(json.dumps(request), encoding="utf-8")
                completed, output = fixture._run(
                    request_path, output_name=f"schema-{analysis_kind}"
                )
                self.assertEqual(completed.returncode, 0, completed.stderr)
                self._validate(
                    "plan-0.4.schema.json", _load_json(output / "plan.json")
                )

            rna_cds_sequences = {
                identifier: sequence.replace("T", "U")
                for identifier, sequence in cds_sequences.items()
            }
            query, candidates, table, translation = fixture._write_nucleotide_bundle(
                molecule_type="coding-rna",
                sequences=rna_cds_sequences,
                translations=translations,
            )
            assert translation is not None
            request = fixture._request_for_bundle(
                "coding-rna", query, candidates, table, translation
            )
            request["molecule"].update(
                {
                    "analysis_kind": "codon",
                    "coding_status": "clean",
                    "genetic_code": 11,
                    "alignment_strategy": "translate-mafft-backtranslate",
                    "trimming_strategy": "protein-mask-to-codons",
                    "source_encoding": "rna-u",
                }
            )
            request_path = fixture.temp_root / "schema-coding-rna.json"
            request_path.write_text(json.dumps(request), encoding="utf-8")
            completed, output = fixture._run(
                request_path, output_name="schema-coding-rna"
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self._validate("plan-0.4.schema.json", _load_json(output / "plan.json"))
        finally:
            fixture.tearDown()

    def test_generated_environment_route_matrix_validates(self) -> None:
        cases: list[tuple[str, dict[str, Any], dict[str, Any], int]] = []
        tools = _snapshot("mafft", "trimal", "iqtree2")
        for name, profile in (
            ("protein", _profile_v02("protein")),
            ("noncoding-dna", _profile_v02("noncoding-dna")),
            ("coding-dna-codon", _profile_v02("coding-dna", genetic_code=11)),
            (
                "coding-dna-nucleotide",
                _profile_v02(
                    "coding-dna", analysis_kind="nucleotide", genetic_code=11
                ),
            ),
            ("coding-rna", _profile_v02("coding-rna", genetic_code=11)),
        ):
            cases.append((name, profile, tools, 0))

        structure_rna = _profile_v02("noncoding-rna")
        structure_rna["requirements"]["rna_structure_aware"] = True
        cases.append(
            (
                "structure-rna",
                structure_rna,
                _snapshot("mafft-qinsi", "trimal", "iqtree2"),
                0,
            )
        )

        disrupted = _profile_v02(
            "coding-dna", coding_status="disrupted", genetic_code=11
        )
        disrupted["requirements"].update(
            {
                "trimming": False,
                "alignment_strategy": "macse",
                "trimming_strategy": "none",
            }
        )
        cases.append(
            ("disrupted-cds", disrupted, _snapshot("macse", "iqtree2"), 2)
        )

        unknown = _profile_v02("unknown")
        cases.append(("unknown", unknown, _snapshot(), 2))

        unknown_visualization = _profile_v02("unknown")
        unknown_visualization["intent"] = "visualization"
        cases.append(
            ("unknown-visualization-without-tree", unknown_visualization, _snapshot(), 2)
        )

        for name, profile, snapshot, returncode in cases:
            with self.subTest(route=name):
                self._validate("environment-profile-0.2.schema.json", profile)
                completed, decision = self._run_route(profile, snapshot, name)
                self.assertEqual(completed.returncode, returncode, completed.stderr)
                self._validate("route-decision-0.2.schema.json", decision)

    def test_legacy_remote_and_mmseqs_search_decisions_validate(self) -> None:
        remote = _profile()
        remote["inputs"].update(
            {"candidates_location": "absent", "candidate_count": None}
        )
        remote["capabilities"]["database_lookup"] = True
        remote["permissions"]["public_database_network"] = True

        local_mmseqs = _profile()
        local_mmseqs["inputs"].update(
            {
                "candidates_location": "absent",
                "sequence_database_location": "compute-target",
                "sequence_database_format": "mmseqs2",
            }
        )

        cases = (
            (
                "legacy-provider",
                remote,
                _snapshot("mafft", "trimal", "iqtree2"),
                "provider-molecule-aware",
            ),
            (
                "legacy-mmseqs",
                local_mmseqs,
                _snapshot("mmseqs2", "mafft", "trimal", "iqtree2"),
                "mmseqs2-protein",
            ),
        )
        for name, profile, snapshot, expected_search in cases:
            with self.subTest(route=name):
                completed, decision = self._run_route(profile, snapshot, name)
                self.assertEqual(completed.returncode, 0, completed.stderr)
                self.assertEqual(decision["source_profile_schema_version"], "0.1")
                self.assertEqual(decision["selected_search_mode"], expected_search)
                self._validate("route-decision-0.2.schema.json", decision)


if __name__ == "__main__":
    unittest.main()
