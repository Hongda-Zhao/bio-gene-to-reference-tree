"""Black-box contracts for molecule-aware nucleotide and codon planning."""

from __future__ import annotations

import csv
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
ASSETS = SKILL_ROOT / "assets"


class NucleotideWorkflowTests(unittest.TestCase):
    maxDiff = None

    def setUp(self) -> None:
        self._temporary_directory = tempfile.TemporaryDirectory(
            prefix="nucleotide gene tree "
        )
        self.temp_root = Path(self._temporary_directory.name)

    def tearDown(self) -> None:
        self._temporary_directory.cleanup()

    def _run(
        self, request_path: Path, *, output_name: str = "out"
    ) -> tuple[subprocess.CompletedProcess[str], Path]:
        output = self.temp_root / output_name
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
                "plan",
                "--request",
                str(request_path),
                "--offline",
                "--dry-run",
                "--out",
                str(output),
            ],
            cwd=self.temp_root,
            env=environment,
            capture_output=True,
            text=True,
            timeout=15,
            check=False,
        )
        return completed, output

    @staticmethod
    def _commands(plan: dict[str, Any]) -> dict[str, dict[str, Any]]:
        return {command["id"]: command for command in plan["planned_commands"]}

    def test_noncoding_dna_example_uses_explicit_nucleotide_commands(self) -> None:
        completed, output = self._run(ASSETS / "request.nucleotide.example.json")
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue((output / "reference_set.fna").is_file())
        self.assertFalse((output / "reference_set.faa").exists())
        plan = json.loads((output / "plan.json").read_text(encoding="utf-8"))
        self.assertEqual(plan["schema_version"], "0.4")
        self.assertEqual(plan["molecule_plan"]["type"], "noncoding-dna")
        self.assertEqual(plan["artifact_plan"]["reference_fasta"]["data_kind"], "nucleotide")
        self.assertEqual(plan["artifact_plan"]["tree_input"]["data_kind"], "nucleotide")
        commands = self._commands(plan)
        self.assertIn("--nuc", commands["align-nucleotides"]["argv"])
        self.assertIn("alignment.raw.fna", commands["trim-balanced"]["argv"])
        iqtree = commands["infer-accurate-tree"]["argv"]
        self.assertEqual(iqtree[iqtree.index("-st") + 1], "DNA")
        with (output / "sequence_metadata.tsv").open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            metadata = list(csv.DictReader(handle, delimiter="\t"))
        query_row = next(row for row in metadata if row["accession"] == "QUERY_NT1")
        self.assertEqual(query_row["coordinates"], "synthetic:1..60")
        self.assertEqual(query_row["analysis_kind"], "nucleotide")

    def test_standard_bootstrap_defaults_bnni_off(self) -> None:
        request = json.loads(
            (ASSETS / "request.nucleotide.example.json").read_text(encoding="utf-8")
        )
        request["query"]["path"] = str(ASSETS / "query.nucleotide.example.fna")
        request["references"].update(
            {
                "candidate_fasta": str(
                    ASSETS / "candidates.nucleotide.example.fna"
                ),
                "candidate_table": str(
                    ASSETS / "candidates.nucleotide.example.tsv"
                ),
            }
        )
        request["tree"]["support"] = {
            "method": "standard-bootstrap",
            "replicates": 1000,
            "sh_alrt": 1000,
        }
        request_path = self.temp_root / "standard-bootstrap-request.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")
        completed, output = self._run(request_path, output_name="standard-bootstrap-out")
        self.assertEqual(completed.returncode, 0, completed.stderr)
        plan = json.loads((output / "plan.json").read_text(encoding="utf-8"))
        self.assertFalse(plan["tree_plan"]["bnni"])
        argv = self._commands(plan)["infer-accurate-tree"]["argv"]
        self.assertIn("-b", argv)
        self.assertNotIn("-bnni", argv)

    def test_positive_target_coverage_threshold_rejects_missing_values(self) -> None:
        request = json.loads(
            (ASSETS / "request.nucleotide.example.json").read_text(encoding="utf-8")
        )
        request["query"]["path"] = str(ASSETS / "query.nucleotide.example.fna")
        request["references"].update(
            {
                "candidate_fasta": str(
                    ASSETS / "candidates.nucleotide.example.fna"
                ),
                "candidate_table": str(
                    ASSETS / "candidates.nucleotide.example.tsv"
                ),
            }
        )
        request["selection"]["min_target_coverage"] = 0.8
        request_path = self.temp_root / "missing-target-coverage-request.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")
        completed, output = self._run(
            request_path, output_name="missing-target-coverage-out"
        )
        self.assertEqual(completed.returncode, 3, completed.stderr)
        rejected = (output / "rejected_references.tsv").read_text(encoding="utf-8")
        self.assertIn("TARGET_COVERAGE_MISSING", rejected)
        plan = json.loads((output / "plan.json").read_text(encoding="utf-8"))
        self.assertEqual(plan["state"], "blocked")

    def _write_nucleotide_bundle(
        self,
        *,
        molecule_type: str,
        sequences: dict[str, str],
        translations: dict[str, str] | None = None,
        genetic_code: int = 11,
    ) -> tuple[Path, Path, Path, Path | None]:
        query_path = self.temp_root / f"query-{molecule_type}.fasta"
        candidate_path = self.temp_root / f"candidates-{molecule_type}.fasta"
        table_path = self.temp_root / f"candidates-{molecule_type}.tsv"
        query_id = "QUERY_NT1"
        query_path.write_text(
            f">{query_id}\n{sequences[query_id]}\n", encoding="utf-8"
        )
        candidate_path.write_text(
            "".join(f">{identifier}\n{sequence}\n" for identifier, sequence in sequences.items()),
            encoding="utf-8",
        )
        fieldnames = [
            "accession",
            "taxon_id",
            "species",
            "role",
            "relation",
            "is_reviewed",
            "is_canonical",
            "is_fragment",
            "query_coverage",
            "sequence_length",
            "bitscore",
            "evalue",
            "source_db",
            "source_release",
            "retrieved_at",
            "clade",
            "molecule_type",
            "feature_name",
            "sequence_region",
            "coordinates",
            "actual_search_database",
            "search_program",
            "search_database_molecule",
            "outgroup_rationale",
            "genetic_code_id",
            "reading_frame_source",
            "strand",
            "complete_cds",
            "internal_stop_count",
            "frameshift_count",
            "translation_matches",
        ]
        metadata = {
            "QUERY_NT1": ("9606", "Homo sapiens", "ingroup", "self", "Mammalia"),
            "REF_NT1": ("10090", "Mus musculus", "ingroup", "one2one_ortholog", "Mammalia"),
            "REF_NT2": ("7955", "Danio rerio", "ingroup", "ortholog", "Actinopterygii"),
            "OUT_NT1": ("8364", "Xenopus tropicalis", "outgroup", "homolog", "Amphibia"),
        }
        coding = molecule_type in {"coding-dna", "coding-rna"}
        with table_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for index, (identifier, sequence) in enumerate(sequences.items()):
                taxid, species, role, relation, clade = metadata[identifier]
                writer.writerow(
                    {
                        "accession": identifier,
                        "taxon_id": taxid,
                        "species": species,
                        "role": role,
                        "relation": relation,
                        "is_reviewed": "true",
                        "is_canonical": "true",
                        "is_fragment": "false",
                        "query_coverage": "1",
                        "sequence_length": str(len(sequence)),
                        "bitscore": str(100 - index),
                        "evalue": "0" if index == 0 else "1e-20",
                        "source_db": "SyntheticDB",
                        "source_release": "2026-08",
                        "retrieved_at": "2026-08-25",
                        "clade": clade,
                        "molecule_type": molecule_type,
                        "feature_name": "synthetic CDS" if coding else "synthetic noncoding locus",
                        "sequence_region": "complete synthetic locus",
                        "coordinates": f"synthetic:1..{len(sequence)}",
                        "actual_search_database": "local-bundle:not-searched",
                        "search_program": "not-searched",
                        "search_database_molecule": "not-applicable",
                        "outgroup_rationale": (
                            "Declared distant homolog for the synthetic test fixture"
                            if role == "outgroup"
                            else ""
                        ),
                        "genetic_code_id": str(genetic_code) if coding else "",
                        "reading_frame_source": "database-CDS" if coding else "",
                        "strand": "plus",
                        "complete_cds": "true" if coding else "",
                        "internal_stop_count": "0" if coding else "",
                        "frameshift_count": "0" if coding else "",
                        "translation_matches": "true" if coding else "",
                    }
                )
        translation_path: Path | None = None
        if translations is not None:
            translation_path = self.temp_root / "translations.faa"
            translation_path.write_text(
                "".join(
                    f">{identifier}\n{sequence}\n"
                    for identifier, sequence in translations.items()
                ),
                encoding="utf-8",
            )
        return query_path, candidate_path, table_path, translation_path

    def _request_for_bundle(
        self,
        molecule_type: str,
        query_path: Path,
        candidate_path: Path,
        table_path: Path,
        translation_path: Path | None = None,
    ) -> dict[str, Any]:
        request = json.loads(
            (ASSETS / "request.nucleotide.example.json").read_text(encoding="utf-8")
        )
        request["query"].update(
            {
                "path": str(query_path),
                "id": "QUERY_NT1",
                "strand": "plus",
                "feature_name": (
                    "synthetic CDS"
                    if molecule_type in {"coding-dna", "coding-rna"}
                    else "synthetic noncoding locus"
                ),
            }
        )
        request["references"].update(
            {
                "candidate_fasta": str(candidate_path),
                "candidate_table": str(table_path),
            }
        )
        request["molecule"]["type"] = molecule_type
        request["molecule"]["source_encoding"] = (
            "rna-u" if molecule_type in {"coding-rna", "noncoding-rna"} else "not-applicable"
        )
        if translation_path is not None:
            request["references"]["translation_fasta"] = str(translation_path)
        return request

    def test_rna_preserves_source_and_materializes_a_t_normalized_copy(self) -> None:
        dna = {
            "QUERY_NT1": "ACGUACGUAA" * 6,
            "REF_NT1": "ACGUACGUAG" * 6,
            "REF_NT2": "ACGUACGCAA" * 6,
            "OUT_NT1": "ACGUUCGUAA" * 6,
        }
        query, candidates, table, _ = self._write_nucleotide_bundle(
            molecule_type="noncoding-rna", sequences=dna
        )
        request = self._request_for_bundle(
            "noncoding-rna", query, candidates, table
        )
        request["molecule"]["rna_structure_aware"] = True
        request["alignment"]["mode"] = "qinsi"
        request_path = self.temp_root / "rna-request.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")
        completed, output = self._run(request_path)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("U", (output / "reference_set.rna.fasta").read_text())
        normalized_sequences = "".join(
            line
            for line in (output / "reference_set.fna").read_text().splitlines()
            if not line.startswith(">")
        )
        self.assertNotIn("U", normalized_sequences)
        commands = self._commands(json.loads((output / "plan.json").read_text()))
        plan = json.loads((output / "plan.json").read_text())
        self.assertEqual(
            plan["artifact_plan"]["source_rna_fasta"]["data_kind"],
            "noncoding-rna",
        )
        with (output / "sequence_metadata.tsv").open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            metadata = list(csv.DictReader(handle, delimiter="\t"))
        query_metadata = next(
            row for row in metadata if row["accession"] == "QUERY_NT1"
        )
        self.assertEqual(query_metadata["source_encoding"], "rna-u")
        self.assertNotEqual(
            query_metadata["source_sequence_sha256"],
            query_metadata["analysis_sequence_sha256"],
        )
        self.assertEqual(commands["align-nucleotides"]["argv"][0], "mafft-qinsi")
        self.assertNotIn("--qinsi", commands["align-nucleotides"]["argv"])
        self.assertIn("--nuc", commands["align-nucleotides"]["argv"])

        t_encoded = {identifier: sequence.replace("U", "T") for identifier, sequence in dna.items()}
        query, candidates, table, _ = self._write_nucleotide_bundle(
            molecule_type="noncoding-rna", sequences=t_encoded
        )
        request = self._request_for_bundle("noncoding-rna", query, candidates, table)
        request["molecule"]["source_encoding"] = "dna-t"
        t_request_path = self.temp_root / "rna-t-request.json"
        t_request_path.write_text(json.dumps(request), encoding="utf-8")
        completed, t_output = self._run(t_request_path, output_name="rna-t-out")
        self.assertEqual(completed.returncode, 0, completed.stderr)
        source_text = (t_output / "reference_set.rna.fasta").read_text()
        analysis_text = (t_output / "reference_set.fna").read_text()
        self.assertEqual(source_text, analysis_text)
        t_plan = json.loads((t_output / "plan.json").read_text())
        self.assertEqual(t_plan["molecule_plan"]["source_encoding"], "dna-t")
        self.assertEqual(
            t_plan["molecule_plan"]["transformations"][0]["operation"],
            "copy-t-encoded-rna-to-analysis-fasta",
        )

    def test_codon_plan_aligns_translations_and_backtranslates_complete_triplets(self) -> None:
        sequences = {
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
        query, candidates, table, translation = self._write_nucleotide_bundle(
            molecule_type="coding-dna",
            sequences=sequences,
            translations=translations,
        )
        assert translation is not None
        request = self._request_for_bundle(
            "coding-dna", query, candidates, table, translation
        )
        request["molecule"].update(
            {
                "analysis_kind": "codon",
                "coding_status": "clean",
                "genetic_code": 11,
                "alignment_strategy": "translate-mafft-backtranslate",
                "trimming_strategy": "protein-mask-to-codons",
                "rna_structure_aware": False,
            }
        )
        request_path = self.temp_root / "codon-request.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")
        completed, output = self._run(request_path)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue((output / "reference_set.translated.faa").is_file())
        plan = json.loads((output / "plan.json").read_text())
        commands = self._commands(plan)
        self.assertEqual(plan["artifact_plan"]["raw_alignment"]["data_kind"], "protein")
        self.assertEqual(plan["artifact_plan"]["tree_input"]["data_kind"], "codon")
        self.assertEqual(plan["artifact_plan"]["unrooted_tree"]["inferred_from"], "codon")
        self.assertEqual(
            plan["artifact_plan"]["backtranslated_alignment"]["path"],
            "alignment.raw.codon.fna",
        )
        self.assertEqual(
            plan["artifact_plan"]["trimmed_alignment"]["path"],
            "alignment.trimmed.balanced.codon.fna",
        )
        self.assertIn("--amino", commands["align-translated-proteins"]["argv"])
        untrimmed = commands["backtranslate-untrimmed-codons"]
        self.assertEqual(untrimmed["outputs"], ["alignment.raw.codon.fna"])
        self.assertEqual(untrimmed["argv"][untrimmed["argv"].index("-gt") + 1], "0")
        self.assertIn("-ignorestopcodon", untrimmed["argv"])
        backtranslate = commands["backtranslate-balanced"]["argv"]
        self.assertIn("-backtrans", backtranslate)
        self.assertIn("reference_set.fna", backtranslate)
        self.assertIn("-ignorestopcodon", backtranslate)
        iqtree = commands["infer-accurate-tree"]["argv"]
        self.assertEqual(iqtree[iqtree.index("-st") + 1], "CODON11")
        with (output / "sequence_metadata.tsv").open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            metadata = list(csv.DictReader(handle, delimiter="\t"))
        query_metadata = next(
            row for row in metadata if row["accession"] == "QUERY_NT1"
        )
        self.assertEqual(
            query_metadata["backtranslation_qc_status"], "pending-execution"
        )
        self.assertEqual(len(query_metadata["translation_sequence_sha256"]), 64)

        code1_sequences = dict(sequences)
        code1_sequences["QUERY_NT1"] = "TTG" + "GCT" * 8 + "TAA"
        code1_sequences["REF_NT1"] = "CTG" + "GCC" * 8 + "TAA"
        query, candidates, table, translation = self._write_nucleotide_bundle(
            molecule_type="coding-dna",
            sequences=code1_sequences,
            translations=translations,
            genetic_code=1,
        )
        assert translation is not None
        request = self._request_for_bundle(
            "coding-dna", query, candidates, table, translation
        )
        request["molecule"].update(
            {
                "analysis_kind": "codon",
                "coding_status": "clean",
                "genetic_code": 1,
                "alignment_strategy": "translate-mafft-backtranslate",
                "trimming_strategy": "protein-mask-to-codons",
            }
        )
        code1_request = self.temp_root / "codon1-request.json"
        code1_request.write_text(json.dumps(request), encoding="utf-8")
        completed, code1_output = self._run(code1_request, output_name="codon1-out")
        self.assertEqual(completed.returncode, 0, completed.stderr)
        commands = self._commands(json.loads((code1_output / "plan.json").read_text()))
        iqtree = commands["infer-accurate-tree"]["argv"]
        self.assertEqual(iqtree[iqtree.index("-st") + 1], "CODON1")

    def test_cds_accepts_reviewed_tblastn_discovery_provenance(self) -> None:
        sequences = {
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
        query, candidates, table, translation = self._write_nucleotide_bundle(
            molecule_type="coding-dna",
            sequences=sequences,
            translations=translations,
        )
        assert translation is not None
        table.write_text(
            table.read_text(encoding="utf-8").replace(
                "local-bundle:not-searched\tnot-searched\tnot-applicable",
                "RefSeq nucleotide release 224\ttblastn\tnucleotide",
            ),
            encoding="utf-8",
        )
        request = self._request_for_bundle(
            "coding-dna", query, candidates, table, translation
        )
        request["molecule"].update(
            {
                "analysis_kind": "codon",
                "coding_status": "clean",
                "genetic_code": 11,
                "alignment_strategy": "translate-mafft-backtranslate",
                "trimming_strategy": "protein-mask-to-codons",
            }
        )
        request_path = self.temp_root / "tblastn-cds-request.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")
        completed, output = self._run(request_path, output_name="tblastn-cds-out")
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue((output / "plan.json").is_file())

    def test_coding_rna_preserves_source_before_translation_and_backtranslation(self) -> None:
        sequences = {
            "QUERY_NT1": "AUG" + "GCU" * 8 + "UAA",
            "REF_NT1": "AUG" + "GCC" * 8 + "UAA",
            "REF_NT2": "AUG" + "GCA" * 8 + "UAA",
            "OUT_NT1": "AUG" + "GUU" * 8 + "UAA",
        }
        translations = {
            "QUERY_NT1": "M" + "A" * 8,
            "REF_NT1": "M" + "A" * 8,
            "REF_NT2": "M" + "A" * 8,
            "OUT_NT1": "M" + "V" * 8,
        }
        query, candidates, table, translation = self._write_nucleotide_bundle(
            molecule_type="coding-rna",
            sequences=sequences,
            translations=translations,
        )
        assert translation is not None
        request = self._request_for_bundle(
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
        request_path = self.temp_root / "coding-rna-request.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")
        completed, output = self._run(request_path)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("U", (output / "reference_set.rna.fasta").read_text())
        analysis_sequences = "".join(
            line
            for line in (output / "reference_set.fna").read_text().splitlines()
            if not line.startswith(">")
        )
        self.assertNotIn("U", analysis_sequences)
        self.assertTrue((output / "reference_set.translated.faa").is_file())
        plan = json.loads((output / "plan.json").read_text())
        self.assertEqual(plan["molecule_plan"]["type"], "coding-rna")
        self.assertEqual(
            plan["artifact_plan"]["backtranslated_alignment"]["path"],
            "alignment.raw.codon.fna",
        )
        operations = {
            item["operation"] for item in plan["molecule_plan"]["transformations"]
        }
        self.assertEqual(
            operations,
            {
                "normalize-u-to-t-derived-copy",
                "use-validated-translations",
                "protein-mask-backtranslation",
            },
        )

    def test_codon_fasttree_and_dna_rna_alphabet_mismatch_fail_closed(self) -> None:
        request = json.loads(
            (ASSETS / "request.nucleotide.example.json").read_text(encoding="utf-8")
        )
        request["tree"].update({"mode": "fast", "tool": "fasttree", "model": "GTR"})
        request["tree"]["support"] = {"method": "sh-like-local", "replicates": 0}
        request["molecule"].update(
            {
                "type": "coding-dna",
                "analysis_kind": "codon",
                "coding_status": "clean",
                "genetic_code": 1,
                "alignment_strategy": "translate-mafft-backtranslate",
                "trimming_strategy": "protein-mask-to-codons",
            }
        )
        request["references"]["translation_fasta"] = "missing.faa"
        request_path = self.temp_root / "fast-codon.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")
        completed, output = self._run(request_path, output_name="fast-codon-out")
        self.assertEqual(completed.returncode, 2)
        self.assertIn("FASTTREE_CODON_MODEL_UNAVAILABLE", completed.stderr)
        self.assertFalse(output.exists())

        bad_query = self.temp_root / "bad-dna.fna"
        bad_query.write_text(">QUERY_NT1\nACGUACGUAA\n", encoding="utf-8")
        dna_request = json.loads(
            (ASSETS / "request.nucleotide.example.json").read_text(encoding="utf-8")
        )
        dna_request["query"]["path"] = str(bad_query)
        dna_request["references"]["candidate_fasta"] = str(
            ASSETS / "candidates.nucleotide.example.fna"
        )
        dna_request["references"]["candidate_table"] = str(
            ASSETS / "candidates.nucleotide.example.tsv"
        )
        dna_request_path = self.temp_root / "bad-dna.json"
        dna_request_path.write_text(json.dumps(dna_request), encoding="utf-8")
        completed, output = self._run(dna_request_path, output_name="bad-dna-out")
        self.assertEqual(completed.returncode, 2)
        self.assertIn("INVALID_NUCLEOTIDE_SEQUENCE", completed.stderr)
        self.assertFalse(output.exists())

    def test_cds_translation_and_nonstandard_code_fail_closed(self) -> None:
        sequences = {
            "QUERY_NT1": "ATG" + "GCT" * 8 + "TAA",
            "REF_NT1": "ATG" + "GCC" * 8 + "TAA",
            "REF_NT2": "ATG" + "GCA" * 8 + "TAA",
            "OUT_NT1": "ATG" + "GTT" * 8 + "TAA",
        }
        translations = {
            "QUERY_NT1": "M" + "A" * 8,
            "REF_NT1": "M" + "A" * 8,
            "REF_NT2": "M" + "A" * 8,
            "OUT_NT1": "M" + "A" * 8,
        }
        query, candidates, table, translation = self._write_nucleotide_bundle(
            molecule_type="coding-dna",
            sequences=sequences,
            translations=translations,
        )
        assert translation is not None
        request = self._request_for_bundle(
            "coding-dna", query, candidates, table, translation
        )
        request["molecule"].update(
            {
                "analysis_kind": "codon",
                "coding_status": "clean",
                "genetic_code": 11,
                "alignment_strategy": "translate-mafft-backtranslate",
                "trimming_strategy": "protein-mask-to-codons",
            }
        )
        request_path = self.temp_root / "translation-mismatch.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")
        completed, output = self._run(request_path, output_name="translation-mismatch-out")
        self.assertEqual(completed.returncode, 2)
        self.assertIn("CDS_TRANSLATION_QC_FAILED", completed.stderr)
        self.assertFalse(output.exists())

        request["molecule"]["genetic_code"] = 5
        request_path = self.temp_root / "nonstandard-code.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")
        completed, output = self._run(request_path, output_name="nonstandard-code-out")
        self.assertEqual(completed.returncode, 2)
        self.assertIn("NONSTANDARD_GENETIC_CODE_ROUTE_REQUIRED", completed.stderr)
        self.assertFalse(output.exists())

    def test_molecule_and_genetic_code_are_decision_bearing(self) -> None:
        first, first_output = self._run(
            ASSETS / "request.nucleotide.example.json", output_name="first"
        )
        self.assertEqual(first.returncode, 0, first.stderr)
        first_plan = json.loads((first_output / "plan.json").read_text())
        changed = json.loads(
            (ASSETS / "request.nucleotide.example.json").read_text(encoding="utf-8")
        )
        changed["molecule"]["type"] = "noncoding-rna"
        changed["molecule"]["source_encoding"] = "rna-u"
        changed["query"]["strand"] = "plus"
        changed["query"]["path"] = str(ASSETS / "query.nucleotide.example.fna")
        changed["references"]["candidate_fasta"] = str(
            ASSETS / "candidates.nucleotide.example.fna"
        )
        changed["references"]["candidate_table"] = str(
            ASSETS / "candidates.nucleotide.example.tsv"
        )
        changed_path = self.temp_root / "changed.json"
        changed_path.write_text(json.dumps(changed), encoding="utf-8")
        completed, _ = self._run(changed_path, output_name="changed-out")
        self.assertEqual(completed.returncode, 2)
        self.assertIn("INVALID_NUCLEOTIDE_SEQUENCE", completed.stderr)
        self.assertTrue(first_plan["molecule_plan"]["no_alphabet_inference"])

    def test_request_03_rejects_unknown_fields_and_invalid_provenance_dates(self) -> None:
        base = json.loads(
            (ASSETS / "request.nucleotide.example.json").read_text(encoding="utf-8")
        )
        mutations = []
        unexpected_root = deepcopy(base)
        unexpected_root["unexpected_top_level"] = True
        mutations.append(("unexpected-root", unexpected_root, "REQUEST_SCHEMA_VIOLATION"))
        unexpected_query = deepcopy(base)
        unexpected_query["query"]["unexpected_query"] = True
        mutations.append(("unexpected-query", unexpected_query, "REQUEST_SCHEMA_VIOLATION"))
        bad_date = deepcopy(base)
        bad_date["query"]["retrieved_at"] = "not-a-date"
        mutations.append(("bad-query-date", bad_date, "INVALID_PROVENANCE_DATE"))
        impossible_date = deepcopy(base)
        impossible_date["query"]["retrieved_at"] = "2026-99-99"
        mutations.append(("impossible-query-date", impossible_date, "INVALID_PROVENANCE_DATE"))
        missing_target_coverage = deepcopy(base)
        del missing_target_coverage["selection"]["min_target_coverage"]
        mutations.append(
            ("missing-target-coverage", missing_target_coverage, "REQUEST_SCHEMA_VIOLATION")
        )
        missing_support = deepcopy(base)
        del missing_support["tree"]["support"]["replicates"]
        del missing_support["tree"]["support"]["sh_alrt"]
        mutations.append(("missing-support", missing_support, "REQUEST_SCHEMA_VIOLATION"))
        duplicate_groups = deepcopy(base)
        duplicate_groups["clustering"]["preserve_analysis_groups"] = [
            "study",
            "study",
            "outgroup",
        ]
        mutations.append(("duplicate-groups", duplicate_groups, "REQUEST_SCHEMA_VIOLATION"))
        incomplete_colors = deepcopy(base)
        incomplete_colors["itol"]["colors"] = {"study": "#E69F00"}
        mutations.append(("incomplete-colors", incomplete_colors, "REQUEST_SCHEMA_VIOLATION"))
        empty_label = deepcopy(base)
        empty_label["itol"]["dataset_label"] = ""
        mutations.append(("empty-itol-label", empty_label, "REQUEST_SCHEMA_VIOLATION"))
        invalid_disabled_taxonomy = deepcopy(base)
        invalid_disabled_taxonomy["taxonomy"]["source"] = "not-ncbi"
        mutations.append(
            ("invalid-disabled-taxonomy", invalid_disabled_taxonomy, "REQUEST_SCHEMA_VIOLATION")
        )
        spaced_original_value = deepcopy(base)
        spaced_original_value["query"]["original_value"] = " QUERY_NT1 "
        mutations.append(
            ("spaced-original-value", spaced_original_value, "REQUEST_SCHEMA_VIOLATION")
        )
        invalid_accurate_model = deepcopy(base)
        invalid_accurate_model["tree"]["model"] = "LG+G4"
        mutations.append(
            ("invalid-accurate-model", invalid_accurate_model, "UNSUPPORTED_TREE_PLAN")
        )
        for optional_object in ("clustering", "trimming", "taxonomy", "itol", "literature"):
            explicit_null = deepcopy(base)
            if optional_object == "trimming":
                explicit_null["molecule"]["trimming_strategy"] = "none"
            explicit_null[optional_object] = None
            mutations.append(
                (
                    f"null-{optional_object}",
                    explicit_null,
                    "REQUEST_SCHEMA_VIOLATION",
                )
            )
        for field_name, mutate in (
            ("schema-version", lambda value: value.update({"schema_version": " 0.3 "})),
            ("objective", lambda value: value.update({"objective": " homolog-context "})),
            (
                "molecule-type",
                lambda value: value["molecule"].update({"type": " noncoding-dna "}),
            ),
            (
                "query-kind",
                lambda value: value["query"].update({"kind": " nucleotide-fasta "}),
            ),
            (
                "reference-strategy",
                lambda value: value["references"].update({"strategy": " local-bundle "}),
            ),
            (
                "alignment-mode",
                lambda value: value["alignment"].update({"mode": " auto "}),
            ),
            (
                "tree-tool",
                lambda value: value["tree"].update({"tool": " iqtree2 "}),
            ),
            (
                "support-method",
                lambda value: value["tree"]["support"].update(
                    {"method": " ultrafast "}
                ),
            ),
        ):
            padded = deepcopy(base)
            mutate(padded)
            mutations.append(
                (f"padded-{field_name}", padded, "REQUEST_SCHEMA_VIOLATION")
            )

        for name, request, reason_code in mutations:
            with self.subTest(mutation=name):
                request_path = self.temp_root / f"{name}.json"
                request_path.write_text(json.dumps(request), encoding="utf-8")
                completed, output = self._run(request_path, output_name=f"out-{name}")
                self.assertEqual(completed.returncode, 2)
                self.assertIn(reason_code, completed.stderr)
                self.assertFalse(output.exists())

    def test_candidate_provenance_and_query_self_row_fail_closed(self) -> None:
        sequences = {
            "QUERY_NT1": "ACGTACGTAA" * 6,
            "REF_NT1": "ACGTACGTAG" * 6,
            "REF_NT2": "ACGTACGCAA" * 6,
            "OUT_NT1": "ACGTTCGTAA" * 6,
        }

        def run_mutation(
            name: str,
            *,
            request_mutator: Any | None = None,
            table_mutator: Any | None = None,
        ) -> subprocess.CompletedProcess[str]:
            query, candidates, table, _ = self._write_nucleotide_bundle(
                molecule_type="noncoding-dna", sequences=sequences
            )
            if table_mutator is not None:
                table.write_text(
                    table_mutator(table.read_text(encoding="utf-8")),
                    encoding="utf-8",
                )
            request = self._request_for_bundle(
                "noncoding-dna", query, candidates, table
            )
            if request_mutator is not None:
                request_mutator(request)
            request_path = self.temp_root / f"{name}.json"
            request_path.write_text(json.dumps(request), encoding="utf-8")
            completed, output = self._run(request_path, output_name=f"out-{name}")
            self.assertFalse(output.exists())
            return completed

        completed = run_mutation(
            "self-organism-mismatch",
            request_mutator=lambda request: request["query"].update(
                {"organism": "Mus musculus"}
            ),
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("QUERY_SELF_PROVENANCE_MISMATCH", completed.stderr)

        completed = run_mutation(
            "invalid-candidate-date",
            table_mutator=lambda text: text.replace(
                "\t2026-08-25\t", "\tnot-a-date\t", 1
            ),
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("INVALID_PROVENANCE_DATE", completed.stderr)

        completed = run_mutation(
            "missing-source-release",
            table_mutator=lambda text: text.replace("\t2026-08\t", "\t\t", 1),
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("MOLECULE_PROVENANCE_REQUIRED", completed.stderr)

        completed = run_mutation(
            "invalid-search-database-pair",
            table_mutator=lambda text: text.replace(
                "\tnot-searched\tnot-applicable", "\tblastp\tnucleotide", 1
            ),
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("SEARCH_DATABASE_PROVENANCE_INVALID", completed.stderr)

        completed = run_mutation(
            "unsupported-candidate-column",
            table_mutator=lambda content: "\n".join(
                f"{line}\tprovider_private_field" for line in content.splitlines()
            )
            + "\n",
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("UNSUPPORTED_CANDIDATE_COLUMNS", completed.stderr)

    def test_request_03_missing_outgroup_rationale_emits_a_blocked_bundle(self) -> None:
        sequences = {
            "QUERY_NT1": "ACGTACGTAA" * 6,
            "REF_NT1": "ACGTACGTAG" * 6,
            "REF_NT2": "ACGTACGCAA" * 6,
            "OUT_NT1": "ACGTTCGTAA" * 6,
        }
        query, candidates, table, _ = self._write_nucleotide_bundle(
            molecule_type="noncoding-dna", sequences=sequences
        )
        table.write_text(
            table.read_text(encoding="utf-8").replace(
                "Declared distant homolog for the synthetic test fixture", ""
            ),
            encoding="utf-8",
        )
        request = self._request_for_bundle(
            "noncoding-dna", query, candidates, table
        )
        request_path = self.temp_root / "missing-outgroup-rationale.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")
        completed, output = self._run(
            request_path, output_name="missing-outgroup-rationale-out"
        )
        self.assertEqual(completed.returncode, 3, completed.stderr)
        self.assertTrue(output.is_dir())
        plan = json.loads((output / "plan.json").read_text(encoding="utf-8"))
        self.assertEqual(plan["state"], "blocked")
        self.assertIn("OUTGROUP_RATIONALE_REQUIRED", plan["hard_stops"])


if __name__ == "__main__":
    unittest.main()
