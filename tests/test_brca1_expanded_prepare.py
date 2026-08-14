"""Contract tests for the manifest-driven expanded BRCA1 preparation helper."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "examples" / "brca1-expanded" / "scripts" / "prepare_brca1_expanded.py"


class Brca1ExpandedPreparationTests(unittest.TestCase):
    def materialize_fixture(self, directory: Path, bad_length: bool = False) -> tuple[Path, Path, Path]:
        fasta = directory / "protein.faa"
        fasta.write_text(
            ">NP_009225.1 BRCA1 [organism=Homo sapiens] [GeneID=672]\n"
            + "A" * 100 + "\n"
            + ">XP_000000001.1 BRCA1 [organism=Species alpha] [GeneID=101] [isoform=X1]\n"
            + "C" * (70 if bad_length else 80) + "\n"
            + ">XP_000000002.1 BRCA1 [organism=Species beta] [GeneID=102] [isoform=X1]\n"
            + "D" * 90 + "\n",
            encoding="utf-8",
        )
        report = directory / "data_report.jsonl"
        records = [
            {"geneId": "672", "taxId": "9606", "taxname": "Homo sapiens", "symbol": "BRCA1", "annotations": []},
            {"geneId": "101", "taxId": "1001", "taxname": "Species alpha", "symbol": "BRCA1", "annotations": []},
            {"geneId": "102", "taxId": "1002", "taxname": "Species beta", "symbol": "BRCA1", "annotations": []},
        ]
        report.write_text("".join(json.dumps(row) + "\n" for row in records), encoding="utf-8")
        manifest = directory / "selection_manifest.tsv"
        manifest.write_text(
            "accession\ttaxon_id\tspecies\tgene_id\tanalysis_group\tclade\tis_reviewed\tselection_rationale\n"
            "NP_009225.1\t9606\tHomo sapiens\t672\tstudy\tPrimates\ttrue\tfocal\n"
            "XP_000000001.1\t1001\tSpecies alpha\t101\toutgroup\tAlpha\tfalse\tfirst outgroup\n"
            "XP_000000002.1\t1002\tSpecies beta\t102\toutgroup\tBeta\tfalse\tsecond outgroup\n",
            encoding="utf-8",
        )
        return fasta, report, manifest

    def run_helper(self, directory: Path, bad_length: bool = False) -> subprocess.CompletedProcess[str]:
        fasta, report, manifest = self.materialize_fixture(directory, bad_length=bad_length)
        return subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--protein-fasta", str(fasta),
                "--data-report", str(report),
                "--selection-manifest", str(manifest),
                "--externally-verified-provider-archive-sha256", "a" * 64,
                "--out", str(directory / "output"),
            ],
            text=True,
            capture_output=True,
            check=False,
        )

    def test_manifest_drives_outputs_and_full_provider_inventory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            completed = self.run_helper(directory)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            output = directory / "output"
            self.assertEqual(
                (output / "candidates.faa").read_text(encoding="utf-8").count(">"),
                3,
            )
            with (output / "candidates.pre-qc.tsv").open(encoding="utf-8", newline="") as handle:
                candidates = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["accession"] for row in candidates], [
                "NP_009225.1", "XP_000000001.1", "XP_000000002.1"
            ])
            self.assertEqual(candidates[1]["relation"], "ortholog")
            self.assertIn("pair-specific evidence", candidates[1]["orthology_evidence"])
            with (output / "provider_species_inventory.tsv").open(encoding="utf-8", newline="") as handle:
                inventory = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(inventory), 3)
            summary = json.loads(
                (output / "provider_scope_summary.json").read_text(encoding="utf-8")
            )
            self.assertEqual(summary["provider_gene_record_count"], 3)
            self.assertEqual(summary["provider_protein_record_count"], 3)
            self.assertEqual(summary["selected_tip_count"], 3)

    def test_length_gate_fails_before_output_materialization(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            completed = self.run_helper(directory, bad_length=True)
            self.assertEqual(completed.returncode, 2)
            self.assertIn("out-of-policy length ratio", completed.stderr)
            self.assertFalse((directory / "output").exists())

    def test_rejection_ledger_requires_a_selected_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            fasta, report, manifest = self.materialize_fixture(directory)
            fasta.write_text(
                fasta.read_text(encoding="utf-8")
                + ">XP_000000003.1 BRCA1 [organism=Species gamma] [GeneID=103]\n"
                + "E" * 85 + "\n",
                encoding="utf-8",
            )
            report.write_text(
                report.read_text(encoding="utf-8")
                + json.dumps({
                    "geneId": "103", "taxId": "1003", "taxname": "Species gamma",
                    "symbol": "BRCA1", "annotations": [],
                }) + "\n",
                encoding="utf-8",
            )
            rejected = directory / "rejected.tsv"
            rejected.write_text(
                "accession\ttaxon_id\tspecies\tgene_id\tstage\tdecision\treason\t"
                "replacement_accession\tevidence\n"
                "XP_000000003.1\t1003\tSpecies gamma\t103\tdomain_qc\trejected\t"
                "missing repeat\tXP_999999999.1\tsynthetic evidence\n",
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    sys.executable, str(SCRIPT),
                    "--protein-fasta", str(fasta),
                    "--data-report", str(report),
                    "--selection-manifest", str(manifest),
                    "--rejected-candidates", str(rejected),
                    "--externally-verified-provider-archive-sha256", "a" * 64,
                    "--out", str(directory / "output"),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn("replacement that is not selected", completed.stderr)
            self.assertFalse((directory / "output").exists())

    def test_committed_review_manifest_is_balanced_and_exact(self) -> None:
        manifest = ROOT / "examples" / "brca1-expanded" / "inputs" / "selection_manifest.tsv"
        with manifest.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 50)
        self.assertEqual(len({row["accession"] for row in rows}), 50)
        self.assertEqual(len({row["taxon_id"] for row in rows}), 50)
        self.assertEqual(sum(row["analysis_group"] == "study" for row in rows), 1)
        self.assertEqual(sum(row["analysis_group"] == "expanded" for row in rows), 44)
        self.assertEqual(sum(row["analysis_group"] == "outgroup" for row in rows), 5)
        self.assertEqual(
            {row["clade"] for row in rows if row["analysis_group"] == "outgroup"},
            {"Anura", "Caudata", "Gymnophiona"},
        )

        example = ROOT / "examples" / "brca1-expanded" / "inputs"
        with (example / "candidates.pre-qc.tsv").open(encoding="utf-8", newline="") as handle:
            candidates = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(
            [row["accession"] for row in candidates],
            [row["accession"] for row in rows],
        )
        self.assertEqual(
            sum(
                1
                for line in (example / "candidates.faa")
                .read_text(encoding="utf-8")
                .splitlines()
                if line.startswith(">")
            ),
            50,
        )

    def test_committed_provider_inventory_and_declared_gaps_are_honest(self) -> None:
        example = ROOT / "examples" / "brca1-expanded" / "inputs"
        with (example / "provider_species_inventory.tsv").open(encoding="utf-8", newline="") as handle:
            inventory = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(inventory), 558)
        self.assertEqual(sum(int(row["provider_protein_count"]) for row in inventory), 3605)
        self.assertEqual(sum(int(row["within_declared_length_gate_count"]) for row in inventory), 2543)
        self.assertEqual(
            Counter(row["selection_status"] for row in inventory),
            Counter({
                "not_selected": 507,
                "selected_for_reference_review": 50,
                "rejected_during_reference_review": 1,
            }),
        )
        no_length = [row for row in inventory if row["reason_code"] == "NO_SEQUENCE_WITHIN_LENGTH_GATE"]
        self.assertEqual(len(no_length), 46)
        self.assertTrue(all("no terminal-domain failure" in row["evidence_note"] for row in no_length))

        with (example / "sampling_gaps.tsv").open(encoding="utf-8", newline="") as handle:
            gaps = {row["species"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(set(gaps), {
            "Ornithorhynchus anatinus", "Tachyglossus aculeatus", "Sphenodon punctatus"
        })
        self.assertEqual(gaps["Sphenodon punctatus"]["provider_status"], "absent")
        self.assertEqual(gaps["Ornithorhynchus anatinus"]["gap_status"], "no_sequence_within_length_gate")

        ascaphus = next(row for row in inventory if row["species"] == "Ascaphus truei")
        self.assertEqual(ascaphus["rejected_accession"], "XP_075436319.1")
        self.assertEqual(ascaphus["reason_code"], "TERMINAL_DOMAIN_QC_FAILED")
        self.assertIn("one distinct C-terminal BRCT", ascaphus["evidence_note"])

    def test_committed_promoted_qc_and_taxonomy_match_the_manifest(self) -> None:
        example = ROOT / "examples" / "brca1-expanded"

        def rows(path: Path) -> list[dict[str, str]]:
            with path.open(encoding="utf-8", newline="") as handle:
                return list(csv.DictReader(handle, delimiter="\t"))

        manifest = rows(example / "inputs" / "selection_manifest.tsv")
        expected = [row["accession"] for row in manifest]
        expected_taxa = {
            row["accession"]: (row["species"], row["taxon_id"])
            for row in manifest
        }
        provenance = rows(example / "inputs" / "candidate_provenance.tsv")
        candidates = rows(example / "inputs" / "candidates.tsv")
        qc = rows(example / "qc" / "candidate_qc.tsv")
        taxonomy = rows(example / "annotation" / "taxonomy_resolution.tsv")
        lineage = rows(example / "annotation" / "taxonomy_lineage.tsv")

        self.assertEqual([row["accession"] for row in provenance], expected)
        self.assertEqual([row["accession"] for row in candidates], expected)
        self.assertEqual([row["accession"] for row in qc], expected)
        self.assertEqual([row["record_id"] for row in taxonomy], expected)
        self.assertEqual({row["accession"] for row in lineage}, set(expected))
        self.assertTrue(all(row["domain_qc_status"] == "pass" for row in qc))
        self.assertTrue(all(int(row["n_terminal_ring_signature_hits"]) >= 1 for row in qc))
        self.assertTrue(all(int(row["brct_repeat_count"]) >= 2 for row in qc))
        for row in taxonomy:
            species, taxon_id = expected_taxa[row["record_id"]]
            self.assertEqual(row["input_name"], species)
            self.assertEqual(row["matched_name"], species)
            self.assertEqual(row["requested_taxon_id"], taxon_id)
            self.assertEqual(row["taxon_id"], taxon_id)
            self.assertEqual(row["status"], "resolved-exact-scientific-name")
        for row in lineage:
            species, taxon_id = expected_taxa[row["accession"]]
            self.assertEqual(row["scientific_name"], species)
            self.assertEqual(row["taxon_id"], taxon_id)
            self.assertEqual(
                row["artifact_role"],
                "post-plan-non-decision-bearing-taxonomy-evidence",
            )
            self.assertGreater(int(row["lineage_node_count"]), 1)

        rejected = rows(example / "inputs" / "rejected_candidates.tsv")
        self.assertEqual(len(rejected), 1)
        self.assertEqual(rejected[0]["accession"], "XP_075436319.1")
        self.assertEqual(rejected[0]["replacement_accession"], "XP_053555561.1")
        self.assertNotIn(rejected[0]["accession"], expected)
        self.assertIn(rejected[0]["replacement_accession"], expected)

        flags = rows(example / "qc" / "manual_review_flags.tsv")
        self.assertEqual({row["accession"] for row in flags}, {
            "XP_019406054.1",
            "XP_034272113.1",
            "NP_001107963.1",
            "XP_053555561.1",
            "XP_026576759.1",
        })
        self.assertTrue(all(row["status"] == "pending-manual-approval" for row in flags))
        bombina = next(row for row in flags if row["accession"] == "XP_053555561.1")
        self.assertIn("SMART-only", bombina["evidence_interpretation"])

        snapshot = json.loads(
            (example / "annotation" / "taxonomy_snapshot.json").read_text(encoding="utf-8")
        )
        self.assertEqual(snapshot["records_resolved"], 50)
        self.assertEqual(
            snapshot["archive_sha256"],
            "e67d56c3e87eac14a28feea8cf710ac612571f00a4f32ba3b85bae0020cd123a",
        )

    def test_committed_qc_checksums_are_portable_and_complete(self) -> None:
        checksum_path = ROOT / "examples" / "brca1-expanded" / "qc" / "checksums.sha256"
        records: dict[str, str] = {}
        for line in checksum_path.read_text(encoding="utf-8").splitlines():
            digest, relative = line.split("  ", 1)
            self.assertRegex(digest, r"^[0-9a-f]{64}$")
            self.assertTrue(relative.startswith("examples/brca1-expanded/"))
            target = ROOT / relative
            self.assertTrue(target.is_file(), relative)
            actual = hashlib.sha256(target.read_bytes()).hexdigest()
            self.assertEqual(actual, digest, relative)
            records[relative] = digest
        self.assertEqual(set(records), {
            "examples/brca1-expanded/inputs/candidates.faa",
            "examples/brca1-expanded/inputs/candidates.pre-qc.tsv",
            "examples/brca1-expanded/inputs/candidates.tsv",
            "examples/brca1-expanded/inputs/query.faa",
            "examples/brca1-expanded/inputs/rejected_candidates.tsv",
            "examples/brca1-expanded/inputs/request.json",
            "examples/brca1-expanded/inputs/selection_manifest.tsv",
            "examples/brca1-expanded/qc/blastp-human-vs-candidates.tsv",
            "examples/brca1-expanded/qc/interproscan-pfam-smart.tsv",
            "examples/brca1-expanded/qc/candidate_qc.tsv",
            "examples/brca1-expanded/qc/manual_review_flags.tsv",
            "examples/brca1-expanded/qc/rejected-attempt/checksums.sha256",
            "examples/brca1-expanded/qc/software_versions.tsv",
        })

        attempt = ROOT / "examples" / "brca1-expanded" / "qc" / "rejected-attempt"
        attempt_records = {}
        for line in (attempt / "checksums.sha256").read_text(encoding="utf-8").splitlines():
            digest, relative = line.split("  ", 1)
            target = ROOT / relative
            self.assertTrue(target.is_file(), relative)
            self.assertEqual(hashlib.sha256(target.read_bytes()).hexdigest(), digest, relative)
            attempt_records[relative] = digest
        self.assertIn("examples/brca1-expanded/qc/rejected-attempt/candidates.faa", attempt_records)
        self.assertIn(
            "examples/brca1-expanded/qc/rejected-attempt/interproscan-pfam-smart.tsv",
            attempt_records,
        )

        raw_interpro = attempt / "interproscan-pfam-smart.tsv"
        ascaphus_intervals = []
        with raw_interpro.open(encoding="utf-8", newline="") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if row[0] != "XP_075436319.1":
                    continue
                if row[4] in {"PF00533", "SM00292"}:
                    ascaphus_intervals.append((int(row[6]), int(row[7])))
        self.assertEqual(len(ascaphus_intervals), 2)
        self.assertLessEqual(max(start for start, _ in ascaphus_intervals), min(
            end for _, end in ascaphus_intervals
        ))

    def test_reference_review_is_pending_and_hash_closed(self) -> None:
        example = ROOT / "examples" / "brca1-expanded"
        review = example / "review"
        plan = json.loads((review / "plan.json").read_text(encoding="utf-8"))
        manifest = json.loads((review / "manifest.json").read_text(encoding="utf-8"))

        self.assertEqual(plan["state"], "pending-reference-approval")
        self.assertEqual(manifest["workflow_state"], "pending-reference-approval")
        self.assertEqual(plan["run_id"], "gtr-9eaf4f8541020b3b")
        self.assertEqual(plan["run_id"], manifest["run_id"])
        self.assertEqual(plan["plan_hash"], manifest["plan_hash"])
        self.assertIsNone(manifest["approved_plan_hash"])
        self.assertFalse((example / "reference-approval.tsv").exists())
        self.assertEqual(
            plan["candidate_and_selection_counts"]["selected_records_including_query"],
            50,
        )
        self.assertEqual(len(plan["outgroup_accessions"]), 5)
        self.assertEqual(plan["clustering_plan"]["status"], "not-needed")
        self.assertFalse(plan["clustering_plan"]["triggered"])
        self.assertTrue(all(not command["executed"] for command in manifest["commands"]))

        for artifact in manifest["input_artifacts"]:
            role = artifact["role"]
            if role in {"taxonomy_names", "taxonomy_nodes"}:
                snapshot = json.loads(
                    (example / "annotation" / "taxonomy_snapshot.json").read_text(
                        encoding="utf-8"
                    )
                )
                key = "names_dmp_sha256" if role == "taxonomy_names" else "nodes_dmp_sha256"
                self.assertEqual(artifact["sha256"], snapshot[key])
                continue
            target = example / artifact["logical_path"]
            self.assertTrue(target.is_file(), artifact["logical_path"])
            self.assertEqual(
                hashlib.sha256(target.read_bytes()).hexdigest(),
                artifact["sha256"],
                artifact["logical_path"],
            )

        for artifact in manifest["output_artifacts"]:
            target = review / artifact["logical_path"]
            self.assertTrue(target.is_file(), artifact["logical_path"])
            self.assertEqual(target.stat().st_size, artifact["bytes"])
            self.assertEqual(
                hashlib.sha256(target.read_bytes()).hexdigest(),
                artifact["sha256"],
                artifact["logical_path"],
            )

    def test_public_reference_review_bundle_is_closed_and_private_path_free(self) -> None:
        example = ROOT / "examples" / "brca1-expanded"
        checksum_path = example / "report" / "checksums.sha256"
        expected_files = {
            path.relative_to(ROOT).as_posix()
            for path in example.rglob("*")
            if path.is_file()
            and "__pycache__" not in path.parts
            and path != checksum_path
        }
        listed_files = set()
        for line in checksum_path.read_text(encoding="utf-8").splitlines():
            digest, relative = line.split("  ", 1)
            target = ROOT / relative
            self.assertTrue(target.is_file(), relative)
            self.assertEqual(hashlib.sha256(target.read_bytes()).hexdigest(), digest, relative)
            listed_files.add(relative)
        self.assertEqual(listed_files, expected_files)

        for markdown in example.rglob("*.md"):
            for destination in re.findall(
                r"\[[^]]+\]\(([^)]+)\)",
                markdown.read_text(encoding="utf-8"),
            ):
                if destination.startswith(("http://", "https://", "mailto:", "#")):
                    continue
                local = destination.split("#", 1)[0]
                self.assertTrue((markdown.parent / local).resolve().exists(), (
                    markdown.relative_to(ROOT), destination
                ))

        private = re.compile(
            r"/(?:Users|home|user[0-9]+|lustre|scratch|gpfs)/|"
            r"(?:private-user|zhaohongda|hongda@)|\b[0-9]+\.fe3-adm\b",
            flags=re.IGNORECASE,
        )
        text_suffixes = {".md", ".json", ".tsv", ".txt", ".faa", ".pbs", ".py"}
        for path in example.rglob("*"):
            if path.is_file() and path.suffix.lower() in text_suffixes:
                self.assertIsNone(private.search(path.read_text(encoding="utf-8")), path)


if __name__ == "__main__":
    unittest.main()
