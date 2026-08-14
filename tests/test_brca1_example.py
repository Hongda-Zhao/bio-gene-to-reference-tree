"""Integrity checks for the executed public BRCA1 worked example."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import unittest
from pathlib import Path
from pathlib import PurePosixPath


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = REPOSITORY_ROOT / "examples" / "brca1"
EXPECTED_TIPS = {
    "NP_009225.1", "NP_001038958.1", "NP_033894.3", "XP_017204702.2",
    "NP_001013434.1", "NP_848668.1", "XP_014595447.2", "XP_003414318.3",
    "XP_058140293.1", "NP_001029141.1", "NP_989500.1", "XP_072775070.1",
    "XP_019406054.1", "XP_023967135.2", "XP_008111382.1", "XP_026576759.1",
    "NP_001107963.1", "XP_029429046.1",
}
EXPECTED_CLADES = {
    "Amphibia", "Mammalia", "Sauropsida", "Primates", "Glires", "Aves",
    "Lepidosauria", "Archelosauria", "Archosauria", "Crocodylia+Testudines",
}

FINAL_UNROOTED_TREE = EXAMPLE / "tree" / "gene-tree.unrooted.nwk"
FINAL_ROOTED_TREE = EXAMPLE / "tree" / "gene-tree.outgroup-rooted.nwk"
FINALIZATION_RECEIPT = EXAMPLE / "rooting-decision.tsv"
FINAL_PUBLICATION_ARTIFACTS = (
    EXAMPLE / "DATA_LICENSES.md",
    FINALIZATION_RECEIPT,
    FINAL_UNROOTED_TREE,
    FINAL_ROOTED_TREE,
    EXAMPLE / "tree" / "gene-tree.topology_sensitivity.tsv",
    EXAMPLE / "figures" / "gene-tree.outgroup-rooted.ggtree.svg",
    EXAMPLE / "figures" / "gene-tree.outgroup-rooted.ggtree.pdf",
    EXAMPLE / "figures" / "gene-tree.outgroup-rooted.ggtree.settings.tsv",
    EXAMPLE / "tree" / "clade_support.tsv",
    EXAMPLE / "tree" / "iqtree" / "brca1-balanced.log",
    EXAMPLE / "tree" / "iqtree" / "brca1-balanced.iqtree",
    EXAMPLE / "tree" / "iqtree" / "brca1-balanced.treefile",
    EXAMPLE / "tree" / "iqtree" / "brca1-balanced.contree",
    EXAMPLE / "tree" / "iqtree" / "brca1-balanced.ufboot",
    EXAMPLE / "tree" / "iqtree" / "redactions.tsv",
    EXAMPLE / "tree" / "iqtree" / "raw-output-receipt.json",
    EXAMPLE / "tree" / "iqtree" / "software_version.txt",
    EXAMPLE / "annotation" / "taxonomy_lineage.tsv",
    EXAMPLE / "report" / "report.md",
    EXAMPLE / "report" / "commands.jsonl",
    EXAMPLE / "report" / "qc_promotion_receipt.json",
    EXAMPLE / "report" / "execution_reconciliation.json",
    EXAMPLE / "report" / "execution_reconciliation.schema.json",
    EXAMPLE / "report" / "software_versions.tsv",
    EXAMPLE / "report" / "checksums.sha256",
)

# Final reports from IQ-TREE and the local post-processing scripts can contain
# the command line used to create them.  Scan every public text/report format
# that could retain an execution path, rather than checking Markdown alone.
PUBLIC_TEXT_SUFFIXES = {
    ".bib", ".bionj", ".contree", ".csv", ".faa", ".iqtree", ".json",
    ".jsonl", ".log", ".md", ".mldist", ".nex", ".nwk", ".pbs", ".py",
    ".r", ".sha256", ".suptree", ".svg", ".treefile", ".tsv", ".txt",
    ".ufboot", ".yaml", ".yml",
}
PRIVATE_EXECUTION_PATTERNS = {
    "private absolute storage/home path": re.compile(
        r"(?<![A-Za-z0-9/])/(?:Users|home|lustre|scratch|user\d+|gpfs)(?:/|\b)",
        re.IGNORECASE,
    ),
    "private file URI": re.compile(
        r"\bfile:/+(?:Users|home|lustre|scratch|user\d+|gpfs)(?:/|\b)",
        re.IGNORECASE,
    ),
    "execution-host login address": re.compile(
        r"(?<![A-Za-z0-9_.-])[A-Za-z][A-Za-z0-9._-]*@(?:gds2|[A-Za-z0-9.-]*-adm)\b",
        re.IGNORECASE,
    ),
    "qualified scheduler job ID": re.compile(
        r"(?<![A-Za-z0-9_.])\d{4,}\.[A-Za-z][A-Za-z0-9_.-]*\b",
    ),
    "contextual scheduler job ID": re.compile(
        r"\b(?:PBS\s+)?job(?:\s+ID)?\s*[:#=]?\s*\d{4,}\b",
        re.IGNORECASE,
    ),
    "scheduler command job ID": re.compile(
        r"\b(?:qstat|qdel|tracejob)\s+(?:-[A-Za-z]+\s+)*\d{4,}\b",
        re.IGNORECASE,
    ),
}
SHA256_LINE = re.compile(
    r"^(?P<digest>[0-9a-f]{64}) (?P<mode>[ *])(?P<path>.+)$"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_records(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier = ""
    chunks: list[str] = []

    def store() -> None:
        if identifier:
            if identifier in records:
                raise AssertionError(f"duplicate FASTA identifier: {path}: {identifier}")
            records[identifier] = "".join(chunks)

    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            store()
            identifier = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line)
    store()
    return records


def receipt(path: Path) -> dict[str, str]:
    pairs = [line.split("\t", 1) for line in path.read_text(encoding="utf-8").splitlines()]
    if any(len(pair) != 2 or not pair[0] for pair in pairs):
        raise AssertionError(f"invalid key/value receipt: {path}")
    return dict(pairs)


def brca1_newick_tips(path: Path) -> list[str]:
    """Extract this fixture's unquoted accession tips without parsing supports."""

    text = path.read_text(encoding="utf-8")
    if text.count(";") != 1 or not text.rstrip().endswith(";"):
        raise AssertionError(f"expected exactly one complete Newick tree: {path}")
    # Newick comments may contain arbitrary provenance text.  They are not tip
    # labels and must not be allowed to satisfy the accession contract.
    without_comments = re.sub(r"\[[^\[\]]*\]", "", text)
    if "[" in without_comments or "]" in without_comments:
        raise AssertionError(f"nested or malformed Newick comment: {path}")
    # All fixed BRCA1 labels are unquoted RefSeq accessions and every final tree
    # carries branch lengths.  A tip therefore occurs only after '(' or ',' and
    # before ':'.  Internal SH-aLRT/UFBoot labels occur after ')' and cannot be
    # mistaken for tips by this fixture-specific expression.
    return re.findall(r"(?<=[(,])\s*([^():,;'\s]+)\s*(?=:)", without_comments)


def public_text_artifacts() -> list[Path]:
    paths = [REPOSITORY_ROOT / "README.md"]
    paths.extend(
        path
        for path in EXAMPLE.rglob("*")
        if path.is_file() and path.suffix.lower() in PUBLIC_TEXT_SUFFIXES
    )
    return sorted(set(paths))


def command_records() -> list[dict[str, object]]:
    """Load the JSONL ledger while rejecting duplicate keys."""

    def reject_duplicate_pairs(pairs: list[tuple[str, object]]) -> dict[str, object]:
        record: dict[str, object] = {}
        for key, value in pairs:
            if key in record:
                raise ValueError(f"duplicate JSON key: {key}")
            record[key] = value
        return record

    lines = (EXAMPLE / "report" / "commands.jsonl").read_text(encoding="utf-8").splitlines()
    if not lines:
        raise AssertionError("commands.jsonl must contain command records")
    records: list[dict[str, object]] = []
    for line_number, line in enumerate(lines, start=1):
        if not line.strip():
            raise AssertionError(f"blank commands.jsonl line: {line_number}")
        record = json.loads(line, object_pairs_hook=reject_duplicate_pairs)
        if not isinstance(record, dict) or not record:
            raise AssertionError(f"invalid command record at line {line_number}")
        records.append(record)
    return records


class Brca1WorkedExampleTests(unittest.TestCase):
    def test_fixed_tip_set_is_identical_across_materialized_stages(self) -> None:
        fasta_paths = (
            EXAMPLE / "inputs" / "candidates.faa",
            EXAMPLE / "review" / "reference_set.faa",
            EXAMPLE / "alignment" / "alignment.raw.faa",
            EXAMPLE / "alignment" / "alignment.trimmed.permissive.faa",
            EXAMPLE / "alignment" / "alignment.trimmed.balanced.faa",
            EXAMPLE / "alignment" / "alignment.trimmed.strict.faa",
        )
        for path in fasta_paths:
            self.assertEqual(set(fasta_records(path)), EXPECTED_TIPS, str(path))

        table_keys = (
            (EXAMPLE / "inputs" / "candidate_provenance.tsv", "accession"),
            (EXAMPLE / "inputs" / "candidates.tsv", "accession"),
            (EXAMPLE / "qc" / "candidate_qc.tsv", "accession"),
            (EXAMPLE / "review" / "taxonomy_resolution.tsv", "record_id"),
        )
        for path, key in table_keys:
            self.assertEqual({row[key] for row in read_tsv(path)}, EXPECTED_TIPS, str(path))

        metadata = read_tsv(EXAMPLE / "review" / "sequence_metadata.tsv")
        selected = {row["tip_id"] for row in metadata if row["inclusion_status"] == "selected"}
        rejected = {row["tip_id"] for row in metadata if row["inclusion_status"] == "rejected"}
        self.assertEqual(selected, EXPECTED_TIPS)
        self.assertEqual(rejected, set())

    def test_taxonomy_terminal_domains_and_trimmed_domains_all_pass(self) -> None:
        taxonomy = read_tsv(EXAMPLE / "review" / "taxonomy_resolution.tsv")
        self.assertEqual(len(taxonomy), 18)
        self.assertTrue(all(row["status"] == "resolved-exact-scientific-name" for row in taxonomy))
        self.assertTrue(all(row["name_class"] == "scientific name" for row in taxonomy))
        self.assertTrue(all(row["requested_taxon_id"] == row["taxon_id"] for row in taxonomy))

        candidate_qc = read_tsv(EXAMPLE / "qc" / "candidate_qc.tsv")
        self.assertEqual(len(candidate_qc), 18)
        self.assertTrue(all(row["domain_qc_status"] == "pass" for row in candidate_qc))
        self.assertTrue(all(int(row["brct_repeat_count"]) >= 2 for row in candidate_qc))

        retention = read_tsv(EXAMPLE / "alignment" / "qc" / "domain_retention.tsv")
        self.assertEqual(len(retention), 18 * 3 * 3)
        self.assertTrue(all(row["status"] == "pass" for row in retention))
        self.assertGreaterEqual(min(float(row["retained_fraction"]) for row in retention), 0.8)

    def test_alignment_metrics_and_approvals_are_hash_bound(self) -> None:
        metrics = {
            row["alignment_id"]: row
            for row in read_tsv(EXAMPLE / "alignment" / "qc" / "alignment_qc.tsv")
        }
        self.assertEqual(
            {name: int(row["columns"]) for name, row in metrics.items()},
            {"raw": 2179, "permissive": 2046, "balanced": 1851, "strict": 1477},
        )
        plan = json.loads((EXAMPLE / "review" / "plan.json").read_text(encoding="utf-8"))
        reference_approval = receipt(EXAMPLE / "reference-approval.tsv")
        alignment_approval = receipt(EXAMPLE / "alignment-approval.tsv")
        self.assertEqual(reference_approval["run_id"], plan["run_id"])
        self.assertEqual(reference_approval["plan_hash"], plan["plan_hash"])

        for path, expected in (
            (EXAMPLE / "review" / "reference_set.faa", reference_approval["reference_set_sha256"]),
            (
                EXAMPLE / "alignment" / "alignment.trimmed.balanced.faa",
                alignment_approval["alignment_sha256"],
            ),
        ):
            observed = hashlib.sha256(path.read_bytes()).hexdigest()
            self.assertEqual(observed, expected)

        manifest = json.loads(
            (EXAMPLE / "review" / "manifest.json").read_text(encoding="utf-8")
        )
        public_inputs = {
            artifact["logical_path"]: artifact
            for artifact in manifest["input_artifacts"]
            if artifact["role"] not in {"taxonomy_names", "taxonomy_nodes"}
        }
        self.assertEqual(
            set(public_inputs),
            {
                "inputs/request.json",
                "inputs/query.faa",
                "inputs/candidates.faa",
                "inputs/candidates.tsv",
            },
        )
        for logical_path, artifact in public_inputs.items():
            path = EXAMPLE / logical_path
            self.assertEqual(
                hashlib.sha256(path.read_bytes()).hexdigest(),
                artifact["sha256"],
                logical_path,
            )

        taxonomy_inputs = {
            artifact["role"]: artifact["sha256"]
            for artifact in manifest["input_artifacts"]
            if artifact["role"] in {"taxonomy_names", "taxonomy_nodes"}
        }
        self.assertEqual(
            taxonomy_inputs,
            {
                "taxonomy_names": plan["taxonomy_plan"]["names_sha256"],
                "taxonomy_nodes": plan["taxonomy_plan"]["nodes_sha256"],
            },
        )

        for artifact in manifest["output_artifacts"]:
            path = EXAMPLE / "review" / artifact["logical_path"]
            data = path.read_bytes()
            self.assertEqual(len(data), artifact["bytes"], artifact["logical_path"])
            self.assertEqual(
                hashlib.sha256(data).hexdigest(),
                artifact["sha256"],
                artifact["logical_path"],
            )

        iqtree_pbs = (EXAMPLE / "scripts" / "run_iqtree.pbs").read_text(encoding="utf-8")
        self.assertIn(alignment_approval["alignment_sha256"], iqtree_pbs)
        domain_approval = receipt(EXAMPLE / "domain-retention-approval.tsv")
        self.assertIn(domain_approval["audit_sha256"], iqtree_pbs)
        self.assertIn("domain-retention-approval.tsv", iqtree_pbs)
        self.assertIn("decision_ok", iqtree_pbs)
        self.assertIn("#PBS -l walltime=12:00:00", iqtree_pbs)

    def test_command_ledger_has_explicit_replay_semantics(self) -> None:
        command_ids: set[str] = set()
        for line_number, record in enumerate(command_records(), start=1):
            self.assertIsInstance(record.get("id"), str)
            self.assertTrue(record["id"])
            self.assertNotIn(record["id"], command_ids, f"duplicate command id: {record['id']}")
            command_ids.add(record["id"])

            argv = record.get("argv")
            self.assertIsInstance(argv, list, f"argv is not an array: {line_number}")
            self.assertTrue(argv, f"argv is empty: {line_number}")
            self.assertTrue(
                all(isinstance(value, str) and value for value in argv),
                f"argv must contain only non-empty strings: {line_number}",
            )
            record_type = record.get("command_record_type")
            self.assertIn(
                record_type,
                {"normalized-replay-argv", "normalized-provenance-argv"},
            )
            cwd = record.get("cwd")
            self.assertIsInstance(cwd, str, f"cwd must be explicit: {line_number}")
            self.assertTrue(cwd, f"cwd is empty: {line_number}")
            tokenized = any("<" in value or ">" in value for value in argv) or (
                "<" in cwd or ">" in cwd
            )
            if record_type == "normalized-provenance-argv":
                self.assertTrue(
                    tokenized,
                    f"provenance argv lacks an explicit path token: {line_number}",
                )
            else:
                self.assertFalse(
                    tokenized,
                    f"replay argv contains an unresolved path token: {line_number}",
                )
                self.assertEqual(
                    cwd,
                    "examples/brca1",
                    f"replay argv has an unexpected base directory: {line_number}",
                )

    def test_example_markdown_local_links_resolve(self) -> None:
        pattern = re.compile(r"\[[^\]]+\]\(([^)]+)\)")
        documents = [EXAMPLE / "README.md", EXAMPLE / "PROVENANCE.md"]
        if receipt(FINALIZATION_RECEIPT).get("decision_status") != "pending-final-tree":
            documents.extend(
                (REPOSITORY_ROOT / "README.md", EXAMPLE / "report" / "report.md")
            )
        for document in documents:
            for target in pattern.findall(document.read_text(encoding="utf-8")):
                if target.startswith(("http://", "https://", "mailto:", "#")):
                    continue
                path_part = target.split("#", 1)[0]
                self.assertTrue((document.parent / path_part).exists(), f"{document}: {target}")

    def test_current_public_text_artifacts_do_not_leak_private_execution_state(self) -> None:
        """Run independently of the asynchronous final-tree publication gate."""

        for path in public_text_artifacts():
            try:
                text = path.read_text(encoding="utf-8")
            except UnicodeDecodeError as error:
                self.fail(f"declared public text artifact is not UTF-8: {path}: {error}")
            for label, pattern in PRIVATE_EXECUTION_PATTERNS.items():
                match = pattern.search(text)
                if match is not None:
                    line_number = text.count("\n", 0, match.start()) + 1
                    self.fail(
                        f"{label} leaked in public artifact "
                        f"{path.relative_to(REPOSITORY_ROOT)}:{line_number}"
                    )

    def test_final_publication_bundle_is_closed_and_portable(self) -> None:
        """Fail closed once the explicit final rooting receipt is finalized.

        The IQ-TREE job and downstream rendering run asynchronously, so the
        neutral receipt is the publication-stage sentinel. While it remains
        pending this one gate is explicitly skipped; once finalized, every
        declared result artifact and audit becomes mandatory in the same run.
        """

        rooting = receipt(FINALIZATION_RECEIPT)
        if rooting.get("decision_status") == "pending-final-tree":
            if os.environ.get("CI"):
                self.fail(
                    "the public BRCA1 homepage example cannot enter CI with a "
                    "pending final-tree receipt"
                )
            self.skipTest(
                "final BRCA1 rooting decision is pending; finalizing "
                "examples/brca1/rooting-decision.tsv activates the complete "
                "publication gate"
            )

        self.assertEqual(rooting.get("decision_status"), "approved-outgroup-rooted")
        self.assertEqual(rooting.get("rooted_tree_created"), "true")
        self.assertRegex(rooting.get("primary_tree_sha256", ""), r"^[0-9a-f]{64}$")

        for path in FINAL_PUBLICATION_ARTIFACTS:
            with self.subTest(artifact=path.relative_to(REPOSITORY_ROOT)):
                self.assertTrue(path.is_file(), f"missing final artifact: {path}")
                self.assertGreater(path.stat().st_size, 0, f"empty final artifact: {path}")

        pdf_path = EXAMPLE / "figures" / "gene-tree.outgroup-rooted.ggtree.pdf"
        pdf_bytes = pdf_path.read_bytes()
        self.assertTrue(pdf_bytes.startswith(b"%PDF-"), "figure is not a PDF file")
        pdf_text = pdf_bytes.decode("latin-1")
        for sensitive_key in ("Author", "Subject", "Keywords"):
            self.assertNotRegex(
                pdf_text,
                rf"/{sensitive_key}\s*(?:\(|<)",
                f"PDF Info dictionary contains /{sensitive_key}",
            )
        for label, pattern in PRIVATE_EXECUTION_PATTERNS.items():
            self.assertIsNone(pattern.search(pdf_text), f"{label} leaked in PDF bytes")

        for document in (
            REPOSITORY_ROOT / "README.md",
            EXAMPLE / "README.md",
            EXAMPLE / "PROVENANCE.md",
            EXAMPLE / "report" / "report.md",
        ):
            with self.subTest(placeholders=document.relative_to(REPOSITORY_ROOT)):
                self.assertNotIn(
                    "<!-- RESULT:",
                    document.read_text(encoding="utf-8"),
                    f"unresolved publication placeholder: {document}",
                )

        for tree in (FINAL_UNROOTED_TREE, FINAL_ROOTED_TREE):
            with self.subTest(tree=tree.relative_to(REPOSITORY_ROOT)):
                tips = brca1_newick_tips(tree)
                self.assertEqual(len(tips), 18, f"tree must contain exactly 18 tips: {tree}")
                self.assertEqual(len(tips), len(set(tips)), f"duplicate tree tip: {tree}")
                self.assertEqual(set(tips), EXPECTED_TIPS, str(tree))

        clade_rows = read_tsv(EXAMPLE / "tree" / "clade_support.tsv")
        self.assertEqual(len(clade_rows), len(EXPECTED_CLADES))
        self.assertEqual({row["clade_id"] for row in clade_rows}, EXPECTED_CLADES)

        for line_number, record in enumerate(command_records(), start=1):
            with self.subTest(command_line=line_number):
                self.assertIs(record.get("executed"), True, f"unexecuted final command: {line_number}")
                self.assertEqual(
                    record.get("status"), "completed", f"incomplete final command: {line_number}"
                )

        checksum_file = EXAMPLE / "report" / "checksums.sha256"
        checksum_lines = checksum_file.read_text(encoding="utf-8").splitlines()
        self.assertTrue(checksum_lines, "checksums.sha256 must not be empty")
        observed_paths: dict[str, str] = {}
        example_resolved = EXAMPLE.resolve(strict=True)
        for line_number, line in enumerate(checksum_lines, start=1):
            match = SHA256_LINE.fullmatch(line)
            self.assertIsNotNone(match, f"invalid checksum line {line_number}")
            assert match is not None  # Narrow the type after unittest's assertion.
            digest = match.group("digest")
            raw_path = match.group("path")
            self.assertNotIn("\\", raw_path, f"non-POSIX checksum path: {raw_path}")
            listed_path = PurePosixPath(raw_path)
            self.assertFalse(listed_path.is_absolute(), f"absolute checksum path: {raw_path}")
            self.assertNotIn("..", listed_path.parts, f"escaping checksum path: {raw_path}")
            self.assertEqual(
                listed_path.as_posix(), raw_path, f"non-normalized checksum path: {raw_path}"
            )
            self.assertGreaterEqual(len(listed_path.parts), 3, raw_path)
            self.assertEqual(
                listed_path.parts[:2], ("examples", "brca1"),
                f"checksum target is outside examples/brca1: {raw_path}",
            )
            self.assertNotIn(raw_path, observed_paths, f"duplicate checksum path: {raw_path}")

            target = (REPOSITORY_ROOT / Path(*listed_path.parts)).resolve(strict=True)
            try:
                target.relative_to(example_resolved)
            except ValueError:
                self.fail(f"checksum target escapes examples/brca1 after resolution: {raw_path}")
            self.assertTrue(target.is_file(), f"checksum target is not a file: {raw_path}")
            observed = hashlib.sha256(target.read_bytes()).hexdigest()
            self.assertEqual(observed, digest, f"checksum mismatch: {raw_path}")
            observed_paths[raw_path] = digest

        required_checksummed = {
            path.relative_to(REPOSITORY_ROOT).as_posix()
            for path in FINAL_PUBLICATION_ARTIFACTS
            if path != checksum_file
        }
        required_checksummed.update(
            {
                "examples/brca1/README.md",
                "examples/brca1/PROVENANCE.md",
            }
        )
        self.assertEqual(
            required_checksummed - set(observed_paths),
            set(),
            "final public artifacts missing from report/checksums.sha256",
        )

if __name__ == "__main__":
    unittest.main()
