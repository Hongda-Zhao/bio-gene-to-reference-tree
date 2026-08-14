"""Focused tests for fail-closed publication of BRCA1 IQ-TREE outputs."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
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
    / "sanitize_iqtree_outputs.py"
)
PRIVATE_ALIGNMENT = "/lustre/scratch/private-run/alignment.trimmed.balanced.faa"
PREFIX = "brca1-balanced"
TREE = b"(NP_009225.1:0.01,(NP_001038958.1:0.02,NP_033894.3:0.03)99/100:0.04);\n"


class Brca1IqtreeSanitizerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory(
            prefix="brca1-iqtree-sanitizer-"
        )
        self.root = Path(self.temporary_directory.name)
        self.raw = self.root / "raw"
        self.raw.mkdir()
        self.public = self.root / "public"
        self.write_minimal_raw()

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def write_minimal_raw(self) -> None:
        (self.raw / f"{PREFIX}.log").write_text(
            "IQ-TREE 3\n"
            f"Command: iqtree3 -s {PRIVATE_ALIGNMENT} -B 1000\n"
            "Host: cs03 (AVX2, 503 GB RAM)\n"
            "Website: https://www.iqtree.org\n",
            encoding="utf-8",
        )
        (self.raw / f"{PREFIX}.iqtree").write_text(
            f"Input file: {PRIVATE_ALIGNMENT}\nHost:\tcs03\n",
            encoding="utf-8",
        )
        for suffix in (".treefile", ".contree", ".ufboot"):
            (self.raw / f"{PREFIX}{suffix}").write_bytes(TREE)
        (self.raw / "software_version.txt").write_text(
            "IQ-TREE multicore version 3.1.3\n", encoding="utf-8"
        )

    def run_script(
        self, *, public: Path | None = None, private_path: str = PRIVATE_ALIGNMENT
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            (
                sys.executable,
                str(SCRIPT),
                "--raw-dir",
                str(self.raw),
                "--public-dir",
                str(public or self.public),
                "--private-alignment-path",
                private_path,
            ),
            capture_output=True,
            text=True,
            check=False,
        )

    def test_publish_redacts_text_omits_model_and_emits_recomputable_receipt(self) -> None:
        model_bytes = gzip.compress(b"binary model payload\x00\xff")
        (self.raw / f"{PREFIX}.model.gz").write_bytes(model_bytes)
        (self.raw / f"{PREFIX}.splits.nex").write_text(
            "#NEXUS\nBEGIN SPLITS;\nEND;\n", encoding="utf-8"
        )
        (self.raw / f"{PREFIX}.ckp.gz").write_bytes(gzip.compress(b"checkpoint"))
        (self.raw / "checksums.sha256").write_text(
            "stale wrapper checksum\n", encoding="utf-8"
        )

        completed = self.run_script()
        self.assertEqual(completed.returncode, 0, msg=completed.stderr)
        expected = {
            f"{PREFIX}{suffix}"
            for suffix in (
                ".log",
                ".iqtree",
                ".treefile",
                ".contree",
                ".ufboot",
                ".splits.nex",
            )
        } | {
            "software_version.txt",
            "redactions.tsv",
            "raw-output-receipt.json",
        }
        self.assertEqual({path.name for path in self.public.iterdir()}, expected)

        public_log = (self.public / f"{PREFIX}.log").read_text(encoding="utf-8")
        self.assertNotIn(PRIVATE_ALIGNMENT, public_log)
        self.assertIn("alignment/alignment.trimmed.balanced.faa", public_log)
        self.assertIn("Host: <batch-compute-node>", public_log)
        self.assertIn("https://www.iqtree.org", public_log)
        self.assertEqual((self.public / f"{PREFIX}.treefile").read_bytes(), TREE)
        self.assertFalse((self.public / f"{PREFIX}.model.gz").exists())
        self.assertIn(
            "Host: <batch-compute-node> (AVX2, 503 GB RAM)", public_log
        )

        with (self.public / "redactions.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            rows = {row["filename"]: row for row in csv.DictReader(handle, delimiter="\t")}
        log_row = rows[f"{PREFIX}.log"]
        self.assertEqual(log_row["redaction_types"], "host_line:1;private_alignment_path:1")
        self.assertEqual(log_row["redaction_count"], "2")
        tree_row = rows[f"{PREFIX}.treefile"]
        tree_sha256 = hashlib.sha256(TREE).hexdigest()
        self.assertEqual(tree_row["redaction_types"], "none")
        self.assertEqual(tree_row["original_sha256"], tree_sha256)
        self.assertEqual(tree_row["public_sha256"], tree_sha256)

        receipt_text = (self.public / "raw-output-receipt.json").read_text(
            encoding="utf-8"
        )
        receipt = json.loads(receipt_text)
        self.assertEqual(receipt["raw_regular_file_count"], 10)
        self.assertEqual(receipt["published_file_count"], 7)
        self.assertEqual(receipt["redacted_file_count"], 2)
        self.assertEqual(receipt["redacted_occurrence_count"], 4)
        self.assertEqual(
            receipt["omitted_file_counts"],
            {
                "iqtree_checkpoint": 1,
                "iqtree_model": 1,
                "wrapper_checksums": 1,
            },
        )
        self.assertRegex(receipt["raw_manifest_sha256"], r"^[0-9a-f]{64}$")
        manifest = receipt["raw_file_manifest"]
        self.assertEqual(
            [record["filename"] for record in manifest],
            sorted(record["filename"] for record in manifest),
        )
        self.assertTrue(
            all(
                set(record) == {"filename", "byte_size", "sha256"}
                for record in manifest
            )
        )
        model_record = next(
            record for record in manifest if record["filename"] == f"{PREFIX}.model.gz"
        )
        self.assertEqual(model_record["byte_size"], len(model_bytes))
        self.assertEqual(model_record["sha256"], hashlib.sha256(model_bytes).hexdigest())
        recomputed = hashlib.sha256()
        for record in manifest:
            name_bytes = record["filename"].encode("utf-8")
            recomputed.update(len(name_bytes).to_bytes(8, byteorder="big"))
            recomputed.update(name_bytes)
            recomputed.update(record["byte_size"].to_bytes(8, byteorder="big"))
            recomputed.update(bytes.fromhex(record["sha256"]))
        self.assertEqual(recomputed.hexdigest(), receipt["raw_manifest_sha256"])
        self.assertNotIn(PRIVATE_ALIGNMENT, receipt_text)
        self.assertNotIn("cs03", receipt_text)

    def test_missing_core_artifact_fails_before_creating_output(self) -> None:
        (self.raw / f"{PREFIX}.ufboot").unlink()
        completed = self.run_script()
        self.assertEqual(completed.returncode, 2)
        self.assertIn("exactly one core .ufboot", completed.stderr)
        self.assertFalse(self.public.exists())

    def test_unknown_file_fails_before_creating_output(self) -> None:
        (self.raw / "unexpected.tmp").write_text("review me\n", encoding="utf-8")
        completed = self.run_script()
        self.assertEqual(completed.returncode, 2)
        self.assertIn("unknown IQ-TREE output", completed.stderr)
        self.assertFalse(self.public.exists())

    def test_symlink_and_nested_directory_are_rejected(self) -> None:
        for name, make_entry in (
            (
                "symlink",
                lambda: (self.raw / "unsafe.log").symlink_to(
                    self.raw / f"{PREFIX}.log"
                ),
            ),
            ("directory", lambda: (self.raw / "nested").mkdir()),
        ):
            with self.subTest(name=name):
                make_entry()
                output = self.root / f"public-{name}"
                completed = self.run_script(public=output)
                self.assertEqual(completed.returncode, 2)
                self.assertFalse(output.exists())
                entry = self.raw / ("unsafe.log" if name == "symlink" else "nested")
                if entry.is_symlink():
                    entry.unlink()
                else:
                    entry.rmdir()

    def test_existing_output_is_never_overwritten(self) -> None:
        self.public.mkdir()
        sentinel = self.public / "sentinel.txt"
        sentinel.write_text("keep\n", encoding="utf-8")
        completed = self.run_script()
        self.assertEqual(completed.returncode, 2)
        self.assertIn("overwrite is forbidden", completed.stderr)
        self.assertEqual(sentinel.read_text(encoding="utf-8"), "keep\n")

    def test_non_utf8_text_fails_before_creating_output(self) -> None:
        (self.raw / f"{PREFIX}.iqtree").write_bytes(b"invalid: \xff\n")
        completed = self.run_script()
        self.assertEqual(completed.returncode, 2)
        self.assertIn("expected UTF-8", completed.stderr)
        self.assertFalse(self.public.exists())

    def test_residual_private_path_login_and_pbs_id_each_fail_closed(self) -> None:
        unsafe_values = (
            "Other path: /home/private/other.faa\n",
            "Login: private-user@gds2\n",
            "Scheduler: 5014284.fe3-adm\n",
        )
        for index, unsafe in enumerate(unsafe_values):
            with self.subTest(unsafe=index):
                path = self.raw / f"{PREFIX}.iqtree"
                baseline = path.read_text(encoding="utf-8")
                path.write_text(baseline + unsafe, encoding="utf-8")
                output = self.root / f"public-unsafe-{index}"
                completed = self.run_script(public=output)
                self.assertEqual(completed.returncode, 2)
                self.assertIn("remaining", completed.stderr)
                self.assertFalse(output.exists())
                path.write_text(baseline, encoding="utf-8")

    def test_mixed_iqtree_prefixes_fail_closed(self) -> None:
        original = self.raw / f"{PREFIX}.contree"
        original.rename(self.raw / "other-run.contree")
        completed = self.run_script()
        self.assertEqual(completed.returncode, 2)
        self.assertIn("does not share the core IQ-TREE prefix", completed.stderr)
        self.assertFalse(self.public.exists())


if __name__ == "__main__":
    unittest.main()
