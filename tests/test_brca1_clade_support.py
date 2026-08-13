"""Static and optional runtime checks for the BRCA1 clade-support audit."""

from __future__ import annotations

import csv
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    REPOSITORY_ROOT
    / "examples"
    / "brca1"
    / "scripts"
    / "summarize_brca1_clade_support.R"
)

AMPHIBIA = ("NP_001107963.1", "XP_029429046.1")
MAMMALIA = (
    "NP_009225.1", "NP_001038958.1", "NP_033894.3", "XP_017204702.2",
    "NP_001013434.1", "NP_848668.1", "XP_014595447.2", "XP_003414318.3",
    "XP_058140293.1", "NP_001029141.1",
)
SAUROPSIDA = (
    "NP_989500.1", "XP_072775070.1", "XP_019406054.1", "XP_023967135.2",
    "XP_008111382.1", "XP_026576759.1",
)
EXPECTED_TIPS = set(AMPHIBIA + MAMMALIA + SAUROPSIDA)


def leaf(accession: str) -> str:
    return f"{accession}:0.1"


def synthetic_tree() -> str:
    primates = f"({leaf('NP_009225.1')},{leaf('NP_001038958.1')})99/100:0.1"
    glires = f"({leaf('NP_033894.3')},{leaf('XP_017204702.2')})98/100:0.1"
    other_mammals = (
        f"({leaf('NP_001013434.1')},"
        f"({leaf('NP_848668.1')},"
        f"({leaf('XP_014595447.2')},"
        f"({leaf('XP_003414318.3')},"
        f"({leaf('XP_058140293.1')},{leaf('NP_001029141.1')})85/90:0.1"
        f")86/91:0.1)87/92:0.1)88/93:0.1)89/94:0.1"
    )
    mammalia = f"({primates},({glires},{other_mammals})90/95:0.1)91/96:0.1"

    aves = f"({leaf('NP_989500.1')},{leaf('XP_072775070.1')})96/100:0.1"
    archosauria = f"({leaf('XP_019406054.1')},{aves})95/99:0.1"
    archelosauria = f"({leaf('XP_023967135.2')},{archosauria})94/98:0.1"
    lepidosauria = f"({leaf('XP_008111382.1')},{leaf('XP_026576759.1')})92/97:0.1"
    sauropsida = f"({lepidosauria},{archelosauria})97/100:0.1"
    amphibia = f"({leaf('NP_001107963.1')},{leaf('XP_029429046.1')})95/100:0.1"
    # The unlabelled amniote node is the one allowed reroot artifact.
    return f"({amphibia},({mammalia},{sauropsida}):0.1)Root;\n"


class Brca1CladeSupportTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory(prefix="brca1-support-test-")
        self.root = Path(self.temporary_directory.name)
        self.tree = self.root / "rooted.nwk"
        self.tree.write_text(synthetic_tree(), encoding="utf-8")
        self.output = self.root / "support.tsv"

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    @staticmethod
    def runtime_available() -> bool:
        rscript = shutil.which("Rscript")
        if not rscript:
            return False
        return subprocess.run(
            [rscript, "-e", 'quit(status=if(requireNamespace("ape", quietly=TRUE)) 0 else 1)'],
            check=False,
            capture_output=True,
        ).returncode == 0

    def run_script(self, *, tree: Path | None = None, output: Path | None = None):
        return subprocess.run(
            [
                shutil.which("Rscript") or "Rscript",
                str(SCRIPT),
                "--tree", str(tree or self.tree),
                "--output", str(output or self.output),
            ],
            text=True,
            capture_output=True,
            check=False,
        )

    def test_script_is_example_specific_local_and_fail_closed(self) -> None:
        content = SCRIPT.read_text(encoding="utf-8")
        for accession in EXPECTED_TIPS:
            self.assertIn(accession, content)
        for required in (
            "ape::read.tree",
            "Tree tips do not exactly equal the fixed 18-accession BRCA1 set",
            "finite, non-negative branch length",
            "numeric SH-aLRT/UFBoot pair",
            "not_applicable_root_artifact",
            "not-recovered",
            "Refusing to overwrite existing output",
            "Structural root split must exactly isolate",
        ):
            self.assertIn(required, content)
        for forbidden in (
            "install.packages", "BiocManager::install", "download.file",
            "system(", "system2(", "http://", "https://",
        ):
            self.assertNotIn(forbidden, content)

    def test_r_syntax_when_rscript_is_available(self) -> None:
        rscript = shutil.which("Rscript")
        if not rscript:
            self.skipTest("Rscript is not installed")
        completed = subprocess.run(
            [rscript, "-e", "parse(file=commandArgs(TRUE)[1])", str(SCRIPT)],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_runtime_reports_expected_clades_and_explicit_alternative(self) -> None:
        if not self.runtime_available():
            self.skipTest("Rscript plus ape is not installed")
        completed = self.run_script()
        self.assertEqual(completed.returncode, 0, completed.stderr)
        with self.output.open(encoding="utf-8", newline="") as handle:
            rows = {row["clade_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(
            set(rows),
            {
                "Amphibia", "Mammalia", "Sauropsida", "Primates", "Glires",
                "Aves", "Lepidosauria", "Archelosauria", "Archosauria",
                "Crocodylia+Testudines",
            },
        )
        for clade_id, row in rows.items():
            if clade_id == "Crocodylia+Testudines":
                self.assertEqual(row["hypothesis_role"], "explicit_alternative")
                self.assertEqual(row["recovery_status"], "not-recovered")
                self.assertEqual(row["support_status"], "not_recovered")
                self.assertEqual(row["node_label"], "")
                self.assertEqual(row["sh_alrt"], "")
                self.assertEqual(row["ufboot"], "")
            else:
                self.assertEqual(row["recovery_status"], "recovered", clade_id)
                self.assertEqual(row["support_status"], "available", clade_id)
                self.assertRegex(row["node_label"], r"^[0-9.]+/[0-9.]+$")
        self.assertEqual(rows["Primates"]["sh_alrt"], "99")
        self.assertEqual(rows["Primates"]["ufboot"], "100")
        self.assertEqual(set(rows["Mammalia"]["target_tips"].split(",")), set(MAMMALIA))

    def test_bad_tip_label_support_root_and_overwrite_fail(self) -> None:
        if not self.runtime_available():
            self.skipTest("Rscript plus ape is not installed")

        bad_tip = self.root / "bad-tip.nwk"
        bad_tip.write_text(synthetic_tree().replace("NP_009225.1", "EXTRA", 1), encoding="utf-8")
        result = self.run_script(tree=bad_tip, output=self.root / "bad-tip.tsv")
        self.assertEqual(result.returncode, 2)
        self.assertIn("fixed 18-accession", result.stderr)

        malformed = self.root / "malformed.nwk"
        malformed.write_text(synthetic_tree().replace("99/100", "99", 1), encoding="utf-8")
        result = self.run_script(tree=malformed, output=self.root / "malformed.tsv")
        self.assertEqual(result.returncode, 2)
        self.assertIn("numeric SH-aLRT/UFBoot pair", result.stderr)

        wrong_root = self.root / "wrong-root.nwk"
        text = synthetic_tree()
        # A rooted binary tree with the right tips but a human+amphibian root side.
        wrong_root.write_text(
            text.replace(leaf("NP_009225.1"), "SWAP_TIP:0.1", 1)
            .replace(leaf("XP_029429046.1"), leaf("NP_009225.1"), 1)
            .replace("SWAP_TIP:0.1", leaf("XP_029429046.1"), 1),
            encoding="utf-8",
        )
        result = self.run_script(tree=wrong_root, output=self.root / "wrong-root.tsv")
        self.assertEqual(result.returncode, 2)
        self.assertIn("Structural root split", result.stderr)

        self.output.write_text("existing\n", encoding="utf-8")
        result = self.run_script()
        self.assertEqual(result.returncode, 2)
        self.assertIn("Refusing to overwrite", result.stderr)
        self.assertEqual(self.output.read_text(encoding="utf-8"), "existing\n")

    def test_negative_and_nonfinite_branch_lengths_fail(self) -> None:
        if not self.runtime_available():
            self.skipTest("Rscript plus ape is not installed")
        for name, replacement in (("negative", ":-0.1"), ("nonfinite", ":NaN")):
            tree = self.root / f"{name}.nwk"
            tree.write_text(synthetic_tree().replace(":0.1", replacement, 1), encoding="utf-8")
            result = self.run_script(tree=tree, output=self.root / f"{name}.tsv")
            self.assertEqual(result.returncode, 2)
            self.assertIn("finite, non-negative branch length", result.stderr)


if __name__ == "__main__":
    unittest.main()
