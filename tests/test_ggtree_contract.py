"""Static and optional runtime checks for the ggtree renderer."""

from __future__ import annotations

import csv
import hashlib
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SKILL_ROOT = REPOSITORY_ROOT / "skills" / "bio-gene-to-reference-tree"
SCRIPT = SKILL_ROOT / "scripts" / "render_tree_ggtree.R"
FIXTURE = REPOSITORY_ROOT / "tests" / "fixtures" / "visualization"


class GgtreeContractTests(unittest.TestCase):
    def test_renderer_is_local_strict_and_uses_vector_outputs(self) -> None:
        content = SCRIPT.read_text(encoding="utf-8")
        for required in (
            "ggtree::ggtree",
            "ggplot2::ggsave",
            "svglite::svglite",
            'ladderize = FALSE',
            '"rectangular", "circular"',
            '"unrooted", "outgroup-rooted"',
            '"fasttree-sh-like", "sh-alrt/ufboot", "sh-alrt/bootstrap"',
            "structural_root_marker_index",
            '!identical(tree$node.label[[label_index]], "Root")',
            "ggplot2::expand_limits(x = max(finite_x) + 0.30 * x_span)",
            '"#E69F00"',
            '"#009E73"',
            '"#999999"',
            '".svg"',
            '".pdf"',
            '".settings.tsv"',
        ):
            self.assertIn(required, content)
        for forbidden in (
            "install.packages",
            "BiocManager::install",
            "download.file",
            "system(",
            "system2(",
            "http://",
            "https://",
        ):
            self.assertNotIn(forbidden, content)

    def test_fixture_tip_sets_and_roles_are_identical(self) -> None:
        with (FIXTURE / "sequence_metadata.tsv").open(encoding="utf-8", newline="") as handle:
            metadata = list(csv.DictReader(handle, delimiter="\t"))
        selected = {row["tip_id"]: row["analysis_role"] for row in metadata if row["inclusion_status"] == "selected"}
        lines = (FIXTURE / "itol_roles.txt").read_text(encoding="utf-8").splitlines()
        data_index = lines.index("DATA")
        itol = {parts[0]: parts[2].lower() for parts in (line.split("\t") for line in lines[data_index + 1 :])}
        self.assertEqual(selected, itol)
        newick = (FIXTURE / "gene-tree.nwk").read_text(encoding="utf-8")
        newick_tips = set(re.findall(r"(?<=[(,])([^():;,]+)(?=:)", newick))
        self.assertEqual(set(selected), newick_tips)

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

    def test_help_and_unknown_option_have_stable_cli_behavior(self) -> None:
        rscript = shutil.which("Rscript")
        if not rscript:
            self.skipTest("Rscript is not installed")
        help_result = subprocess.run(
            [rscript, str(SCRIPT), "--help"],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(help_result.returncode, 0, help_result.stderr)
        self.assertIn("Usage: render_tree_ggtree.R", help_result.stdout)
        unknown = subprocess.run(
            [rscript, str(SCRIPT), "--unknown", "value"],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(unknown.returncode, 2)
        self.assertIn("Unknown argument: --unknown", unknown.stderr)

    def test_vector_smoke_render_when_packages_are_available(self) -> None:
        rscript = shutil.which("Rscript")
        if not rscript:
            self.skipTest("Rscript is not installed")
        dependency_check = subprocess.run(
            [
                rscript,
                "-e",
                'quit(status=if(all(vapply(c("ape","ggplot2","ggtree","openssl","svglite"), requireNamespace, quietly=TRUE, FUN.VALUE=logical(1)))) 0 else 1)',
            ],
            check=False,
        )
        if dependency_check.returncode != 0:
            self.skipTest("Local ggtree runtime packages are not all installed")
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary) / "gene-tree.unrooted.ggtree"
            command = [
                rscript,
                str(SCRIPT),
                "--tree",
                str(FIXTURE / "gene-tree.nwk"),
                "--metadata",
                str(FIXTURE / "sequence_metadata.tsv"),
                "--itol-roles",
                str(FIXTURE / "itol_roles.txt"),
                "--out-prefix",
                str(prefix),
                "--root-state",
                "unrooted",
                "--layout",
                "rectangular",
                "--branch-length",
                "auto",
                "--support-format",
                "sh-alrt/ufboot",
                "--width",
                "8",
                "--height",
                "6",
            ]
            completed = subprocess.run(command, text=True, capture_output=True, check=False)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            for suffix in (".svg", ".pdf", ".settings.tsv"):
                output = Path(f"{prefix}{suffix}")
                self.assertTrue(output.is_file())
                self.assertGreater(output.stat().st_size, 50)
            settings_path = Path(f"{prefix}.settings.tsv")
            settings = settings_path.read_text(encoding="utf-8")
            self.assertIn("rendered_branch_mode\tphylogram", settings)
            self.assertIn("support_format\tsh-alrt/ufboot", settings)
            self.assertIn("ladderized\tfalse", settings)
            with settings_path.open(encoding="utf-8", newline="") as handle:
                settings_by_name = {
                    row["setting"]: row["value"]
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            expected_tree_hash = hashlib.sha256(
                (FIXTURE / "gene-tree.nwk").read_bytes()
            ).hexdigest()
            self.assertEqual(settings_by_name["tree_sha256"], expected_tree_hash)
            self.assertEqual(settings_by_name["study_color"], "#E69F00")
            self.assertEqual(settings_by_name["expanded_color"], "#009E73")
            self.assertEqual(settings_by_name["outgroup_color"], "#999999")

            repeated = subprocess.run(command, text=True, capture_output=True, check=False)
            self.assertEqual(repeated.returncode, 2)
            self.assertIn("Refusing to overwrite", repeated.stderr)

    def test_runtime_join_is_by_named_tip_id_not_first_tsv_column(self) -> None:
        rscript = shutil.which("Rscript")
        if not rscript:
            self.skipTest("Rscript is not installed")
        dependency_check = subprocess.run(
            [
                rscript,
                "-e",
                'quit(status=if(all(vapply(c("ape","ggplot2","ggtree","openssl","svglite"), requireNamespace, quietly=TRUE, FUN.VALUE=logical(1)))) 0 else 1)',
            ],
            check=False,
        )
        if dependency_check.returncode != 0:
            self.skipTest("Local ggtree runtime packages are not all installed")
        with tempfile.TemporaryDirectory() as temporary:
            temporary_path = Path(temporary)
            rows: list[dict[str, str]]
            with (FIXTURE / "sequence_metadata.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            reordered = temporary_path / "metadata-reordered.tsv"
            columns = ["species", "accession", "analysis_role", "inclusion_status", "tip_id"]
            with reordered.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            prefix = temporary_path / "reordered"
            completed = subprocess.run(
                [
                    rscript,
                    str(SCRIPT),
                    "--tree",
                    str(FIXTURE / "gene-tree.nwk"),
                    "--metadata",
                    str(reordered),
                    "--out-prefix",
                    str(prefix),
                    "--root-state",
                    "unrooted",
                    "--layout",
                    "rectangular",
                    "--support-format",
                    "sh-alrt/ufboot",
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            svg = Path(f"{prefix}.svg").read_text(encoding="utf-8")
            for tip_id in ("QUERY_001", "MOUSE_CAN", "CHICKEN_OK", "CIONA_OUT"):
                self.assertIn(tip_id, svg)

    def test_outgroup_rooted_declaration_rejects_unrooted_tree(self) -> None:
        rscript = shutil.which("Rscript")
        if not rscript:
            self.skipTest("Rscript is not installed")
        dependency_check = subprocess.run(
            [
                rscript,
                "-e",
                'quit(status=if(all(vapply(c("ape","ggplot2","ggtree","openssl","svglite"), requireNamespace, quietly=TRUE, FUN.VALUE=logical(1)))) 0 else 1)',
            ],
            check=False,
        )
        if dependency_check.returncode != 0:
            self.skipTest("Local ggtree runtime packages are not all installed")
        with tempfile.TemporaryDirectory() as temporary:
            completed = subprocess.run(
                [
                    rscript,
                    str(SCRIPT),
                    "--tree",
                    str(FIXTURE / "gene-tree.nwk"),
                    "--metadata",
                    str(FIXTURE / "sequence_metadata.tsv"),
                    "--out-prefix",
                    str(Path(temporary) / "invalid-root"),
                    "--root-state",
                    "outgroup-rooted",
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn("structurally rooted", completed.stderr)

    def test_outgroup_rooted_render_requires_and_accepts_root_split(self) -> None:
        rscript = shutil.which("Rscript")
        if not rscript:
            self.skipTest("Rscript is not installed")
        dependency_check = subprocess.run(
            [
                rscript,
                "-e",
                'quit(status=if(all(vapply(c("ape","ggplot2","ggtree","openssl","svglite"), requireNamespace, quietly=TRUE, FUN.VALUE=logical(1)))) 0 else 1)',
            ],
            check=False,
        )
        if dependency_check.returncode != 0:
            self.skipTest("Local ggtree runtime packages are not all installed")
        with tempfile.TemporaryDirectory() as temporary:
            completed = subprocess.run(
                [
                    rscript,
                    str(SCRIPT),
                    "--tree",
                    str(FIXTURE / "gene-tree-rooted.nwk"),
                    "--metadata",
                    str(FIXTURE / "sequence_metadata.tsv"),
                    "--out-prefix",
                    str(Path(temporary) / "rooted"),
                    "--root-state",
                    "outgroup-rooted",
                    "--support-format",
                    "sh-alrt/ufboot",
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            rooted_svg = Path(temporary, "rooted.svg")
            self.assertTrue(rooted_svg.is_file())
            svg_text = rooted_svg.read_text(encoding="utf-8")
            self.assertIn(">95/99<", svg_text)
            self.assertNotIn(">Root<", svg_text)

    def test_root_support_exception_is_rejected_away_from_structural_root(self) -> None:
        rscript = shutil.which("Rscript")
        if not rscript:
            self.skipTest("Rscript is not installed")
        dependency_check = subprocess.run(
            [
                rscript,
                "-e",
                'quit(status=if(all(vapply(c("ape","ggplot2","ggtree","openssl","svglite"), requireNamespace, quietly=TRUE, FUN.VALUE=logical(1)))) 0 else 1)',
            ],
            check=False,
        )
        if dependency_check.returncode != 0:
            self.skipTest("Local ggtree runtime packages are not all installed")
        with tempfile.TemporaryDirectory() as temporary:
            temporary_path = Path(temporary)
            non_root_marker = temporary_path / "non-root-marker.nwk"
            non_root_marker.write_text(
                "((QUERY_001:0.1,(MOUSE_CAN:0.1,CHICKEN_OK:0.1)Root:0.1)95/99:0.1,"
                "CIONA_OUT:0.3)98/100;\n",
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    rscript,
                    str(SCRIPT),
                    "--tree",
                    str(non_root_marker),
                    "--metadata",
                    str(FIXTURE / "sequence_metadata.tsv"),
                    "--out-prefix",
                    str(temporary_path / "non-root-marker"),
                    "--root-state",
                    "outgroup-rooted",
                    "--support-format",
                    "sh-alrt/ufboot",
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn("explicit value/value form", completed.stderr)

            unrooted_marker = temporary_path / "unrooted-marker.nwk"
            unrooted_marker.write_text(
                "(QUERY_001:0.1,(MOUSE_CAN:0.1,CHICKEN_OK:0.1)95/99:0.1,"
                "CIONA_OUT:0.3)Root;\n",
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    rscript,
                    str(SCRIPT),
                    "--tree",
                    str(unrooted_marker),
                    "--metadata",
                    str(FIXTURE / "sequence_metadata.tsv"),
                    "--out-prefix",
                    str(temporary_path / "unrooted-marker"),
                    "--root-state",
                    "unrooted",
                    "--support-format",
                    "sh-alrt/ufboot",
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn("explicit value/value form", completed.stderr)


if __name__ == "__main__":
    unittest.main()
