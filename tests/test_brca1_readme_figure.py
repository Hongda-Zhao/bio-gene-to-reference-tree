"""Contract and optional runtime tests for the BRCA1 README hero figure."""

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
EXAMPLE = REPOSITORY_ROOT / "examples" / "brca1"
SCRIPT = EXAMPLE / "scripts" / "render_brca1_readme.R"
FIGURE_PREFIX = EXAMPLE / "figures" / "brca1-readme"
TREE = EXAMPLE / "tree" / "gene-tree.outgroup-rooted.nwk"
METADATA = EXAMPLE / "annotation" / "sequence_metadata.tsv"
CLADE_SUPPORT = EXAMPLE / "tree" / "clade_support.tsv"
CHECKSUMS = EXAMPLE / "report" / "checksums.sha256"

EXPECTED_WEAK_SUPPORTS = {"95.6/89", "80.3/77", "55.2/63", "88.9/88", "77.3/79"}
EXPECTED_BANDS = {
    "Mammalia": "#E8F1F8",
    "Sauropsida": "#F7EAD7",
    "Amphibia": "#EEEAF4",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_settings(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    if any(set(row) != {"setting", "value"} for row in rows):
        raise AssertionError(f"invalid settings rows: {path}")
    settings = {row["setting"]: row["value"] for row in rows}
    if len(settings) != len(rows):
        raise AssertionError(f"duplicate settings keys: {path}")
    return settings


class Brca1ReadmeFigureTests(unittest.TestCase):
    def test_committed_hero_artifacts_are_self_contained_and_legible(self) -> None:
        svg_path = FIGURE_PREFIX.with_suffix(".svg")
        pdf_path = FIGURE_PREFIX.with_suffix(".pdf")
        settings_path = Path(f"{FIGURE_PREFIX}.settings.tsv")
        for path in (SCRIPT, svg_path, pdf_path, settings_path):
            self.assertTrue(path.is_file(), str(path))
            self.assertGreater(path.stat().st_size, 100, str(path))

        svg = svg_path.read_text(encoding="utf-8")
        self.assertRegex(svg, r"width='828\.00pt' height='504\.00pt'")
        self.assertIn("BRCA1 protein gene tree", svg)
        self.assertIn("provisional amphibian-outgroup display root", svg)
        self.assertIn("SH-aLRT &gt;= 80 and UFBoot &gt;= 95", svg)
        self.assertIn("0.5 substitutions/site", svg)
        self.assertIn("Tip order was not ladderized", svg)
        self.assertIn('font-family: "Arial"', svg)
        for clade, color in EXPECTED_BANDS.items():
            self.assertIn(clade, svg)
            self.assertIn(color.upper(), svg.upper())
        for color in ("#E69F00", "#009E73", "#999999"):
            self.assertIn(color, svg.upper())

        metadata = [
            row for row in read_tsv(METADATA) if row["inclusion_status"] == "selected"
        ]
        self.assertEqual(len(metadata), 18)
        for row in metadata:
            self.assertIn(row["species"], svg)
            self.assertIn(row["accession"], svg)
        for support in EXPECTED_WEAK_SUPPORTS:
            self.assertIn(f">{support}<", svg)
        self.assertNotIn(">100/100<", svg)

        species_font_sizes = []
        for row in metadata:
            match = re.search(
                rf"<text[^>]*style='[^']*font-size: ([0-9.]+)px;[^']*'[^>]*>"
                rf"{re.escape(row['species'])}</text>",
                svg,
            )
            self.assertIsNotNone(match, row["species"])
            assert match is not None
            species_font_sizes.append(float(match.group(1)))
        self.assertGreaterEqual(min(species_font_sizes), 9.0)

        pdf_bytes = pdf_path.read_bytes()
        self.assertTrue(pdf_bytes.startswith(b"%PDF-"))
        self.assertNotIn(b"/Author", pdf_bytes)

    def test_settings_bind_inputs_and_visual_semantics(self) -> None:
        settings = read_settings(Path(f"{FIGURE_PREFIX}.settings.tsv"))
        expected = {
            "display_style": "compact-rectangular-phylogram",
            "root_state": "provisional-outgroup-rooted-display-derivative",
            "branch_length": "substitutions-per-site",
            "support_format": "sh-alrt/ufboot",
            "sh_alrt_joint_threshold": "80",
            "ufboot_joint_threshold": "95",
            "joint_strong_node_count": "10",
            "below_joint_threshold_node_count": "5",
            "undisplayed_internal_label_slot_count": "2",
            "ladderized": "false",
            "tip_count": "18",
            "study_tip_count": "1",
            "expanded_tip_count": "15",
            "outgroup_tip_count": "2",
            "width_inches": "11.5",
            "height_inches": "7",
            "font_family": "Arial",
            "broad_clade_bands": "Mammalia,Sauropsida,Amphibia",
            "study_color": "#E69F00",
            "expanded_color": "#009E73",
            "outgroup_color": "#999999",
        }
        for key, value in expected.items():
            self.assertEqual(settings.get(key), value, key)
        for key, path in (
            ("tree_sha256", TREE),
            ("metadata_sha256", METADATA),
            ("clade_support_sha256", CLADE_SUPPORT),
        ):
            self.assertEqual(settings[key], hashlib.sha256(path.read_bytes()).hexdigest())
        for clade, color in EXPECTED_BANDS.items():
            self.assertEqual(settings[f"{clade}_band_color"], color)

    def test_readme_uses_hero_and_preserves_detailed_audit_links(self) -> None:
        readme = (REPOSITORY_ROOT / "README.md").read_text(encoding="utf-8")
        self.assertIn(
            "![Outgroup-rooted BRCA1 protein gene tree for 18 vertebrates]"
            "(examples/brca1/figures/brca1-readme.svg)",
            readme,
        )
        self.assertIn(
            "[detailed audit SVG]"
            "(examples/brca1/figures/gene-tree.outgroup-rooted.ggtree.svg)",
            readme,
        )
        self.assertIn(
            "[detailed audit PDF]"
            "(examples/brca1/figures/gene-tree.outgroup-rooted.ggtree.pdf)",
            readme,
        )

    def test_example_checksum_manifest_covers_new_artifacts(self) -> None:
        lines = CHECKSUMS.read_text(encoding="utf-8").splitlines()
        manifest: dict[str, str] = {}
        for line in lines:
            digest, path = line.split("  ", 1)
            manifest[path] = digest
        required = (
            SCRIPT,
            FIGURE_PREFIX.with_suffix(".svg"),
            FIGURE_PREFIX.with_suffix(".pdf"),
            Path(f"{FIGURE_PREFIX}.settings.tsv"),
        )
        for path in required:
            relative = path.relative_to(REPOSITORY_ROOT).as_posix()
            self.assertIn(relative, manifest)
            self.assertEqual(manifest[relative], hashlib.sha256(path.read_bytes()).hexdigest())

    def test_script_is_offline_fail_closed_and_parses_when_r_is_available(self) -> None:
        script = SCRIPT.read_text(encoding="utf-8")
        for required in (
            "ladderize = FALSE",
            "validate_root_split",
            "read_broad_clades",
            "sh_alrt >= 80 & ufboot >= 95",
            "geom_hilight",
            "aligned display metadata",
            "ggplot2::geom_segment",
            "svglite::svglite",
            "grDevices::cairo_pdf",
            "Refusing to overwrite",
        ):
            self.assertIn(required, script)
        for forbidden in (
            "install.packages",
            "BiocManager::install",
            "download.file",
            "system(",
            "system2(",
            "http://",
            "https://",
        ):
            self.assertNotIn(forbidden, script)

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

    def test_runtime_render_and_refuse_overwrite_when_packages_are_available(self) -> None:
        rscript = shutil.which("Rscript")
        if not rscript:
            self.skipTest("Rscript is not installed")
        dependency_check = subprocess.run(
            [
                rscript,
                "-e",
                'quit(status=if(capabilities("cairo") && '
                'all(vapply(c("ape","ggplot2","ggtree","openssl","svglite"), '
                "requireNamespace, quietly=TRUE, FUN.VALUE=logical(1)))) 0 else 1)",
            ],
            check=False,
        )
        if dependency_check.returncode != 0:
            self.skipTest("Local BRCA1 README figure runtime is incomplete")

        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary) / "brca1-readme"
            command = [
                rscript,
                str(SCRIPT),
                "--tree",
                str(TREE),
                "--metadata",
                str(METADATA),
                "--clade-support",
                str(CLADE_SUPPORT),
                "--out-prefix",
                str(prefix),
            ]
            completed = subprocess.run(command, text=True, capture_output=True, check=False)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            for suffix in (".svg", ".pdf", ".settings.tsv"):
                self.assertTrue(Path(f"{prefix}{suffix}").is_file())
            runtime_settings = read_settings(Path(f"{prefix}.settings.tsv"))
            self.assertEqual(runtime_settings["joint_strong_node_count"], "10")
            self.assertEqual(runtime_settings["below_joint_threshold_node_count"], "5")

            repeated = subprocess.run(command, text=True, capture_output=True, check=False)
            self.assertEqual(repeated.returncode, 2)
            self.assertIn("Refusing to overwrite", repeated.stderr)


if __name__ == "__main__":
    unittest.main()
