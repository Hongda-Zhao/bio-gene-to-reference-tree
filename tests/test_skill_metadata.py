"""Dependency-free checks for the portable Agent Skill package."""

from __future__ import annotations

import csv
import json
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from datetime import date
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SKILL_ROOT = REPOSITORY_ROOT / "skills" / "bio-gene-to-reference-tree"
STEP_FILENAMES = (
    "01-resolve-query.md",
    "02-define-objective.md",
    "03-discover-candidates.md",
    "04-select-references-and-outgroups.md",
    "05-cluster-expanded-candidates.md",
    "06-align-and-assess-conservation.md",
    "07-trim-and-test-sensitivity.md",
    "08-infer-root-and-check-tree.md",
    "09-annotate-and-visualize.md",
    "10-compare-evidence-and-report.md",
)
STEP_REQUIRED_HEADINGS = (
    "## When to read",
    "## Required inputs",
    "## Procedure",
    "## Required outputs",
    "## Review gate and stop conditions",
    "## Supporting references",
)


class SkillPackageTests(unittest.TestCase):
    def test_frontmatter_uses_the_portable_common_subset(self) -> None:
        content = (SKILL_ROOT / "SKILL.md").read_text(encoding="utf-8")
        lines = content.splitlines()
        self.assertEqual(lines[0], "---")
        closing_index = lines[1:].index("---") + 1
        frontmatter = "\n".join(lines[1:closing_index])
        body = "\n".join(lines[closing_index + 1 :])
        keys = re.findall(r"^([a-z_][a-z0-9_-]*):", frontmatter, flags=re.MULTILINE)

        self.assertEqual(keys, ["name", "description"])
        self.assertEqual(frontmatter.count("name:"), 1)
        self.assertEqual(frontmatter.count("description:"), 1)
        self.assertNotIn("TODO", content)
        self.assertGreater(len(body.strip()), 0)

        name_match = re.search(r"^name:\s*(.+)$", frontmatter, flags=re.MULTILINE)
        self.assertIsNotNone(name_match)
        name = name_match.group(1).strip()
        self.assertLessEqual(len(name), 64)
        self.assertRegex(name, r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
        self.assertEqual(name, SKILL_ROOT.name)

        description_match = re.search(r"^description:\s*(.+)$", frontmatter, flags=re.MULTILINE)
        self.assertIsNotNone(description_match)
        description = description_match.group(1).strip()
        self.assertGreater(len(description), 0)
        self.assertLessEqual(len(description), 1024)
        self.assertNotRegex(description, r"[<>]")
        self.assertIn("gene-family conservation", description)
        self.assertIn("explicit taxonomic scope", description)

    def test_progressive_disclosure_and_local_links(self) -> None:
        skill_document = SKILL_ROOT / "SKILL.md"
        skill_lines = skill_document.read_text(encoding="utf-8").splitlines()
        self.assertLessEqual(len(skill_lines), 130)

        local_link_pattern = re.compile(r"\[[^\]]+\]\(([^)]+)\)")
        reference_documents = sorted((SKILL_ROOT / "references").rglob("*.md"))
        repository_documents = [
            REPOSITORY_ROOT / "README.md",
            *sorted((REPOSITORY_ROOT / "docs").rglob("*.md")),
        ]
        documents = [
            skill_document,
            *repository_documents,
            *reference_documents,
        ]
        link_graph: dict[Path, set[Path]] = {document.resolve(): set() for document in documents}

        def heading_anchors(document: Path) -> set[str]:
            anchors: set[str] = set()
            occurrences: dict[str, int] = {}
            for heading in re.findall(
                r"^#{1,6}\s+(.+?)\s*$",
                document.read_text(encoding="utf-8"),
                flags=re.MULTILINE,
            ):
                plain = re.sub(r"[`*_~]", "", heading).lower()
                slug = re.sub(r"[^\w\-\s]", "", plain)
                slug = re.sub(r"\s+", "-", slug.strip())
                count = occurrences.get(slug, 0)
                occurrences[slug] = count + 1
                anchors.add(slug if count == 0 else f"{slug}-{count}")
            return anchors

        for document in documents:
            content = document.read_text(encoding="utf-8")
            for raw_target in local_link_pattern.findall(content):
                if raw_target.startswith(("http://", "https://", "mailto:")):
                    continue
                relative_target, separator, fragment = raw_target.partition("#")
                resolved_target = (
                    (document.parent / relative_target).resolve()
                    if relative_target
                    else document.resolve()
                )
                allowed_root = (
                    SKILL_ROOT
                    if document.resolve().is_relative_to(SKILL_ROOT.resolve())
                    else REPOSITORY_ROOT
                )
                self.assertTrue(
                    resolved_target.is_relative_to(allowed_root.resolve()),
                    f"Local link escapes its package: {document}: {raw_target}",
                )
                self.assertTrue(
                    resolved_target.exists(),
                    f"Broken local link: {document}: {raw_target}",
                )
                if resolved_target.suffix.lower() == ".md":
                    link_graph.setdefault(document.resolve(), set()).add(resolved_target)
                    if separator and fragment:
                        self.assertIn(
                            fragment,
                            heading_anchors(resolved_target),
                            f"Broken Markdown anchor: {document}: {raw_target}",
                        )

        for reference in reference_documents:
            lines = reference.read_text(encoding="utf-8").splitlines()
            if len(lines) > 100:
                self.assertIn(
                    "## Contents",
                    lines[:40],
                    f"Long reference needs an early table of contents: {reference}",
                )

        reachable = {skill_document.resolve()}
        frontier = [skill_document.resolve()]
        while frontier:
            source = frontier.pop()
            for target in link_graph.get(source, set()):
                if target not in reachable:
                    reachable.add(target)
                    frontier.append(target)
        orphaned = sorted(
            str(path.relative_to(SKILL_ROOT))
            for path in reference_documents
            if path.resolve() not in reachable
        )
        self.assertEqual(orphaned, [], "Every reference must be reachable from SKILL.md")

    def test_step_router_is_complete_ordered_and_uniform(self) -> None:
        steps_directory = SKILL_ROOT / "references" / "steps"
        actual = tuple(path.name for path in sorted(steps_directory.glob("*.md")))
        self.assertEqual(actual, STEP_FILENAMES)
        self.assertEqual(len(list(SKILL_ROOT.rglob("SKILL.md"))), 1)

        skill = (SKILL_ROOT / "SKILL.md").read_text(encoding="utf-8")
        routed_steps = tuple(
            re.findall(r"\(references/steps/([^)#]+\.md)(?:#[^)]*)?\)", skill)
        )
        self.assertEqual(routed_steps, STEP_FILENAMES)

        readme = (REPOSITORY_ROOT / "README.md").read_text(encoding="utf-8")
        self.assertIn("## Progressive Skill", readme)
        readme_steps = tuple(
            re.findall(
                r"\(skills/bio-gene-to-reference-tree/references/steps/([^)#]+\.md)(?:#[^)]*)?\)",
                readme,
            )
        )
        self.assertEqual(readme_steps, STEP_FILENAMES)

        for step_number, filename in enumerate(STEP_FILENAMES, start=1):
            with self.subTest(step=filename):
                path = steps_directory / filename
                content = path.read_text(encoding="utf-8")
                self.assertEqual(len(re.findall(r"^#\s+", content, flags=re.MULTILINE)), 1)
                self.assertRegex(content, rf"(?m)^# Step {step_number} — .+$")
                self.assertLessEqual(len(content.splitlines()), 100)
                h2_headings = tuple(
                    re.findall(r"^##\s+.+$", content, flags=re.MULTILINE)
                )
                self.assertEqual(h2_headings, STEP_REQUIRED_HEADINGS)
                section_offsets = [content.index(heading) for heading in h2_headings]
                for index, heading in enumerate(h2_headings):
                    section_end = (
                        section_offsets[index + 1]
                        if index + 1 < len(section_offsets)
                        else len(content)
                    )
                    section_body = content[
                        section_offsets[index] + len(heading) : section_end
                    ].strip()
                    self.assertTrue(section_body, f"Empty task section: {path}: {heading}")
                self.assertEqual(
                    skill.count(f"(references/steps/{filename})"),
                    1,
                    "Each task must have one authoritative router entry",
                )

        retired_flat_references = (
            "query-resolution.md",
            "reference-selection.md",
            "alignment-and-tree.md",
            "itol-and-literature.md",
        )
        for filename in retired_flat_references:
            self.assertFalse((SKILL_ROOT / "references" / filename).exists())

        for exact_task_command in (
            "mmseqs easy-linclust",
            "iqtree2 -s",
            '["trimal", "-in"',
        ):
            self.assertNotIn(exact_task_command, skill)
            self.assertNotIn(exact_task_command, readme)

    def test_readme_is_a_concise_navigation_page(self) -> None:
        readme = (REPOSITORY_ROOT / "README.md").read_text(encoding="utf-8")
        self.assertEqual(
            tuple(re.findall(r"^## .+$", readme, flags=re.MULTILINE)),
            (
                "## What it does",
                "## Install",
                "## Progressive Skill",
                "## Nucleotide quick start",
                "## BRCA1 example",
                "## Quick links",
            ),
        )
        self.assertLessEqual(len(readme.splitlines()), 150)
        self.assertLessEqual(len(readme.encode("utf-8")), 12_000)

        for required_link in (
            "(skills/bio-gene-to-reference-tree/SKILL.md)",
            "(docs/installation.md)",
            "(examples/brca1/README.md)",
            "(examples/brca1/figures/brca1-readme.svg)",
            "(examples/brca1/report/execution_reconciliation.json)",
        ):
            self.assertIn(required_link, readme)

        for retired_heading in (
            "## Why this project exists",
            "## Capability matrix",
            "## Scientific guardrails",
            "## Privacy and safe use",
            "## Optional local tools",
        ):
            self.assertNotIn(retired_heading, readme)
        self.assertNotIn("**Protocol-deviation notice.**", readme)

    def test_declared_resources_exist(self) -> None:
        expected = {
            "scripts/gene_to_tree.py",
            "scripts/ncbi_taxonomy.py",
            "scripts/render_tree_ggtree.R",
            "assets/request.example.json",
            "assets/query.example.faa",
            "assets/candidates.example.faa",
            "assets/candidates.example.tsv",
            "assets/conservation-assessment.example.tsv",
            "assets/environment-profile.example.json",
            "assets/request.nucleotide.example.json",
            "assets/query.nucleotide.example.fna",
            "assets/candidates.nucleotide.example.fna",
            "assets/candidates.nucleotide.example.tsv",
            "references/workflow.md",
            "references/output-contract.md",
            "references/tool-routing.md",
            "references/environment-routing.md",
            "references/sequence-type-routing.md",
            "references/environment-profile-0.1.schema.json",
            "references/environment-profile-0.2.schema.json",
            "references/environment-snapshot-0.1.schema.json",
            "references/route-decision-0.1.schema.json",
            "references/route-decision-0.2.schema.json",
            "references/taxonomy-resolution.md",
            "references/recent-msa-trimming-evidence.md",
            "references/recent-msa-trimming-evidence.tsv",
            "references/ggtree-visualization.md",
            "references/request-0.2.schema.json",
            "references/request-0.3.schema.json",
            "references/plan-0.2.schema.json",
            "references/plan-0.3.schema.json",
            "references/plan-0.4.schema.json",
            "agents/openai.yaml",
            "LICENSE",
        }
        expected.update(f"references/steps/{filename}" for filename in STEP_FILENAMES)
        missing = sorted(path for path in expected if not (SKILL_ROOT / path).is_file())
        self.assertEqual(missing, [])
        expected_references = {
            Path(path) for path in expected if path.startswith("references/")
        }
        actual_references = {
            path.relative_to(SKILL_ROOT)
            for path in (SKILL_ROOT / "references").rglob("*")
            if path.is_file()
        }
        self.assertEqual(actual_references, expected_references)

    def test_conservation_assessment_template_is_real_tsv(self) -> None:
        template = SKILL_ROOT / "assets" / "conservation-assessment.example.tsv"
        raw = template.read_text(encoding="utf-8")
        self.assertTrue(raw.endswith("\n"))
        self.assertNotIn("\r", raw)
        rows = list(csv.DictReader(raw.splitlines(), delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(
            set(rows[0]),
            {
                "analysis_unit",
                "conservation_class",
                "conservation_scope",
                "conservation_basis",
                "evidence_ids",
                "assessment_status",
                "limitations",
            },
        )
        row = rows[0]
        self.assertEqual(row["conservation_class"], "broadly-conserved-gene-family")
        self.assertEqual(row["assessment_status"], "provisional")
        self.assertIn("Template only", row["limitations"])
        self.assertIn("NP_009225.1", row["evidence_ids"])

    def test_recent_msa_trimming_evidence_catalog_is_auditable(self) -> None:
        catalog = SKILL_ROOT / "references" / "recent-msa-trimming-evidence.tsv"
        raw = catalog.read_text(encoding="utf-8")
        self.assertNotIn("\r", raw)
        self.assertNotIn("\x00", raw)

        rows = list(csv.DictReader(raw.splitlines(), delimiter="\t"))
        self.assertTrue(rows)
        required_columns = {
            "analysis_id",
            "citation_id",
            "publication_date",
            "title",
            "journal",
            "doi",
            "pmid",
            "pmcid",
            "broad_group",
            "taxon_scope",
            "gene_or_markers",
            "conservation_class",
            "conservation_scope",
            "conservation_basis",
            "molecule_type",
            "dataset_scale",
            "msa_tool",
            "msa_version",
            "msa_parameters",
            "msa_reporting_status",
            "trimming_method",
            "trimming_version",
            "trimming_parameters",
            "trimming_status",
            "tree_method",
            "evidence_location",
            "source_url",
            "retrieved_at",
            "curator_note",
        }
        self.assertEqual(set(rows[0]), required_columns)
        self.assertGreaterEqual(len(rows), 50)
        self.assertGreaterEqual(len({row["citation_id"] for row in rows}), 40)

        analysis_ids = [row["analysis_id"] for row in rows]
        self.assertEqual(len(analysis_ids), len(set(analysis_ids)))
        for analysis_id in analysis_ids:
            self.assertRegex(analysis_id, r"^[a-z0-9]+(?:-[a-z0-9]+)*$")

        window_start = date(2023, 8, 24)
        window_end = date(2026, 8, 24)
        required_values = {
            "citation_id",
            "title",
            "journal",
            "taxon_scope",
            "gene_or_markers",
            "conservation_class",
            "conservation_scope",
            "conservation_basis",
            "molecule_type",
            "dataset_scale",
            "msa_tool",
            "msa_version",
            "msa_parameters",
            "trimming_method",
            "trimming_version",
            "trimming_parameters",
            "tree_method",
            "evidence_location",
            "curator_note",
        }
        bibliographic_fields = (
            "publication_date",
            "title",
            "journal",
            "doi",
            "pmid",
            "pmcid",
            "source_url",
        )
        citations: dict[str, tuple[str, ...]] = {}

        for row in rows:
            self.assertNotIn(None, row)
            for field in required_values:
                self.assertTrue(row[field].strip(), f"Empty {field}: {row['analysis_id']}")

            publication_date = date.fromisoformat(row["publication_date"])
            self.assertGreaterEqual(publication_date, window_start)
            self.assertLessEqual(publication_date, window_end)
            self.assertEqual(date.fromisoformat(row["retrieved_at"]), window_end)
            self.assertIn(row["msa_reporting_status"], {"exact", "partial"})
            self.assertIn(
                row["trimming_status"],
                {"exact", "partial", "explicit-none", "not-reported"},
            )
            if row["msa_reporting_status"] == "exact":
                msa_report = f"{row['msa_version']} {row['msa_parameters']}".lower()
                self.assertNotIn("not reported", msa_report)
                self.assertNotIn("not restated", msa_report)
                self.assertNotEqual(row["msa_version"], "NA")
            if row["trimming_status"] == "exact":
                trimming_report = (
                    f"{row['trimming_method']} {row['trimming_version']} "
                    f"{row['trimming_parameters']}"
                ).lower()
                self.assertNotIn("not reported", trimming_report)
                self.assertNotIn("not restated", trimming_report)
                self.assertNotEqual(row["trimming_version"], "NA")
            elif row["trimming_status"] == "explicit-none":
                self.assertEqual(row["trimming_method"], "None")
                self.assertEqual(row["trimming_version"], "NA")
                self.assertIn("explicit", row["trimming_parameters"].lower())
            elif row["trimming_status"] == "not-reported":
                self.assertEqual(row["trimming_method"], "not reported")
                self.assertEqual(row["trimming_version"], "NA")
                self.assertRegex(
                    row["trimming_parameters"].lower(),
                    r"(?:not reported|no .* reported)",
                )
            self.assertTrue(row["source_url"].startswith("https://"))
            self.assertTrue(
                any(row[field] != "NA" for field in ("doi", "pmid", "pmcid")),
                f"No stable identifier: {row['analysis_id']}",
            )
            if row["doi"] != "NA":
                self.assertTrue(row["doi"].startswith("10."))
            if row["pmid"] != "NA":
                self.assertRegex(row["pmid"], r"^[0-9]+$")
            if row["pmcid"] != "NA":
                self.assertRegex(row["pmcid"], r"^PMC[0-9]+$")

            citation = tuple(row[field] for field in bibliographic_fields)
            previous = citations.setdefault(row["citation_id"], citation)
            self.assertEqual(previous, citation)

        broad_groups = {row["broad_group"] for row in rows}
        self.assertTrue(
            {"animals", "plants", "fungi", "protists", "bacteria", "archaea", "viruses"}
            <= broad_groups
        )
        msa_tools = " ".join(row["msa_tool"] for row in rows)
        trimming_tools = " ".join(row["trimming_method"] for row in rows)
        for tool in ("MAFFT", "MUSCLE", "Clustal Omega"):
            self.assertIn(tool, msa_tools)
        for tool in ("trimAl", "ClipKIT", "BMGE", "Gblocks"):
            self.assertIn(tool, trimming_tools)
        trimming_states = {row["trimming_status"] for row in rows}
        self.assertIn("explicit-none", trimming_states)
        self.assertIn("not-reported", trimming_states)

        conservation_classes = {
            "universal-or-deep-core",
            "clade-conserved-marker",
            "broadly-conserved-gene-family",
            "variable-multigene-family",
            "lineage-specific-or-rapidly-evolving",
            "mixed-conservation-panel",
            "not-applicable",
        }
        self.assertEqual(
            {row["conservation_class"] for row in rows}, conservation_classes
        )
        for row in rows:
            self.assertNotIn(
                row["conservation_class"],
                {"conserved", "non-conserved", "nonconserved"},
            )
            self.assertGreaterEqual(len(row["conservation_scope"].split()), 2)
            self.assertGreaterEqual(len(row["conservation_basis"].split()), 4)

        by_analysis = {row["analysis_id"]: row for row in rows}
        expected_routing = {
            "myrmicini-uce": "clade-conserved-marker",
            "eukaryote-ccm-enzyme-trees": "broadly-conserved-gene-family",
            "early-hexapod-receptor-sensitivity": "variable-multigene-family",
            "ciliate-lgt-proteins": "lineage-specific-or-rapidly-evolving",
            "ska2-gubbins-recomb": "not-applicable",
            "alifilter-cyanobacteria-benchmark": "not-applicable",
            "cloak-mammal-ortholog-benchmark": "not-applicable",
            "asgard-phylome-ensemble": "mixed-conservation-panel",
        }
        for analysis_id, expected_class in expected_routing.items():
            self.assertEqual(
                by_analysis[analysis_id]["conservation_class"], expected_class
            )

        pseudoalignments = [
            row for row in rows if "pseudoalignment" in row["molecule_type"].lower()
        ]
        self.assertTrue(pseudoalignments)
        for row in pseudoalignments:
            self.assertNotIn("MAFFT", row["msa_tool"])
            self.assertNotIn("trimAl", row["trimming_method"])
            self.assertIn("pseudoalignment", row["curator_note"].lower())
            self.assertEqual(row["conservation_class"], "not-applicable")

    def test_codex_metadata_matches_the_skill(self) -> None:
        content = (SKILL_ROOT / "agents" / "openai.yaml").read_text(encoding="utf-8")
        self.assertIn('display_name: "Gene-to-Reference Tree"', content)
        self.assertRegex(content, r'short_description: "[^"\n]{25,64}"')
        self.assertIn("$bio-gene-to-reference-tree", content)
        self.assertIn("conservation", content)
        self.assertIn("explicit taxonomic scope", content)

    def test_conservation_assessment_is_a_hash_bound_review_artifact(self) -> None:
        alignment_step = (
            SKILL_ROOT
            / "references"
            / "steps"
            / "06-align-and-assess-conservation.md"
        ).read_text(encoding="utf-8")
        workflow = (SKILL_ROOT / "references" / "workflow.md").read_text(encoding="utf-8")
        contract = (SKILL_ROOT / "references" / "output-contract.md").read_text(
            encoding="utf-8"
        )
        for content in (alignment_step, workflow, contract):
            self.assertIn("conservation_assessment.tsv", content)
        for field in (
            "conservation_class",
            "conservation_scope",
            "conservation_basis",
            "evidence_ids",
            "assessment_status",
            "limitations",
        ):
            self.assertIn(field, contract)
        self.assertIn("SHA-256", contract)
        self.assertIn("reopen", contract)
        expected_header = "\t".join(
            (
                "analysis_unit",
                "conservation_class",
                "conservation_scope",
                "conservation_basis",
                "evidence_ids",
                "assessment_status",
                "limitations",
            )
        )
        self.assertIn(expected_header, contract)
        self.assertNotIn("analysis_unit, conservation_class", contract)

    def test_v04_workflow_and_current_schema_surfaces_are_synchronized(self) -> None:
        script = (SKILL_ROOT / "scripts" / "gene_to_tree.py").read_text(encoding="utf-8")
        protein_request = json.loads(
            (SKILL_ROOT / "assets" / "request.example.json").read_text(encoding="utf-8")
        )
        nucleotide_request = json.loads(
            (SKILL_ROOT / "assets" / "request.nucleotide.example.json").read_text(
                encoding="utf-8"
            )
        )

        self.assertIn('VERSION = "0.4.0"', script)
        self.assertIn('OUTPUT_SCHEMA_VERSION = "0.4"', script)
        self.assertIn('ENVIRONMENT_PROFILE_SCHEMA_VERSION = "0.2"', script)
        self.assertIn('ROUTE_DECISION_SCHEMA_VERSION = "0.2"', script)
        self.assertEqual(protein_request["schema_version"], "0.2")
        self.assertEqual(nucleotide_request["schema_version"], "0.3")

        legacy_request_schema = json.loads(
            (SKILL_ROOT / "references" / "request-0.2.schema.json").read_text(encoding="utf-8")
        )
        current_request_schema = json.loads(
            (SKILL_ROOT / "references" / "request-0.3.schema.json").read_text(encoding="utf-8")
        )
        current_plan_schema = json.loads(
            (SKILL_ROOT / "references" / "plan-0.4.schema.json").read_text(encoding="utf-8")
        )
        current_environment_schema = json.loads(
            (SKILL_ROOT / "references" / "environment-profile-0.2.schema.json").read_text(
                encoding="utf-8"
            )
        )
        current_route_schema = json.loads(
            (SKILL_ROOT / "references" / "route-decision-0.2.schema.json").read_text(
                encoding="utf-8"
            )
        )
        legacy_plan_schema = json.loads(
            (SKILL_ROOT / "references" / "plan-0.3.schema.json").read_text(encoding="utf-8")
        )
        self.assertIn("taxonomy", legacy_request_schema["properties"])
        self.assertIn("molecule", current_request_schema["required"])
        self.assertIn("molecule_plan", current_plan_schema["required"])
        self.assertEqual(
            current_environment_schema["properties"]["schema_version"]["const"],
            "0.2",
        )
        self.assertEqual(
            current_route_schema["properties"]["schema_version"]["const"], "0.2"
        )
        self.assertIn("taxonomy_plan", legacy_plan_schema["required"])

    def test_portable_bundle_copies_to_codex_cursor_and_claude_paths(self) -> None:
        source_files = {
            path.relative_to(SKILL_ROOT)
            for path in SKILL_ROOT.rglob("*")
            if path.is_file()
            and "__pycache__" not in path.parts
            and path.suffix != ".pyc"
        }
        self.assertTrue(source_files)
        self.assertFalse(any(path.is_symlink() for path in SKILL_ROOT.rglob("*")))

        vendor_skill_copies: list[Path] = []
        for vendor_root_name in (".agents", ".codex", ".cursor", ".claude"):
            vendor_root = REPOSITORY_ROOT / vendor_root_name
            if vendor_root.exists():
                vendor_skill_copies.extend(vendor_root.rglob("SKILL.md"))
        self.assertEqual(vendor_skill_copies, [])

        destinations = (
            Path(".agents/skills/bio-gene-to-reference-tree"),
            Path(".codex/skills/bio-gene-to-reference-tree"),
            Path(".cursor/skills/bio-gene-to-reference-tree"),
            Path(".claude/skills/bio-gene-to-reference-tree"),
        )
        local_link_pattern = re.compile(r"\[[^\]]+\]\(([^)]+)\)")

        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary_root = Path(temporary_directory)
            for relative_destination in destinations:
                with self.subTest(destination=str(relative_destination)):
                    destination = temporary_root / relative_destination
                    shutil.copytree(
                        SKILL_ROOT,
                        destination,
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
                    )
                    copied_files = {
                        path.relative_to(destination)
                        for path in destination.rglob("*")
                        if path.is_file()
                    }
                    self.assertEqual(copied_files, source_files)
                    for relative_path in source_files:
                        self.assertEqual(
                            (destination / relative_path).read_bytes(),
                            (SKILL_ROOT / relative_path).read_bytes(),
                        )

                    for document in destination.rglob("*.md"):
                        content = document.read_text(encoding="utf-8")
                        for raw_target in local_link_pattern.findall(content):
                            if raw_target.startswith(("http://", "https://", "mailto:", "#")):
                                continue
                            relative_target = raw_target.split("#", 1)[0]
                            if not relative_target:
                                continue
                            resolved_target = (document.parent / relative_target).resolve()
                            self.assertTrue(resolved_target.is_relative_to(destination.resolve()))
                            self.assertTrue(resolved_target.exists())

                    help_result = subprocess.run(
                        [
                            sys.executable,
                            str(destination / "scripts" / "gene_to_tree.py"),
                            "--help",
                        ],
                        cwd=destination,
                        capture_output=True,
                        text=True,
                        timeout=10,
                        check=False,
                    )
                    self.assertEqual(help_result.returncode, 0, help_result.stderr)
                    self.assertIn("usage:", help_result.stdout.lower())

    def test_claude_marketplace_wraps_only_the_canonical_skill(self) -> None:
        marketplace_path = REPOSITORY_ROOT / ".claude-plugin" / "marketplace.json"
        marketplace = json.loads(marketplace_path.read_text(encoding="utf-8"))

        self.assertEqual(marketplace["name"], "hongda-zhao-bio-skills")
        self.assertEqual(marketplace["owner"]["name"], "Hongda Zhao")
        self.assertEqual(len(marketplace["plugins"]), 1)
        plugin = marketplace["plugins"][0]
        self.assertEqual(plugin["name"], "bio-gene-to-reference-tree")
        self.assertEqual(plugin["source"], "./")
        self.assertIs(plugin["strict"], False)
        self.assertEqual(
            plugin["skills"],
            ["./skills/bio-gene-to-reference-tree"],
        )

        referenced_skill = (REPOSITORY_ROOT / plugin["skills"][0]).resolve()
        self.assertEqual(referenced_skill, SKILL_ROOT.resolve())
        self.assertTrue((referenced_skill / "SKILL.md").is_file())
        self.assertFalse((REPOSITORY_ROOT / ".claude-plugin" / "skills").exists())

    def test_public_discovery_surfaces_are_documented(self) -> None:
        readme = (REPOSITORY_ROOT / "README.md").read_text(encoding="utf-8")
        installation = (REPOSITORY_ROOT / "docs" / "installation.md").read_text(
            encoding="utf-8"
        )
        self.assertIn(
            "https://skills.sh/hongda-zhao/bio-gene-to-reference-tree/bio-gene-to-reference-tree",
            readme,
        )
        self.assertNotIn("https://skills.sh/b/", readme)
        self.assertIn("actions/workflows/validate.yml/badge.svg", readme)
        self.assertIn("npx skills add Hongda-Zhao/bio-gene-to-reference-tree", readme)
        self.assertIn("https://agentskills.io/specification", readme)
        self.assertIn("$bio-gene-to-reference-tree", readme)
        self.assertIn("/bio-gene-to-reference-tree", readme)
        self.assertIn("(docs/installation.md)", readme)

        self.assertIn("--agent codex", installation)
        self.assertIn("--agent cursor", installation)
        self.assertIn("--agent claude-code", installation)
        self.assertIn("https://cursor.com/docs/skills", installation)
        self.assertIn("https://code.claude.com/docs/en/skills", installation)
        self.assertIn("https://code.claude.com/docs/en/plugin-marketplaces", installation)
        self.assertIn(".codex/skills/bio-gene-to-reference-tree/", installation)
        self.assertIn("~/.codex/skills/bio-gene-to-reference-tree/", installation)
        self.assertIn(".cursor/skills/bio-gene-to-reference-tree/", installation)
        self.assertIn("~/.cursor/skills/bio-gene-to-reference-tree/", installation)
        self.assertIn(".claude/skills/bio-gene-to-reference-tree/", installation)
        self.assertIn("~/.claude/skills/bio-gene-to-reference-tree/", installation)
        self.assertIn("Remote Rule (Github)", installation)
        self.assertIn("Cursor 2.4+", installation)
        self.assertIn("Custom Mode", installation)
        self.assertIn("complete `skills/bio-gene-to-reference-tree/` package", installation)
        self.assertIn("a Claude Code project installation belongs in", installation)
        self.assertIn(
            "/plugin marketplace add Hongda-Zhao/bio-gene-to-reference-tree",
            installation,
        )
        self.assertIn(
            "/bio-gene-to-reference-tree:bio-gene-to-reference-tree", installation
        )
        self.assertIn(
            "no forked Cursor- or Claude-specific prompt is required", installation
        )

        tool_routing = (
            SKILL_ROOT / "references" / "tool-routing.md"
        ).read_text(encoding="utf-8")
        self.assertIn("Codex, Cursor, Claude Code", tool_routing)
        self.assertFalse((SKILL_ROOT / "CLAUDE.md").exists())
        self.assertFalse((SKILL_ROOT / ".cursor").exists())
        self.assertFalse((SKILL_ROOT / ".cursorrules").exists())


if __name__ == "__main__":
    unittest.main()
