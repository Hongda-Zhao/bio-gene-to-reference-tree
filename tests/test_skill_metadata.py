"""Dependency-free checks for the portable Agent Skill package."""

from __future__ import annotations

import re
import unittest
import json
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SKILL_ROOT = REPOSITORY_ROOT / "skills" / "bio-gene-to-reference-tree"


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

    def test_progressive_disclosure_and_local_links(self) -> None:
        skill_document = SKILL_ROOT / "SKILL.md"
        skill_lines = skill_document.read_text(encoding="utf-8").splitlines()
        self.assertLessEqual(len(skill_lines), 500)

        local_link_pattern = re.compile(r"\[[^\]]+\]\(([^)]+)\)")
        documents = [skill_document, *sorted((SKILL_ROOT / "references").glob("*.md"))]
        for document in documents:
            content = document.read_text(encoding="utf-8")
            for raw_target in local_link_pattern.findall(content):
                if raw_target.startswith(("http://", "https://", "mailto:", "#")):
                    continue
                relative_target = raw_target.split("#", 1)[0]
                if not relative_target:
                    continue
                resolved_target = (document.parent / relative_target).resolve()
                self.assertTrue(
                    resolved_target.is_relative_to(SKILL_ROOT.resolve()),
                    f"Local link escapes the Skill directory: {document}: {raw_target}",
                )
                self.assertTrue(
                    resolved_target.exists(),
                    f"Broken local link: {document}: {raw_target}",
                )

        for reference in sorted((SKILL_ROOT / "references").glob("*.md")):
            lines = reference.read_text(encoding="utf-8").splitlines()
            if len(lines) > 100:
                self.assertIn(
                    "## Contents",
                    lines[:40],
                    f"Long reference needs an early table of contents: {reference}",
                )

    def test_declared_resources_exist(self) -> None:
        expected = {
            "scripts/gene_to_tree.py",
            "scripts/ncbi_taxonomy.py",
            "scripts/render_tree_ggtree.R",
            "assets/request.example.json",
            "assets/query.example.faa",
            "assets/candidates.example.faa",
            "assets/candidates.example.tsv",
            "references/workflow.md",
            "references/reference-selection.md",
            "references/output-contract.md",
            "references/tool-routing.md",
            "references/query-resolution.md",
            "references/taxonomy-resolution.md",
            "references/alignment-and-tree.md",
            "references/itol-and-literature.md",
            "references/ggtree-visualization.md",
            "references/request-0.2.schema.json",
            "references/plan-0.2.schema.json",
            "references/plan-0.3.schema.json",
            "agents/openai.yaml",
            "LICENSE",
        }
        missing = sorted(path for path in expected if not (SKILL_ROOT / path).is_file())
        self.assertEqual(missing, [])

    def test_codex_metadata_matches_the_skill(self) -> None:
        content = (SKILL_ROOT / "agents" / "openai.yaml").read_text(encoding="utf-8")
        self.assertIn('display_name: "Gene-to-Reference Tree"', content)
        self.assertRegex(content, r'short_description: "[^"\n]{25,64}"')
        self.assertIn("$bio-gene-to-reference-tree", content)

    def test_v03_workflow_and_v02_schema_surfaces_are_synchronized(self) -> None:
        script = (SKILL_ROOT / "scripts" / "gene_to_tree.py").read_text(encoding="utf-8")
        request = json.loads((SKILL_ROOT / "assets" / "request.example.json").read_text(encoding="utf-8"))
        readme = (REPOSITORY_ROOT / "README.md").read_text(encoding="utf-8")

        self.assertIn('VERSION = "0.3.0"', script)
        self.assertIn('OUTPUT_SCHEMA_VERSION = "0.3"', script)
        self.assertEqual(request["schema_version"], "0.2")
        self.assertIn("v0.3 review candidate", readme)

        request_schema = json.loads(
            (SKILL_ROOT / "references" / "request-0.2.schema.json").read_text(encoding="utf-8")
        )
        plan_schema = json.loads(
            (SKILL_ROOT / "references" / "plan-0.3.schema.json").read_text(encoding="utf-8")
        )
        self.assertIn("taxonomy", request_schema["properties"])
        self.assertIn("taxonomy_plan", plan_schema["required"])

    def test_public_discovery_surfaces_are_documented(self) -> None:
        readme = (REPOSITORY_ROOT / "README.md").read_text(encoding="utf-8")
        self.assertIn(
            "https://skills.sh/hongda-zhao/bio-gene-to-reference-tree/bio-gene-to-reference-tree",
            readme,
        )
        self.assertNotIn("https://skills.sh/b/", readme)
        self.assertIn("actions/workflows/validate.yml/badge.svg", readme)
        self.assertIn("npx skills add Hongda-Zhao/bio-gene-to-reference-tree", readme)
        self.assertIn("https://agentskills.io/specification", readme)


if __name__ == "__main__":
    unittest.main()
