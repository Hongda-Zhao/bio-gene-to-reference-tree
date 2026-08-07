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
        self.assertTrue(content.startswith("---\n"))
        _opening, frontmatter, body = content.split("---", 2)
        keys = re.findall(r"^([a-z_][a-z0-9_-]*):", frontmatter, flags=re.MULTILINE)

        self.assertEqual(keys, ["name", "description"])
        self.assertIn("name: bio-gene-to-reference-tree", frontmatter)
        self.assertNotIn("TODO", content)
        self.assertGreater(len(body.strip()), 0)
        description_match = re.search(r"^description:\s*(.+)$", frontmatter, flags=re.MULTILINE)
        self.assertIsNotNone(description_match)
        description = description_match.group(1).strip()
        self.assertLessEqual(len(description), 1024)
        self.assertNotRegex(description, r"[<>]")

    def test_declared_resources_exist(self) -> None:
        expected = {
            "scripts/gene_to_tree.py",
            "assets/request.example.json",
            "assets/query.example.faa",
            "assets/candidates.example.faa",
            "assets/candidates.example.tsv",
            "references/workflow.md",
            "references/reference-selection.md",
            "references/output-contract.md",
            "references/tool-routing.md",
            "references/query-resolution.md",
            "references/alignment-and-tree.md",
            "references/itol-and-literature.md",
            "references/request-0.2.schema.json",
            "references/plan-0.2.schema.json",
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

    def test_v02_version_and_schema_surfaces_are_synchronized(self) -> None:
        script = (SKILL_ROOT / "scripts" / "gene_to_tree.py").read_text(encoding="utf-8")
        request = json.loads((SKILL_ROOT / "assets" / "request.example.json").read_text(encoding="utf-8"))
        readme = (REPOSITORY_ROOT / "README.md").read_text(encoding="utf-8")

        self.assertIn('VERSION = "0.2.0"', script)
        self.assertIn('OUTPUT_SCHEMA_VERSION = "0.2"', script)
        self.assertEqual(request["schema_version"], "0.2")
        self.assertIn("v0.2 review candidate", readme)


if __name__ == "__main__":
    unittest.main()
