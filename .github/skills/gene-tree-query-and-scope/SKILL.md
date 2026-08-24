---
name: gene-tree-query-and-scope
description: Resolve a protein, nucleotide, or CDS accession, local sequence, or name-plus-organism query and define its gene-tree objective. Use for molecule identity, provenance, exact taxonomy, analysis space, ingroup, or deliverable decisions before homolog discovery.
---

# Gene-Tree Query and Scope

This repository-level adapter delegates its scientific rules to the canonical Skill.

1. Read the [core claims](../../../skills/bio-gene-to-reference-tree/SKILL.md) and [workflow gates](../../../skills/bio-gene-to-reference-tree/references/workflow.md).
2. Follow [Step 1 — Resolve the query](../../../skills/bio-gene-to-reference-tree/references/steps/01-resolve-query.md).
3. After Gate 1 passes, follow [Step 2 — Define the objective](../../../skills/bio-gene-to-reference-tree/references/steps/02-define-objective.md).
4. Load [sequence-type routing](../../../skills/bio-gene-to-reference-tree/references/sequence-type-routing.md) or [taxonomy resolution](../../../skills/bio-gene-to-reference-tree/references/taxonomy-resolution.md) only when its canonical trigger applies.

When both linked exit conditions are satisfied, hand off to the reference-curation adapter.
