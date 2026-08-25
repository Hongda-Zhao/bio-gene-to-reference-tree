---
name: gene-tree-reference-curation
description: Discover and review homolog or ortholog candidates, references, clustering, and outgroups for a protein or nucleotide gene tree. Use after query identity and scope are fixed and before sequence alignment begins.
---

# Gene-Tree Reference Curation

This repository-level adapter delegates its scientific rules to the canonical Skill.

1. Read the [core claims](../../../skills/bio-gene-to-reference-tree/SKILL.md) and [workflow gates](../../../skills/bio-gene-to-reference-tree/references/workflow.md).
2. Follow [Step 3 — Discover candidates](../../../skills/bio-gene-to-reference-tree/references/steps/03-discover-candidates.md).
3. Follow [Step 4 — Select references and outgroups](../../../skills/bio-gene-to-reference-tree/references/steps/04-select-references-and-outgroups.md).
4. Follow conditional [Step 5 — Cluster expanded candidates](../../../skills/bio-gene-to-reference-tree/references/steps/05-cluster-expanded-candidates.md) only when the canonical trigger applies, then return to Step 4.

After the linked review gate passes, hand off to the alignment-and-inference adapter.
