---
name: bio-gene-to-reference-tree
description: Orchestrate this repository's end-to-end auditable protein or nucleotide gene-tree workflow. Use when a request spans several stages or the active stage is unclear; select one task adapter, then follow the canonical Skill and its review gates.
---

# Gene-to-Reference Tree Orchestrator

This is a GitHub Copilot discovery adapter, not a separate protocol or installable package.

1. Read the [canonical Skill](../../../skills/bio-gene-to-reference-tree/SKILL.md) and [workflow gates](../../../skills/bio-gene-to-reference-tree/references/workflow.md).
2. Determine the earliest incomplete state.
3. Continue through exactly one matching adapter: [query and scope](../gene-tree-query-and-scope/SKILL.md), [reference curation](../gene-tree-reference-curation/SKILL.md), [alignment and inference](../gene-tree-alignment-and-inference/SKILL.md), or [visualization and reporting](../gene-tree-visualization-reporting/SKILL.md).
4. Use [environment routing](../gene-tree-environment-routing/SKILL.md) as a sidecar whenever software or compute placement must be selected.

Keep all scientific decisions, commands, and artifacts in the canonical package. Stop at its review gates.
