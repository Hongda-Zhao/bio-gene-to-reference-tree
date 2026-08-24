# Gene-to-Reference Tree

[![Validation](https://github.com/Hongda-Zhao/bio-gene-to-reference-tree/actions/workflows/validate.yml/badge.svg)](https://github.com/Hongda-Zhao/bio-gene-to-reference-tree/actions/workflows/validate.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**Build the reference set before you build the tree.**

[What it does](#what-it-does) · [Install](#install) · [Nucleotide support](#nucleotide-support) · [BRCA1 example](#brca1-example) · [Progressive Skill](#progressive-skill) · [Quick links](#quick-links)

## What it does

An open, portable Agent Skill for building an auditable **protein or nucleotide gene tree** from an accession, declared local sequence, or gene/protein/feature name plus organism.

It resolves the query, finds and reviews homolog or ortholog candidates, selects references and outgroups, plans molecule-specific MAFFT/trimAl and FastTree/IQ-TREE2 commands through available host tools, and produces iTOL or local ggtree/ggplot2 outputs. It avoids the unsafe shortcut of treating the top similarity hits as a ready-made reference set.

> **Execution model:** deterministic `plan`/`route` modes do not run workflow tools. `doctor` defaults to passive `PATH` discovery; live retrieval, active version probes, and bioinformatics execution require separately approved host or local capabilities.

| Stage | Main result |
|---|---|
| Query | Versioned record or hash-bound, explicitly declared local molecule with provenance |
| References | Selected and rejected candidates, reason codes, taxonomic balance, and outgroup rationale |
| Alignment | Raw MSA, conservation assessment, QC, and trimming sensitivity |
| Tree | Unrooted gene tree, optional approved rooted copy, model, and explicitly documented support method |
| Reporting | Sequence metadata, iTOL roles, optional SVG/PDF, and current-literature comparison |

The Skill reports a gene tree, not automatically a species tree; similarity is not treated as proof of orthology.

## Install

Browse the rendered Skill on [skills.sh](https://skills.sh/hongda-zhao/bio-gene-to-reference-tree/bio-gene-to-reference-tree), or install it with the third-party `skills` CLI:

```bash
npx skills add Hongda-Zhao/bio-gene-to-reference-tree \
  --skill bio-gene-to-reference-tree
```

For an explicit global client installation, add `--agent codex`, `--agent cursor`, or `--agent claude-code` together with `--global`.

| Client | Entry point |
|---|---|
| Codex | `$bio-gene-to-reference-tree ...` |
| Cursor | `/bio-gene-to-reference-tree ...` |
| Claude Code | `/bio-gene-to-reference-tree ...` |
| GitHub Copilot | Ask naturally; the matching project Skill is selected automatically |

See [detailed installation](docs/installation.md) for global commands, manual directories, project-level GitHub Copilot routing, Cursor discovery, Claude Marketplace installation, and telemetry settings.

Minimal example:

```text
$bio-gene-to-reference-tree Build an auditable ortholog protein tree for human BRCA1 NP_009225.1, retain the unrooted result, and evaluate amphibian outgroups.
```

## Nucleotide support

Protein, noncoding DNA/RNA, and clean CDS are routed explicitly; molecule type is never inferred from sequence letters. Noncoding loci use nucleotide alignment and DNA models, clean CDS uses translation-guided codon-preserving alignment, and disrupted CDS stops for reviewed MACSE handling.

```text
$bio-gene-to-reference-tree Build an auditable noncoding-dna gene tree from my resolved local FASTA, use only comparable homologous loci and an actual nucleotide search database, retain the unrooted tree, and report MAFFT --nuc plus IQ-TREE -st DNA commands.
```

See the [minimal nucleotide request](skills/bio-gene-to-reference-tree/assets/request.nucleotide.example.json) and [sequence-type router](skills/bio-gene-to-reference-tree/references/sequence-type-routing.md) for RNA/CDS details and hard stops.

## BRCA1 example

![Outgroup-rooted BRCA1 protein gene tree for 18 vertebrates](examples/brca1/figures/brca1-readme.svg)

The executed example starts from human RefSeq BRCA1 [`NP_009225.1`](https://www.ncbi.nlm.nih.gov/protein/NP_009225.1) and reviews a fixed 18-protein vertebrate set. It validates taxonomy and terminal RING/BRCT architecture, compares three trim profiles, and infers the primary tree with IQ-TREE2 ModelFinder plus 1,000 UFBoot2 and 1,000 SH-aLRT replicates.

| Result | Executed value |
|---|---|
| Sampling | 16 amniote ingroup proteins + 2 amphibian outgroups |
| Taxonomy/QC | 18/18 exact scientific-name/TaxID matches; 162/162 terminal-domain checks passed |
| Alignment | MAFFT E-INS-i; 2,179 raw columns; balanced trim retained 1,851 (84.95%) |
| Primary tree | IQ-TREE 2.4.0; `Q.bird+F+I+R3`; 1,000 SH-aLRT + 1,000 UFBoot2 |
| Rooting | Amphibian-rooted copy retained as a provisional display hypothesis |

Orange marks the focal human sequence, green marks added references, and gray marks outgroups. This remains a protein gene tree; the complete audit preserves commands, versions, hashes, warnings, non-promoted attempts, and execution reconciliation outside the homepage.

[full executed audit record](examples/brca1/README.md) ·
[detailed audit SVG](examples/brca1/figures/gene-tree.outgroup-rooted.ggtree.svg) ·
[detailed audit PDF](examples/brca1/figures/gene-tree.outgroup-rooted.ggtree.pdf) ·
[unrooted Newick](examples/brca1/tree/gene-tree.unrooted.nwk) ·
[rooted derivative](examples/brca1/tree/gene-tree.outgroup-rooted.nwk) ·
[metadata](examples/brca1/annotation/sequence_metadata.tsv) ·
[execution reconciliation](examples/brca1/report/execution_reconciliation.json) ·
[50-candidate expanded review](examples/brca1-expanded/README.md)

## Progressive Skill

The Agent layer translates a biological request into a reviewable, environment-matched plan. It loads only the active task guidance, preserves provenance and approval gates, and can hand approved work to local, browser, SSH, or HPC capabilities without changing the scientific plan.

```mermaid
flowchart TD
    U[Biological request] --> D{Client discovery}
    D -->|Codex, Cursor, Claude| C[Canonical SKILL.md]
    D -->|GitHub Copilot| A[Thin task adapter]
    A --> C
    C --> T[One active step]
    T --> P[Reviewable plan]
    P --> E[Approved local / browser / SSH / HPC]
    E --> O[Auditable tree and report]
```

The installable package follows the open [Agent Skills specification](https://agentskills.io/specification). GitHub Copilot gets thin project-level discovery adapters; every client then delegates to the same canonical router and progressively loaded task files:

```text
.
├── .github/
│   ├── copilot-instructions.md
│   └── skills/                          # GitHub Copilot task discovery
│       ├── bio-gene-to-reference-tree/SKILL.md
│       ├── gene-tree-query-and-scope/SKILL.md
│       ├── gene-tree-reference-curation/SKILL.md
│       ├── gene-tree-alignment-and-inference/SKILL.md
│       ├── gene-tree-visualization-reporting/SKILL.md
│       └── gene-tree-environment-routing/SKILL.md
└── skills/
    └── bio-gene-to-reference-tree/
        ├── SKILL.md                     # canonical router and core claims
        ├── agents/openai.yaml           # discovery metadata
        ├── references/
        │   ├── workflow.md              # states and review gates
        │   ├── steps/                   # load only the active task
        │   │   ├── 01-resolve-query.md
        │   │   ├── ...
        │   │   └── 10-compare-evidence-and-report.md
        │   ├── sequence-type-routing.md
        │   ├── environment-routing.md
        │   └── output-contract.md
        ├── scripts/                     # deterministic helpers
        └── assets/                      # minimal requests and data
```

The `.github/skills/` files are repository-local adapters, not separate installable packages. Task modules: [intake](skills/bio-gene-to-reference-tree/references/steps/01-resolve-query.md) · [sampling](skills/bio-gene-to-reference-tree/references/steps/04-select-references-and-outgroups.md) · [alignment](skills/bio-gene-to-reference-tree/references/steps/06-align-and-assess-conservation.md) · [tree inference](skills/bio-gene-to-reference-tree/references/steps/08-infer-root-and-check-tree.md) · [visualization](skills/bio-gene-to-reference-tree/references/steps/09-annotate-and-visualize.md) · [reporting](skills/bio-gene-to-reference-tree/references/steps/10-compare-evidence-and-report.md).

## Quick links

| Need | Go to |
|---|---|
| Install on Codex, Cursor, Claude Code, or GitHub Copilot | [Installation guide](docs/installation.md) |
| Read the agent entrypoint | [Canonical `SKILL.md`](skills/bio-gene-to-reference-tree/SKILL.md) |
| Use the repository task layer | [GitHub Copilot orchestrator](.github/skills/bio-gene-to-reference-tree/SKILL.md) |
| Follow states and approval gates | [Workflow](skills/bio-gene-to-reference-tree/references/workflow.md) |
| Choose protein, noncoding nucleotide, or clean CDS analysis | [Sequence-type router](skills/bio-gene-to-reference-tree/references/sequence-type-routing.md) |
| Match tasks to local, browser, or HPC capabilities | [Environment and software routing](skills/bio-gene-to-reference-tree/references/environment-routing.md) |
| Inspect artifacts and schemas | [Output contract](skills/bio-gene-to-reference-tree/references/output-contract.md) |
| Resolve exact NCBI scientific names/TaxIDs | [Taxonomy policy](skills/bio-gene-to-reference-tree/references/taxonomy-resolution.md) |
| Review recent MSA/trimming evidence | [Evidence guide](skills/bio-gene-to-reference-tree/references/recent-msa-trimming-evidence.md) |
| Render a local tree figure | [ggtree/ggplot2 guide](skills/bio-gene-to-reference-tree/references/ggtree-visualization.md) |
| Review privacy boundaries | [Privacy policy](PRIVACY.md) |
| Reuse or inspect the worked example | [BRCA1 audit](examples/brca1/README.md) |

MIT licensed. See [LICENSE](LICENSE).
