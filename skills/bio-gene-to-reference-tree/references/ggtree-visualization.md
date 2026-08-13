# Local ggtree and ggplot2 visualization

Use the bundled R renderer for reproducible local figures and retain iTOL as the interactive/web annotation path. `ggtree` supplies tree geometry and metadata attachment; `ggplot2` supplies the composable graphics and export system.

## Contents

- [Choose this route](#choose-this-route)
- [Required inputs](#required-inputs)
- [Safe rendering defaults](#safe-rendering-defaults)
- [Run the renderer](#run-the-renderer)
- [Outputs and provenance](#outputs-and-provenance)

## Choose this route

Use `render_tree_ggtree.R` when the user wants a local scriptable figure, vector publication output, or further ggplot2 layers. Use iTOL when the user wants interactive exploration, very large-tree navigation, or template-based online annotation. Both routes consume the same stable tip IDs and role metadata.

Required local R packages are `ape`, `ggplot2`, `ggtree`, `openssl`, and `svglite`. Install `ggtree` from the Bioconductor release compatible with the active R version. The renderer checks packages but never installs anything or contacts a network endpoint.

Official package documentation:

- ggtree: <https://bioconductor.org/packages/release/bioc/html/ggtree.html>
- ggplot2 export: <https://ggplot2.tidyverse.org/reference/ggsave.html>

## Required inputs

- a reviewed Newick tree;
- `sequence_metadata.tsv` from the same approved reference set;
- optional `itol_roles.txt` to reuse its orange/green/gray palette exactly.

Join only by exact `tip_id`. The Newick tip set and selected metadata tip set must be identical. If iTOL roles are supplied, its DATA tip set and role labels must also be identical. Never strip accession versions, join on species names, fuzzy-match labels, or silently discard unmatched tips.

Rejected metadata rows remain in the audit table but are not drawn. Display labels may combine species and accession; accessions remain the stable tree identifiers.

## Safe rendering defaults

- layout: rectangular is the default and circular is an explicit alternative; for an unrooted tree both are display conveniences with an arbitrary root position/orientation, not an inferred root;
- ladderization: disabled, because visual ordering has no evolutionary direction;
- rooting: declared as `unrooted` or `outgroup-rooted`; the renderer never reroots, and an outgroup-rooted declaration requires a structurally rooted tree whose root split isolates exactly the selected outgroup tips;
- branches: `auto` uses a phylogram only when every edge length is finite, non-negative, and not all zero; otherwise it draws a cladogram;
- support: hidden unless its semantics are explicitly declared;
- colors: focal study orange `#E69F00`, expanded references green `#009E73`, outgroup gray `#999999`;
- encoding: color points at tips rather than entire branches, avoiding an unsupported ancestral-role claim;
- export: SVG and PDF vectors plus a settings TSV; no lossy publication PNG by default.

Support modes are deliberately separate:

- `fasttree-sh-like`: internal values must be numeric 0–1;
- `sh-alrt/ufboot`: internal labels must be `value/value`, each 0–100;
- `sh-alrt/bootstrap`: internal labels must be `value/value`, each 0–100;
- `none`: draw no node labels.

The script never guesses a support scale or converts 0.95 to 95.

## Run the renderer

For an unrooted IQ-TREE result:

```text
Rscript <skill-root>/scripts/render_tree_ggtree.R \
  --tree gene-tree.unrooted.treefile \
  --metadata sequence_metadata.tsv \
  --itol-roles itol_roles.txt \
  --out-prefix figures/gene-tree.unrooted.ggtree \
  --root-state unrooted \
  --layout rectangular \
  --branch-length auto \
  --support-format sh-alrt/ufboot \
  --show-tip-labels true \
  --width 10 --height 8
```

For a separately approved rooted copy, pass that file and `--root-state outgroup-rooted`. Do not declare an unrooted file rooted merely because it is displayed left-to-right.

Use `--layout circular` only when space or broad group structure benefits from it. With a rooted input it is a rooted circular tree; with `--root-state unrooted`, its root position and orientation remain arbitrary display choices and must not be interpreted biologically.

## Outputs and provenance

The prefix above creates:

```text
figures/gene-tree.unrooted.ggtree.svg
figures/gene-tree.unrooted.ggtree.pdf
figures/gene-tree.unrooted.ggtree.settings.tsv
```

The settings table records root declaration, layout, requested and actual branch mode, support semantics, ladderization, tip/node counts, canvas dimensions, palette source, R version, and package versions. The renderer writes through temporary files, refuses to overwrite any existing output, and removes incomplete temporary outputs on failure.

State in the figure caption whether the tree is rooted, whether branches represent substitutions per site or only topology, which support method is shown, that tip order was not ladderized, and what colors mean. A gene-tree figure remains a gene-tree claim, not a species-tree claim.
