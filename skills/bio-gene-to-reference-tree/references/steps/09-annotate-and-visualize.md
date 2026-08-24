# Step 9 — Annotate and visualize

## When to read

Read after tree inference when preparing iTOL files, full evolutionary metadata, or a local ggtree/ggplot2 figure.

## Required inputs

- approved unrooted tree and optional separately rooted copy;
- selected and rejected sequence metadata with stable `tip_id` mappings and exact molecule/analysis-space provenance;
- declared root state, branch-length meaning, and support format;
- permission decision for any remote iTOL upload.

## Procedure

Assign exactly one display/provenance role to each selected tip:

| Role | Color | Meaning |
|---|---|---|
| `study` | `#E69F00` | User-supplied or focal research sequence |
| `expanded` | `#009E73` | Reference added by the workflow |
| `outgroup` | `#999999` | Approved rooting/context sequence |

Keep display role separate from biological ingroup/outgroup scope. Generate an official iTOL `DATASET_COLORSTRIP` in `itol_roles.txt`; it remains valid even when role tips are scattered. Use `DATASET_RANGE` only after the final topology shows a requested set is a meaningful contiguous clade. Never use a range to manufacture monophyly.

Keep complete information in `sequence_metadata.tsv`, one row per selected or rejected candidate. Preserve accession/version, molecule and analysis space, feature/region/coordinates/strand, inclusion state, reason codes, taxonomy, gene/protein/feature labels, relationship and orthology evidence, review/canonical/fragment flags, similarity and coverage values, applicable domain architecture, actual search-database provenance, clustering, outgroup rationale, sequence hash, and notes. The planner records RNA source encoding plus source/analysis-copy hashes and, for CDS, code/frame/completeness/stop/frameshift/translation evidence and translation hashes. Backtranslation QC is truthfully marked `pending-execution`; the execution/report stage must replace it with reviewed pass/fail evidence from the real command before claiming a completed analysis. Missing values remain empty, never invented.

For local figures, run `scripts/render_tree_ggtree.R` only on approved local files. Require exact equality among Newick tips, selected metadata `tip_id` values, and optional iTOL DATA tips. The renderer must not reroot, ladderize, guess support scale, install packages, or contact the network. Prefer SVG/PDF and preserve renderer settings TSV.

Treat iTOL upload as a separate remote action. Obtain explicit permission before sending unpublished trees or metadata.

## Required outputs

- official `itol_roles.txt` color strip and optional topology-validated range dataset;
- complete molecule-aware `sequence_metadata.tsv` and accession-to-tip mapping;
- requested SVG/PDF figure plus settings TSV, or a reproducible rendering plan;
- root-state, branch-length, support-format, and tip-set validation records;
- upload permission/provenance if a remote service was used.

## Review gate and stop conditions

Stop on any tip-set mismatch, duplicate tip ID, mixed/undeclared molecule space, unknown role/color, undeclared root/support semantics, non-contiguous range request, or missing permission for upload. Visualization must not change the scientific tree or obscure RNA/CDS derivation.

## Supporting references

- [Complete ggtree/ggplot2 rendering contract](../ggtree-visualization.md)
- [Sequence-type provenance contract](../sequence-type-routing.md)
- [Annotation and metadata contract](../output-contract.md)
- [Workflow annotation/evidence gate](../workflow.md#annotation-and-evidence-gate)

Official iTOL templates: color strip <https://itol.embl.de/help/dataset_color_strip_template.txt>, ranges <https://itol.embl.de/help/dataset_ranges_template.txt>, and text labels <https://itol.embl.de/help/dataset_text_template.txt>.
