# Step 10 — Compare evidence and report

## When to read

Read after tree annotation to compare the gene tree with current evolutionary evidence and compile the final auditable report.

## Required inputs

- final unrooted tree, optional approved rooted copy, and support/model semantics;
- complete molecule-aware sequence metadata, comparable-region definition, and declared taxonomic scope;
- all approved hashes, commands, versions, warnings, and review decisions;
- access to current primary literature/taxonomy, or a recorded limitation.

## Procedure

Search evidence in this order:

1. phylogenomic or formal taxonomic studies directly covering selected species;
2. genus-level studies;
3. family-level studies;
4. order or broader-clade studies;
5. recent reviews and foundational analyses;
6. recognized taxonomy and synthetic-tree resources for context.

Prefer direct relevance, adequate data scale, method transparency, and current recognized taxonomy over journal prestige alone. Label evidence as `exact-species`, `genus`, `family`, `order`, `broader`, or `taxonomy-only`. For viruses, include current ICTV taxonomy and gene/segment-specific research; qualify recombinant or reassorted segments.

Create `literature_evidence.tsv` only from retrieved sources. Record DOI/PMID, year, journal, taxa, molecule/data architecture, evidence type, method/model, topology claim, directness, conflicts, limitations, URL, and retrieval time. When a source informs alignment/trimming, also retain version, parameters, `exact`/`partial` reporting state, `explicit-none`/`not-reported` trimming state, and retained sites/fraction where reported. A protein, direct-nucleotide, or codon precedent is comparable only to the same declared analysis architecture.

Compare the gene tree qualitatively with accepted species relationships. Record agreements, unsupported nodes, and conflicts. Discordance can reflect duplication/loss, incomplete lineage sorting, introgression, horizontal transfer, recombination, alignment error, model misspecification, or incorrect orthology; it is not automatic pipeline failure.

Compile the completion bundle defined in [output-contract.md](../output-contract.md). Call the result a gene tree and distinguish observed results from planned, inferred, indirect, or unavailable evidence.

## Required outputs

- `literature_evidence.tsv` with real citations and directness labels, or an explicit search limitation;
- gene-tree/species-relationship comparison with conflicts and alternative explanations;
- final report with molecule/analysis space, comparable-region and actual search-database provenance, query resolution, selection/rejections, clustering, MSA/QC, reviewed conservation assessment, inference, rooting, annotations, literature, decisions, versions, commands, and checksums;
- RNA source-encoding/analysis-copy receipt or CDS frame/code/translation/backtranslation/MACSE-review status when applicable;
- all native artifacts required by the completion contract.

## Review gate and stop conditions

Mark complete only when every requested artifact exists, molecule/region/derivation hashes match approved decisions, support and root semantics are explicit, citations were actually retrieved, and limitations are visible. Do not fabricate a citation, taxonomy claim, tool result, database release, translation, or cross-analysis result when a capability is unavailable.

## Supporting references

- [Recent MSA/trimming evidence catalog](../recent-msa-trimming-evidence.md)
- [Sequence-type and provenance routing](../sequence-type-routing.md)
- [Workflow completion contract](../workflow.md#completion-contract)
- [Final executed-report and artifact contract](../output-contract.md)
- [Tool, literature, and privacy boundaries](../tool-routing.md)

Context resources: Open Tree of Life <https://tree.opentreeoflife.org/about/open-tree-of-life>, NCBI Taxonomy Common Tree <https://www.ncbi.nlm.nih.gov/books/NBK54428/>, and ICTV <https://ictv.global/>. Use them according to the evidence limitations above.
