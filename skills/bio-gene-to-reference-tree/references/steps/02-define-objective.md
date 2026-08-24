# Step 2 — Define the biological objective

## When to read

Read after the query is resolved and before choosing databases, sampling breadth, paralog policy, or an outgroup.

## Required inputs

- the resolved molecule-matched query, comparable-region metadata, and any RNA/CDS QC provenance;
- the user's biological question;
- intended ingroup, taxonomic breadth, and comparison unit;
- speed/accuracy preference and whether rooting is required;
- current local/host/HPC capabilities, materialized inputs, and the environment route status;
- viral, within-species, noncoding-locus, RNA, or codon-specific context when applicable.

## Procedure

Choose exactly one primary objective:

| Objective | Intended claim | Default relationship policy |
|---|---|---|
| `ortholog-tree` | Compare corresponding genes across species | Prefer curated orthology; exclude paralogs unless explicitly allowed |
| `homolog-context` | Place a sequence in a broader family or subfamily | Retain labelled paralogs when they answer the question |
| `within-species` | Compare alleles, strains, isolates, or close copies | Preserve biologically meaningful copies without requiring multiple TaxIDs |

Declare the analysis unit and one analysis space: protein, direct noncoding nucleotide, or clean protein-guided codon. Do not change space to accommodate missing software. If the requested claim is a species tree, a duplication/loss model, HGT test, divergence date, positive-selection analysis, recombination analysis, RNA-structure analysis, disrupted-CDS analysis, or genome-scale phylogeny, route it to a dedicated reviewed workflow rather than silently broadening this Skill.

Record:

- ingroup definition and desired sampling depth;
- exact molecule type, homologous region, strand/orientation policy, and protein/nucleotide/codon interpretation;
- RNA source encoding and analysis-copy transformation, or CDS genetic-code/frame/translation requirements when applicable;
- whether the query is a focal `study` sequence or one member of a broader panel;
- whether paralogs, co-orthologs, fragments, or domain-only matches may be retained;
- target reference count, per-taxon cap, minimum coverage, and clustering trigger as project policies rather than universal constants;
- whether an unrooted result is acceptable and what would constitute a defensible outgroup;
- exploratory FastTree versus primary IQ-TREE2 intent;
- expected annotations, local figure, and literature comparison.

For `sequence_context: viral`, define the homologous gene/segment explicitly. Review recombination, reassortment, segmentation, and mosaic ancestry before assuming that one bifurcating tree can answer the question.

## Required outputs

- a reviewed objective and scope in the request/plan record;
- explicit relationship, paralog, fragment, coverage, taxon-balance, and clustering policies;
- declared rooting requirement and outgroup eligibility concept;
- declared execution mode, deliverables, and privacy constraints;
- selected molecule-compatible workflow tier, actual search database, compute target, `route_hash`, blockers/limitations, and any external handoff requiring later authorization;
- a routing note for any out-of-scope request.

## Review gate and stop conditions

Do not begin discovery while the objective, molecule/analysis space, comparable region, ingroup, or relationship policy is ambiguous. Stop if the requested inference is incompatible with a single gene tree, if viral recombination/segment identity has not been addressed, or if a protein↔nucleotide↔codon fallback is being proposed.

## Supporting references

- [Workflow states and decision gates](../workflow.md)
- [Environment-aware software and workflow routing](../environment-routing.md)
- [Sequence-type and analysis-space routing](../sequence-type-routing.md)
- [Tool and privacy boundaries](../tool-routing.md)
- [Request and output contract](../output-contract.md)
