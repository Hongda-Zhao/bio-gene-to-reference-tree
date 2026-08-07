# iTOL annotation and literature context

Keep tree-tip styling, full sequence metadata, and external phylogenetic evidence as separate reproducible artifacts.

## Analysis roles

Assign exactly one display role to every selected tip:

| Role | Default color | Meaning |
|---|---|---|
| `study` | `#E69F00` | User-supplied or focal research sequence |
| `expanded` | `#009E73` | Reference sequence added by the workflow |
| `outgroup` | `#999999` | Approved sequence used only for rooting/context |

Keep biological `ingroup/outgroup` scope separate from provenance/display role. Preserve accession-to-tip mappings and use stable tool-safe tip IDs.

## Default iTOL file

Generate `itol_roles.txt` as `DATASET_COLORSTRIP` because it classifies individual tips even when roles are scattered across the topology. Use this shape:

```text
DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL	Sequence roles
COLOR	#000000
LEGEND_TITLE	Sequence role
LEGEND_COLORS	#E69F00	#009E73	#999999
LEGEND_LABELS	Study	Expanded	Outgroup
DATA
tip_1	#E69F00	Study
tip_2	#009E73	Expanded
tip_3	#999999	Outgroup
```

Generate `DATASET_RANGE` only after the final topology demonstrates that the requested leaves form a meaningful contiguous clade. Do not draw a range around scattered tips or use one to imply monophyly. Uploading to the iTOL service is a separate remote action that requires permission for unpublished trees or metadata.

Official templates:

- Color strip: <https://itol.embl.de/help/dataset_color_strip_template.txt>
- Range: <https://itol.embl.de/help/dataset_ranges_template.txt>
- Text labels: <https://itol.embl.de/help/dataset_text_template.txt>

## Sequence metadata table

Keep full information in `sequence_metadata.tsv`. Include one row per candidate, selected or rejected, with these fields when available:

```text
tip_id, accession, accession_version, analysis_role,
inclusion_status, selection_order, reason_codes,
taxon_id, species, lineage, clade,
gene_name, protein_name, relation,
orthology_source, orthology_evidence,
is_reviewed, is_canonical, is_fragment,
query_coverage, target_coverage, percent_identity, alignment_length,
sequence_length, bitscore, evalue, domain_architecture,
source_db, source_release, retrieved_at, retrieval_query_id,
cluster_id, cluster_representative, outgroup_rationale,
sequence_sha256, notes
```

Use empty values for unavailable metadata. Never invent a release, lineage, score, orthology call, or citation to fill the table.

## Literature evidence search

Search in this order:

1. phylogenomic or formal taxonomic studies directly covering the selected species;
2. studies covering the corresponding genus;
3. family-level studies;
4. order or broader-clade studies;
5. recent reviews and foundational analyses;
6. authoritative taxonomy and synthetic tree resources for context.

Prioritize direct taxonomic relevance, data scale, method transparency, and current recognized taxonomy before journal prestige. A recent high-impact article that does not sample the relevant lineage is weaker evidence than a directly relevant rigorous study.

For viruses, include current ICTV taxonomy and gene/segment-specific studies. Do not infer a virus species relationship from one recombinant or reassorted segment without qualification.

## Literature evidence table

Create `literature_evidence.tsv` only after retrieving real evidence. Use:

```text
citation_id, title, year, journal, doi, pmid,
taxon_rank, taxa_covered, evidence_type, data_type,
inference_method, model, topology_claim, directness,
conflicts, limitations, source_url, retrieved_at
```

Label `directness` as exact-species, genus, family, order, broader, or taxonomy-only. Include stable DOI/PMID links where available. Keep paraphrased topology claims concise and distinguish author conclusions from workflow inference.

## Gene-tree/species-tree comparison

Compare the inferred gene tree qualitatively with directly relevant phylogenomic relationships. Record agreements, unsupported relationships, and conflicts. Treat discordance as a biological or methodological signal that can arise from duplication/loss, incomplete lineage sorting, introgression, horizontal transfer, recombination, alignment error, model misspecification, or incorrect orthology.

Open Tree of Life synthesizes published phylogenies and taxonomy and is useful for discovery/context: <https://tree.opentreeoflife.org/about/open-tree-of-life>. NCBI Taxonomy Common Tree represents classification and is not itself a statistically inferred phylogenetic tree: <https://www.ncbi.nlm.nih.gov/books/NBK54428/>. Use current ICTV resources for formal virus taxonomy: <https://ictv.global/>.
