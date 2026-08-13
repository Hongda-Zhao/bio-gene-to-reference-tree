---
name: bio-gene-to-reference-tree
description: Build an auditable protein gene tree from an accession, raw amino-acid sequence, or protein or gene name plus source organism. Use when an agent must resolve query metadata, discover and curate ortholog or homolog references and an outgroup, cluster redundant candidates, align and trim proteins, run FastTree or IQ-TREE2, generate iTOL annotations and metadata tables, or compare a gene tree with current phylogenetic literature. Require review gates before reference selection and tree inference.
---

# Gene to Reference Tree

Build a protein gene tree through explicit, reviewable decisions. Treat database acquisition, reference sampling, alignment, rooting, and literature comparison as scientific analyses rather than mechanical top-hit processing.

## Preserve the core claim

- Call the result a **gene tree**, not a species tree.
- Treat similarity as evidence of homology, never as proof of orthology or identical function.
- Treat support as repeatability under a method, not proof that the topology is correct.
- Preserve the raw candidate pool, untrimmed alignment, unrooted tree, commands, versions, database snapshots, and every exclusion.
- Never upload an unpublished sequence, candidate set, tree, or metadata to a remote service without explicit permission.

## Run the workflow

Read [workflow.md](references/workflow.md) before the first run. Revisit it whenever a state transition, approval gate, invalidated hash, or failure condition affects what may run next.

### 1. Classify and resolve the query

Classify the input as:

- a versioned NCBI, RefSeq, UniProt, Ensembl, or other accession;
- a raw amino-acid sequence;
- a protein or gene name plus source organism or TaxID.

Resolve an accession against its authoritative namespace and record status, version, sequence, organism, TaxID, lineage, gene/protein names, source release, retrieval time, and sequence SHA-256. For a name, require an organism or TaxID and use literature only to disambiguate identity and context; retrieve the sequence from an authoritative database record. For a raw sequence, validate the alphabet, length, low complexity, and likely molecule type before searching.

Before any live lookup, including a public accession or symbol lookup, record whether remote queries are allowed. If they are not, stop and request a local resolved FASTA/metadata handoff. Obtain separate explicit permission before uploading an unpublished sequence.

Read [query-resolution.md](references/query-resolution.md) before resolving an accession, name, raw sequence, CDS, or viral query.

### 2. Define the biological objective

Choose one objective explicitly:

- `ortholog-tree`: compare corresponding genes across species; prefer curated orthology and exclude paralogs by default;
- `homolog-context`: place a sequence in a broader family; retain labelled paralogs when they answer the question;
- `within-species`: compare alleles, strains, isolates, or closely related copies without requiring multiple TaxIDs.

Record `sequence_context: viral` separately when applicable. For viral data, require a segment/gene definition and check recombination, reassortment, segmentation, and mosaic ancestry before interpreting a single tree.

### 3. Discover and annotate candidates

Use the objective-specific route in [tool-routing.md](references/tool-routing.md):

1. Prefer curated ortholog resources for an ortholog tree.
2. Search RefSeq protein first for a raw protein.
3. Escalate unresolved cases through reviewed UniProt/Swiss-Prot, broader UniProtKB or `nr`, then profile/domain or structure-aware methods.
4. Retrieve substantially more candidates than the final tree needs.

Record query and target coverage, percent identity, alignment span, E-value, bit score, domain architecture, orthology evidence, accession version, taxonomy, database release, retrieval time, and the raw provider response or its checksum. Never select only BLAST top N.

### 4. Select references and outgroup candidates

Apply [reference-selection.md](references/reference-selection.md). Balance references across the declared taxonomic scope and keep explicit inclusion/exclusion reason codes. Prefer a nearby homologous sister lineage outside the ingroup as an outgroup; never choose the lowest-scoring or most distant hit automatically. Keep two or more outgroup candidates when feasible and retain an unrooted interpretation if no defensible outgroup exists.

Stop for **reference approval** after showing selected and rejected accessions, taxa, paralogs, fragments, domain warnings, sampling gaps, and outgroup rationales.

### 5. Cluster only when needed

If the expanded candidate pool crosses the declared trigger, cluster only the `expanded` group with MMseqs2. Preserve every `study` and `outgroup` sequence outside clustering. Specify both sequence identity and coverage; `-c` alone is not a redundancy threshold.

Use this full-length starting profile unless the objective justifies another value:

```text
mmseqs easy-linclust expanded_candidates.faa clusters mmseqs_tmp \
  --min-seq-id 0.95 -c 0.8 --cov-mode 0 --threads <fixed>
```

Keep the representative/member mapping, annotate `cluster_id`, and re-run reference planning. Do not silently erase accessions or distinct taxa.

### 6. Align and inspect

Use MAFFT and choose the mode from sequence count and architecture, not divergence alone:

- `auto` for general routing;
- L-INS-i for a small set with one alignable domain and difficult flanks;
- G-INS-i for globally alignable full-length proteins;
- E-INS-i for conserved motifs separated by long insertions, when motif order is shared.

Inspect coverage, gap fraction, occupancy, conserved motifs, mixed domains, fragments, fusions, duplicate tip IDs, and suspicious long branches. Preserve `alignment.raw.faa`. Read [alignment-and-tree.md](references/alignment-and-tree.md) before choosing or running MAFFT, trimAl, FastTree, or IQ-TREE2.

### 7. Treat trimming as a sensitivity analysis

Use trimAl profiles with explicit `-gt` semantics. A value of `0.98` is extremely strict; `0.10` or `0.05` is extremely permissive. Never infer a threshold solely from “close,” “distant,” or “viral.”

Retain each profile, report columns removed and retained fraction, verify conserved regions, and compare key topology when profiles differ. Stop if trimming removes too much information or changes the biological conclusion. Obtain **alignment/trimming approval** before tree inference.

### 8. Infer and label support correctly

- Use FastTree only for an exploratory approximate-ML tree. Label its default node values as SH-like local support, not global bootstrap.
- Use IQ-TREE2 for the primary accurate workflow. Default to `-m MFP -B 1000 -bnni -alrt 1000` with fixed threads and seed.
- Use `-b 1000`, not `-B 1000`, only when the user explicitly requests standard nonparametric bootstrap.
- Consider C60/PMSF or other richer models for deep, heterogeneous, long-branch-prone protein data.

Retain the unrooted tree. Produce a separate rooted copy only from an approved outgroup. Never interpret high support as immunity to alignment error, model misspecification, or long-branch attraction.

### 9. Generate iTOL and metadata outputs

Generate an official `DATASET_COLORSTRIP` file by default:

- `study`: orange `#E69F00`;
- `expanded`: green `#009E73`;
- `outgroup`: gray `#999999`.

Generate `DATASET_RANGE` only after the final topology shows that a requested group is a meaningful contiguous clade; do not use a range to imply monophyly. Keep full evolutionary metadata in TSV rather than overloading tree labels. Read [itol-and-literature.md](references/itol-and-literature.md) before generating iTOL/metadata outputs or beginning the literature comparison.

### 10. Compare with current phylogenetic evidence

Search directly relevant phylogenomic and taxonomic literature first, then recent reviews, foundational studies, and recognized taxonomy. Escalate exact species to genus, family, then order when direct evidence is absent, and label the evidence as indirect. For viruses, use current ICTV taxonomy in addition to primary literature.

Record DOI/PMID, year, taxon coverage, data type, inference method/model, topology claim, directness, conflicts, and limitations. Compare the gene tree qualitatively with accepted species relationships; do not treat discordance as automatic pipeline failure.

## Compile the deterministic review bundle

Read [output-contract.md](references/output-contract.md) before creating a request, running the planner, interpreting its review bundle, or assembling the final executed report.

After an authorized host agent has materialized a resolved protein and candidate TSV/FASTA bundle, locate this `SKILL.md`, treat its directory as the skill root, and run:

```text
python3 <skill-root>/scripts/gene_to_tree.py plan \
  --request <request.json> --offline --dry-run --out <new-output-directory>
```

Review `selected_references.tsv`, `rejected_references.tsv`, `reference_set.faa`, `sequence_metadata.tsv`, `itol_roles.txt`, `plan.json`, and `manifest.json`. The helper launches no network request or external executable in plan mode, refuses overwrite, stores commands as argument arrays, and invalidates approval when a decision-bearing input changes.

Inspect optional local executables with:

```text
python3 <skill-root>/scripts/gene_to_tree.py doctor --json
```

## Keep capability claims honest

Use host-provided, authorized database, literature, browser, and shell capabilities for live acquisition and execution. The bundled helper validates a materialized local bundle and compiles deterministic plans and annotations; it does not itself query NCBI, UniProt, Ensembl, OMA, OrthoDB, literature indexes, Open Tree, ICTV, or iTOL, and it does not run MMseqs2, MAFFT, trimAl, FastTree, or IQ-TREE2 during `plan`.

Do not fabricate a lookup result, sequence, TaxID, orthology call, citation, tool version, or database release when a capability is unavailable.

## Bundled resources

- [query-resolution.md](references/query-resolution.md): accession, name, raw-sequence, CDS, fallback, privacy, and viral routing.
- [tool-routing.md](references/tool-routing.md): authoritative databases, search tiers, and executable boundaries.
- [reference-selection.md](references/reference-selection.md): selection, clustering, taxonomic balance, outgroups, and reason codes.
- [alignment-and-tree.md](references/alignment-and-tree.md): MAFFT, trimAl, FastTree, IQ-TREE2, QC, and support semantics.
- [itol-and-literature.md](references/itol-and-literature.md): iTOL files, metadata, evidence search, and gene-tree/species-tree comparison.
- [workflow.md](references/workflow.md): states, gates, failure conditions, and viral branch.
- [output-contract.md](references/output-contract.md): request, artifact, plan, manifest, and final-report contracts.
- `references/request-0.2.schema.json` and `references/plan-0.2.schema.json`: portable JSON schemas.
- `scripts/gene_to_tree.py`: standard-library offline review-bundle compiler and tool doctor.
