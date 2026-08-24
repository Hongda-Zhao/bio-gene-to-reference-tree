# Recent MSA and trimming evidence

Use this reference when recent literature is being used to choose, challenge, or justify an alignment and post-alignment filtering strategy. The machine-readable companion is [recent-msa-trimming-evidence.tsv](recent-msa-trimming-evidence.tsv).

## Scope and provenance

This is a curated snapshot of 59 analysis workflows from 45 primary evolutionary or methods papers published from 2023-08-24 through 2026-08-24 and checked on 2026-08-24. Searches routed through Europe PMC, PubMed, PubMed Central, and official publisher pages using combinations of phylogenomics, gene-family evolution, MAFFT, multiple sequence alignment, trimAl, ClipKIT, BMGE, Gblocks, and alignment filtering. Methods were verified against accessible full text rather than inferred from abstracts.

The catalog is deliberately diverse rather than systematic or exhaustive. It contains animal, plant, fungal, protist, bacterial, archaeal, and viral examples; single genes, multigene families, UCEs, rDNA, mitogenomes, BUSCO matrices, phylogenomic marker sets, one reference-mapped SNP pseudoalignment boundary case, and a recent filtering benchmark. It includes explicit trimming, manual curation, untrimmed controls, and incompletely reported workflows. Each TSV row is one analysis workflow, so one citation can have several rows.

All method descriptions are paraphrased factual metadata. Follow each `source_url` and check the article's license before reusing figures, tables, supplements, or substantial text.

## Reporting states

Interpret `msa_reporting_status` and `trimming_status` literally. The former reports completeness; the latter also preserves the two special actions that must not be collapsed into missing data:

| State | Meaning |
|---|---|
| `exact` | Tool, version, and the material mode or threshold were reported for that stage. |
| `partial` | The tool was named, but at least one version, mode, threshold, or dataset-specific detail was missing. |
| `explicit-none` | The authors explicitly analyzed an untrimmed alignment or stated that no trimming was performed. |
| `not-reported` | No trimming or masking step was found in the reported workflow. This is not evidence that none occurred. |

Manual masking or terminal cleanup is recorded as the trimming method and receives `exact` or `partial` according to how reproducibly it was described. Never convert prose such as “sites with gaps in more than 75% of sequences were removed” into a trimAl `-gt` value without independently checking the program's semantics and the original command.

`exact` measures reporting completeness, not whether the authors interpreted a flag correctly. Before execution, verify every decision-bearing flag against the documentation for the reported software version. If article prose conflicts with the reported argv or official semantics, preserve both as a documented conflict, do not silently “repair” the command, and require explicit approval for either exact reproduction or a corrected implementation.

## What the evidence supports

- Classify the alignment architecture before selecting a precedent. A conventional protein MSA, protein-guided codon alignment, target-capture or concatenated locus matrix, and reference-mapped SNP pseudoalignment are not interchangeable. The pseudoalignment row is a routing boundary, not permission to apply a genome pipeline inside this protein-gene-tree Skill.
- Very large curated protein families can use a stable template plus MAFFT `--add`; preserve the template and inspect every addition.
- For coding-gene analyses, align proteins and back-translate. If filtering is needed, trim the protein MSA before back-translation or prove that nucleotide-column removal preserves complete triplets; the cataloged NPC2 workflow reports direct trimAl filtering after PAL2NAL without documenting triplet preservation and is a caution, not a default.
- Conserved nucleotide loci and rDNA may still need Gblocks, explicit manual endpoint cleanup, or removal of genuinely unalignable ITS regions.
- Deep protein phylogenomics commonly uses L-INS-i or E-INS-i with trimAl `gappyout`, ClipKIT, BMGE, or several competing filters; mixture-model and outgroup sensitivity remain separate requirements.
- Large phylogenomic matrices also need locus completeness, composition, homology, saturation, and long-branch checks. Column trimming cannot substitute for those checks.
- Viral workflows range from explicit no-trim placement of short RdRp fragments to E-INS-i plus gap filtering for divergent giant-virus proteins. Gene or segment definition and recombination screening must decide the route.
- When topology may depend on alignment treatment, compare aligners or trimming profiles and include a defensible untrimmed control. A recent plastid study compared seven filters with an untrimmed matrix; a fungal study compared two aligner-trimmer pairs; and the 2026 AliFilter benchmark compared six filters with unfiltered alignments.
- The AliFilter benchmark recovered a dataset-specific mix of identical and different topologies among filters. The CLOAK benchmark found that gentle filtering of focal single-gene alignments and stricter filtering of large substitution-model training sets served different purposes. Together these support task-specific sensitivity testing, not a universal ranking of AliFilter, CLOAK, ClipKIT, trimAl, BMGE, or Gblocks.
- Recent high-impact papers still omit point versions, modes, thresholds, or support replicate counts. Journal prestige and recency do not convert `partial` into `exact`; never reconstruct missing commands from cited software papers or current defaults.
- For a single-gene tree, filtering should be no stronger than the biological and alignment evidence justify. Retain the raw alignment and include an untrimmed or gentler-filter control when a conclusion changes with column removal.

These are empirical precedents, not defaults. Select a row only after matching molecule type, sequence architecture, dataset size, taxonomic depth, missing-data pattern, and analysis objective. Do not choose a method because the focal organism shares only the row's broad taxon label.

## Agent use

1. Filter the TSV by `molecule_type`, `gene_or_markers`, `dataset_scale`, and `taxon_scope`; use `broad_group` only as a secondary cue.
2. Prefer rows with `exact` reporting, but retain directly relevant `partial` rows as limitations rather than filling missing parameters.
3. Retrieve the source again before treating a row as decisive. Record the retrieval date and whether the current article or correction differs from this snapshot.
4. Propose at least one biologically justified primary alignment and one reasonable sensitivity analysis when the evidence is method-sensitive.
5. Preserve raw and filtered alignments, exact argv arrays, versions, retained-site statistics, and topology comparisons.
6. Cite the paper that informed the decision and explain why its data architecture is comparable to the current dataset.

If `molecule_type` contains `pseudoalignment`, stop this Skill's alignment route. A reference-mapped SNP pseudoalignment needs its own reference, callable-core, repeat/recombination-mask, invariant-site, and ascertainment records; it must not be passed through MAFFT or trimAl as though it were a protein-family MSA.

Do not let the catalog bypass the reference-approval or alignment/trimming-approval gates. If no row is genuinely comparable, follow primary tool documentation and say that recent empirical precedent was not identified.

## Refresh protocol

For a later update, define a new fixed three-year window and retrieval date, search primary studies, verify full Methods text, add one row per distinct workflow, and preserve stable DOI/PMID/PMCID identifiers. Keep `not-reported` separate from `explicit-none`, retain exact prose thresholds when a CLI mapping is uncertain, and run the repository metadata tests before release.
