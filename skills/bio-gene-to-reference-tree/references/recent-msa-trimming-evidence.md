# Recent MSA and trimming evidence

Use this reference when recent literature is being used to choose, challenge, or justify an alignment and post-alignment filtering strategy. The machine-readable companion is [recent-msa-trimming-evidence.tsv](recent-msa-trimming-evidence.tsv).

Literature does not choose molecule type. First apply [sequence-type-routing.md](sequence-type-routing.md): legacy requests remain protein-only, while request 0.3 explicitly separates protein, noncoding DNA/RNA, and clean CDS with an approved nucleotide-site or codon analysis. A paper using another analysis space is not permission to fall back into it.

## Scope and provenance

This is a curated snapshot of 59 analysis workflows from 45 primary evolutionary or methods papers published from 2023-08-24 through 2026-08-24 and checked on 2026-08-24. Searches routed through Europe PMC, PubMed, PubMed Central, and official publisher pages using combinations of phylogenomics, gene-family evolution, MAFFT, multiple sequence alignment, trimAl, ClipKIT, BMGE, Gblocks, and alignment filtering. Methods were verified against accessible full text rather than inferred from abstracts.

The catalog is deliberately diverse rather than systematic or exhaustive. It contains animal, plant, fungal, protist, bacterial, archaeal, and viral examples; single genes, multigene families, UCEs, rDNA, mitogenomes, BUSCO matrices, phylogenomic marker sets, one reference-mapped SNP pseudoalignment boundary case, and recent filtering benchmarks. It includes explicit trimming, manual curation, untrimmed controls, and incompletely reported workflows. Each TSV row is one analysis workflow, so one citation can have several rows.

All method descriptions are paraphrased factual metadata. Follow each `source_url` and check the article's license before reusing figures, tables, supplements, or substantial text.

## Conservation is scale-dependent

“Conserved gene” and “non-conserved gene” are incomplete labels unless the comparison scale and biological unit are stated. A family may have an ancient conserved catalytic core, rapidly evolving loops, lineage-specific duplications, and patchy gene loss at the same time. Conservation also does not establish orthology, single-copy status, function, or a suitable outgroup.

The TSV therefore uses three routing fields:

- `conservation_class`: one controlled, coarse workflow-level category;
- `conservation_scope`: the taxonomic scale at which that category is useful;
- `conservation_basis`: the marker design or reported distribution, copy-number, orthology, or family evidence supporting the assignment.

Use the classes as a first-pass map, not as measured evolutionary rates:

| Class | Approximate use | Typical examples |
|---|---|---|
| `universal-or-deep-core` | Ancient, high-retention core or mostly single-copy marker panels; the particular application may still sample a shallow clade | ribosomal proteins, translation or RNA-polymerase markers, validated prokaryotic core sets, BUSCO matrices, deep phylogenomic panels |
| `clade-conserved-marker` | Markers designed or retained within a named lineage rather than universally | 16S/18S/28S rDNA, lineage single-copy orthologues, UCEs, Angiosperms353, COI/cytb or rbcL/matK at a suitable scale |
| `broadly-conserved-gene-family` | A homologous family spans broad taxa but is not itself a universal core panel | central-metabolism enzymes, hydrogenases, aminoacyl-tRNA synthetases, cross-plant DFR homologues, some polymerase families |
| `variable-multigene-family` | A homologous family has important duplication, loss, domain, or insertion variation | odorant/gustatory receptors, GT1, NPC2, other expanded paralogous families |
| `lineage-specific-or-rapidly-evolving` | Distribution is patchy or the focal family is recent, transferred, effector-like, or unusually divergent | fungal accessory/cluster genes, candidate LGT families, lineage-restricted paralogues, highly divergent viral proteins |
| `mixed-conservation-panel` | One biological analysis deliberately combines loci or families with different retention or evolutionary properties | heterogeneous marker supermatrices, phylomes, and broad multi-family collections |
| `not-applicable` | The analysis unit is not one biological gene or marker panel, so one conservation label would misdescribe it | method benchmarks over unrelated datasets and reference-mapped genome-wide SNP pseudoalignments |

These categories are intentionally more informative than a binary conserved/non-conserved flag. Reassign them when the focal sampling scale changes, and record uncertainty when distribution or copy-number evidence is incomplete. Do not derive an MSA mode or trim threshold from the class alone; sequence architecture, observed alignment quality, missingness, domains, and the analysis objective remain decisive.

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

- Classify the alignment architecture before selecting a precedent. A conventional protein MSA, direct noncoding nucleotide MSA, protein-guided codon alignment, target-capture or concatenated locus matrix, and reference-mapped SNP pseudoalignment are not interchangeable. The pseudoalignment row is a routing boundary, not permission to apply a genome pipeline inside this single-gene-tree Skill.
- Very large curated protein families can use a stable template plus MAFFT `--add`; preserve the template and inspect every addition.
- For clean coding-gene analyses, verify frame/code/translation, align proteins with explicit MAFFT `--amino`, and use trimAl `-backtrans <matched-CDS.fna>` so filtering decisions are projected as complete codons. This deterministic trimAl route is limited to NCBI codes 1/11; other codes need a reviewed MACSE/PAL2NAL/precomputed handoff. Preserve both protein and codon outputs and re-check translation equality. The cataloged NPC2 workflow reports direct trimAl filtering after PAL2NAL without documenting triplet preservation and is a caution, not a default.
- Conserved noncoding nucleotide loci and rDNA use explicit MAFFT `--nuc` and may still need Gblocks, explicit manual endpoint cleanup, or removal of genuinely unalignable ITS regions. Require comparable locus boundaries and orientation; RNA needs declared `rna-u`/`dna-t` source encoding and source/analysis-copy provenance. Requested Q-INS-i uses the dedicated `mafft-qinsi` executable.
- Frameshifted, disrupted, internal-stop-containing, pseudogene, or uncertain-frame CDS requires reviewed MACSE handling. A documented terminal stop in an otherwise clean, translation-validated CDS is not an internal stop. Disrupted CDS must not be silently aligned as direct nucleotide or promoted into the clean codon route; only a reviewed/exported, typed precomputed codon alignment may be re-profiled for inference.
- Deep protein phylogenomics commonly uses L-INS-i or E-INS-i with trimAl `gappyout`, ClipKIT, BMGE, or several competing filters; mixture-model and outgroup sensitivity remain separate requirements.
- Large phylogenomic matrices also need locus completeness, composition, homology, saturation, and long-branch checks. Column trimming cannot substitute for those checks.
- Viral workflows range from explicit no-trim placement of short RdRp fragments to E-INS-i plus gap filtering for divergent giant-virus proteins. Gene or segment definition and recombination screening must decide the route.
- When topology may depend on alignment treatment, compare aligners or trimming profiles and include a defensible untrimmed control. A recent plastid study compared seven filters with an untrimmed matrix; a fungal study compared two aligner-trimmer pairs; and the 2026 AliFilter benchmark compared six filters with unfiltered alignments.
- The AliFilter benchmark recovered a dataset-specific mix of identical and different topologies among filters. The CLOAK benchmark found that gentle filtering of focal single-gene alignments and stricter filtering of large substitution-model training sets served different purposes. Together these support task-specific sensitivity testing, not a universal ranking of AliFilter, CLOAK, ClipKIT, trimAl, BMGE, or Gblocks.
- Recent high-impact papers still omit point versions, modes, thresholds, or support replicate counts. Journal prestige and recency do not convert `partial` into `exact`; never reconstruct missing commands from cited software papers or current defaults.
- For a single-gene tree, filtering should be no stronger than the biological and alignment evidence justify. Retain the raw alignment and include an untrimmed or gentler-filter control when a conclusion changes with column removal.

These are empirical precedents, not defaults. Select a row only after exactly matching molecule type and analysis architecture, then comparing conservation class and scope, dataset size, homologous region, taxonomic depth, missing-data pattern, and analysis objective. Do not choose a method because the focal organism shares only the row's broad taxon label.

## Agent use

1. State the exact request molecule/analysis space and the focal dataset's provisional `conservation_class`, `conservation_scope`, and `conservation_basis` in the hash-bound assessment. Filter the TSV first by compatible `molecule_type`/architecture, then by `gene_or_markers`, class/scope, region, scale, and `taxon_scope`; use `broad_group` only as a secondary cue.
2. Prefer rows with `exact` reporting, but retain directly relevant `partial` rows as limitations rather than filling missing parameters.
3. Retrieve the source again before treating a row as decisive. Record the retrieval date and whether the current article or correction differs from this snapshot.
4. Propose at least one biologically justified primary alignment and one reasonable sensitivity analysis when the evidence is method-sensitive.
5. Preserve raw and filtered alignments, exact argv arrays, versions, retained-site statistics, and topology comparisons.
6. Cite the paper that informed the decision and explain why its data architecture is comparable to the current dataset.

If `molecule_type` contains `pseudoalignment`, stop this Skill's alignment route. A reference-mapped SNP pseudoalignment needs its own reference, callable-core, repeat/recombination-mask, invariant-site, and ascertainment records; it must not be passed through MAFFT or trimAl as though it were a protein-family MSA.

When the catalog lacks a comparable route, follow the primary [MAFFT](https://mafft.cbrc.jp/alignment/software/manual/manual.html), [trimAl](https://trimal.readthedocs.io/en/latest/usage.html), [IQ-TREE](https://iqtree.github.io/doc/Command-Reference), and [FastTree](https://morgannprice.github.io/fasttree/) documentation. Do not import a protein threshold/model into DNA or codon analysis, or use FastTree nucleotide mode as a codon model.

Do not let the catalog bypass the reference-approval or alignment/trimming-approval gates. If no row is genuinely comparable, follow primary tool documentation and say that recent empirical precedent was not identified.

## Refresh protocol

For a later update, define a new fixed three-year window and retrieval date, search primary studies, verify full Methods text, add one row per distinct workflow, and preserve stable DOI/PMID/PMCID identifiers. Assign conservation only with an explicit comparison scope and evidence basis; revise it when sampling or marker architecture changes. Keep `not-reported` separate from `explicit-none`, retain exact prose thresholds when a CLI mapping is uncertain, and run the repository metadata tests before release.
