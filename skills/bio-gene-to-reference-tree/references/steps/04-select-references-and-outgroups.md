# Step 4 — Select references and outgroups

## When to read

Read after the broad candidate pool is materialized. Finish with a reviewable selection proposal or route the expanded pool through Step 5.

## Required inputs

- immutable molecule-matched candidate FASTA/TSV, resolved query, and any required translation/RNA-provenance files;
- objective, taxonomic scope, and relationship/paralog policy;
- comparable-region, coverage, length, domain/feature, per-taxon, and maximum-reference policies;
- orthology evidence and proposed outgroup taxonomy.

## Procedure

Apply decisions in a deterministic order:

1. Require one metadata row and one sequence per accession; include the query exactly once as `relation=self`.
2. Reject molecule/region/orientation mismatches, invalid sequences, flagged fragments, inadequate query/target coverage, incompatible length, incompatible domain/feature architecture, and failed RNA/CDS QC.
3. Apply the objective-specific relationship policy. Treat one-to-many, many-to-many, co-ortholog, and conflicting calls as review conditions.
4. Rank eligible records within each TaxID by relationship evidence, reviewed status, canonical status, coverage, bit score, then accession.
5. Apply the declared per-taxon cap, then sample round-robin across declared clades until the reference cap is reached.
6. Report unsampled major clades and all caps; do not manufacture breadth from weak or fragmentary hits.

For `within-species`, preserve meaningful strains, alleles, and copies even when highly similar. In request schemas 0.2 and 0.3, `max_references` counts every retained ingroup, outgroup, and additional `study` reference but not the query. Protected study records consume capacity before expanded ingroup sampling. If the protected studies and required outgroups alone exceed the cap, retain them for audit but block with `REFERENCE_CAP_EXCEEDED_BY_PROTECTED_STUDIES` instead of silently exceeding the approved design.

Use stable rejection codes: `FRAGMENT_FLAG`, `LOW_QUERY_COVERAGE`, `TARGET_COVERAGE_MISSING`, `LOW_TARGET_COVERAGE`, `LENGTH_RATIO_OUT_OF_RANGE`, `DOMAIN_ARCHITECTURE_MISMATCH`, `RELATION_NOT_ALLOWED`, `MMSEQS_CLUSTER_REDUNDANT`, `PER_TAXON_LIMIT`, `OUTGROUP_LIMIT`, and `MAX_REFERENCE_LIMIT`. A positive target-coverage threshold rejects a reference whose target coverage is unknown; protected study records remain visible review conditions. Join independent reasons in fixed evaluation order and write each rejected record once.

Treat rooting as a separate biological decision. An outgroup must be homologous, outside the ingroup, and preferably from a nearby sister lineage. Retain two or more candidates when feasible, with taxonomy and prose rationales. Never choose the weakest, most distant, or longest-branched hit automatically. Preserve an unrooted interpretation when no defensible candidate exists.

If the expanded pool crosses the approved clustering trigger, do not finalize the reference set; hand it to Step 5. Records with `analysis_group=study` or `analysis_group=outgroup` remain protected.

## Required outputs

- `selected_references.tsv`, `rejected_references.tsv`, and molecule-matched `reference_set.faa` or `reference_set.fna` proposal;
- exact CDS↔translation selection mapping or RNA source↔derivative mapping when applicable;
- counts before/after each rule, selected taxa, unsampled clades, warnings, and stable reason codes;
- candidate outgroups with taxonomic evidence and explicit rationales;
- planned clustering/alignment/trimming/tree argv arrays;
- current plan hash and approval status.

## Review gate and stop conditions

Pause for reference/outgroup approval. Stop on ambiguous identifiers, mixed molecules/non-comparable regions, failed translation/normalization mapping, conflicting orthology that changes membership, fragment/fusion/domain-only matches, inadequate breadth, an outgroup inside the ingroup, a non-homologous or objective-changing paralog, or fewer than four usable taxa/sequences. Never recover by changing analysis space.

## Supporting references

- [Workflow Gate 2](../workflow.md#gate-2-reference-and-outgroup-approval)
- [Taxonomy exact-match policy](../taxonomy-resolution.md)
- [Molecule, region, and CDS/RNA gates](../sequence-type-routing.md)
- [Selection artifacts and reason-code contract](../output-contract.md)
