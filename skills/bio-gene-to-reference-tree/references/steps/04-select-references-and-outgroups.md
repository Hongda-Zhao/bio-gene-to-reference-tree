# Step 4 — Select references and outgroups

## When to read

Read after the broad candidate pool is materialized. Finish with a reviewable selection proposal or route the expanded pool through Step 5.

## Required inputs

- immutable candidate FASTA/TSV and resolved query;
- objective, taxonomic scope, and relationship/paralog policy;
- coverage, length, domain, per-taxon, and maximum-reference policies;
- orthology evidence and proposed outgroup taxonomy.

## Procedure

Apply decisions in a deterministic order:

1. Require one metadata row and one sequence per accession; include the query exactly once as `relation=self`.
2. Reject invalid sequences, flagged fragments, inadequate query/target coverage, incompatible length, and incompatible domain architecture.
3. Apply the objective-specific relationship policy. Treat one-to-many, many-to-many, co-ortholog, and conflicting calls as review conditions.
4. Rank eligible records within each TaxID by relationship evidence, reviewed status, canonical status, coverage, bit score, then accession.
5. Apply the declared per-taxon cap, then sample round-robin across declared clades until the reference cap is reached.
6. Report unsampled major clades and all caps; do not manufacture breadth from weak or fragmentary hits.

For `within-species`, preserve meaningful strains, alleles, and copies even when highly similar. In request schema 0.2, `max_references` counts retained ingroup and outgroup references but not the query.

Use stable rejection codes: `FRAGMENT_FLAG`, `LOW_QUERY_COVERAGE`, `LOW_TARGET_COVERAGE`, `LENGTH_RATIO_OUT_OF_RANGE`, `DOMAIN_ARCHITECTURE_MISMATCH`, `RELATION_NOT_ALLOWED`, `MMSEQS_CLUSTER_REDUNDANT`, `PER_TAXON_LIMIT`, `OUTGROUP_LIMIT`, and `MAX_REFERENCE_LIMIT`. Join independent reasons in fixed evaluation order and write each rejected record once.

Treat rooting as a separate biological decision. An outgroup must be homologous, outside the ingroup, and preferably from a nearby sister lineage. Retain two or more candidates when feasible, with taxonomy and prose rationales. Never choose the weakest, most distant, or longest-branched hit automatically. Preserve an unrooted interpretation when no defensible candidate exists.

If the expanded pool crosses the approved clustering trigger, do not finalize the reference set; hand it to Step 5. Records with `analysis_group=study` or `analysis_group=outgroup` remain protected.

## Required outputs

- `selected_references.tsv`, `rejected_references.tsv`, and `reference_set.faa` proposal;
- counts before/after each rule, selected taxa, unsampled clades, warnings, and stable reason codes;
- candidate outgroups with taxonomic evidence and explicit rationales;
- planned clustering/alignment/trimming/tree argv arrays;
- current plan hash and approval status.

## Review gate and stop conditions

Pause for reference/outgroup approval. Stop on ambiguous identifiers, conflicting orthology that changes membership, fragment/fusion/domain-only matches, inadequate breadth, an outgroup inside the ingroup, a non-homologous or objective-changing paralog, or fewer than four usable taxa/sequences.

## Supporting references

- [Workflow Gate 2](../workflow.md#gate-2-reference-and-outgroup-approval)
- [Taxonomy exact-match policy](../taxonomy-resolution.md)
- [Selection artifacts and reason-code contract](../output-contract.md)
