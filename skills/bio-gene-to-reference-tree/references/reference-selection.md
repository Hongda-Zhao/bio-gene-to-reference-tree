# Reference-selection policy

Select references to answer the declared biological objective, not to maximize similarity. Keep the complete candidate pool and make every decision reproducible.

## Contents

- [Objective-specific policy](#objective-specific-policy)
- [Candidate handoff schema](#candidate-handoff-schema)
- [Deterministic decision order](#deterministic-decision-order)
- [Stable reason codes](#stable-reason-codes)
- [MMseqs2 clustering policy](#mmseqs2-clustering-policy)
- [Coverage, similarity, and domains](#coverage-similarity-and-domains)
- [Taxonomic balance](#taxonomic-balance)
- [Outgroup policy](#outgroup-policy)
- [Mandatory review conditions](#mandatory-review-conditions)

## Objective-specific policy

### `ortholog-tree`

- Prefer curated one-to-one orthology from an appropriate taxonomic resource.
- Exclude labelled paralogs when `allow_paralogs` is false.
- Treat one-to-many, many-to-many, co-ortholog, and conflicting calls as review conditions.
- Use similarity, coverage, length, and domain architecture as quality checks, not proof of orthology.
- Do not call the resulting topology a species tree.

### `homolog-context`

- Retain labelled paralogs when they provide requested family/subfamily context.
- Balance close matches, reviewed representatives, taxonomic breadth, and relevant paralog groups.
- Do not choose only the highest bit score from each species; that can hide duplications or select the wrong copy.
- Do not infer identical function from homology alone.

### `within-species`

- Preserve strain, isolate, allele, and copy identifiers.
- Count usable sequences rather than distinct TaxIDs for minimum-size gates.
- Do not collapse biologically meaningful alleles or copies merely because they are highly similar.
- Require a rooted interpretation only when a defensible external homolog is available.

## Candidate handoff schema

Require the following TSV columns and exactly one matching FASTA record per `accession`:

```text
accession, taxon_id, species, role, relation,
is_reviewed, is_canonical, is_fragment,
query_coverage, sequence_length, bitscore, evalue,
source_db, source_release, retrieved_at, clade
```

Accept and preserve these optional fields:

```text
accession_version, gene_name, protein_name, lineage, analysis_group,
target_coverage, percent_identity, alignment_length,
orthology_source, orthology_evidence, domain_architecture,
cluster_id, cluster_representative, outgroup_rationale,
retrieval_query_id, notes
```

Use real tabs. Represent unavailable optional values as empty. Never invent confidence, taxonomy, release, retrieval metadata, or an outgroup rationale.

When local NCBI taxdump validation is enabled, every `species` value must be the intended character-for-character NCBI `scientific name`, and every `taxon_id` must agree with its single exact match in the supplied snapshot. Do not rewrite candidate rows from a case-insensitive, fuzzy, substring, synonym, common-name, or first-result match. Preserve the independent resolution evidence in `taxonomy_resolution.tsv`; stop before selection if any name is unresolved, ambiguous, missing from `nodes.dmp`, or inconsistent with its supplied TaxID. This check does not derive or validate the free-text `lineage` or `clade` columns; curate those separately and do not label them dump-validated.

Keep `role=ingroup|outgroup` as the biological scope. Derive the iTOL/display `analysis_group` as:

- query/self record → `study`;
- retained ingroup reference → `expanded`;
- outgroup candidate → `outgroup`.

Permit additional `study` records only when explicitly labelled and validated.

## Deterministic decision order

Apply rules in this order and use accession as the final tie-breaker:

1. Validate one unique metadata row and one sequence per accession.
2. Include the query exactly once as `relation=self`.
3. Reject invalid sequences, flagged fragments, insufficient query/target coverage, incompatible length, and incompatible domain architecture.
4. Apply objective-specific relationship rules.
5. If MMseqs2 clustering is triggered, protect study/outgroup sequences and retain the reviewed, canonical, best-supported expanded representative per cluster; preserve every member in metadata.
6. Rank within each TaxID by relationship evidence, reviewed status, canonical status, query coverage, bit score, then accession.
7. Apply `max_per_taxon`.
8. Round-robin across declared clades until `max_references` is reached; do not remove required outgroups.
9. Verify ingroup breadth, total sequence count, outgroup rationale, and unresolved warnings.

In schema 0.2, `max_references` counts retained ingroup and outgroup references but not the query.

## Stable reason codes

| Code | Meaning |
|---|---|
| `FRAGMENT_FLAG` | Metadata labels the record as a fragment. |
| `LOW_QUERY_COVERAGE` | Query coverage is below policy. |
| `LOW_TARGET_COVERAGE` | Target coverage is below policy when configured. |
| `LENGTH_RATIO_OUT_OF_RANGE` | Candidate/query length ratio is outside policy. |
| `DOMAIN_ARCHITECTURE_MISMATCH` | Domain composition does not answer the full-length objective. |
| `RELATION_NOT_ALLOWED` | Relationship conflicts with objective or paralog policy. |
| `MMSEQS_CLUSTER_REDUNDANT` | A better expanded representative occupies the precomputed MMseqs2 cluster. |
| `PER_TAXON_LIMIT` | A better record already fills that TaxID quota. |
| `OUTGROUP_LIMIT` | Better approved outgroup candidates fill the quota. |
| `MAX_REFERENCE_LIMIT` | Candidate passed filters but fell beyond the balanced cap. |

Join multiple independent failures with semicolons in fixed evaluation order. Write every non-selected record exactly once. A query/self record is included once and never appears in the rejection table.

## MMseqs2 clustering policy

Use clustering only when the expanded pool crosses the configured trigger or when the user requests it. Require all of:

```text
--min-seq-id <identity>
-c <coverage>
--cov-mode <mode>
```

Use `--cov-mode 0` as a full-length starting point because coverage is measured against the longer sequence. Consider target-oriented modes only for a documented fragment use case. Do not describe `-c 0.7` as “70% identity.”

Cluster `expanded_candidates.faa`, never the mixed study/outgroup set. Preserve `cluster_id`, representative/member mapping, role, taxon, and accession. Re-run selection after clustering; do not use a pre-clustering selection as the final approved set.

## Coverage, similarity, and domains

- Treat coverage and length bounds as explicit project policy, not universal constants.
- Prefer query and target coverage together. A high-scoring local domain match can be unsuitable for a full-length tree.
- Compare E-values only within the same database snapshot/search; use bit score only after eligibility checks.
- Treat reviewed and canonical flags as ranking evidence, not proof of biological relevance.
- Flag fusions, domain losses, long insertions, low complexity, and extreme length rather than forcing them into the MSA.

## Taxonomic balance

Sample the declared ingroup rather than allowing database-rich species or one dense clade to dominate. Apply validated TaxID quotas, report unsampled major clades, and record every cap. For one-to-one orthologs, one record per species is reasonable only after the relationship is established. For duplication-rich families, retain multiple copies and label them explicitly. Treat NCBI parent/rank data as classification evidence for sampling, not as a statistically inferred species tree.

Do not manufacture breadth using extremely weak or fragmentary hits. Stop when the requested scope cannot be represented credibly.

## Outgroup policy

Treat rooting as a separate biological decision.

- Require an outgroup to be homologous and outside the declared ingroup.
- Prefer a nearby sister lineage over the most distant available hit.
- Keep two or more candidates when feasible and record a prose rationale plus taxonomic evidence.
- Check domain compatibility, extreme branch length, paralogy, and stability across candidates.
- Preserve the unrooted result and create a separate rooted copy only after approval.

Stop when the proposed outgroup is non-homologous, inside the ingroup, an objective-changing paralog, excessively divergent, or unsupported by available taxonomy.

## Mandatory review conditions

Pause when identifier resolution is ambiguous, orthology sources conflict, one-to-many calls affect selection, a candidate is a fragment/fusion/domain-only hit, clustering removes a unique taxon, taxonomic scope is inadequately sampled, no defensible outgroup remains, or fewer than four usable taxa/sequences survive.

The bundled helper evaluates supplied metadata only. It must not infer missing TaxIDs, upgrade a generic ortholog call, or validate an outgroup from sequence distance alone.
