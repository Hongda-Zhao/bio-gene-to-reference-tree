# Step 3 — Discover and annotate candidates

## When to read

Read after the objective is fixed. This step builds a broad, auditable candidate pool; it does not choose the final tree tips.

## Required inputs

- resolved query protein and objective;
- ingroup scope, relationship policy, and allowed database/search capabilities;
- remote-access decision and database release/snapshot requirements;
- coverage, domain, and candidate-pool targets.

## Procedure

Use the narrowest suitable evidence route:

1. For `ortholog-tree`, prefer curated one-to-one/one-to-many calls from an appropriate orthology resource.
2. For a raw or unresolved protein, search RefSeq protein first when suitable.
3. Escalate through reviewed Swiss-Prot/UniProtKB, broader UniProtKB or `nr`, then domain/profile searches such as InterPro/Pfam, jackhmmer, HHsearch, or iterative MMseqs2.
4. Use structure-aware search only when sequence/profile evidence is insufficient. A shared fold supports a hypothesis; it does not prove homology or identical function.
5. Retrieve substantially more candidates than the final tree requires and retain every raw response or its checksum.

Keep evidence levels separate: authoritative accession resolution; curated orthology/reviewed annotation; compatible full-length similarity; domain/profile evidence; structure evidence; and name/literature context. Do not average them into one confidence score or select only BLAST top N.

For every candidate, record accession/version, sequence, TaxID/species, lineage/clade, gene/protein names, relation/orthology evidence, reviewed/canonical/fragment status, query and target coverage, percent identity, alignment span, E-value, bit score, sequence length, domain architecture, source database/release, retrieval time, retrieval-query ID, and notes. Use empty fields for unavailable values; never invent them.

Require these TSV fields (real tabs):

```text
accession, taxon_id, species, role, relation,
is_reviewed, is_canonical, is_fragment,
query_coverage, sequence_length, bitscore, evalue,
source_db, source_release, retrieved_at, clade
```

Preserve optional accession version, names, lineage, `analysis_group`, target coverage, identity/alignment length, orthology/domain evidence, cluster fields, outgroup rationale, retrieval-query ID, and notes. Keep biological `role=ingroup|outgroup` separate from display provenance: query → `study`, retained ingroup reference → `expanded`, and outgroup candidate → `outgroup`.

Materialize exactly one FASTA sequence per candidate TSV accession. Preserve stable raw provider identifiers and tool-safe future tip IDs separately. Validate intended species/TaxID pairs with the same exact NCBI taxdump snapshot when that policy is enabled.

## Required outputs

- immutable raw candidate FASTA and metadata TSV;
- search queries/parameters, provider response or checksum, database release, and retrieval time;
- candidate counts by acquisition tier, taxon, relationship, and warning class;
- explicit unresolved annotations or evidence conflicts;
- a handoff that satisfies the TSV schema above and the request/provenance rules in [output-contract.md](../output-contract.md).

## Review gate and stop conditions

Stop if only a generic domain/fold is supported, profile drift is evident, candidate records cannot be mapped one-to-one to sequences, taxonomy is unresolved, permissions forbid required discovery, or the pool cannot represent the declared scope. Do not label similarity hits as orthologs without separate evidence.

## Supporting references

- [Tool and database routing](../tool-routing.md)
- [Taxonomy exact-match policy](../taxonomy-resolution.md)
- [Output and provenance contract](../output-contract.md)

Primary search documentation: NCBI BLAST database descriptions <https://www.ncbi.nlm.nih.gov/books/NBK62345/>, UniProt REST API <https://www.uniprot.org/help/api>, and InterProScan <https://interproscan-docs.readthedocs.io/en/v5/Introduction.html>.
