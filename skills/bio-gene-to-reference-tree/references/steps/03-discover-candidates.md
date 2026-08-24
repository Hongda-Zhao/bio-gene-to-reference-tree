# Step 3 — Discover and annotate candidates

## When to read

Read after the objective is fixed. This step builds a broad, auditable candidate pool; it does not choose the final tree tips.

## Required inputs

- resolved molecule-matched query, comparable-region definition, and objective;
- ingroup scope, relationship policy, and allowed database/search capabilities;
- remote-access decision and actual compatible search-database/index requirements;
- coverage, domain, and candidate-pool targets.

## Procedure

Use the narrowest suitable evidence route:

1. For `ortholog-tree`, prefer curated one-to-one/one-to-many calls from an appropriate molecule/feature-aware resource.
2. For protein, search RefSeq protein or another actual protein database with `blastp` when suitable; protein-only domain/profile escalation remains available.
3. For noncoding DNA/RNA, search a verified nucleotide database with `blastn`, then retrieve the same homologous locus/feature and orientation.
4. For clean CDS, prefer curated CDS plus matched translation; reviewed nucleotide or translated searches may discover records, but every retained record still requires comparable CDS boundaries and translation QC.
5. Use structure-aware search only for proteins when sequence/profile evidence is insufficient. A shared fold supports a hypothesis; it does not prove homology or identical function.
6. Retrieve substantially more candidates than the final tree requires and retain every raw response or its checksum.

Keep evidence levels separate: authoritative accession/feature resolution; curated orthology/reviewed annotation; compatible full-region similarity; protein domain/profile evidence; structure evidence; and name/literature context. Do not average them into one confidence score or select only BLAST top N.

For every candidate, record accession/version, exact molecule type, sequence/feature region, coordinates/strand/orientation, TaxID/species, lineage/clade, gene/protein/feature names, relation/orthology evidence, reviewed/canonical/fragment status, query and target coverage, percent identity, alignment span, E-value, bit score, sequence length, protein domain architecture where applicable, source database/release, retrieval time, retrieval-query ID, actual search database, and notes. For CDS also record genetic code, frame source, completeness, stop/frameshift counts, and translation equality. Use empty fields for unavailable values; never invent them.

Require these TSV fields (real tabs):

```text
accession, taxon_id, species, role, relation,
is_reviewed, is_canonical, is_fragment,
query_coverage, sequence_length, bitscore, evalue,
source_db, source_release, retrieved_at, clade
```

Preserve optional accession version, names, lineage, `coordinates`, `analysis_group`, target coverage, identity/alignment length, orthology/domain evidence, cluster fields, outgroup rationale, retrieval-query ID, and notes. Unsupported columns are rejected instead of silently discarded; use `notes` or a separate acquisition record for provider-specific fields. Keep biological `role=ingroup|outgroup` separate from display provenance: query → `study`, retained ingroup reference → `expanded`, and outgroup candidate → `outgroup`.

For request 0.3, `molecule_type`, `sequence_region`, `strand`, `actual_search_database`, `search_program`, and `search_database_molecule` are mandatory decision evidence even when the legacy table core leaves them optional. The supported program/database-molecule pairs are defined in [sequence-type-routing.md](../sequence-type-routing.md); use `search_program=not-searched`, `search_database_molecule=not-applicable`, and `actual_search_database=local-bundle:not-searched` only for a record that was not acquired through a search. Coding rows additionally require `genetic_code_id`, `reading_frame_source`, `complete_cds`, `internal_stop_count`, `frameshift_count`, and `translation_matches`.

Materialize exactly one molecule-matched FASTA sequence per candidate TSV accession. For coding nucleotide, materialize an exact-ID translation FASTA; for RNA, preserve the declared source encoding and separately named DNA-alphabet analysis copy. Preserve stable raw provider identifiers and tool-safe future tip IDs separately. Validate intended species/TaxID pairs with the same exact NCBI taxdump snapshot when that policy is enabled.

The helper validates that the three required search fields are present and form a supported molecule-compatible pair; it does not execute the search or verify the database label. Record query molecule and any actual database release, physical location, index format, or build/checksum evidence in the separate host-side acquisition record. A visible BLAST/MMseqs2 executable or a prose label such as “RefSeq” is not proof that the compatible database was queried.

## Required outputs

- immutable molecule-matched raw candidate FASTA, metadata TSV, and required translation/RNA provenance;
- search queries/parameters, provider response or checksum, database release, and retrieval time;
- candidate counts by acquisition tier, taxon, relationship, and warning class;
- explicit unresolved annotations or evidence conflicts;
- a handoff that satisfies the TSV schema above and the request/provenance rules in [output-contract.md](../output-contract.md).

## Review gate and stop conditions

Stop if molecule types/regions are mixed, only a generic domain/fold is supported, profile drift is evident, candidate records cannot be mapped one-to-one to sequences/translations, the actual compatible search database is unknown, CDS or RNA provenance fails, taxonomy is unresolved, permissions forbid required discovery, or the pool cannot represent the declared scope. Do not label similarity hits as orthologs without separate evidence or change analysis space as a fallback.

## Supporting references

- [Tool and database routing](../tool-routing.md)
- [Sequence-type and search compatibility](../sequence-type-routing.md)
- [Taxonomy exact-match policy](../taxonomy-resolution.md)
- [Output and provenance contract](../output-contract.md)

Primary search documentation: NCBI BLAST database descriptions <https://www.ncbi.nlm.nih.gov/books/NBK62345/>, UniProt REST API <https://www.uniprot.org/help/api>, and InterProScan <https://interproscan-docs.readthedocs.io/en/v5/Introduction.html>.
