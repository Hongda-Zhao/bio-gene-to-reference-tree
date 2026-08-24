# Step 1 — Resolve the query

## When to read

Read this file at intake, before candidate discovery. Finish with one stable local protein record, versioned when an authoritative record exists, or stop with an explicit unresolved state.

## Required inputs

- the unchanged user input: accession, raw sequence, protein/gene name, or CDS;
- source organism or TaxID when the input is a name or the namespace is ambiguous, and before a raw sequence enters taxon-aware reference planning;
- the intended molecule and any isoform, segment, ORF, or genetic-code context;
- the recorded decision on remote lookups and unpublished-data submission.

## Procedure

Classify the input before looking it up:

| Input | Resolution target | Mandatory stop |
|---|---|---|
| Accession | Authoritative accession.version and protein sequence | Withdrawn/suppressed record, ambiguous namespace, or sequence mismatch |
| Raw protein | Validated local query record; authoritative match when established | Invalid alphabet, extreme low complexity, likely nucleotide, or unresolved molecule type |
| Protein/gene name | Stable gene/protein ID and versioned protein | Missing taxon, several distinct genes, or unresolved isoform |
| CDS | Valid protein translation plus reversible CDS mapping | Internal stop, invalid frame, uncertain genetic code, or probable non-coding input |

For an accession, confirm the namespace with its provider. Retrieve record status, accession.version, sequence, organism, TaxID, lineage, gene/protein names, reviewed/canonical/isoform status, release or snapshot, and retrieval time. Preserve obsolete or replaced identifiers and compare any user-supplied sequence with the authoritative sequence. Never silently replace a requested isoform.

For a name, require organism or TaxID. Search stable identifiers using symbols, synonyms, locus tags, and family terminology. Use literature to disambiguate identity, never as the sequence source. Show biologically plausible transcript or isoform alternatives and require review if they alter length, domains, localization, or interpretation. Never default a symbol to human or another model organism.

For a raw protein, validate alphabet, length, ambiguous residues, stop characters, low complexity, repeats, and likely signal/transmembrane regions, then assign a stable local identifier bound to its sequence hash. Candidate/database searching belongs to Step 3. A nearest hit proposes family membership; it does not prove orthology.

An unknown-source raw protein may be validated and searched locally, but its nearest homolog cannot establish the query's source organism. Stop before the current taxon-aware planner/reference-approval path until the organism is supplied or independently resolved from sample provenance; never copy a hit's TaxID onto the query.

For CDS, verify the genetic code, frame, length modulo three, and internal stops. Align translated proteins and back-translate with a codon-aware method if nucleotide analysis is required. Route genomic, non-coding, or RNA inputs elsewhere.

For viral inputs, record host, genome type, segment, ORF coordinates, polyprotein processing, and mature protein. Do not combine different segments or mature proteins.

When assigning a TaxID from an organism string, apply the exact-match policy in [taxonomy-resolution.md](../taxonomy-resolution.md); do not normalize, fuzzy-match, accept an alias, or select the first result.

## Required outputs

- original input and complete transformation chain;
- a local resolved protein FASTA with a stable tip-safe identifier and a versioned accession when one was established;
- authoritative metadata and sequence SHA-256;
- provider request/release/snapshot, retrieval time, and raw-response checksum when available;
- explicit unresolved alternatives and limitations;
- taxonomy-resolution evidence when a name was mapped to a TaxID.

## Review gate and stop conditions

Pass Gate 1 only when protein identity, organism, sequence, and provenance are resolved. Stop on ambiguous identity, unresolved isoforms, unknown source taxon for taxon-aware planning, conflicting TaxIDs, invalid translation, generic fold/domain-only evidence, or missing permission for a required remote operation. Never log credentials, cookies, unpublished sequences, or private paths.

## Supporting references

- [Workflow states and decision gates](../workflow.md)
- [Taxonomy exact-match policy](../taxonomy-resolution.md)
- [Tool and privacy boundaries](../tool-routing.md)
- [Output and provenance contract](../output-contract.md)

Primary provider documentation: NCBI Datasets gene metadata <https://www.ncbi.nlm.nih.gov/datasets/docs/how-tos/genes/get-gene-metadata/>, NCBI Taxonomy dumps <https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/>, and UniProt entry retrieval <https://www.uniprot.org/help/api_retrieve_entries>.
