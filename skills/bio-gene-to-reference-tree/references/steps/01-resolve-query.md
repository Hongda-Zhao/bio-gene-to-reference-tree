# Step 1 — Resolve the query

## When to read

Read this file at intake, before candidate discovery. Finish with one stable local record matching the explicitly declared molecule and homologous region, versioned when an authoritative record exists, or stop with an explicit unresolved state.

## Required inputs

- the unchanged user input: accession, raw sequence, protein/gene/feature name, or CDS;
- source organism or TaxID when the input is a name or the namespace is ambiguous, and before a raw sequence enters taxon-aware reference planning;
- for request 0.3, exact `protein`, `noncoding-dna`, `noncoding-rna`, `coding-dna`, or `coding-rna`; requests 0.1/0.2 remain protein-only;
- any isoform, feature/locus, comparable region, coordinates, strand, segment, ORF, or genetic-code context;
- the recorded decision on remote lookups and unpublished-data submission.

## Procedure

Classify the input before looking it up:

| Input | Resolution target | Mandatory stop |
|---|---|---|
| Accession | Authoritative accession.version and molecule-matched sequence/feature | Withdrawn/suppressed record, ambiguous namespace/molecule, or sequence mismatch |
| Raw protein | Validated local query record; authoritative match when established | Protein-alphabet failure, extreme low complexity, or declaration conflict |
| Protein/gene/feature name | Stable ID and versioned molecule-matched record | Missing taxon, several distinct records, or unresolved isoform/feature |
| Raw noncoding DNA/RNA | Validated local locus/feature with comparable-region definition | Molecule inferred only from letters, unknown orientation/boundaries, or mixed features |
| Clean CDS | Versioned CDS, exact translation, and reversible mapping | Internal stop, frameshift, invalid/uncertain frame or genetic code, translation mismatch |

For an accession, confirm the namespace and molecule with its provider. Retrieve record status, accession.version, sequence, organism, TaxID, lineage, gene/protein/feature names, coordinates/strand, reviewed/canonical/isoform status, release or snapshot, and retrieval time. Preserve obsolete or replaced identifiers and compare any user-supplied sequence with the authoritative sequence. Never silently replace a requested isoform or feature.

For a name, require organism or TaxID. Search stable identifiers using symbols, synonyms, locus tags, feature names, and family terminology. Use literature to disambiguate identity, never as the sequence source. Show biologically plausible transcript, isoform, locus, or feature alternatives and require review if they alter length, boundaries, domains, localization, or interpretation. Never default a symbol to human or another model organism.

For a raw sequence, validate against the declared molecule rather than guessing from its alphabet. Check length, ambiguity, gap/stop symbols, complexity, orientation, and comparable-region boundaries; for proteins also inspect repeats and likely signal/transmembrane regions. Assign a stable local identifier bound to its sequence hash. Candidate/database searching belongs to Step 3. A nearest hit proposes family membership; it does not prove molecule type, source organism, or orthology.

An unknown-source raw sequence may be validated and searched locally, but its nearest homolog cannot establish the query's source organism. Stop before the current taxon-aware planner/reference-approval path until the organism is supplied or independently resolved from sample provenance; never copy a hit's TaxID onto the query.

For noncoding DNA/RNA, require the same homologous feature/region and orientation across candidates. Declare RNA `source_encoding` as `rna-u` or `dna-t`, preserve it unchanged, and create a separately named DNA-alphabet analysis copy with source/copy hashes. Convert U→T only for `rna-u`; keep the `dna-t` copy byte-equivalent.

For CDS, verify strand, frame source, completeness, genetic-code ID, length modulo three, internal stops/frameshifts, exact CDS↔translation IDs, and translation equality. Clean code-1/11 CDS may proceed by protein-guided alignment and deterministic trimAl `-backtrans`; a documented terminal stop is permitted only after that exact translation gate. Another genetic code requires a reviewed MACSE/PAL2NAL/precomputed handoff; disrupted, internal-stop-containing, pseudogene, or uncertain CDS emits `MACSE_ROUTE_REQUIRED`. MACSE output must be reviewed/exported under an explicit policy and reintroduced only as a typed precomputed codon alignment. Never silently treat either case as direct noncoding nucleotide.

For viral inputs, record host, genome type, segment, ORF coordinates, polyprotein processing, and mature protein. Do not combine different segments or mature proteins.

When assigning a TaxID from an organism string, apply the exact-match policy in [taxonomy-resolution.md](../taxonomy-resolution.md); do not normalize, fuzzy-match, accept an alias, or select the first result.

## Required outputs

- original input and complete transformation chain;
- a local resolved molecule-matched FASTA with a stable tip-safe identifier and a versioned accession when one was established;
- comparable-region/orientation evidence and any RNA source→derived receipt;
- for coding nucleotide, exact translation FASTA plus frame/code/translation QC;
- authoritative metadata and sequence SHA-256;
- provider request/release/snapshot, retrieval time, and raw-response checksum when available;
- explicit unresolved alternatives and limitations;
- taxonomy-resolution evidence when a name was mapped to a TaxID.

## Review gate and stop conditions

Pass Gate 1 only when molecule, identity/feature, organism, comparable sequence region, and provenance are resolved. Stop on alphabet-based molecule inference, mixed/non-comparable regions, ambiguous identity, unresolved isoforms/features, unknown source taxon for taxon-aware planning, conflicting TaxIDs, failed RNA/CDS provenance or translation QC, an unreviewed MACSE requirement, generic fold/domain-only evidence, or missing permission for a required remote operation. Never log credentials, cookies, unpublished sequences, or private paths.

## Supporting references

- [Workflow states and decision gates](../workflow.md)
- [Sequence-type and analysis-space routing](../sequence-type-routing.md)
- [Taxonomy exact-match policy](../taxonomy-resolution.md)
- [Tool and privacy boundaries](../tool-routing.md)
- [Output and provenance contract](../output-contract.md)

Primary provider documentation: NCBI Datasets gene metadata <https://www.ncbi.nlm.nih.gov/datasets/docs/how-tos/genes/get-gene-metadata/>, NCBI Taxonomy dumps <https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/>, and UniProt entry retrieval <https://www.uniprot.org/help/api_retrieve_entries>.
