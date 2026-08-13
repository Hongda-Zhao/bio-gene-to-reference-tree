# Query resolution

Resolve every user input to a local, versioned protein record before candidate selection. Keep the original input unchanged and record the transformation chain.

## Input classifier

| Input | Required context | Resolution target | Mandatory stop |
|---|---|---|---|
| Accession | Supplied token; organism if namespace is ambiguous | Authoritative versioned record | Withdrawn, suppressed, ambiguous namespace, or sequence mismatch |
| Raw protein | Amino-acid sequence; remote-search permission | Stable database matches plus local query record | Invalid alphabet, extreme low complexity, likely nucleotide, or unresolved molecule type |
| Protein/gene name | Organism name or TaxID | Stable gene/protein ID and versioned protein | Missing taxon, multiple biologically distinct genes, or unresolved isoform |
| CDS | Genetic code, frame, and biological objective | Valid translated protein plus reversible CDS mapping | Internal stop, invalid frame, uncertain code, or probable genomic/non-coding input |

Do not decide whether a protein “belongs to a database.” Instead, identify the input type, resolve it against an authoritative record when possible, and record uncertainty.

## Accession route

1. Recognize a possible namespace only as a routing hint; confirm it with the provider.
2. Retrieve record status, accession.version, sequence, organism, TaxID, lineage, gene/protein names, reviewed/canonical/isoform status, source release or snapshot, and retrieval time.
3. Compare any user-supplied sequence with the authoritative sequence and report exact or version-specific differences.
4. Preserve replaced, merged, obsolete, or suppressed identifiers in the provenance chain.
5. Never silently substitute a canonical isoform when the requested isoform changes the biological question.

## Name route

Require a source organism or TaxID. Search stable identifiers using the supplied name, synonyms, locus tags, and gene-family terminology. Use literature to disambiguate nomenclature, function, and historical synonyms; obtain the actual sequence from an authoritative database record.

When several transcripts or isoforms resolve:

- prefer a reviewed canonical record only as a proposed default;
- show all biologically plausible alternatives;
- require approval if the alternatives change length, domain architecture, cellular localization, or tree interpretation.

Never default a symbol to human or another model organism.

## Organism-name and TaxID validation

When a TaxID must be assigned from an organism name, prefer already-extracted `names.dmp` and `nodes.dmp` files from one verified official NCBI Taxonomy snapshot. Match the supplied text character for character to `names.dmp` `name_txt`, restrict accepted rows to `name class = scientific name`, require exactly one unique TaxID, and confirm that TaxID exists in the same snapshot's `nodes.dmp`.

Do not trim, case-fold, normalize punctuation or Unicode, remove taxonomic qualifiers, accept substrings, or silently promote a synonym, equivalent name, common name, or misspelling. Secondary-name matches may be reported only as diagnostics. Stop when there is no exact scientific-name match, more than one TaxID matches, the node is absent, or an existing TaxID disagrees with the resolved TaxID.

Keep name resolution separate from phylogenetic inference. `names.dmp` maps names to TaxIDs; `nodes.dmp` validates current nodes and supplies parent/rank classification. Neither proves an evolutionary topology or makes an outgroup suitable. Record the snapshot label, exact official archive URL, retrieval time, and SHA-256 of both files in `taxonomy_resolution.tsv` and the manifest. Never invent the snapshot date from file modification times.

## Raw-protein route

1. Validate the amino-acid alphabet, sequence length, internal stop characters, ambiguous residues, low complexity, signal peptides, repeats, and probable transmembrane regions.
2. Search RefSeq protein first when that database is appropriate and available.
3. Escalate unresolved cases through reviewed Swiss-Prot/UniProt records, broader UniProtKB or `nr`, then InterPro/Pfam or iterative profile search.
4. Use structure-aware search only after sequence/profile evidence is insufficient, and treat fold similarity as supporting evidence rather than automatic proof of homology.
5. Record the search database snapshot, parameters, E-value, bit score, query and target coverage, percent identity, and aligned coordinates.

Do not infer orthology from the nearest search hit. Use the search to identify candidate family membership, then apply orthology evidence and tree-aware review.

## Evidence ladder

Keep evidence types separate:

1. exact authoritative accession resolution;
2. curated orthology or reviewed protein annotation;
3. strong full-length similarity with compatible domain architecture;
4. profile/domain evidence;
5. structure-aware evidence;
6. name/literature context without a resolved record.

Do not collapse these levels into one confidence number. Stop before tree planning when the evidence only supports a generic fold, domain, or broad superfamily that does not answer the declared objective.

## CDS and codon analyses

Translate CDS with the correct genetic code, verify length modulo three and internal stops, align the proteins, and back-translate with a codon-aware method when nucleotide/codon analysis is required. Do not directly align divergent CDS as ordinary DNA. Route genomic or non-coding inputs to an annotation or RNA-specific workflow.

## Viral route

Record virus taxon, host, genome type, segment, ORF coordinates, polyprotein processing, and gene/protein name. Do not combine non-homologous segments or different mature proteins from a polyprotein. Before interpreting a viral gene tree, search for known recombination, reassortment, segmentation, and mosaic ancestry; route to a recombination-aware workflow when these processes invalidate a single-tree assumption.

## Privacy and provenance

Obtain explicit permission before transmitting unpublished sequences or unpublished metadata. When remote access is disallowed, require local FASTA and candidate metadata. Record sequence SHA-256, provider request identifier, release/snapshot, retrieval time, and raw-response checksum. Never log API keys, authorization headers, cookies, private filesystem paths, or unpublished sequence content in command logs.

## Primary documentation

- NCBI Datasets gene metadata: <https://www.ncbi.nlm.nih.gov/datasets/docs/how-tos/genes/get-gene-metadata/>
- NCBI Taxonomy dump files: <https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/>
- NCBI taxdump field definitions: <https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/taxdump_readme.txt>
- NCBI BLAST database descriptions: <https://www.ncbi.nlm.nih.gov/books/NBK62345/>
- UniProt REST entry retrieval: <https://www.uniprot.org/help/api_retrieve_entries>
- InterProScan introduction: <https://interproscan-docs.readthedocs.io/en/v5/Introduction.html>
