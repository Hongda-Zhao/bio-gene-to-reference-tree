# Sequence-type and analysis-space routing

Read this reference before query resolution whenever a request may contain nucleotide sequence. Sequence letters cannot identify molecule type reliably: the nucleotide IUPAC alphabet overlaps valid protein symbols. Require an explicit declaration and keep the analysis in that declared space.

## Contents

- [Version boundary](#version-boundary)
- [Mandatory intake record](#mandatory-intake-record)
- [Supported routes](#supported-routes)
- [Comparable-region gate](#comparable-region-gate)
- [RNA normalization provenance](#rna-normalization-provenance)
- [Coding-sequence QC and MACSE boundary](#coding-sequence-qc-and-macse-boundary)
- [Search and database compatibility](#search-and-database-compatibility)
- [Clustering, alignment, trimming, and inference](#clustering-alignment-trimming-and-inference)
- [Approval and failure rules](#approval-and-failure-rules)
- [Primary documentation](#primary-documentation)

## Version boundary

Requests using schema 0.1 or 0.2 retain the legacy protein-only contract and outputs. Do not reinterpret an existing legacy request as nucleotide because its letters resemble DNA or RNA.

Request schema 0.3 requires one explicit `molecule.type`: `protein`, `noncoding-dna`, `noncoding-rna`, `coding-dna`, or `coding-rna`. It supports three analysis spaces: protein, direct nucleotide, and codon. A change of molecule type or analysis space is a new decision-bearing request, not a fallback.

## Mandatory intake record

Before database search or local validation, record:

- the exact molecule type and intended analysis space;
- coding status, feature/locus identity, coordinates, strand/orientation, and comparable region;
- source organism/TaxID and accession/version when known;
- original file path, identifier set, byte hash, and per-sequence hash;
- for CDS, the genetic-code ID, frame source, completeness, translation source, and translation FASTA;
- for RNA, `source_encoding: rna-u` or `source_encoding: dna-t`, and the separately named DNA-alphabet analysis copy;
- the approved search program, database molecule, and actual database label. Physical index/build evidence belongs to the host-side acquisition record when available.

Never infer molecule type from alphabet, filename, accession shape, gene name, or the closest hit. Every candidate row must declare the same molecule type as the request.

## Supported routes

| Declared input | Analysis route | Alignment and trimming | Tree model space |
|---|---|---|---|
| `protein` | Protein | MAFFT `--amino`; optional trimAl columns | FastTree protein model or IQ-TREE `-st AA` |
| `noncoding-dna` | Direct nucleotide | MAFFT `--nuc`; optional trimAl columns | FastTree `-nt -gtr` or IQ-TREE `-st DNA` |
| `noncoding-rna` | Direct nucleotide on a provenance-bound DNA-alphabet analysis copy | MAFFT `--nuc`, or the dedicated `mafft-qinsi --nuc` executable when `rna_structure_aware: true`; optional trimAl columns | FastTree `-nt -gtr` or IQ-TREE `-st DNA` |
| clean `coding-dna` or `coding-rna`, code 1/11, `analysis_kind: codon` | Protein-guided, codon-preserving | MAFFT `--amino` on verified translations; trimAl `-backtrans` with matched CDS | IQ-TREE `-st CODON<n>` |
| clean `coding-dna` or `coding-rna`, code 1/11, `analysis_kind: nucleotide` | Protein-guided, codon-preserving alignment analyzed as nucleotide sites | Same MAFFT/backtranslation route; triplet QC is still mandatory | FastTree `-nt -gtr` or IQ-TREE `-st DNA`; this is not a codon model |
| other-code clean CDS, or frameshifted/disrupted/uncertain CDS | Reviewed MACSE/PAL2NAL/precomputed route | Review outputs remain separate artifacts | No automatic tree plan |

Do not substitute between rows. In particular, direct nucleotide MAFFT is not a rescue for failed CDS translation, a protein tree is not a fast fallback for a codon request, and FastTree nucleotide mode is not a codon model. To analyze another space, create and approve a separate request with its own hashes and interpretation.

## Comparable-region gate

Every selected sequence must represent the same homologous biological region. Record the feature type and region definition, including exon/intron/UTR status, locus or viral-segment boundaries, coordinates, strand, and orientation where applicable.

Do not mix genomic loci, spliced transcripts, CDS, mature RNA, amplicons, partial domains, or different viral segments merely because a local similarity exists. Extract a documented comparable region first, preserve the source-to-derived coordinate map, and re-run reference approval. Stop when homology or boundaries cannot be established.

## RNA normalization provenance

Require one explicit source encoding. `rna-u` accepts RNA-IUPAC input with `U`; `dna-t` accepts an RNA record archived with DNA-style `T`. Preserve that source encoding unchanged in `reference_set.rna.fasta`, then create a separately named `reference_set.fna` analysis copy. For `rna-u`, convert only `U`→`T`; for `dna-t`, make a byte-equivalent T-preserving copy. Record:

```text
source_path, source_sha256, derived_path, derived_sha256,
source_encoding, transformation, transformed_record_count,
identifier_map, tool_or_procedure, created_at
```

For `rna-u`, the transformation may change only `U`→`T`; for `dna-t`, source and analysis-copy sequence bytes remain equal. Identifiers and sequence lengths remain identical in both cases. Label the computation as derived from RNA in the plan, metadata, manifest, and report. Never overwrite the source or call the analysis copy raw DNA.

For coding RNA, materialize the DNA-alphabet CDS copy before translation/backtranslation and retain the same provenance. For noncoding RNA with `rna_structure_aware: true`, use the dedicated `mafft-qinsi` executable; `mafft --qinsi` is not a valid substitute. Covariance-model RNA phylogenetics beyond Q-INS-i remains a separate reviewed workflow.

## Coding-sequence QC and MACSE boundary

Accept a CDS for the clean codon route only after all of the following are explicit and pass for every sequence:

- one approved NCBI genetic-code ID and recorded strand/frame source; the deterministic trimAl backtranslation planner accepts only codes 1 and 11;
- comparable complete CDS boundaries and length divisible by three;
- unique, identical first-token IDs in CDS and translation FASTAs;
- no frameshift and no unexpected internal stop;
- translation produced with the declared code matches the supplied protein exactly, apart from a documented terminal stop convention;
- CDS contains no alignment gaps before backtranslation.

The deterministic code-1/11 commands use trimAl `-ignorestopcodon` only after the planner has independently established `internal_stop_count=0` and exact translation equality; this permits the common terminal-stop convention without weakening the input gate. Treat warnings about truncation, padding, unmatched or duplicate IDs, triplet failure, or translation mismatch as QC failures. Because trimAl's stop handling uses the universal stop set, every other genetic code must use a reviewed MACSE, PAL2NAL, or precomputed codon-alignment handoff rather than the deterministic `-backtrans` planner.

If any sequence has a frameshift, internal stop, disrupted ORF, pseudogene annotation, uncertain frame, or failed translation equality, emit `MACSE_ROUTE_REQUIRED`. Preserve the original CDS and review [MACSE](https://www.agap-ge2pop.org/alignsequences/) output, frameshift symbols, stop codons, and corrections separately. The MACSE route is fail-closed: after review, apply an explicit export policy, verify exact IDs and complete codon triplets, and rerun routing with the result declared as a typed precomputed codon alignment. Do not silently promote MACSE output into the clean codon route or align the disrupted CDS as ordinary noncoding nucleotide.

## Search and database compatibility

A PATH-visible BLAST executable is not a database. Request 0.3 requires every candidate row to declare `actual_search_database`, `search_program`, and `search_database_molecule`. The helper checks that fields are present and that the program/database-molecule pair is supported:

- `blastp`→`protein`, `blastn`→`nucleotide`, `blastx`→`protein`, `tblastn`/`tblastx`→`nucleotide`;
- `mmseqs2` or `curated-provider`→`protein` or `nucleotide`, subject to the declared input route;
- `not-searched`→`not-applicable`, with `actual_search_database: local-bundle:not-searched`.

Protein, noncoding, and coding inputs allow only their molecule-compatible subset of those pairs. The helper does not execute the search or independently verify the database label, release, physical index, or checksum. Preserve name/release/build/index/checksum details in the host-side acquisition provenance when known.

- Protein discovery normally uses `blastp` against a verified protein database.
- Direct-nucleotide discovery normally uses `blastn` against a verified nucleotide database.
- A reviewed `blastx`, protein-translation-driven `tblastn`, or translated MMseqs2 search may support CDS discovery, but final candidates still require molecule-matched retrieval and CDS/translation QC.

Never run `blastp` against a nucleotide index, `blastn` against a protein index, or claim RefSeq search when the actual indexed database is unknown. If the required compatible database is absent, block that route; do not switch analysis space.

## Clustering, alignment, trimming, and inference

MMseqs2 may cluster protein or nucleotide expanded candidates only when its input type and index are verified. Record identity, coverage, coverage mode, and the molecule type; protect `study` and `outgroup` sequences. `-c` remains coverage, not identity.

Use explicit molecule flags and record argv arrays:

```text
mafft --amino --auto reference_set.faa
mafft --nuc --auto reference_set.fna
mafft-qinsi --nuc --thread <n> reference_set.fna
trimal -in alignment.raw.translated.faa -backtrans reference_set.fna -ignorestopcodon \
  -out alignment.trimmed.<profile>.codon.fna -gt <reviewed> -fasta
FastTree -nt -gtr alignment.<approved>.fna
iqtree2 -s alignment.<approved>.fna -st DNA -m MFP ...
iqtree2 -s alignment.<approved>.codon.fna -st DNA -m MFP ...
iqtree2 -s alignment.<approved>.codon.fna -st CODON<n> -m MFP ...
```

For the codon route, save the raw protein MSA, the untrimmed backtranslation, any trimmed protein MSA, and every corresponding codon alignment. The same trim profile must define both the reviewed protein columns and its backtranslated codon output. Verify identical ID sets, codon-alignment length divisible by three, intact triplets, and translation equality after every backtranslation.

`CODON<n>` must use the approved genetic-code ID; code 1 emits explicit `CODON1`, never bare `CODON`. FastTree has no codon model, so absence of IQ-TREE blocks codon inference. Direct DNA/RNA and an explicitly requested nucleotide-site CDS analysis use `-st DNA`, never `-st CODON<n>`; the latter still uses a codon-preserving backtranslation and must not be described as a codon substitution model.

## Approval and failure rules

Invalidate reference or alignment approval when molecule declaration, analysis space, region boundaries, orientation, RNA normalization, CDS frame/code/translation, search database, identifier mapping, or backtranslation changes.

Stop before planning or execution on mixed molecule declarations, non-comparable regions, missing RNA source-encoding/analysis-copy provenance, incompatible or unidentified search databases, failed CDS QC, a non-1/11 code without a reviewed alignment handoff, an unreviewed MACSE requirement, or any attempted cross-analysis fallback. Preserve the blocker, original inputs, and review evidence in the audit bundle.

## Primary documentation

- [MAFFT input-type options](https://mafft.cbrc.jp/alignment/software/manual/manual.html)
- [trimAl usage, including `-backtrans`](https://trimal.readthedocs.io/en/latest/usage.html)
- [FastTree nucleotide models](https://morgannprice.github.io/fasttree/)
- [IQ-TREE command reference](https://iqtree.github.io/doc/Command-Reference) and [substitution models](https://iqtree.github.io/doc/Substitution-Models)
- [NCBI BLAST program/database molecule matrix](https://www.ncbi.nlm.nih.gov/BLAST/about/)
- [NCBI BLAST database descriptions](https://www.ncbi.nlm.nih.gov/books/NBK62345/)
- [MMseqs2 user guide](https://github.com/soedinglab/MMseqs2/wiki)
- [MACSE documentation](https://www.agap-ge2pop.org/macse/macse-documentation/)
