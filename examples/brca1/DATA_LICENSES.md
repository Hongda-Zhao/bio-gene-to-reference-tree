# BRCA1 example data and license notice

The repository's MIT license covers the authors' original code and text. It
does **not** automatically relicense records obtained from NCBI, third-party
software, external service output, or the articles cited by this example.

NCBI states that it places no restrictions on use or distribution of the
molecular data it makes available, while also warning that submitters or other
parties may assert patent, copyright, or other rights and that NCBI cannot
transfer or grant those rights. Users remain responsible for checking the
status of particular records and their intended use. See the official
[NCBI data-use and copyright policy](https://www.ncbi.nlm.nih.gov/home/about/policies/).

| Material in this example | Source or creation route | Redistribution note |
| --- | --- | --- |
| Protein FASTA, candidate metadata, and the scoped source inventory | NCBI Gene, RefSeq, and the NCBI Datasets ortholog package retrieved 2026-08-13 | Provider-derived molecular records are redistributed for reproducibility under the NCBI policy above; the repository MIT license does not create additional rights in those records. |
| Scientific names, TaxIDs, ranks, and lineage evidence | NCBI Taxonomy snapshot `new_taxdump_2026-08-01` | Derived from the hash-recorded official taxonomy snapshot; use and attribution remain subject to the NCBI policy above. The full taxdump is not redistributed. |
| BLAST, InterProScan, MAFFT, trimAl, FastTree, IQ-TREE, and R-derived tables, trees, logs, and figures | Computed from the provider-derived sequences with the versions and commands recorded in `PROVENANCE.md` and `report/` | These are analysis outputs, but their inclusion here does not relicense embedded or underlying database records or third-party software. Consult each tool's own license before redistributing the software itself. |
| Literature evidence table and bibliography | Bibliographic facts and original summaries of linked publications | Article text is not redistributed. Each linked publication remains under its publisher or author license. |
| Scripts, tests, documentation, and original explanatory graphics | Created for this repository | Covered by the repository MIT license unless a file says otherwise. |

The large NCBI download package, complete taxdump, private scheduler logs, and
private native execution directory are not part of the public example. Their
public receipts retain only the provenance needed to identify and audit the
sources. This notice is a provenance aid, not legal advice.
