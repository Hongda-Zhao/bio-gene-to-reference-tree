# BRCA1 example provenance

This receipt separates immutable facts about the public inputs from results
created by later analysis stages. It intentionally omits usernames, login
addresses, credentials, unpublished data, and absolute execution paths.

## Source acquisition

| Field | Value |
| --- | --- |
| Retrieval date | 2026-08-13 |
| Query | Human BRCA1, NCBI GeneID 672 |
| Query protein | RefSeq `NP_009225.1`, 1,863 aa |
| Query annotation | `GCF_000001405.40-RS_2025_08`, annotation release date 2025-08-01 |
| Client | NCBI Datasets CLI 18.6.0 |
| Command | `datasets download gene gene-id 672 --ortholog all --include protein,rna,gene --filename brca1-orthologs.zip` |
| Download result | 558 gene records; 3,605 protein FASTA records before review |
| Package SHA-256 | `30f1ef86dbe2989981bed005c837e7538226cae24b800aef44a32a589bfe04c0` |
| Provider documentation | [Retrieve ortholog data and metadata](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/tutorials/download-ortholog-dataset/) |
| Ortholog method | [How are orthologs calculated?](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/ortholog-calculation/) |

The download command is a retrieval operation, not a phylogenetic inference.
NCBI ortholog membership was used to discover candidates. It did not constrain
the BRCA1 alignment or tree topology.

The downloaded ZIP and raw provider payload are not committed. The fixed
selected sequences, per-sequence hashes, GeneIDs, TaxIDs, annotation labels,
and selection notes are retained in
[`inputs/candidate_provenance.tsv`](inputs/candidate_provenance.tsv).

## NCBI taxonomy snapshot

All 18 candidate organism names and TaxIDs were checked against two files from
one dated NCBI taxonomy archive.

| Field | Value |
| --- | --- |
| Snapshot label | `new_taxdump_2026-08-01` |
| Archive URL | <https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/new_taxdump_2026-08-01.zip> |
| Retrieved | 2026-08-13 |
| Archive size | 156,769,236 bytes |
| Archive SHA-256 | `e67d56c3e87eac14a28feea8cf710ac612571f00a4f32ba3b85bae0020cd123a` |
| `names.dmp` SHA-256 | `c3aac85da1bf473792e95bf8fe2d6933dcb31f2d649a3f8579e21011d772dca9` |
| `nodes.dmp` SHA-256 | `502ef4e7f1a8383e89d759b056524e8bcc16fe6950965f972d13f44411430f17` |
| Match policy | Character-for-character `scientific name` match plus exact TaxID agreement |

The dated archive directory did not provide a sidecar MD5 file for this ZIP.
No provider checksum is invented; the table records SHA-256 values calculated
for the downloaded archive and extracted files. The archive itself and full
taxdump are not redistributed in this example.

The successful per-tip resolution table is emitted as
`review/taxonomy_resolution.tsv` and copied to
`annotation/taxonomy_resolution.tsv` for the final bundle. Validation is
all-or-nothing: a synonym-only hit, alias, fuzzy match, missing node, malformed
dump, mixed snapshot, or TaxID mismatch blocks planning.

## Fixed sequence set

The materialized analysis set contains 18 proteins chosen as a manual
fixed-tip demonstration:

- one study sequence: human `NP_009225.1`;
- 15 expanded amniote ortholog representatives;
- two proposed amphibian outgroups, *Xenopus tropicalis*
  `NP_001107963.1` and *Rhinatrema bivittatum* `XP_029429046.1`.

The set is frozen by accession version and amino-acid SHA-256 in
[`inputs/candidate_provenance.tsv`](inputs/candidate_provenance.tsv). It spans
16 amniote tips total, not 18 amniotes plus two extra outgroups.

The acquisition audit
[`inputs/source_protein_inventory.tsv`](inputs/source_protein_inventory.tsv)
records every provider protein for all 18 selected species plus platypus and
tuatara. It contains 518 present records: 18 selected proteins, 498
not-promoted alternatives from selected species, and two platypus proteins.
Those alternatives are inventory records, not inferred QC failures. The two
platypus records (`XP_028930515.1` and `XP_028930516.1`, both 1,346 aa) were
excluded solely because their length ratio 0.722491 was below the configured
0.75 minimum; no domain failure is claimed. *Sphenodon punctatus* TaxID 8508
was resolvable in the frozen taxdump but absent from both the frozen gene report
and protein FASTA. The shorter *Pseudonaja textilis* model was eligible only if
its N-terminal RING and tandem C-terminal BRCT evidence passed the same
executable QC as every other tip.

The fixed post-review set is below the request's 200-sequence MMseqs2 trigger.
The 3,605 proteins in the provider package are raw archive breadth, not 3,605
eligible analysis candidates.

This inventory is intentionally scope-limited to the selected species plus two
explicitly screened taxa. It does not enumerate all taxa in the 3,605-protein
archive, and the planner's 18 selected, zero rejected, and zero unsampled-clade
counts describe only its supplied 18-row candidate table. They are not
archive-wide discovery results.

Two metadata fields need narrow interpretation. `relation=one2one_ortholog`
is the analyst's fixed-set selection label for members returned through the
NCBI Ortholog group; the retained provider payload does not expose a
protein-accession-specific pairwise evidence record, so the field is not proof
that every method component was observed for that exact pair. Likewise,
`source_release` is the first gene-record annotation context reported for the
taxon, not a demonstrated accession-specific product-to-assembly linkage.
Neither field was used as a topological constraint.

The plan-bound taxonomy resolver records exact species-name/TaxID matches,
direct parent TaxIDs, and ranks; its legacy `lineage` column remains empty and
`clade` remains an analyst-assigned sampling label. A post-plan, explicitly
non-decision-bearing traversal of the same hash-verified `names.dmp` and
`nodes.dmp` produced
[`annotation/taxonomy_lineage.tsv`](annotation/taxonomy_lineage.tsv), which
contains the complete named/ranked root-to-tip lineage for all 18 accessions.
Its SHA-256 is
`c547db0a97c4b9f5db3951495e671bbc490ff90b88a2830f7e9c49714badd62b`.

The deterministic QC summarizer originally wrote its tables under
`qc/summary/`; the historical utility and time that copied them to
`qc/candidate_qc.tsv` and `inputs/candidates.tsv` were not retained. The
post-hoc [`report/qc_promotion_receipt.json`](report/qc_promotion_receipt.json)
regenerates both outputs from hash-recorded inputs and proves byte equality to
the promoted destinations. It does not invent an exact replay of the original
copy operation.

## Planner and execution boundary

The bundled `gene_to_tree.py plan` command is a deterministic, offline plan
compiler. For this example it validates the local query, candidate FASTA and
metadata, exact NCBI taxonomy evidence, selection thresholds, outgroup roles,
and planned command arrays. It writes a `pending-reference-approval` bundle and
does not contact NCBI, run an aligner, submit batch jobs, infer a tree, or render
a figure.

Actual sequence QC, alignment, trimming, exploratory topology screening, and
final tree inference run as host-side jobs. Compute-heavy commands were staged
on Lustre and submitted to PBS compute nodes on the gds2 environment; they were
not run interactively on the login node. Public scripts accept a caller-owned
run directory through `BRCA1_RUN_DIR`, while this receipt deliberately omits
the actual path and account identity.

The committed PBS files are gds2 execution receipts and site templates. Their
queue, module names, module-init path, Lustre assumptions, and GNU
`sha256sum` use are not a portability claim; another cluster must adapt them.
Raw scheduler output belongs in the ignored `pbs-logs/` directory and is not a
public artifact.

Prospective approval must be tied to both the final `plan_hash` and the selected
alignment SHA-256. A changed sequence, candidate field, taxonomy dump,
threshold, role, outgroup, command, MSA, trim profile, model, support method, or
seed invalidates exact conformance with that approval.

The historical receipts in this example approve the reviewed reference set,
scientific policy, and selected alignment, but the executed host argv differed
from the planned argv. They are not reinterpreted as prospective authorization
for those command changes. Instead,
[`report/execution_reconciliation.json`](report/execution_reconciliation.json)
binds the old plan hash, exact planned and normalized executed argv, input and
output hashes, timing, and a post-hoc reviewer acknowledgement. Its final status
must state that the protocol deviation was accepted for publication, not that
the launch was retroactively approved. The receipt is not retroactive approval.

### Execution deviations retained in the audit

- The first alignment submission stopped before MAFFT because its trimAl
  version probe used an option not supported by the pinned module. It produced
  no alignment or tree promoted into this example. The public script records
  the verified module version without that probe, and the replacement run
  completed normally.
- The login and batch environments did not expose the same durable home path,
  so materialized inputs were staged in shared Lustre storage for compute and
  copied back after validation. Public records keep content hashes and
  repository-relative paths, not account-specific storage paths.
- The planned MAFFT and trimAl paths were relocated in the host run. trimAl also
  added `-fasta`, an output-serialization option that did not change any gap
  threshold. These are exact-command deviations even though the input/output
  bytes and scientific thresholds are hash-bound.
- The approved plan named `iqtree2 ... --prefix gene-tree`; the successful host
  run used IQ-TREE 2.4.0 as `iqtree2 ... --prefix brca1-balanced` with the
  approved alignment at a private hash-bound path. These path and naming
  differences are non-identical argv, so no blanket scientific equivalence is
  claimed, even though every explicit model, support, thread, and seed setting
  was unchanged.
- An earlier IQ-TREE 3.1.3 (`iqtree3`) attempt ran single-threaded despite the
  requested setting and was canceled after about four hours with exit 271. It
  produced no promoted output; the command ledger and reconciliation retain it
  only as a failed-attempt audit record.

These deviations did not change the frozen accessions, sequences, taxdump
snapshot, approved alignment bytes, or explicit inference settings. That fact
narrows the deviation; it does not convert post-hoc review into prospective
command approval.

The canceled IQ-TREE 3 attempt began at `2026-08-13T09:12:23Z`; the
domain-retention receipt was created at `2026-08-13T09:17:42Z`, so it did not
gate that attempt. All 162 checks passed without changing the MSA or inference
arguments. The successful IQ-TREE 2 run began later, at
`2026-08-13T13:10:04Z`, and the committed reusable PBS script now makes this
receipt a prerequisite.

## Verified software snapshot

Availability and versions were inspected on the execution environment before
submission. The final report and command ledger record the versions actually
used by each completed job.

| Component | Inspected version or module |
| --- | --- |
| NCBI Datasets CLI | 18.6.0 |
| BLASTP | BLAST+ 2.17.0 |
| InterProScan | 5.72-103.0; Pfam and SMART applications |
| MAFFT | 7.526 |
| trimAl | 1.5.1 |
| FastTree | 2.2.0 |
| IQ-TREE | 2.4.0, executable `iqtree2`; canceled non-promoted attempt used 3.1.3 `iqtree3` |
| R | 4.5.1 available on gds2; the separately recorded local rendering runtime is authoritative for the figure |

The promoted primary tree was produced by “IQ-TREE multicore version 2.4.0 for
Linux x86 64-bit built Feb 7 2025”; its public version receipt has SHA-256
`4ec85e058a26f2eb7dc3c0c674e2be03d7dc0698812a87dadf2dfed115a969ed`.
The figure was rendered under R 4.4.1 with ape 5.8.1, ggplot2 4.0.1, ggtree
3.14.0, openssl 2.3.4, and svglite 2.2.2, as recorded in its settings receipt.

## Integrity links filled by execution

These markers are replaced only after the named outputs exist and their
checksums have been added to `report/checksums.sha256`.

Planner run ID: `gtr-51d4906dff86327a`.

Plan hash: `8b03bd3fae2882a07351fc83f1b8879ee616ed61f60570d029c5ef8471f65da8`.

Approved `reference_set.faa` SHA-256:
`cea86b4416bcde376e6d8177a2b51d5cd2a501e458af6da9405e9acd6ed69dd8`.

Raw E-INS-i alignment SHA-256:
`9fd8ffaf87c8400175ddb2f42138e0d63bc7d6598d1d2fc3a44cb8f57c304b10`.

Approved profile: balanced trimAl `-gt 0.5`; SHA-256
`b3ea66343c7b4ad4129459715d5de7fdf2f053d9d1eb8b77a0a5637b5ef84dd6`.

Preserved unrooted ML tree SHA-256:
`5df0dd2067f52c008637d81d944fc532ef2e3b67fa242b5eb40ec74b7d494380`.

Validated amphibian-rooted display derivative SHA-256:
`8738dfc3e1a18a4d58a200042d644d7557b6806e8be75d2058b805de2f5d5202`.

Sequence metadata SHA-256:
`e3f304c50a36d6383311fc11622a6124c03ec595d74a1c2193ba18dda5350d9f`;
iTOL roles SHA-256:
`7c9c2f4ef51cf911cb0ce28ddb130b98d4d665c33ec619d71cafed5c65cb72a8`;
full taxonomic lineages SHA-256:
`c547db0a97c4b9f5db3951495e671bbc490ff90b88a2830f7e9c49714badd62b`;
SVG/PDF/settings SHA-256 values:
`817a49b9a50aada987bc01209aeb005f98aa26f99592b661d4786fd6f632dad4`,
`640efa23e318bb81db80860ad4ae3f9da9a5046b3f36c7c090c77046ef96be6f`,
and `b4c3196aa034237f0bd7828a40321cd4cca2a0dfd674e73845382a68832cefc6`.

Final normalized command ledger SHA-256:
`4d14bbb9fa5c6daf35f916a4a70fc67026c17852ec34f41a90b458d90f00be7a`;
executed-workflow report SHA-256:
`c34061264585fd066fd010ff8911c33aed8aa4b77c806e9b4683072c3701fbeb`.

## Root-state provenance

IQ-TREE infers the preserved unrooted ML result. A rooted file is permitted
only after checking that both approved amphibian outgroups form the exclusive
split opposite all 16 amniote tips and that reasonable trim-profile screening
does not make that placement unstable. If the test fails, the absence of a
rooted artifact is an expected audited result, not a missing file.

The implemented gate checks that isolating split in the final tree and four
FastTree trim-profile trees. It does not implement an explicit long-branch
test, outgroup-removal analysis, alternative-outgroup comparison, or
richer-model rooting sensitivity. A passing rooted derivative is therefore a
topology-based display hypothesis, not evidence that long-branch attraction
has been excluded or that the biological root is resolved.

The ggtree/ggplot2 renderer cannot reroot, ladderize, or guess support
semantics. Its settings TSV must record the input Newick hash, metadata hash,
optional iTOL role hash, declared root state, branch-length mode, support
format, layout, canvas dimensions, palette source, R version, and package
versions.

Raw IQ-TREE outputs remain in durable execution storage. Public text reports
are exact-content copies except for documented private path and compute-node
redactions; `tree/iqtree/redactions.tsv` retains raw and public SHA-256 values.
The safe manifest in `tree/iqtree/raw-output-receipt.json` binds all raw regular
files by filename, byte size, and SHA-256, including omitted checkpoints,
wrapper checksums, and compressed model internals.

## Database mutability and interpretation

NCBI Gene, RefSeq annotations, NCBI Orthologs, UniProt, and literature indexes
are dynamic. The retrieval date, accession versions, release labels, dated
taxdump, and hashes make this run reviewable; they do not imply that a future
query will return the same 558 records or the same models.

Finally, this is BRCA1 locus history. Its long variable central region,
lineage-specific selection, annotation errors, paralogy mistakes, alignment
uncertainty, and stochastic gene-tree discordance can all produce differences
from accepted species relationships. External phylogenomics is therefore
recorded as comparison evidence in
[`evidence/literature_evidence.tsv`](evidence/literature_evidence.tsv), never as
a constraint that forces the BRCA1 tree to agree.

Provider-derived data and original repository material have distinct license
boundaries documented in [`DATA_LICENSES.md`](DATA_LICENSES.md).
