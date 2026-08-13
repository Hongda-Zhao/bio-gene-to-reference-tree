# BRCA1 executed-workflow report

This report summarizes the public, fixed 18-tip BRCA1 protein gene-tree run.
It is a single-locus protein gene tree, not a species-tree estimate.
Machine-readable argument arrays are retained in
[commands.jsonl](commands.jsonl). Each record states whether it is replayable
or provenance-only; caller-owned directories in provenance records are
replaced with declared public path tokens. Content hashes are retained in
`checksums.sha256` after finalization.

## Query, references, and provenance

- Query: human BRCA1 RefSeq `NP_009225.1` (GeneID 672; TaxID 9606).
- Retrieval: NCBI Datasets CLI on 2026-08-13; 558 ortholog gene records and
  3,605 protein records before fixed-tip and isoform review.
- Analysis set: 18 proteins—16 amniotes and two amphibian outgroup candidates.
- Taxonomy: all 18 scientific names and TaxIDs matched exactly against the
  dated NCBI `new_taxdump_2026-08-01` snapshot; a separate hash-bound table
  records every full named/ranked root-to-tip lineage from that snapshot.
- Planner: run `gtr-51d4906dff86327a`, plan hash
  `8b03bd3fae2882a07351fc83f1b8879ee616ed61f60570d029c5ef8471f65da8`.

The raw provider archive and full taxonomy dump are intentionally omitted;
their source URLs and hashes are recorded in [../PROVENANCE.md](../PROVENANCE.md).
The 18 tips were manually fixed before planning. The scoped 518-protein
inventory covers the targeted species plus two screened taxa, not the whole
3,605-protein archive; planner counts of 18 selected and zero rejected apply
only to the supplied candidate table.

## Candidate and alignment QC

All 18 proteins passed the same N-terminal RING plus tandem C-terminal BRCT
architecture gate. The balanced trimAl `-gt 0.5` profile was approved before
final inference. It retains 1,851 of 2,179 raw E-INS-i columns (84.95%), has
mean occupancy 94.45%, and contains 1,418 parsimony-informative columns. Every
audited RING/BRCT interval retained at least 80% of its residues under all
three trim profiles; the observed minimum was 93.68%. That projection audit
finished five minutes after final-tree submission, so it validates the same
unchanged, hash-bound alignment but is not presented as a launch gate.

The summarizer first wrote `qc/summary/candidate_qc.tsv` and
`qc/summary/candidates.tsv`. The historical copy utility and time were not
retained. The post-hoc
[`qc_promotion_receipt.json`](qc_promotion_receipt.json) deterministically
regenerates both files and verifies byte identity with `qc/candidate_qc.tsv`
and `inputs/candidates.tsv`; it does not invent an exact historical replay.

Exploratory FastTree screening found normalized unrooted RF distance 0.066667
(one split) between the balanced tree and each of the raw and permissive trees;
the strict and balanced topologies were identical. The amphibian pair formed
an isolating outgroup split under all four exploratory profiles. FastTree
SH-like support is not interpreted as bootstrap support.

That screen is the complete implemented rooting gate: it tests the isolating
split across trim profiles but includes no explicit long-branch diagnostic,
outgroup-removal analysis, alternative-outgroup comparison, or richer-model
rooting sensitivity. Any rooted derivative is therefore a topology-based
display hypothesis; the biological root remains provisional.

Similarity coverage is an identity check, not an alignment-adequacy guarantee.
For *Crocodylus porosus*, BLAST interval-union coverage was only 0.112721 of
the human query and 0.113821 of the target (210 HSP-summed columns), even
though the balanced MSA retained 1,737 of 1,845 crocodile residues. NCBI
ortholog-group membership plus terminal domains supports the BRCA1 identity;
occupancy trimming alone does not establish positional homology across the
long variable center. The full-length deep topology remains provisional
without a conserved-block or terminal-domain-only sensitivity analysis.

## Primary maximum-likelihood result

IQ-TREE 2.4.0 selected `Q.bird+F+I+R3` according to BIC. The final ML tree
had log-likelihood −45,445.3350, standard-error field 432.3746, and total tree
length 10.1282 substitutions per site. The 8-thread run consumed 789.735 CPU
seconds and 99.242 wall-clock seconds and completed 1,000 UFBoot2 replicates
with `-bnni` plus 1,000 SH-aLRT replicates. The report defines internal labels
as SH-aLRT/UFBoot percentages. No bootstrap-correlation coefficient and no
non-convergence warning were emitted, so numerical UFBoot convergence is not
claimed. Snake `XP_026576759.1` was the sole composition chi-square failure
(1/18; *p* < 0.05, df = 19). The only emitted warning class was the same
state-frequency-normalization warning repeated 1,025 times across model testing
and the subsequent UFBoot/tree search; it is recorded neutrally, and no
`ERROR` was present.

The preserved unrooted tree SHA-256 is
`5df0dd2067f52c008637d81d944fc532ef2e3b67fa242b5eb40ec74b7d494380`.
The amphibian isolating split passed in the IQ-TREE topology and all four
FastTree trim screens, permitting a separately written rooted derivative with
SHA-256 `8738dfc3e1a18a4d58a200042d644d7557b6806e8be75d2058b805de2f5d5202`.
Each FastTree topology was RF 8/30 from the final IQ-TREE tree. Root approval
is topology-based and does not include a long-branch, outgroup-removal,
alternative-outgroup, or richer-model sensitivity analysis.

Recovered prespecified clades and SH-aLRT/UFBoot support were: Amphibia
99.6/100, Mammalia 100/100, Sauropsida 97.1/97, Primates 100/100, Glires
80.3/77, Aves 100/100, Lepidosauria 100/100, Archelosauria 99.9/100, and
Archosauria 77.3/79. The explicit Crocodylia+Testudines alternative was not
recovered. Glires and Archosauria therefore lack the joint ≥80/≥95 support
combination used here as a conservative descriptive threshold. This protein
gene tree remains provisional for deep relationships because snake composition
is heterogeneous, crocodile full-length similarity coverage is low, and no
conserved-block or terminal-domain-only sensitivity tree was run.

## Figure and annotations

The verified 14 × 10 inch vector outputs are
[`gene-tree.outgroup-rooted.ggtree.svg`](../figures/gene-tree.outgroup-rooted.ggtree.svg)
(SHA-256 `817a49b9a50aada987bc01209aeb005f98aa26f99592b661d4786fd6f632dad4`)
and [`gene-tree.outgroup-rooted.ggtree.pdf`](../figures/gene-tree.outgroup-rooted.ggtree.pdf)
(SHA-256 `640efa23e318bb81db80860ad4ae3f9da9a5046b3f36c7c090c77046ef96be6f`).
The settings receipt SHA-256 is
`b4c3196aa034237f0bd7828a40321cd4cca2a0dfd674e73845382a68832cefc6`.
Visual inspection confirmed all 18 labels, the three role colors, support
labels, scale axis, and unclipped tips; the PDF carries no author, subject, or
keyword field and no private execution residue.

The figure role colors are fixed: study orange `#E69F00`, expanded green
`#009E73`, and outgroup gray `#999999`. Branch lengths are substitutions per
site. Internal labels are shown only after confirming IQ-TREE's emitted
SH-aLRT/UFBoot ordering.

## Reproducibility boundary

The bundled skill helper compiled and audited the offline plan; an authorized
host executed the declared external tools on batch compute. Public command
records normalize caller-owned directories to repository-relative paths and
do not expose usernames, hostnames, credentials, scheduler identifiers, or
private absolute paths. Database records may change after the frozen retrieval
date, so reruns must create new provenance rather than overwrite this example.
The command ledger consequently marks the immutable historical plan compiler
entry as normalized provenance with an explicit output token, not as a replay
command against the committed `review/` directory.

This run has an explicit protocol deviation. The prospective plan bound the
reference set and scientific settings but named different exact MAFFT/trimAl
paths, omitted trimAl's output-format-only `-fasta` flag, and specified
`iqtree2 --prefix gene-tree` rather than the successful IQ-TREE 2.4.0 command
with prefix `brca1-balanced` and a private hash-bound input path. No scientific equivalence
is asserted for non-identical argv. The post-hoc
[`execution_reconciliation.json`](execution_reconciliation.json) binds both
argv sets and their inputs, but does not retroactively authorize the launch.
It also records a prior IQ-TREE 3.1.3 (`iqtree3`) attempt that was canceled
after running single-threaded; no output from that attempt was promoted.

IQ-TREE text reports in the public bundle are path/host-redacted publication
copies, with raw/public hashes in `../tree/iqtree/redactions.tsv`. A safe
per-file raw manifest in `../tree/iqtree/raw-output-receipt.json` binds files
that were not published, including checkpoints and compressed model internals.

One initial alignment submission failed before MAFFT because a trimAl version
probe was incompatible with the pinned module. No result from that attempt was
promoted. The corrected public script and its replacement run are the source
of every alignment and exploratory tree in this bundle.

See [../DATA_LICENSES.md](../DATA_LICENSES.md) for the separate license and
redistribution boundaries of provider-derived records, analysis outputs,
literature references, and original repository material.
