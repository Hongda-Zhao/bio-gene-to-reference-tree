# BRCA1 expanded reference review

**Status: 50 candidate proteins are ready for human reference/outgroup review;
no expanded alignment or tree has been inferred.** This directory is a new,
independent review candidate. It does not overwrite or retroactively change the
executed 18-tip example in [`../brca1/`](../brca1/README.md).

The expansion asks a narrower and more defensible question than “add the top
BLAST hits”: can a taxonomically balanced set of amniote BRCA1 ortholog-group
members, plus multiple amphibian outgroup candidates, pass exact taxonomy,
length, local-similarity, and terminal-domain review before alignment?

## Review snapshot

| Item | Auditable result |
| --- | --- |
| Frozen provider scope | 558 NCBI Ortholog gene records and 3,605 protein records |
| Length screen | 2,543 proteins from 512 provider species fell within 0.75–1.10× the 1,863-aa human query length; this is not an orthology test |
| Proposed set | 50 unique accessions and TaxIDs: 1 study sequence, 44 expanded amniote references, and 5 amphibian outgroup candidates |
| Taxonomy | 50/50 exact, case-sensitive scientific-name and TaxID matches in one official NCBI `new_taxdump` snapshot |
| Terminal architecture | 50/50 promoted candidates passed the declared N-terminal RING + two distinct C-terminal BRCT interval gate |
| Pre-plan exclusion | `XP_075436319.1` (*Ascaphus truei*) was rejected for having only one distinct BRCT interval and replaced by `XP_053555561.1` (*Bombina bombina*) |
| Clustering | Not triggered: 44 expanded references are below the declared 200-sequence MMseqs2 threshold |
| Deterministic plan | Run `gtr-9eaf4f8541020b3b`; plan hash `0fb0a32e613304af10e696c9ebd40f0c2fc0d6b2c64965945748446385186bc6`; state `pending-reference-approval` |
| Inference | Not started; MAFFT, trimAl, FastTree, IQ-TREE, rooting, and the expanded figure remain behind review gates |

“Passed” above means the executable rule was satisfied. It does not mean every
candidate is already approved, complete, canonical, compositionally adequate,
or safe to use as an outgroup.

## Sampling design

The proposed set contains 45 amniote ingroup tips including the human query,
plus five amphibian outgroup candidates.

| Broad sample | Tips | Constituent sampling labels |
| --- | ---: | --- |
| Mammalia | 25 | Primates 5; Glires 5; Laurasiatheria 7; Afrotheria 3; Xenarthra 2; Marsupialia 3 |
| Sauropsida | 20 | Aves 9; Crocodylia 3; Testudines 3; Lepidosauria 5 |
| Amphibia candidate outgroups | 5 | Anura 2; Caudata 2; Gymnophiona 1 |

The exact accession, TaxID, species, GeneID, role, sampling label, and rationale
are in [`selection_manifest.tsv`](inputs/selection_manifest.tsv). Selection is
one record per TaxID and is intentionally balanced; it is not a claim that the
other provider records failed QC. The archive-wide species inventory is
[`provider_species_inventory.tsv`](inputs/provider_species_inventory.tsv).

Three desired sampling gaps remain explicit:

- the two frozen platypus models were 1,346 aa (0.722× the human length) and
  fell below the declared gate;
- the frozen echidna model was 1,370 aa (0.735×) and also fell below it;
- the frozen provider package contained no *Sphenodon punctatus* BRCA1 record.

See [`sampling_gaps.tsv`](inputs/sampling_gaps.tsv). No terminal-domain failure
is inferred for records that were never promoted to executable domain QC.

## What the QC does—and does not prove

The gds2 screen used BLASTP 2.17.0 and InterProScan 5.72-103.0 with Pfam 37.1
and SMART 9.0. [`candidate_qc.tsv`](qc/candidate_qc.tsv) reports union query and
subject coverage across local HSP intervals, plus the raw HSP-weighted identity
and terminal-domain evidence. Overlapping HSP lengths can be counted more than
once in `alignment_length`; it is not a global unique-alignment length.

Human-query BLAST coverage is deliberately treated as a **local-similarity
diagnostic**, not as proof of sequence completeness or orthology. Deeply
divergent BRCA1 proteins can retain terminal domains while the long variable
centre aligns poorly to human. Five records with query coverage below 0.50 are
therefore listed in [`manual_review_flags.tsv`](qc/manual_review_flags.tsv) and
remain pending manual approval. In particular:

- *Crocodylus porosus* has only 0.113/0.114 query/target union coverage;
- *Bombina bombina* has 0.433/0.461 coverage, a RING interval beginning at
  residue 1 with one signature, and a second BRCT interval supported only by
  SMART at E-value 0.36;
- the amphibian outgroups still require branch-length, leave-one-out,
  outgroup-removal, and conserved-domain sensitivity checks before any root is
  accepted.

The `is_reviewed` flag is a conservative RefSeq accession-class rule
(`NP_` = true, `XP_` = false), not a UniProt-style manual-review claim.
`is_canonical=true` is reserved for the chosen human focal isoform and says
nothing about canonical isoforms in other species. `is_fragment=false` means
the exported header did not explicitly say fragment/partial and the sequence
passed the length gate; it is not independent proof of model completeness.

## Rejection evidence

The first 50-candidate attempt used *Ascaphus truei*. The exact historical
manifest, FASTA, pre-QC table, raw BLAST/InterProScan outputs, software receipt,
and hashes are retained in [`qc/rejected-attempt/`](qc/rejected-attempt/README.md).
This makes the one-BRCT exclusion auditable instead of hiding it behind a prose
note. The promoted set contains *Bombina bombina* instead; it is still a
candidate, not an approved outgroup.

## Frozen NCBI taxonomy evidence

All 50 names were matched by exact `name_txt` equality to `scientific name` in
`names.dmp`; `nodes.dmp` verified the requested node and supplied its direct
parent and rank. No synonym, case-folded, substring, or fuzzy fallback was
used. The machine-readable snapshot receipt is
[`taxonomy_snapshot.json`](annotation/taxonomy_snapshot.json), the exact-match
table is [`taxonomy_resolution.tsv`](annotation/taxonomy_resolution.tsv), and
the complete root-to-tip lineages are published separately as
`annotation/taxonomy_lineage.tsv`.

The archive SHA-256 is
`e67d56c3e87eac14a28feea8cf710ac612571f00a4f32ba3b85bae0020cd123a`;
the extracted `names.dmp` and `nodes.dmp` hashes are embedded in every
resolution row. The large taxdump files themselves are intentionally not
committed.

## Files to inspect before approval

1. [`selection_manifest.tsv`](inputs/selection_manifest.tsv) — the proposed
   one-per-TaxID sampling and role decisions.
2. [`candidates.tsv`](inputs/candidates.tsv) and
   [`candidate_provenance.tsv`](inputs/candidate_provenance.tsv) — promoted
   evidence, sequence hashes, and provider annotation context.
3. [`candidate_qc.tsv`](qc/candidate_qc.tsv) and
   [`manual_review_flags.tsv`](qc/manual_review_flags.tsv) — executable results
   and the five candidates still requiring judgement.
4. [`rejected_candidates.tsv`](inputs/rejected_candidates.tsv) and the
   [raw rejected-attempt receipt](qc/rejected-attempt/README.md).
5. [`taxonomy_resolution.tsv`](review/taxonomy_resolution.tsv) — exact NCBI
   name–TaxID evidence used by the planner.
6. [`selected_references.tsv`](review/selected_references.tsv),
   [`plan.json`](review/plan.json), and [`manifest.json`](review/manifest.json)
   — the deterministic 50-tip handoff and its hashes.

The planner selected all 50 promoted records and recorded zero blockers, but
that is not human approval. There is deliberately no approval receipt in this
directory.

## Planned work after approval

If the reference and outgroup choices are accepted, the next run should:

1. align all 50 proteins with MAFFT E-INS-i and retain the raw MSA;
2. compare trimAl gap-threshold profiles 0.1, 0.5, and 0.9, audit RING/BRCT
   retention, and inspect flagged sequences in the full MSA;
3. add conserved-block or RING+BRCT-only sensitivity analyses because BRCA1's
   variable centre is difficult to align across deep vertebrate splits;
4. screen trim-profile topologies with FastTree, then run IQ-TREE2 ModelFinder
   with 1,000 UFBoot2 and 1,000 SH-aLRT replicates only after alignment approval;
5. retain the unrooted ML tree and test amphibian combinations, leave-one-tip
   outgroups, branch lengths, and composition before producing a rooted copy;
6. render a detailed vector tree plus a simplified README overview. A 50-tip
   overview will use selected labels and taxonomy/role bands rather than
   shrinking every accession into an unreadable image.

## Reproducibility boundary

The public scripts are deterministic, but the frozen NCBI archive and extracted
taxdump are not committed. An authorized host materialized those inputs and ran
the executable QC and offline planner under PBS; private paths and scheduler
logs are excluded from GitHub. Representative commands are:

```bash
python3 scripts/prepare_brca1_expanded.py \
  --protein-fasta <frozen-package>/protein.faa \
  --data-report <frozen-package>/data_report.jsonl \
  --selection-manifest inputs/selection_manifest.tsv \
  --rejected-candidates inputs/rejected_candidates.tsv \
  --externally-verified-provider-archive-sha256 \
    30f1ef86dbe2989981bed005c837e7538226cae24b800aef44a32a589bfe04c0 \
  --out <new-staging-directory>

qsub -v BRCA1_RUN_DIR=<absolute-run-directory> ../brca1/scripts/run_brca1_qc.pbs
qsub -v BRCA1_RUN_DIR=<absolute-run-directory>,TAXDUMP_DIR=<frozen-taxdump> \
  scripts/run_taxonomy_resolution.pbs
qsub -v BRCA1_RUN_DIR=<absolute-run-directory> scripts/run_reference_plan.pbs
```

All commands refuse to overwrite their primary output. See
[`qc/checksums.sha256`](qc/checksums.sha256) for the promoted QC input/output
binding and [`report/checksums.sha256`](report/checksums.sha256) for the full
public reference-review bundle. Data-source terms and attribution follow the existing
[`DATA_LICENSES.md`](../brca1/DATA_LICENSES.md).
