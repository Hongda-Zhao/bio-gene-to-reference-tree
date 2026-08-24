# Step 7 — Trim and test sensitivity

## When to read

Read after raw-MSA QC when filtering is proposed, when published post-alignment processing must be reproduced, or when an explicit untrimmed decision must be recorded.

## Required inputs

- immutable route-specific raw alignment(s), exact molecule/analysis space, and QC;
- for codon analysis, the approved CDS FASTA, exact translation mapping, genetic code, and raw protein/codon alignments;
- reviewed conservation assessment and biological motifs/domains to preserve;
- proposed trimAl profiles and exact tool version;
- reporting provenance from any decision-bearing article.

## Procedure

Interpret `trimal -gt x` as the minimum non-gap occupancy retained per column. Thus `0.98` is extremely strict, `0.90` is strict for heterogeneous samples, and `0.10`/`0.05` are extremely permissive. Never choose a threshold solely because sequences were described as close, distant, or viral.

Treat values such as `0.98/0.95/0.90` and `0.10/0.05` as named sensitivity profiles, not biological truths. Preserve every output and its exact argument array. Protein and direct-nucleotide routes trim their own raw MSA, for example:

```text
["trimal", "-in", "alignment.raw.faa", "-out", "alignment.trimmed.balanced.faa", "-gt", "0.9"]
```

For a clean NCBI-code-1/11 codon route, apply the reviewed profile to the protein MSA and project exactly those columns through `-backtrans`; also save the trimmed protein counterpart. Add `-ignorestopcodon` only after the independent input gate proves zero internal stops and exact translation equality, allowing a documented terminal stop without weakening post-backtranslation triplet/translation QC. Other codes require a reviewed MACSE/PAL2NAL/precomputed handoff because current trimAl stop handling uses the universal stop set.

```text
["trimal", "-in", "alignment.raw.translated.faa", "-out", "alignment.trimmed.balanced.faa", "-gt", "0.9", "-fasta"]
["trimal", "-in", "alignment.raw.translated.faa", "-backtrans", "reference_set.fna", "-ignorestopcodon", "-out", "alignment.trimmed.balanced.codon.fna", "-gt", "0.9", "-fasta"]
```

Require unique identical first-token IDs in protein and CDS FASTAs. Treat trimAl truncation, padding, unmatched-ID, duplicate-ID, triplet, or translation warnings as failures. Verify output length divisible by three, complete triplets, and translation equality. A disrupted/frameshift/internal-stop-containing CDS needs reviewed MACSE handling; MACSE output is not automatically tree-ready. Review frameshift/stop handling, apply an explicit export policy, verify IDs/triplets, and rerun routing with a typed precomputed codon alignment. Never trim it as ordinary nucleotide.

For each profile, report input/retained columns, retained fraction, occupancy distribution, conserved-motif/domain/region retention, and affected sequences. Codon reports include retained amino-acid and nucleotide columns plus phase/translation checks. Compare key topology with molecule-compatible profile trees when reasonable profiles differ materially; FastTree is not a codon-model screen. Permissive trimming cannot repair mixed domains/regions, CDS disruption, or replace viral recombination analysis.

Preserve publication reporting states separately: `exact`, `partial`, `explicit-none`, and `not-reported`. “Not reported” never means no trimming. An exact reported argv does not prove the paper interpreted a flag correctly; expose any prose/argv/version conflict and keep reproduction versus correction as separate reviewed choices.

Select one primary alignment explicitly: raw, a named trim profile, or another justified output. Bind approval to exact molecule/analysis space, alignment/backtranslation, CDS/translation or RNA-derivative, and conservation-assessment hashes.

## Required outputs

- all route-specific trim profile FASTAs, backtranslations where applicable, versioned argv arrays, logs, and hashes;
- retained-column metrics, motif checks, and topology-sensitivity summary;
- reporting-status evidence and any semantics conflict;
- selected primary alignment, rationale, hash, and approval record.

## Review gate and stop conditions

Pass the alignment/trimming gate before tree inference. Stop if filtering deletes key homologous regions, leaves too little signal, hides domain/region incompatibility, breaks triplets/translation, or changes the biological conclusion across reasonable profiles. Any changed molecule space, MSA/backtranslation, threshold, profile, CDS/RNA mapping, or conservation assessment invalidates approval.

## Supporting references

- [Recent MSA/trimming evidence and reporting states](../recent-msa-trimming-evidence.md)
- [Sequence-type, codon, and MACSE boundaries](../sequence-type-routing.md)
- [Workflow Gate 3](../workflow.md#gate-3-alignment-and-trimming-approval)
- [Plan, approval, and checksum contract](../output-contract.md)

Primary software documentation: trimAl usage, including `-backtrans`, <https://trimal.readthedocs.io/en/latest/usage.html>.
