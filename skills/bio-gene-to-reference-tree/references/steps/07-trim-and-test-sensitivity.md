# Step 7 — Trim and test sensitivity

## When to read

Read after raw-MSA QC when filtering is proposed, when published post-alignment processing must be reproduced, or when an explicit untrimmed decision must be recorded.

## Required inputs

- immutable `alignment.raw.faa` and its QC;
- reviewed conservation assessment and biological motifs/domains to preserve;
- proposed trimAl profiles and exact tool version;
- reporting provenance from any decision-bearing article.

## Procedure

Interpret `trimal -gt x` as the minimum non-gap occupancy retained per column. Thus `0.98` is extremely strict, `0.90` is strict for heterogeneous samples, and `0.10`/`0.05` are extremely permissive. Never choose a threshold solely because sequences were described as close, distant, or viral.

Treat values such as `0.98/0.95/0.90` and `0.10/0.05` as named sensitivity profiles, not biological truths. Preserve every output and its exact argument array, for example:

```text
["trimal", "-in", "alignment.raw.faa", "-out", "alignment.trimmed.balanced.faa", "-gt", "0.9"]
```

For each profile, report input/retained columns, retained fraction, occupancy distribution, conserved-motif/domain retention, and affected sequences. Compare key topology with fast profile trees when reasonable profiles differ materially. Permissive trimming cannot repair mixed domains or replace viral recombination analysis.

Preserve publication reporting states separately: `exact`, `partial`, `explicit-none`, and `not-reported`. “Not reported” never means no trimming. An exact reported argv does not prove the paper interpreted a flag correctly; expose any prose/argv/version conflict and keep reproduction versus correction as separate reviewed choices.

Select one primary alignment explicitly: raw, a named trim profile, or another justified output. Bind approval to the exact alignment and conservation-assessment hashes.

## Required outputs

- all trim profile FASTAs, versioned argv arrays, logs, and hashes;
- retained-column metrics, motif checks, and topology-sensitivity summary;
- reporting-status evidence and any semantics conflict;
- selected primary alignment, rationale, hash, and approval record.

## Review gate and stop conditions

Pass the alignment/trimming gate before tree inference. Stop if filtering deletes key homologous regions, leaves too little signal, hides domain incompatibility, or changes the biological conclusion across reasonable profiles. Any changed MSA, threshold, profile, or conservation assessment invalidates approval.

## Supporting references

- [Recent MSA/trimming evidence and reporting states](../recent-msa-trimming-evidence.md)
- [Workflow Gate 3](../workflow.md#gate-3-alignment-and-trimming-approval)
- [Plan, approval, and checksum contract](../output-contract.md)

Primary software documentation: trimAl guide <https://vicfero.github.io/trimal/whatcanido.html>.
