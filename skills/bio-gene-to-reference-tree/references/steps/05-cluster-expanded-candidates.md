# Step 5 — Cluster expanded candidates

## When to read

Read only when the approved size trigger is crossed or the user explicitly requests redundancy reduction. Otherwise skip this step without loading or simulating it.

## Required inputs

- the full candidate pool and pre-clustering selection proposal;
- `analysis_group` values that distinguish `study`, `expanded`, and `outgroup`, kept separate from biological `role=ingroup|outgroup`;
- approved identity, coverage, coverage-mode, thread, and representative-ranking policies;
- a fresh output location and preserved pre-clustering plan hash.

## Procedure

Cluster only records with `analysis_group=expanded` in `expanded_candidates.faa`. Never place `analysis_group=study` or `analysis_group=outgroup` records into the clustering input, regardless of biological `role`. Specify both identity and coverage; `-c 0.7` means coverage, not 70% identity.

Use this full-length starting profile unless the objective justifies another value:

```text
mmseqs easy-linclust expanded_candidates.faa clusters mmseqs_tmp \
  --min-seq-id 0.95 -c 0.8 --cov-mode 0 --threads <fixed>
```

`--cov-mode 0` measures coverage against the longer sequence and is a reasonable full-length default. Use target-oriented coverage only for a documented fragment use case. Record the exact MMseqs2 version and argv.

Preserve every representative/member mapping, accession, TaxID, biological `role`, `analysis_group`, and `cluster_id`. Within each cluster, choose the representative by the same reviewed/canonical/relationship/coverage ordering used for selection. Do not silently erase unique taxa, alleles, biologically meaningful copies, or rejected members.

Re-import cluster membership, apply `MMSEQS_CLUSTER_REDUNDANT` where appropriate, and run Step 4 again. A pre-clustering approval cannot authorize the changed set.

## Required outputs

- immutable MMseqs2 input FASTA, native cluster output, and representative/member mapping;
- version, exact argv, threads, exit status, and output hashes;
- candidate TSV with preserved `cluster_id` and representative fields;
- re-planned selected/rejected tables, reference FASTA, and new plan hash.

## Review gate and stop conditions

Return to the reference/outgroup approval gate after clustering. Stop if clustering touches a protected `analysis_group`, removes the only member of a required taxon, cannot reproduce membership, or lacks an explicit identity/coverage/cov-mode triplet.

## Supporting references

- [Tool and executable boundaries](../tool-routing.md)
- [Workflow Gate 2 and approval invalidation](../workflow.md#gate-2-reference-and-outgroup-approval)
- [Plan and provenance contract](../output-contract.md)

Primary software documentation: MMseqs2 user guide <https://github.com/soedinglab/MMseqs2/wiki>.
