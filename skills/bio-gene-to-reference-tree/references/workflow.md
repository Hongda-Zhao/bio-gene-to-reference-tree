# Workflow states and decision gates

Use a staged workflow so that automated acquisition cannot silently become an approved biological analysis.

## State model

```text
intake
  -> query-resolved
  -> candidates-materialized
  -> pending-clustering       # only when the size trigger is met
  -> pending-reference-approval
  -> references-approved
  -> alignment-qc
  -> pending-alignment-approval
  -> alignment-approved
  -> tree-inference
  -> annotation-and-context
  -> complete
```

Any state may enter `blocked` or `failed`. Invalidate prior approval when a query sequence, candidate record, threshold, cluster mapping, outgroup, MSA, trim profile, model, support method, command, or decision-bearing color changes.

## Gate 1: query resolution

Require a local resolved protein record with a stable ID, sequence, organism, and provenance. For accession/name routes, retain the original input and resolution evidence. Stop on ambiguity, unresolved isoforms, non-protein input, invalid CDS translation, or missing organism/TaxID for a name.

When NCBI taxdump validation is enabled, validate every candidate `species`/`taxon_id` pair against already-extracted `names.dmp` and `nodes.dmp` from the same recorded snapshot before selection. Accept only a unique character-for-character `scientific name` match and an exact TaxID agreement. Stop on aliases, fuzzy or normalized matches, ambiguity, missing nodes, or mixed/unrecorded snapshots. Bind approval to the dump hashes and `taxonomy_resolution.tsv`.

For unpublished material, stop until the user approves each remote submission class. Permission to query one database does not automatically authorize iTOL upload or another external service.

## Gate 2: reference and outgroup approval

Present:

- acquisition tier and database provenance;
- exact-name taxonomy evidence and taxdump hashes when enabled;
- counts before/after every filter and cluster;
- retained taxa and unsampled clades;
- study, expanded, and outgroup roles;
- one-to-many, paralog, fragment, fusion, domain, and low-complexity warnings;
- selected and rejected accessions with stable reason codes;
- each proposed outgroup and taxonomic rationale;
- planned MMseqs2, MAFFT, trimAl, and tree command arrays.

Tie approval to the current `plan_hash`. If MMseqs2 is required, execute it only on expanded candidates, import cluster membership, re-plan, and request approval on the new hash.

## Gate 3: alignment and trimming approval

After MAFFT, present raw-MSA length, per-tip gap/coverage statistics, column occupancy, conserved motif checks, unusual insertions, excluded sequences, and any domain conflict. Never silently remove a sequence.

When trimming is enabled, present every trimAl profile, threshold semantics, retained length/fraction, removed-column record, motif retention, and topology sensitivity if fast profile trees were compared. Require an explicit choice of the primary alignment before IQ-TREE2.

Stop when the MSA has mixed molecule types/domains, duplicate tip IDs, severe coverage failure, unresolvable homology, fewer than four usable taxa/sequences, or a key conclusion that is unstable across reasonable MSA/trimming decisions.

## Tree-inference gate

Verify the approved MSA hash, exact executable/version, resource limits, thread count, seed, model plan, support method, output directory, and current plan hash. Distinguish:

- FastTree approximate ML with SH-like local support;
- IQ-TREE2 UFBoot2 using `-B` and `-bnni`;
- standard bootstrap using `-b`;
- SH-aLRT using `-alrt`.

Preserve the unrooted tree and all native logs. Create a rooted derivative only from approved outgroup tips. Re-open the decision gate if outgroup behavior, long-branch attraction, or model sensitivity makes the root unreliable.

When the rooted copy retains node support, map labels by canonical unrooted bipartition rather than by internal-node number, then verify that every original split retains its exact label once. Rerooting software may otherwise shift support labels along the reroot path.

## Annotation and evidence gate

Generate the iTOL color strip and metadata locally. Generate a range dataset only after checking contiguity/monophyly. Obtain separate permission before uploading unpublished material to iTOL.

When a local figure is requested, run the bundled ggtree/ggplot2 renderer only on an approved Newick tree and the corresponding metadata. Require exact equality between the tree tip set and selected `tip_id` values, declare root state, branch-length mode, and support format, and preserve SVG, PDF, and renderer settings TSV outputs. The renderer must never install packages, contact the network, reroot, ladderize, or guess support semantics.

Search current phylogenetic evidence, label direct versus broader-taxonomic sources, and compare it with the gene tree without forcing agreement. Record conflicts and plausible causes. For viral analyses, require recombination/reassortment/segment review before completion.

## Completion contract

Mark the workflow complete only when the final report includes:

- resolved query and acquisition provenance;
- selected/rejected references and cluster mapping;
- raw and approved MSA plus QC;
- unrooted tree and optional separately rooted tree;
- correctly named support measures and model;
- iTOL roles and full sequence metadata;
- optional NCBI taxonomy resolution evidence with snapshot and dump hashes;
- requested local ggtree SVG/PDF figures and their settings TSV;
- real literature/taxonomy evidence or an explicit evidence-search limitation;
- exact commands, versions, hashes, warnings, manual decisions, and approved plan hashes.

## Bundled helper behavior

The planner compiles only the local pre-execution review bundle. It may validate supplied local taxdump files but never downloads them. It uses `pending-clustering`, `blocked`, or `pending-reference-approval`; tree rendering and later states belong to separately invoked local or host-agent execution and must not be claimed by the planner.

The helper accepts request schema 0.1 for migration, emits a deprecation warning, and never mutates the source request. It emits plan/output schema 0.3; use request schema 0.2 for all new work.

## Out-of-scope routing

Route species-tree inference, gene-tree/species-tree reconciliation, duplication/loss modeling, HGT analysis, divergence dating, positive-selection tests, recombination-aware inference, and publication-grade figure design to dedicated workflows. This skill may identify the need; it must not silently expand the analysis.
