# Workflow states and decision gates

Use a staged workflow so that automated acquisition cannot silently become an approved biological analysis.

## Contents

- [State model](#state-model)
- [Environment preflight](#environment-preflight)
- [Sequence-type preflight](#sequence-type-preflight)
- [Query resolution](#gate-1-query-resolution)
- [Reference and outgroup approval](#gate-2-reference-and-outgroup-approval)
- [Alignment and trimming approval](#gate-3-alignment-and-trimming-approval)
- [Tree inference](#tree-inference-gate)
- [Annotation and evidence](#annotation-and-evidence-gate)
- [Completion contract](#completion-contract)
- [Bundled helper behavior](#bundled-helper-behavior)
- [Out-of-scope routing](#out-of-scope-routing)

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

Any state may enter `blocked` or `failed`. Invalidate prior approval when a query sequence, molecule declaration, analysis space, comparable-region definition, orientation, RNA normalization, CDS frame/code/translation, candidate record, search database, threshold, cluster mapping, outgroup, MSA/backtranslation, trim profile, model, support method, command, or decision-bearing color changes.

## Environment preflight

Environment feasibility is a sidecar to the scientific state model, not another biological approval state. Before choosing execution mode, create an explicit profile, record every input's host/compute location, run `doctor` on the actual compute target, add a separate host snapshot whenever the route depends on host-local planning, search, rooting, annotation, or rendering executables, and compile a route using [environment-routing.md](environment-routing.md). Re-run after any profile or used snapshot change; declared host API/tree-I/O capabilities remain explicit profile assertions rather than inferred binaries.

Bind the decision to the profile and all used snapshots with `route_hash`, but do not include machine state in `plan_hash`. A newly installed tool cannot approve a changed reference set or MSA, and an approved scientific plan cannot prove that a machine is capable or authorized to execute it.

## Sequence-type preflight

Before resolving any request, read [sequence-type-routing.md](sequence-type-routing.md). Requests 0.1 and 0.2 retain their protein-only behavior. Request 0.3 requires an exact molecule declaration and routes to protein, direct noncoding nucleotide, or clean CDS with protein-guided codon-preserving alignment followed by an explicit nucleotide-site or codon model. Never infer molecule type from letters or silently cross analysis spaces.

Require comparable homologous regions. Declare RNA `source_encoding` as `rna-u` or `dna-t`, preserve that source, and hash a separately named DNA-alphabet analysis copy (U→T for `rna-u`, byte-equivalent for `dna-t`). For CDS, require frame, strand, genetic-code, completeness, translation, and stop/frameshift QC. The deterministic trimAl backtranslation route accepts only NCBI codes 1/11; other codes, disrupted CDS, and uncertain CDS stop for reviewed MACSE/PAL2NAL/precomputed handling and do not continue as ordinary nucleotide input. A MACSE review can continue only after an explicit export policy and rerouting with a typed precomputed codon alignment.

## Gate 1: query resolution

Require a local molecule-matched record with a stable ID, sequence, organism, region definition, and provenance. For accession/name routes, retain the original input and resolution evidence. Stop on ambiguity, unresolved isoforms/features, missing molecule declaration, non-comparable regions, failed RNA source-encoding/analysis-copy provenance, invalid CDS translation, or missing organism/TaxID for a name.

When NCBI taxdump validation is enabled, validate every candidate `species`/`taxon_id` pair against already-extracted `names.dmp` and `nodes.dmp` from the same recorded snapshot before selection. Accept only a unique character-for-character `scientific name` match and an exact TaxID agreement. Stop on aliases, fuzzy or normalized matches, ambiguity, missing nodes, or mixed/unrecorded snapshots. Bind approval to the dump hashes and `taxonomy_resolution.tsv`.

For unpublished material, stop until the user approves each remote submission class. Permission to query one database does not automatically authorize iTOL upload or another external service.

## Gate 2: reference and outgroup approval

Present:

- acquisition tier and database provenance;
- declared molecule/analysis space, region compatibility, the three request-0.3 search fields, and any separately established host-side index evidence;
- RNA normalization or CDS translation/backtranslation evidence when applicable;
- exact-name taxonomy evidence and taxdump hashes when enabled;
- counts before/after every filter and cluster;
- retained taxa and unsampled clades;
- `study`, `expanded`, and `outgroup` analysis groups, separate from biological ingroup/outgroup roles;
- one-to-many, paralog, fragment, fusion, domain, and low-complexity warnings;
- selected and rejected accessions with stable reason codes;
- each proposed outgroup and taxonomic rationale;
- planned MMseqs2, MAFFT, trimAl, and tree command arrays.

Tie approval to the current `plan_hash`. If MMseqs2 is required, execute it only on expanded candidates, import cluster membership, re-plan, and request approval on the new hash.

## Gate 3: alignment and trimming approval

Before choosing the primary alignment, review `evidence/conservation_assessment.tsv`: one focal analysis unit, one controlled class, an explicit comparison scope, an evidence basis, stable evidence IDs, and limitations. Accept or revise the provisional row, mark it `reviewed`, and bind its SHA-256 to the alignment approval. A changed assessment reopens this gate.

After MAFFT, present the explicit `--amino` or `--nuc` route, raw-MSA length, per-tip gap/coverage statistics, column occupancy, conserved motif/region checks, unusual insertions, excluded sequences, and any domain conflict. For a codon route, also present exact CDS↔translation IDs, genetic code, raw protein MSA, trimAl `-backtrans` argv, triplet integrity, and translation-equality QC. Never silently remove a sequence.

When trimming is enabled, present every trimAl profile, threshold semantics, retained length/fraction, removed-column record, motif retention, and topology sensitivity if fast profile trees were compared. Require an explicit choice of the primary alignment before IQ-TREE2.

Stop when the MSA has mixed molecule types/domains, non-comparable regions, duplicate tip IDs, severe coverage failure, unresolvable homology, invalid codon/backtranslation QC, an unreviewed MACSE requirement, fewer than four usable taxa/sequences, or a key conclusion that is unstable across reasonable MSA/trimming decisions.

## Tree-inference gate

Verify the approved MSA hash, exact executable/version, resource limits, thread count, seed, model plan, support method, output directory, and current plan hash. Distinguish:

- FastTree approximate ML with SH-like local support;
- IQ-TREE2 UFBoot2 using `-B` and `-bnni`;
- standard bootstrap using `-b`;
- SH-aLRT using `-alrt`.

Bind the inference command to analysis space: protein uses an explicit protein model; noncoding DNA/RNA and an explicitly requested nucleotide-site CDS analysis use FastTree `-nt -gtr` or IQ-TREE `-st DNA`; codon analysis uses IQ-TREE `-st CODON<n>` with the approved genetic code, including explicit `CODON1`. Both CDS analysis spaces retain the reviewed codon-preserving backtranslation, but only the latter is a codon substitution model. FastTree cannot replace a codon model, and a missing analysis-specific executable blocks rather than changing molecule space.

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
- exact molecule/analysis-space declaration, comparable-region evidence, and actual search database provenance;
- RNA source-encoding/analysis-copy receipt or CDS frame/code/translation/backtranslation QC when applicable;
- reviewed, hash-bound conservation assessment;
- unrooted tree and optional separately rooted tree;
- correctly named support measures and model;
- iTOL roles and full sequence metadata;
- optional NCBI taxonomy resolution evidence with snapshot and dump hashes;
- requested local ggtree SVG/PDF figures and their settings TSV;
- real literature/taxonomy evidence or an explicit evidence-search limitation;
- exact commands, versions, hashes, warnings, manual decisions, and approved plan hashes.

## Bundled helper behavior

The planner compiles only the local pre-execution review bundle. It may validate supplied local taxdump files but never downloads them. It uses `pending-clustering`, `blocked`, or `pending-reference-approval`; tree rendering and later states belong to separately invoked local or host-agent execution and must not be claimed by the planner.

The helper keeps request schemas 0.1 and 0.2 as legacy protein-only inputs and never mutates them. Use request schema 0.3 for new protein or nucleotide work; it emits the molecule-aware plan/output schema 0.4. A failed nucleotide or codon route never falls back to the legacy protein plan.

## Out-of-scope routing

Route species-tree inference, gene-tree/species-tree reconciliation, duplication/loss modeling, HGT analysis, divergence dating, positive-selection tests, recombination-aware inference, and publication-grade figure design to dedicated workflows. This skill may identify the need; it must not silently expand the analysis.
