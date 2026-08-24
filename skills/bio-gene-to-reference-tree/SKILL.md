---
name: bio-gene-to-reference-tree
description: Build an auditable protein or nucleotide gene tree from an accession, explicitly declared local sequence, or gene, protein, or feature name plus source organism. Use when an agent must choose an environment-appropriate route; resolve metadata and exact NCBI TaxIDs; classify gene-family conservation at an explicit taxonomic scope; curate references and outgroups; align proteins, noncoding DNA or RNA, or clean CDS through a protein-guided codon-preserving route; run FastTree or IQ-TREE2; generate iTOL or ggtree/ggplot2 outputs; or compare a gene tree with current literature. Require review gates before reference selection, alignment choice, and tree inference.
---

# Gene to Reference Tree

Build a single-gene protein or nucleotide tree through explicit, reviewable decisions. Treat molecule/region definition, acquisition, sampling, alignment, rooting, and literature comparison as scientific analyses rather than mechanical top-hit processing.

## Preserve the core claim

- Call the result a **gene tree**, not a species tree.
- Require an exact molecule and analysis-space declaration; never infer it from sequence letters or cross protein, nucleotide, and codon spaces as a fallback.
- Treat similarity as evidence of homology, never proof of orthology or identical function.
- Treat support as repeatability under a method, not proof that a topology is correct.
- Preserve the raw candidate pool, untrimmed alignment, unrooted tree, commands, versions, database snapshots, decisions, and exclusions.
- Never upload an unpublished sequence, candidate set, tree, or metadata without explicit permission for that remote action.
- Never fabricate a lookup, sequence, TaxID, orthology call, citation, version, release, or executed result.

## Load progressively

1. Read [workflow.md](references/workflow.md) once at the start to establish states, approval gates, invalidation rules, and completion criteria.
2. Determine the current state and load **only the matching step file** from the table below. Do not preload all ten steps.
3. Treat supporting links as discovery pointers, not automatic loads. Open a shared reference only when the current procedure requires it or the explicit trigger below applies.
4. Produce the step's required outputs and satisfy its stop/review conditions before advancing.
5. On resume, verify hashes and approvals, then continue from the earliest invalid or incomplete step.

## Route the current task

| Step | Load when | Task file | Exit condition |
|---:|---|---|---|
| 1 | Intake is an accession, raw protein/nucleotide/CDS, or name plus organism | [Resolve the query](references/steps/01-resolve-query.md) | One stable molecule-matched local record and provenance pass Gate 1 |
| 2 | Query identity is resolved but the biological question/scope is not fixed | [Define the objective](references/steps/02-define-objective.md) | Objective, ingroup, relationship policy, and deliverables are explicit |
| 3 | Objective is fixed and a broad candidate pool is needed | [Discover candidates](references/steps/03-discover-candidates.md) | Candidate FASTA/metadata and acquisition provenance are materialized |
| 4 | Candidate pool exists and tree tips/outgroups must be proposed | [Select references and outgroups](references/steps/04-select-references-and-outgroups.md) | A hash-bound proposal reaches clustering or Gate 2 review |
| 5 | Expanded candidates cross the declared clustering trigger | [Cluster expanded candidates](references/steps/05-cluster-expanded-candidates.md) | Cluster mapping is imported and Step 4 is rerun; otherwise skip |
| 6 | Reference/outgroup set is approved | [Align and assess conservation](references/steps/06-align-and-assess-conservation.md) | Route-specific raw MSA, conservation assessment, and QC are reviewable |
| 7 | Trimming/untrimmed sensitivity and primary MSA choice are required | [Trim and test sensitivity](references/steps/07-trim-and-test-sensitivity.md) | One exact alignment hash passes Gate 3 |
| 8 | Alignment and inference plan are approved | [Infer, root, and check the tree](references/steps/08-infer-root-and-check-tree.md) | Unrooted tree and any separately approved rooted copy are validated |
| 9 | Tree tips require iTOL, metadata, or local ggtree/ggplot2 output | [Annotate and visualize](references/steps/09-annotate-and-visualize.md) | Tip sets, semantics, annotations, and requested figures agree exactly |
| 10 | Current evolutionary evidence and final reporting are required | [Compare evidence and report](references/steps/10-compare-evidence-and-report.md) | Completion contract is satisfied with evidence or explicit limitations |

Step 5 is conditional. Never cluster records whose `analysis_group` is `study` or `outgroup`, and never advance on a pre-clustering approval after membership changes.

## Apply shared rules only when triggered

- Read [sequence-type-routing.md](references/sequence-type-routing.md) before resolving request 0.3 or any potentially nucleotide input. Requests 0.1/0.2 remain protein-only. Require comparable regions, explicit RNA source encoding (`rna-u` or `dna-t`) with a separately named analysis copy, and CDS frame/code/translation QC. Deterministic trimAl backtranslation is limited to NCBI codes 1/11; other codes and disrupted CDS require a reviewed handoff.
- Read [taxonomy-resolution.md](references/taxonomy-resolution.md) before deriving or validating a TaxID from an organism name. Use one verified NCBI taxdump snapshot and accept only one character-for-character `scientific name` match whose node exists.
- Read [environment-routing.md](references/environment-routing.md) before choosing an execution mode or handing work to HPC/SSH. Snapshot the actual compute target, route every required step, and keep `route_hash` independent from the scientific `plan_hash`.
- Read [tool-routing.md](references/tool-routing.md) before the first live lookup, unpublished-data submission, external executable, or search-database choice.
- Read [recent-msa-trimming-evidence.md](references/recent-msa-trimming-evidence.md) when conservation classification or recent MSA/trimming precedent affects a decision; filter its [TSV catalog](references/recent-msa-trimming-evidence.tsv) by data architecture and scope.
- Read [ggtree-visualization.md](references/ggtree-visualization.md) before local publication-oriented rendering.
- Read [output-contract.md](references/output-contract.md) before creating a request, running the planner, approving hashes, or assembling an executed report.

## Compile the deterministic review bundle

After an authorized host agent has materialized a resolved molecule-matched query and candidate TSV/FASTA bundle plus any required CDS translation or RNA provenance, locate this file, treat its directory as `<skill-root>`, and run:

```text
python3 <skill-root>/scripts/gene_to_tree.py plan \
  --request <request.json> --offline --dry-run --out <new-output-directory>
```

Review the selected/rejected tables, reference FASTA, sequence metadata, iTOL roles, plan, manifest, hashes, warnings, and planned argv arrays. The helper performs no network request or external bioinformatics execution in `plan` mode, refuses overwrite, and invalidates approval after decision-bearing changes.

Inspect optional local executables with:

```text
python3 <skill-root>/scripts/gene_to_tree.py doctor --json
```

This default inspection only resolves command names on `PATH`; it launches no discovered executable. Use `--run-version-probes` only after approving local version/help commands and R namespace loading.

For an environment-matched route, materialize a profile and snapshot, then run:

```text
python3 <skill-root>/scripts/gene_to_tree.py route \
  --profile <environment-profile.json> \
  --snapshot <environment-snapshot.json> \
  --host-snapshot <optional-host-snapshot.json> \
  --out <route-decision.json>
```

Neither environment command executes a workflow stage. Default `doctor` is passive `PATH` discovery; its opt-in active mode invokes only local version/help commands and R namespace checks, without claiming those programs are side-effect-free. A route decision is evidence of capability fit, not authorization. The host supplies separately authorized database, literature, browser, shell, file-transfer, and scheduler capabilities. The bundled helper validates materialized handoffs and compiles deterministic plans; it does not download databases, search literature, run workflow stages, connect by SSH, submit a scheduler job, or upload to iTOL during `plan` or `route`.

## Route unsupported analyses

Route species-tree inference, gene-tree/species-tree reconciliation, duplication/loss modeling, HGT analysis, divergence dating, positive selection, recombination-aware inference, covariance-model RNA phylogenetics beyond the optional dedicated `mafft-qinsi` route, automated disrupted/frameshift CDS analysis, genome-scale phylogeny, and publication figure design beyond the bundled renderer to dedicated workflows. The clean-CDS route does not silently promote MACSE review output. Identify the need without expanding the claim.
