# Output contract

Write a new output directory for every planning run. Refuse overwrite. Use UTF-8, LF endings, real tab delimiters, stable ordering, relative artifact paths, and SHA-256 hashes.

## Contents

- [Request contract](#request-contract)
- [Environment preflight artifacts](#environment-preflight-artifacts)
- [Deterministic planning artifacts](#deterministic-planning-artifacts)
- [Final host-agent execution artifacts](#final-host-agent-execution-artifacts)
- [Literature evidence schema](#literature-evidence-schema)
- [Failure semantics](#failure-semantics)

## Request contract

Use `schema_version: "0.3"` for new requests and validate the structure against [request-0.3.schema.json](request-0.3.schema.json). Requests 0.1 and 0.2 remain supported as legacy protein-only contracts; the helper must not infer or migrate them into nucleotide analysis.

Every live accession or name route must be materialized before planning as:

- one resolved local FASTA matching the declared molecule type;
- one molecule-matched candidate FASTA containing the query/self record;
- one candidate metadata TSV;
- for clean coding DNA/RNA, one exact-ID translation FASTA plus CDS/frame/code/translation QC; deterministic trimAl backtranslation accepts only NCBI codes 1/11;
- for RNA computation, explicit `molecule.source_encoding`: `rna-u` for an RNA-IUPAC source or `dna-t` for a T-encoded archival source, plus a separately named and hashed DNA-alphabet analysis copy;
- recorded provenance in the request and candidate table.

The helper never performs the live resolution itself.

Request 0.3 must declare `protein`, `noncoding-dna`, `noncoding-rna`, `coding-dna`, or `coding-rna` and the corresponding protein, nucleotide-site, or codon analysis space. Every candidate row must match it. Clean CDS may explicitly choose nucleotide-site or codon inference after the same codon-preserving protein-guided alignment; that choice remains decision-bearing. Apply the complete [sequence-type routing contract](sequence-type-routing.md); a failed declared route cannot fall back across analysis spaces.

`privacy.remote_search_allowed` is a coarse capability gate, not consent to submit an unpublished sequence, tree, or metadata. Record each such submission/upload permission separately in host-side decision provenance, default it to denied when absent, and never treat a planned bundle as authorization for a remote action.

## Environment preflight artifacts

Before external execution, preserve an environment profile, a `doctor` snapshot from the actual compute target, any separate snapshot used to prove host-local executables, and the compiled route decision described in [environment-routing.md](environment-routing.md). Validate them against `environment-profile-0.2.schema.json`, `environment-snapshot-0.1.schema.json`, and `route-decision-0.2.schema.json`. The profile records locations for the query, atomic candidate bundle, cluster mapping, existing alignment, tree, metadata, taxdump, and sequence database; its bundle kind covers a required translation FASTA. The scientific plan and manifest separately describe RNA derivatives and alignment/backtranslation outputs. A reachable audited cluster mapping can satisfy the Step 5 software requirement, but it must be imported and Step 4 rerun before reference approval. Keep declared host API/tree-I/O capabilities explicit in the profile.

The route decision hashes the normalized profile plus compute/host snapshots into an independent `route_hash`. Do not add machine state to the scientific `plan_hash`. Do not store hostnames, usernames, private paths, raw command output, module commands, credentials, or scheduler job IDs in the portable snapshot. `execution_authorized` remains false because capability selection cannot authorize network, transfer, SSH, scheduler, upload, or external-tool execution.

## Deterministic planning artifacts

A successful plan emits [plan/output schema 0.4](plan-0.4.schema.json). Legacy requests retain protein behavior; request 0.3 uses molecule-specific artifact names and produces:

```text
selected_references.tsv
rejected_references.tsv
reference_set.faa            # protein
reference_set.fna            # nucleotide or CDS; mutually exclusive with .faa
reference_set.rna.fasta      # preserved RNA source encoding, RNA routes only
reference_set.translated.faa # verified translations, coding routes only
sequence_metadata.tsv
taxonomy_resolution.tsv     # when local NCBI taxdump validation is enabled
itol_roles.txt              # when iTOL annotation is enabled
plan.json
manifest.json
```

When automatic clustering is triggered but no cluster mapping is supplied, also produce `expanded_candidates.faa` for protein or `expanded_candidates.fna` for nucleotide/CDS, set state `pending-clustering`, and require re-planning after molecule-matched MMseqs2 clustering.

Do not write an empty `literature_evidence.tsv`; create it only after a real evidence search.

When taxonomy validation is enabled, treat `names.dmp` and `nodes.dmp` as hashed inputs from one recorded NCBI snapshot. The emitted `taxonomy_resolution.tsv` must contain one row per candidate and retain record ID, raw input name, requested TaxID, matched scientific name, resolved TaxID, parent TaxID, rank, snapshot provenance, and both dump hashes. Do not emit a partial successful table after any unresolved, ambiguous, alias-only, missing-node, or TaxID-mismatch result.

### `selected_references.tsv`

Include the query exactly once plus every retained ingroup and outgroup. Preserve candidate metadata and add `sequence_sha256`, `selection_order`, and `decision_reason`. Sort query first, selected ingroup deterministically, and outgroups last. Use accession, not a display label, as the unique key.

### `rejected_references.tsv`

Include every non-selected candidate exactly once. Preserve identifying/provenance columns and add deterministic semicolon-delimited `reason_codes`.

### `reference_set.faa` or `reference_set.fna`

Write exactly one ungapped sequence per selected row in the route-matched reference set: `.faa` for protein and `.fna` for DNA or the RNA analysis copy. Use accession as the first FASTA token and require unique tool-safe identifiers. Keep descriptive metadata in TSV. Verify molecule declaration, identifier-set equality, and sequence length. For DNA/protein, verify SHA-256 against the selected table. For RNA, preserve the declared source encoding in `reference_set.rna.fasta`. With `source_encoding: rna-u`, create `reference_set.fna` by U→T normalization and record both hashes; with `source_encoding: dna-t`, create a separately named, byte-equivalent `.fna` analysis copy and record both hashes. Never overwrite the source or relabel its biological molecule. For coding input, write the exact-ID validated translations to `reference_set.translated.faa`; do not replace the selected CDS FASTA with translations.

### `sequence_metadata.tsv`

Include every candidate exactly once, selected or rejected. Add:

```text
tip_id, analysis_role, inclusion_status,
selection_order, reason_codes
```

Then preserve all base and optional candidate metadata. For request 0.3, each row must include non-empty `actual_search_database`, `search_program`, and `search_database_molecule`; the helper validates their declared combination but does not execute or independently verify the search. Include molecule type, feature/region, strand/orientation, and all applicable CDS/frame/code/translation fields. The selected tip set must equal the molecule-matched reference FASTA identifiers.

### `itol_roles.txt`

Use the official `DATASET_COLORSTRIP` format and include selected tips only. Default colors are study orange `#E69F00`, expanded green `#009E73`, and outgroup gray `#999999`. Do not generate a range dataset before topology review establishes a meaningful contiguous group.

### `plan.json`

Make this the reviewable, decision-bearing artifact. Include:

```text
schema_version, request_schema_version, project_id, run_id, state,
objective, sequence_context, molecule_plan, query, query_resolution,
taxon_scope, taxonomy_plan, privacy, reference_discovery, selection_parameters,
clustering_plan, alignment_plan, trimming_plan, tree_plan,
annotation_plan, literature_plan, selection_summary,
selected_accessions, rejected_accessions_and_reasons,
decision_gates, artifact_plan, planned_commands, rooting_plan,
warnings, hard_stops, candidate_semantic_hash, plan_hash, approval
```

Use `pending-clustering`, `blocked`, or `pending-reference-approval`. Keep `approval` null in the helper.

The `artifact_plan` map makes planned data flow explicit with `query_fasta`, `reference_fasta`, nullable `source_rna_fasta` and `translation_fasta`, `raw_alignment`, nullable `backtranslated_alignment` and `trimmed_alignment`, `tree_input`, and `unrooted_tree`. It does not contain a `candidate_fasta` descriptor. Each non-null descriptor records its `path`, `artifact_type`, and analysis `data_kind` (or `inferred_from` for a tree); coding requests require `translation_fasta`. Do not label a protein MSA, nucleotide MSA, or codon MSA interchangeably.

Store each command as an argv array with stable ID, stage, tool, logical inputs/outputs, status, and `executed: false`. Never store a shell string. Label support semantics explicitly for FastTree, UFBoot2, standard bootstrap, and SH-aLRT.

Compute `plan_hash` over canonical decision-bearing content and semantic candidate/sequence hashes. Exclude volatile presentation time, approval, and raw byte-order hashes. A changed molecule/analysis space, comparable region, orientation, RNA transformation, CDS frame/code/translation, actual search database, sequence, metadata value, taxonomy snapshot or dump hash, threshold, cluster, outgroup, MSA/backtranslation, trim profile, model, support method, seed, command, or iTOL color must change it. Reordering semantically identical candidates must not.

### `manifest.json`

Record:

```text
schema_version, workflow_version, run_id, workflow_state,
offline, mode, project_id, request_path_and_hash,
input_artifacts, output_artifacts, query, database_provenance,
policy, decisions, plan_hash, approved_plan_hash,
tool_versions, commands, execution, warnings, errors
```

For each artifact, record logical relative path, media type, byte size, and SHA-256. Do not store the manifest's own digest inside itself. Mark unavailable tool versions `not-inspected`; never guess. Record zero network calls and external processes in plan mode.

The planning manifest's `database_provenance` is derived from candidate `source_db`, `source_release`, and `retrieved_at` values; enabled taxonomy validation adds the recorded taxdump source and hashes. The request-0.3 candidate tables separately retain `actual_search_database`, `search_program`, and `search_database_molecule`. The helper validates required declarations and supported program/database-molecule pairs, but it does not query a database, inspect a physical index, or verify a search-index build/checksum. If the host acquisition workflow records physical location, index format, release/build, or checksum, preserve those facts in a separate `search-database-provenance.tsv` or raw acquisition record without claiming that the planning manifest verified them. For coding input, preserve CDS and translation hashes plus frame/code/QC evidence. For RNA, preserve source and analysis-copy hashes and the encoding-specific transformation receipt.

When taxonomy validation is enabled, include logical `taxonomy_names` and `taxonomy_nodes` inputs, their SHA-256 values, exact archive URL, snapshot label, retrieval time, match mode, and `taxonomy_resolution.tsv` output. When it is disabled, record `taxonomy_plan.status = not-requested`; do not imply that unvalidated names were checked.

Never write credentials, headers, cookies, unpublished sequence content, home-directory paths, or absolute host paths into manifests or planned commands.

## Final host-agent execution artifacts

After approval and actual execution, preserve this structure or an equivalent manifest-linked layout:

```text
acquisition/
  raw-provider-responses/
  candidate-provenance.tsv
  search-database-provenance.tsv
  rna-analysis-copy-provenance.tsv # RNA only: U→T or T-preserving copy
  cds-translation-qc.tsv       # coding nucleotide only
  cluster-membership.tsv       # when clustering ran
alignment/
  alignment.raw.faa            # protein route
  alignment.raw.fna            # direct nucleotide route
  alignment.raw.translated.faa # codon route: protein-guided MSA
  alignment.raw.codon.fna      # codon route: untrimmed backtranslation
  alignment.trimmed.<profile>.faa  # protein route or codon protein mask
  alignment.trimmed.<profile>.fna
  alignment.trimmed.<profile>.codon.fna
  alignment-qc.tsv
  trimming-qc.tsv
  id-map.tsv
tree/
  gene-tree.fast.unrooted.nwk  # optional exploratory result
  gene-tree.unrooted.treefile
  gene-tree.rooted.nwk         # only with approved outgroup
  iqtree/
annotation/
  itol_roles.txt
  itol_ranges.txt              # only after topology review
  sequence_metadata.tsv
  taxonomy_resolution.tsv      # when local taxdump validation ran
figures/
  gene-tree.<root-state>.ggtree.svg
  gene-tree.<root-state>.ggtree.pdf
  gene-tree.<root-state>.ggtree.settings.tsv
evidence/
  conservation_assessment.tsv
  literature_evidence.tsv
  references.bib
report/
  report.md
  commands.jsonl
  checksums.sha256
```

Emit only the files applicable to the declared route; mutually exclusive examples above are not a request to duplicate one analysis in several spaces. For codon work, preserve raw and trimmed protein and backtranslated codon alignments, exact trimAl `-backtrans` argv, one-to-one IDs, triplet integrity, and translation-equality results.

Link every executed artifact to exact input hashes, the approved reference plan hash, the approved MSA hash, molecule/analysis space, and exact tool versions. For a local ggtree/ggplot2 figure, additionally record the Newick hash, metadata hash, optional iTOL-role hash, declared root state, branch-length mode, support format, layout, canvas dimensions, palette source, R version, and package versions. The renderer must refuse missing/extra/duplicate tip IDs and must not install packages or contact the network.

### `conservation_assessment.tsv`

Write one row for the focal analysis before choosing the primary alignment:

```text
analysis_unit	conservation_class	conservation_scope	conservation_basis	evidence_ids	assessment_status	limitations
```

Encode the file as UTF-8, LF-terminated, tab-delimited text; do not emit comma-separated content under a `.tsv` name. Start from `assets/conservation-assessment.example.tsv`, replace its illustrative row, and keep exactly one row per focal analysis unit. Use only the seven classes defined in `recent-msa-trimming-evidence.md`; `assessment_status` is `provisional` until the alignment/trimming review accepts or revises it, then `reviewed`. Use `not-applicable` for a method benchmark or genome-level pseudoalignment whose analysis unit cannot support one biological gene-conservation label. Keep `evidence_ids` as stable accession, orthology, domain, catalog-analysis, or citation identifiers rather than unsupported prose. Bind the reviewed file's SHA-256 to the alignment approval and final checksum manifest. If class, scope, basis, or evidence changes, reopen the alignment decision instead of silently replacing the row. The bundled planner does not infer or emit this host-side biological judgment.

## Literature evidence schema

Create rows only from real retrieved records:

```text
citation_id, title, year, journal, doi, pmid,
taxon_rank, taxa_covered, evidence_type, data_type,
inference_method, model, topology_claim, directness,
conflicts, limitations, source_url, retrieved_at
```

When a retrieved paper informs MSA or trimming choices, append these optional backward-compatible fields:

```text
alignment_method, alignment_version, alignment_parameters, alignment_reporting_status,
trimming_method, trimming_version, trimming_parameters, trimming_status,
retained_sites_or_fraction
```

Preserve exact, partial, explicit-none, and not-reported states. Never turn an unreported method into an explicit no-trim decision or infer an exact CLI flag from prose without recording that derivation.

## Failure semantics

- Exit non-zero on malformed JSON/TSV/FASTA, missing sequence matches, duplicate IDs, invalid values, unresolved query handoff, insufficient taxa/sequences, missing required outgroup/rationale, pending required clustering, or an existing output directory.
- Exit non-zero on missing/mixed molecule declarations, non-comparable regions, missing required search-database provenance, invalid RNA source/analysis-copy provenance, CDS/frame/code/translation failure, a non-1/11 code submitted to deterministic trimAl planning, mismatched CDS/translation IDs, invalid backtranslation/triplets, `MACSE_ROUTE_REQUIRED`, or any cross-analysis fallback attempt. MACSE output is review evidence only: inspect frameshift/stop handling, export a valid codon alignment under an explicit policy, verify IDs and triplets, then rerun environment routing with a typed precomputed codon alignment. Never promote MACSE output directly to tree input.
- Exit non-zero before selection when enabled NCBI taxonomy validation finds no unique character-for-character `scientific name` match, an alias-only match, ambiguity, a missing node, malformed dump evidence, a non-NCBI archive URL, or a TaxID mismatch. The host-side archive verifier, not the two-file resolver, must reject mixed snapshots and checksum failures before invocation.
- Exit non-zero before figure output when tree tips and selected metadata IDs differ, packages are unavailable, support semantics are invalid or undeclared, or any SVG/PDF/settings target already exists.
- For biological blockers discovered after valid selection, emit a blocked audit bundle and return status 3.
- Emit validation diagnostics to stderr and keep structured stdout machine-readable.
- Never convert failed validation into an empty successful plan.
- Never execute downstream work when approval hashes are absent, stale, or mismatched.
