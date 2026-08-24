# Environment-aware software and workflow routing

Choose a route from observed capabilities, explicit permissions, materialized inputs, and the requested scientific claim. A Codex, Cursor, Claude Code, browser, SSH client, or scheduler name does not prove that a bioinformatics executable is installed or that a remote action is authorized.

## Contents

- [Preflight artifacts](#preflight-artifacts)
- [Task and software matrix](#task-and-software-matrix)
- [Supported environment classes](#supported-environment-classes)
- [Official installation references](#official-installation-references)
- [Deterministic selection policy](#deterministic-selection-policy)
- [Status and fallback semantics](#status-and-fallback-semantics)
- [Inspect and route](#inspect-and-route)
- [HPC and SSH handoff](#hpc-and-ssh-handoff)
- [Review checklist](#review-checklist)

## Preflight artifacts

Keep environment feasibility separate from the scientific request and `plan_hash`:

1. `environment-profile.json` records user-visible inputs, intended accuracy, required deliverables, declared host/compute capabilities, permissions, and one proposed compute class. Validate it against [environment-profile-0.1.schema.json](environment-profile-0.1.schema.json).
2. `compute-snapshot.json` records the actual compute target. An optional `host-snapshot.json` separately records host-local executables used for planning, local search, rooting, annotation, or rendering; declared host API/tree-I/O capabilities remain explicit profile assertions. Validate each snapshot against [environment-snapshot-0.1.schema.json](environment-snapshot-0.1.schema.json).
3. `route-decision.json` records every task as `ready`, `conditional`, `blocked`, or `skipped`, plus considered routes, stable reason codes, the compute snapshot hash, an optional/effective host hash, and a separate `route_hash`. On one-machine routes, the effective host hash equals the compute hash; on remote routes without host-local executables it is null. Validate the decision against [route-decision-0.1.schema.json](route-decision-0.1.schema.json).

For every query, candidate bundle, cluster mapping, alignment, tree, metadata table, taxdump, and sequence database, declare `absent`, `host`, `compute-target`, or `both`. Also declare an existing tree as `unrooted` or `rooted`. These fields prevent a host-side download from being mistaken for a file already visible to HPC, and prevent an unrooted tree from satisfying a rooted-tree request.

For `public-sequence` or `unpublished-sequence` with `query_resolved: false`, `query_location` is the raw sequence file's location; a local or remote similarity search cannot start when it is `absent`. For a resolved query it is the resolved record's location. Pair `sequence_database_location` with its verified `blast`, `mmseqs2`, or `both` format; a visible search executable is not compatible with an arbitrary database layout.

Set `compute_tree_io` or `host_tree_io` only when that site has an explicitly reviewed procedure that reroots a copy, maps support by canonical unrooted bipartition, and verifies every original split/label exactly once. `Rscript` plus `ape` alone is not sufficient evidence; the bundled ggtree renderer never reroots.

The route decision is not execution authority. Re-check authorization immediately before any database request, unpublished-sequence submission, file transfer, SSH connection, scheduler submission, package installation, or iTOL upload.

## Task and software matrix

`Required when` is conditional on the selected objective. Do not require MMseqs2 below the clustering trigger, trimAl when a reviewed untrimmed route was explicitly chosen, or a figure stack when no figure was requested.

| Step | Task | Required when | Preferred software or capability | Legitimate alternative | Network |
|---:|---|---|---|---|---|
| 1 | Resolve query | Query is not materialized | Python 3.10+ for local FASTA validation; NCBI E-utilities/Datasets, UniProt REST, or Ensembl for public IDs | Hash-bound local protein; local BLAST+/MMseqs2 for sequence context | Public lookup only; raw sequence submission is a separate permission |
| 2 | Define objective | Always | Reviewed biological decision | None needed | No |
| 3 | Discover candidates | Candidate bundle is absent | Curated orthology API, RefSeq/UniProt search, or BLAST+ against a local sequence database | MMseqs2 local search; jackhmmer/HHsearch/InterPro/Foldseek for justified distant-homology escalation | Remote discovery only when permitted |
| 4 | Select references/outgroups | Planning or later | Python 3.10+ and bundled planner | None | No; exact taxonomy uses local official `names.dmp` + `nodes.dmp` |
| 5 | Cluster expanded pool | Candidate count reaches the declared trigger | MMseqs2 with explicit identity, coverage, and coverage mode | Audited precomputed cluster mapping | No |
| 6 | Align proteins | No approved alignment exists | MAFFT | No silent algorithm substitution | No |
| 7 | Trim/test sensitivity | Trimming is enabled | trimAl; optional FastTree profile screens | Explicitly reviewed untrimmed primary MSA | No |
| 8 | Infer/check tree | No approved tree exists | FastTree for quick exploration; IQ-TREE2 for accurate ML | Existing provenance-checked tree | No |
| 8 | Root derivative | Rooted copy requested | Explicit annotation-preserving tree-I/O capability with split-label remapping and verification | R/ape only when wrapped by that verified procedure; otherwise keep the unrooted tree | No |
| 9 | Annotate | Tree and metadata exist | Bundled local iTOL writer | None | No |
| 9 | Visualize | Figure requested | Rscript + ape, ggplot2, ggtree, openssl, svglite | Separately authorized iTOL upload | iTOL only |
| 10 | Current evidence/report | Current comparison requested | Host scholarly search with DOI/PMID-linked records | Dated cached evidence plus explicit search limitation | Usually yes |

The deterministic router automatically selects only the local close-homology routes whose database location is declared (`blastp`, then MMseqs2). `ncbi-datasets`, InterProScan, jackhmmer, HHsearch, and Foldseek are inventory signals for a reviewed Step 1/3 escalation; their presence alone does not prove that the required database, profile, or structural index exists. Materialize and document that discovery result before continuing rather than claiming that the router executed it. Likewise, a visible `ssh` command never establishes a remote target.

Important semantics remain unchanged:

- MMseqs2 `-c` is coverage, not sequence identity.
- trimAl `-gt` is the minimum non-gap occupancy retained per column.
- FastTree SH-like local support is not bootstrap.
- IQ-TREE2 UFBoot2 `-B` and standard bootstrap `-b` are different methods.
- Exact TaxID matching validates nomenclature, not orthology, topology, or outgroup suitability.

## Supported environment classes

Support is capability-based. `Supported` never means the repository installs third-party software automatically.

| Environment | Support level | Best use | Boundary |
|---|---|---|---|
| Linux workstation | Full candidate | Small/medium end-to-end execution | Probe every executable and package first |
| macOS workstation | Capability-based | Planning, small/medium analysis, local R figures | No macOS CI claim; external packages vary by installation |
| Windows native | Planner plus tool-by-tool support | MAFFT, FastTree, IQ-TREE, and a separately verified trimAl build can run natively | The mixed toolchain still requires one successful probe per command; MMseqs2 recommends WSL for stability/performance |
| WSL2 | Preferred unified Windows route | Linux-compatible local workflow, especially MMseqs2 and Conda-based stacks | Probe inside WSL, not from Windows PowerShell; keep I/O-heavy databases in the Linux filesystem |
| PBS HPC, including gds2-like systems | Example-backed handoff | MMseqs2, MAFFT, trimAl, IQ-TREE2, larger jobs | Run `doctor` after modules are loaded on the actual compute environment; BRCA1 PBS scripts are example-specific |
| Slurm/LSF/other HPC | Planned handoff | Larger jobs when tools are verified | No scheduler submission adapter is bundled |
| SSH compute host | Snapshot-based handoff | Remote analysis with a verified target | Presence of an `ssh` binary is not authorization or proof of a usable target |
| Host agent with database/browser tools | Acquisition/reporting partner | Steps 1, 3, and 10 | Client support does not imply local binaries |
| Browser only | Partial | Public record review, objective design, literature, existing-tree iTOL use | Cannot execute the auditable local planner, MSA, trimming, or tree inference |
| Offline workstation/HPC | Local-bundle route | Steps 2 and 4–9 from frozen inputs | Steps 1, 3, and current Step 10 need cached/materialized evidence |

Common workflow tiers:

- `planning-offline`: Python 3.10+ plus materialized query/candidate files; emits review artifacts and commands without claiming execution.
- `workstation-online`: authorized database/literature capabilities plus a locally verified analysis stack.
- `hybrid-agent-hpc`: host agent performs acquisition and current-literature work; a verified HPC target performs compute-heavy Steps 5–8; local R or authorized iTOL handles Step 9.
- `hpc-offline`: frozen local databases, taxdump, inputs, and analysis tools on HPC; current literature must be cached or reported as unavailable.
- `browser-or-host-only`: public lookup/review only; hand off materialized files before planning or computation.

## Official installation references

Platform support changes independently of this Skill. Check the upstream pages at installation time:

- [MMseqs2 system requirements and installation](https://github.com/soedinglab/MMseqs2/wiki#system-requirements) covers Linux/macOS builds, a Windows preview, and the upstream WSL recommendation.
- [MAFFT downloads](https://mafft.ddbj.nig.ac.jp/alignment/software/) lists Linux, macOS, Windows, and source distributions.
- [trimAl upstream build matrix](https://github.com/inab/trimal/blob/trimAl/.github/workflows/build.yml) exercises Linux x86_64/ARM64, macOS Intel/ARM, and Windows MSYS2/MinGW builds.
- [FastTree installation](https://morgannprice.github.io/fasttree/) provides Linux/Windows binaries and source-build guidance for macOS and other systems.
- [IQ-TREE downloads](https://iqtree.github.io/) and [Quickstart](https://iqtree.github.io/doc/Quickstart) cover current Linux, macOS, and Windows binaries. Current IQ-TREE 3 commonly uses `iqtree3`; this Skill inventories `iqtree3`, `iqtree2`, and `iqtree` separately and never treats them as silent substitutes.
- [Bioconductor ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html) records current R/Bioconductor dependencies and platform binaries.
- [Microsoft's WSL filesystem guidance](https://learn.microsoft.com/en-us/windows/wsl/filesystems) recommends storing Linux-tool workloads in the WSL filesystem rather than `/mnt/c` when Linux command-line performance matters.

## Deterministic selection policy

Apply this order; do not use an opaque score:

1. Derive required tasks from the intent, materialized inputs, clustering trigger, trimming choice, rooting requirement, visualization request, and report request.
2. Apply privacy and data-egress gates before considering a remote lookup or upload.
3. Treat `missing`, `unknown`, `incompatible`, and `probe-failed` as unavailable. Default `path-only` discovery proves only that a command name resolves; it adds a review limitation because neither version nor runtime compatibility was exercised.
4. Honor explicit intent:
   - `accurate` requires IQ-TREE2 and must not fall back to FastTree;
   - `quick` requires FastTree;
   - `planning` requires no downstream executable;
   - `visualization` starts from an existing tree;
   - `auto` tries accurate, then quick, then planning-only.
5. Track where each input is created and whether authorized transfer makes it reachable by the next task. Host acquisition, compute-target inference, host rooting/rendering, and iTOL upload are distinct transitions.
6. Prefer materialized local inputs over repeated network retrieval when provenance is adequate.
7. Prefer a local workstation when it fully satisfies the request. Use an HPC/SSH handoff when the selected target is remote or local resources are inadequate. The router reports the handoff; it never performs it.
8. Prefer host or compute-target ggtree rendering over iTOL when both satisfy `visualization: either`, avoiding an unnecessary upload.
9. Emit one decision record for all ten steps. Re-run routing after any profile field or used snapshot changes, including intent, requirements, query state, artifact locations/root state, capabilities, permissions, candidate count/trigger, scheduler, resources, or compute target.

One invocation evaluates one proposed compute target. When local, HPC, and SSH targets are all candidates, compile one profile/decision per target with the same scientific requirements, then choose deterministically: discard blocked routes; preserve explicit `accurate`/`quick` intent; prefer a fully satisfying local route over handoff at equal quality; then prefer fewer limitations; use stable `profile_id` ordering only as a final tie-break. Preserve every rejected decision rather than claiming that a single snapshot compared all machines.

The environment `route_hash` changes with the normalized profile, snapshot, or decision. It is deliberately independent of the scientific `plan_hash`; a tool installation does not silently approve a new reference set, alignment, model, or topology claim.

## Status and fallback semantics

Overall route statuses:

| Status | Meaning |
|---|---|
| `ready` | Every required stage is feasible from the declared permissions, capabilities, artifact locations, and snapshot evidence, with no recorded limitation |
| `ready-with-limitations` | The route is feasible, but conditional acquisition/evidence, PATH-only discovery, or resource uncertainty/overcommit remains in `limitations` |
| `planning-only` | The best allowed route can prepare/review inputs and commands but cannot claim downstream execution |
| `handoff-required` | An HPC/SSH target is selected from its snapshot; transfer/submission remains external, separately authorized, and subject to every recorded limitation |
| `blocked` | No permitted route satisfies the explicit objective |

Never silently:

- change `accurate` to FastTree;
- disable requested trimming because trimAl is missing;
- skip triggered MMseqs2 clustering;
- upload unpublished data because public metadata lookup was allowed;
- treat a missing local figure stack as permission to use iTOL;
- infer an organism or orthology from the closest local hit;
- call scheduler/SSH availability an executed handoff.

## Inspect and route

Start from [the example profile](../assets/environment-profile.example.json), replace every illustrative assertion with the user's actual situation, preserve denied permissions as `false`, and record every input's actual location.

Create a snapshot on the environment where commands would actually run:

```text
python3 <skill-root>/scripts/gene_to_tree.py doctor \
  --environment-id local --kind local \
  --out environment-snapshot.json
```

Default `doctor` performs `PATH` discovery only and launches no discovered command. It records OS/architecture, Python, CPU/memory/scratch hints, visible command basenames, and `unknown` R-package status; `network` remains `not-probed`. This is the safest portable inventory, but it does not verify versions or runtime compatibility.

After approving local probes, opt in to version/help commands and R namespace checks:

```text
python3 <skill-root>/scripts/gene_to_tree.py doctor \
  --environment-id local --kind local \
  --run-version-probes \
  --out environment-snapshot.active.json
```

Active mode executes local programs and can trigger wrapper or package-load hooks; do not describe it as side-effect-free or guaranteed offline. The snapshot retains only executable basenames and parsed version tokens, never raw output, hostnames, or paths. The router honors an imported `incompatible` status but does not infer compatibility from a version number; review every selected tool against the project and platform requirements.

Compile the route:

```text
python3 <skill-root>/scripts/gene_to_tree.py route \
  --profile environment-profile.json \
  --snapshot compute-snapshot.json \
  --host-snapshot host-snapshot.json \
  --out route-decision.json
```

Omit `--host-snapshot` when compute and host are the same machine or no host-side executable is proposed. Omit `--snapshot` only when routing the current local machine; the command then performs default `PATH` discovery. Both commands refuse to overwrite an existing output file. A blocked decision exits non-zero while still emitting its audit record.

## HPC and SSH handoff

For PBS, Slurm, LSF, or SSH targets:

1. Load or activate the exact intended tool environment on the target.
2. Run `doctor` there with a non-sensitive environment ID and the correct `--kind`; opt into active probes only if approved.
3. If the route depends on host-local planning, search, rooting, annotation, or rendering executables, create a separate `--kind local` host snapshot.
4. Transfer only the JSON compute snapshot back through an authorized channel.
5. Route using the compute snapshot, the host snapshot when applicable, and a profile that gives each input its current location.
6. Review `handoff-required`, probe modes, software versions, scheduler visibility, artifact transfers, root state, resource hints, and every blocked/conditional stage.
7. Obtain authorization separately, then transfer data and submit commands outside the helper.

Do not store hostnames, usernames, private paths, module-load commands, tokens, scheduler job IDs, or SSH configuration in the portable snapshot. The presence of `qsub`, `sbatch`, `bsub`, or `ssh` shows only that a launcher is visible.

## Review checklist

- Profile values describe the user's current request, not an assumed default.
- Compute snapshot was generated on the actual target after the intended environment was activated; any host snapshot represents a different, explicit role.
- Every materialized artifact has a declared location, and every host/compute transition has file-transfer capability.
- Existing-tree root state is explicit; an unrooted tree never bypasses the rooting requirement.
- Unpublished-sequence submission and iTOL upload permissions are separate and default denied.
- Exact-taxonomy requests have one verified local NCBI taxdump snapshot or an explicit acquisition prerequisite.
- `auto` fallback is visible; explicit quick/accurate choices are not changed.
- Every required step is `ready`, `conditional`, or deliberately handed off; blockers are shown to the user.
- Remote execution is separately authorized after route review.
- Scientific approvals and `plan_hash` remain independent from `route_hash`.
