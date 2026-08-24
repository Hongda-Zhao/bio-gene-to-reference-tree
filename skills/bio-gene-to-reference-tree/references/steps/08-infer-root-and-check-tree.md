# Step 8 — Infer, root, and check the tree

## When to read

Read only after an alignment hash and inference plan have passed the alignment/trimming approval gate.

## Required inputs

- approved primary alignment and hash;
- exact executable/version, resource limits, fixed threads, fixed seed, and fresh output prefix;
- exploratory versus primary-inference intent;
- approved outgroup accessions or an explicit unrooted-only decision.

## Procedure

For rapid exploration, FastTree may provide an approximate-ML topology:

```text
FastTree -wag -gamma alignment.raw.faa
```

Label its default internal values as SH-like local support, never as 1,000-replicate bootstrap. Re-run the curated final alignment with IQ-TREE2 before making a primary biological claim.

For primary accurate inference, use ModelFinder and dual support unless the approved design says otherwise:

```text
iqtree2 -s alignment.trimmed.balanced.faa -m MFP \
  -B 1000 -bnni -alrt 1000 -T <fixed> -seed <fixed> --prefix gene-tree
```

Keep support methods distinct: UFBoot2 uses `-B` (often interpreted with a 95 threshold); SH-aLRT uses `-alrt` (often 80); standard nonparametric bootstrap uses `-b` (often 70) and must not be combined with `-bnni`. These are method-specific repeatability measures, not proof of correctness.

For deep, compositionally heterogeneous, or long-branch-prone proteins, test whether site-homogeneous ModelFinder candidates are adequate; consider C10–C60/PMSF, faster-site sensitivity, and outgroup sensitivity. High support under one inadequate model may be wrong.

Always retain the native unrooted tree and logs. Create a separately named rooted derivative only from approved homologous outgroups outside the ingroup. Never midpoint-root automatically. If support is retained after rerooting, map exact labels by canonical unrooted bipartition—not internal-node number—and verify one-to-one preservation; otherwise show support only on the unrooted tree.

## Required outputs

- unrooted Newick/treefile, model-selection results, native reports, stdout/stderr, and hashes;
- executable version, argv array, threads, seed, support semantics, exit status, and alignment/plan hashes;
- separately named rooted copy plus outgroup list and root rationale when approved;
- model, fast-site, trim-profile, and outgroup sensitivity notes where applicable.

## Review gate and stop conditions

Stop if the alignment hash differs from approval, output already exists, the outgroup is non-homologous/inside the ingroup/a distant paralog, rooting shifts support without split-safe remapping, or topology is unstable under required sensitivity checks. Re-open the relevant gate after any decision-bearing change.

## Supporting references

- [Workflow tree-inference gate](../workflow.md#tree-inference-gate)
- [Command, manifest, and executed-report contract](../output-contract.md)
- [Tool and executable boundaries](../tool-routing.md)

Primary software documentation: FastTree <https://morgannprice.github.io/fasttree/> and IQ-TREE tutorial <https://iqtree.github.io/doc/Tutorial>.
