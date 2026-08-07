# Alignment, trimming, and tree inference

Treat the alignment as a biological homology hypothesis. A tree inference program will fit misaligned columns confidently, so stop before inference when the MSA is not defensible.

## MAFFT routing

| Mode | Arguments | Use |
|---|---|---|
| Auto | `--auto` | General routing by dataset size |
| L-INS-i | `--localpair --maxiterate 1000` | Small set, one alignable domain, difficult flanks |
| G-INS-i | `--globalpair --maxiterate 1000` | Small globally alignable full-length proteins |
| E-INS-i | `--genafpair --maxiterate 1000` | Shared motif order with long insertions or long unalignable regions |

Record a fixed thread count and exact MAFFT version. Route primarily by sequence count and domain architecture, not a close/distant label alone.

## Raw-alignment QC

Preserve `alignment.raw.faa` and inspect:

- unique and reversible tip IDs;
- equal aligned length and at least four usable sequences;
- per-sequence ungapped length, query coverage, gap fraction, and terminal gaps;
- per-column occupancy and retained conserved motifs;
- fragments, fusions, mixed domains, low complexity, and unusually long insertions;
- unexpected near-identical duplicates or long isolated branches.

Stop when molecule types or domains are incompatible, the alignment is dominated by gaps, homologous regions cannot be identified, or a biologically necessary sequence is only a weak local-domain match.

## trimAl semantics

Use `trimal -gt x` as a minimum non-gap occupancy threshold. For example, `-gt 0.9` retains a column only when roughly 90% of sequences contain a residue. Therefore:

- `0.98` is extremely strict and can delete informative indel-rich regions;
- `0.90` is strict for heterogeneous sampling;
- `0.10` and `0.05` are extremely permissive and may retain noise.

Treat proposed values such as `0.98/0.95/0.90` or `0.10/0.05` as named sensitivity profiles, not automatic biological truths. Preserve every profile, record retained columns and retained fraction, and verify key motifs. Require review when the primary topology changes across reasonable profiles.

Example argument array:

```text
["trimal", "-in", "alignment.raw.faa", "-out", "alignment.trimmed.balanced.faa", "-gt", "0.9"]
```

Do not use permissive trimming as a substitute for viral recombination analysis or for fixing mixed domain architecture.

## Fast exploratory inference

Use FastTree only for a rapid first-pass approximate-ML topology, especially for a very large reference pool. A common protein command is:

```text
FastTree -wag -gamma alignment.raw.faa
```

Label its default internal-node values as SH-like local support. Do not call them 1,000-replicate bootstrap values. Re-run the curated final alignment with IQ-TREE2 when the result will support a biological claim.

## IQ-TREE2 accurate inference

Use ModelFinder and dual support by default:

```text
iqtree2 -s alignment.trimmed.balanced.faa -m MFP \
  -B 1000 -bnni -alrt 1000 -T <fixed> -seed <fixed> --prefix gene-tree
```

Interpret the methods separately:

| Measure | Flag | Common strong-support rule | Meaning |
|---|---|---|---|
| UFBoot2 | `-B 1000` | at least 95 | Fast resampling-based repeatability |
| SH-aLRT | `-alrt 1000` | at least 80 | Local likelihood support against NNI alternatives |
| Standard bootstrap | `-b 1000` | often at least 70 | Traditional nonparametric site-resampling frequency |

Use `-bnni` with UFBoot2 under model violation. Do not use `-bnni` with standard bootstrap. Never apply a 70% standard-bootstrap rule to UFBoot values.

For deep protein relationships, test whether site-homogeneous ModelFinder candidates are adequate. Consider C10–C60 profile mixtures and PMSF for site-compositional heterogeneity or long-branch attraction. Test sensitivity to faster sites and outgroup choice. High support under one inadequate model can be confidently wrong.

## Rooting and topology checks

Infer and retain an unrooted tree. Create a rooted derivative only with approved homologous outgroup accessions that lie outside but near the ingroup. Prefer multiple nearby candidates when feasible. Re-check the topology after removing a long-branched outgroup or switching to a richer model.

Stop when the proposed outgroup is a distant paralog, non-homologous, excessively long-branched, inside the ingroup, or unstable across reasonable analyses. Do not midpoint-root automatically.

## Required records

Record executable version, argv array, working-directory-neutral inputs and outputs, threads, seed, model-selection result, support method, alignment hash, plan hash, exit status, stdout/stderr files, and native reports. Never represent an argv array as a shell string or interpolate untrusted FASTA headers into a command.

## Primary documentation

- MAFFT manual: <https://mafft.cbrc.jp/alignment/software/manual/manual.html>
- MAFFT algorithm guide: <https://mafft.cbrc.jp/alignment/software/algorithms/>
- trimAl documentation: <https://vicfero.github.io/trimal/whatcanido.html>
- FastTree documentation: <https://morgannprice.github.io/fasttree/>
- IQ-TREE tutorial: <https://iqtree.github.io/doc/Tutorial>
