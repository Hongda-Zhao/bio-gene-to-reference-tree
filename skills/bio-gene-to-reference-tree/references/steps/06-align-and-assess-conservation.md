# Step 6 — Align and assess conservation

## When to read

Read after the reference set is approved and before any trimming choice. This step produces and reviews the untrimmed protein MSA.

## Required inputs

- approved homologous protein reference FASTA and its hash;
- objective, domain architecture, taxonomic scope, and copy-number/orthology evidence;
- approved MAFFT executable/version and fixed resource settings;
- relevant recent-method evidence when literature will affect the route.

## Procedure

Confirm the data architecture first. This path accepts homologous protein sequences for one gene/family. Codon alignments, concatenated loci, target-capture matrices, pangenome core alignments, and reference-mapped SNP pseudoalignments require architecture-specific workflows; do not add MAFFT/trimAl merely to resemble this pipeline.

Assign a provisional scale-aware `conservation_class`, explicit `conservation_scope`, and evidence-based `conservation_basis`. Use taxonomic distribution, orthology/copy number, marker design, and domains—not a familiar gene name or similarity alone. Record and review `evidence/conservation_assessment.tsv`; bind its SHA-256 to alignment approval.

Choose MAFFT by sequence count and architecture:

| Mode | Arguments | Suitable use |
|---|---|---|
| Auto | `--auto` | General routing by dataset size |
| L-INS-i | `--localpair --maxiterate 1000` | Small set, one alignable domain, difficult flanks |
| G-INS-i | `--globalpair --maxiterate 1000` | Small, globally alignable full-length proteins |
| E-INS-i | `--genafpair --maxiterate 1000` | Shared motif order separated by long insertions |

Record fixed threads and the exact version. Preserve `alignment.raw.faa`. Inspect unique/reversible tip IDs, equal aligned length, ungapped length, coverage, per-tip gap fraction, terminal gaps, column occupancy, conserved motifs, fragments, fusions, mixed domains, low complexity, long insertions, near-duplicates, and suspicious isolated branches.

When recent literature is decision-bearing, filter the evidence catalog by molecule, gene/marker architecture, conservation class/scope, dataset scale, and taxonomic depth. Re-check the original source; an analogous paper is precedent, not an automatic default.

## Required outputs

- reviewed conservation assessment and hash;
- `alignment.raw.faa` with exact MAFFT argv/version/threads and input hash;
- per-tip and per-column QC, motif/domain review, warnings, and excluded-sequence proposals;
- explicit proposed trimming profiles or an untrimmed-analysis rationale.

## Review gate and stop conditions

Do not silently remove a sequence. Stop on mixed molecules/domains, duplicate IDs, fewer than four usable sequences, severe coverage/gap failure, unidentifiable homologous regions, or a necessary record supported only by a weak local-domain match. A changed reference set or conservation assessment invalidates this stage.

## Supporting references

- [Recent MSA/trimming evidence and conservation vocabulary](../recent-msa-trimming-evidence.md)
- [Workflow Gate 3](../workflow.md#gate-3-alignment-and-trimming-approval)
- [Alignment artifact contract](../output-contract.md)

Primary software documentation: MAFFT manual <https://mafft.cbrc.jp/alignment/software/manual/manual.html> and algorithm guide <https://mafft.cbrc.jp/alignment/software/algorithms/>.
