# Step 6 — Align and assess conservation

## When to read

Read after the reference set is approved and before any trimming choice. This step produces and reviews the untrimmed MSA in the declared protein, direct-nucleotide, or protein-guided codon route.

## Required inputs

- approved molecule-matched homologous reference FASTA and its hash;
- exact molecule/analysis space, comparable region, objective, domain/feature architecture, taxonomic scope, and copy-number/orthology evidence;
- RNA source-encoding/analysis-copy receipt or exact CDS/translation mapping and frame/code/translation QC when applicable;
- approved MAFFT executable/version and fixed resource settings;
- relevant recent-method evidence when literature will affect the route.

## Procedure

Confirm the data architecture first. This path accepts one homologous protein family, one comparable noncoding locus/feature, or clean comparable CDS analyzed through verified translations. Q-INS-i is available only for explicitly requested structure-aware noncoding RNA through the dedicated `mafft-qinsi` executable; more specialized covariance-model RNA analysis is outside this route. Concatenated loci, target-capture matrices, pangenome core alignments, disrupted CDS, and reference-mapped SNP pseudoalignments require architecture-specific workflows; do not add MAFFT/trimAl merely to resemble this pipeline.

Choose exactly one route:

| Route | MAFFT input and required flag | Raw outputs |
|---|---|---|
| Protein | `reference_set.faa`, `--amino` | `alignment.raw.faa` |
| Noncoding DNA/RNA | comparable `.fna`; RNA uses its provenance-bound DNA-alphabet analysis copy; `mafft --nuc` | `alignment.raw.fna` |
| Structure-aware noncoding RNA | same `.fna` analysis copy; `mafft-qinsi --nuc` (dedicated executable) | `alignment.raw.fna` |
| Clean code-1/11 CDS | `reference_set.translated.faa`, `--amino`; then trimAl `-backtrans reference_set.fna` | `alignment.raw.translated.faa`, `alignment.raw.codon.fna` |

Do not rely on MAFFT alphabet auto-detection. For clean CDS, verify unique identical IDs, length divisible by three, complete triplets, and translation equality after untrimmed backtranslation. Deterministic trimAl planning rejects codes other than 1/11 for reviewed MACSE/PAL2NAL/precomputed handling. A frame/stop/translation failure emits `MACSE_ROUTE_REQUIRED`, not a direct-nucleotide fallback.

Assign a provisional scale-aware `conservation_class`, explicit `conservation_scope`, and evidence-based `conservation_basis`. Use taxonomic distribution, orthology/copy number, marker design, and domains—not a familiar gene name or similarity alone. Record and review `evidence/conservation_assessment.tsv`; bind its SHA-256 to alignment approval.

After fixing `--amino` or `--nuc`, choose MAFFT mode by sequence count and architecture. Q-INS-i is invoked as `mafft-qinsi`, not as a `mafft --qinsi` flag:

| Mode | Arguments | Suitable use |
|---|---|---|
| Auto | `--auto` | General routing by dataset size |
| L-INS-i | `--localpair --maxiterate 1000` | Small set, one alignable domain, difficult flanks |
| G-INS-i | `--globalpair --maxiterate 1000` | Small, globally alignable full-length proteins |
| E-INS-i | `--genafpair --maxiterate 1000` | Shared motif order separated by long insertions |
| Q-INS-i | dedicated `mafft-qinsi --nuc` | Explicitly requested structure-aware noncoding RNA |

Record fixed threads and the exact version. Preserve the route-specific raw outputs. Inspect unique/reversible tip IDs, equal aligned length, ungapped length, coverage, per-tip gap fraction, terminal gaps, column occupancy, conserved motifs/regions, fragments, fusions, mixed domains/features, low complexity, long insertions, near-duplicates, and suspicious isolated branches. For nucleotide routes also inspect orientation and region boundaries; for codons inspect triplet phase and translation equality.

When recent literature is decision-bearing, filter the evidence catalog by molecule, gene/marker architecture, conservation class/scope, dataset scale, and taxonomic depth. Re-check the original source; an analogous paper is precedent, not an automatic default.

## Required outputs

- reviewed conservation assessment and hash;
- route-specific raw protein/nucleotide/codon alignments with exact MAFFT and trimAl `-backtrans` argv/version/threads and input hashes;
- per-tip and per-column QC, motif/domain review, warnings, and excluded-sequence proposals;
- explicit proposed trimming profiles or an untrimmed-analysis rationale.

## Review gate and stop conditions

Do not silently remove a sequence. Stop on mixed molecules/domains/features, non-comparable regions/orientations, duplicate or CDS/translation-mismatched IDs, failed triplet/translation QC, an unreviewed MACSE requirement, fewer than four usable sequences, severe coverage/gap failure, unidentifiable homologous regions, or a necessary record supported only by a weak local-domain match. A changed reference set, molecule route, or conservation assessment invalidates this stage.

## Supporting references

- [Recent MSA/trimming evidence and conservation vocabulary](../recent-msa-trimming-evidence.md)
- [Sequence-type alignment and CDS/RNA QC](../sequence-type-routing.md)
- [Workflow Gate 3](../workflow.md#gate-3-alignment-and-trimming-approval)
- [Alignment artifact contract](../output-contract.md)

Primary software documentation: MAFFT manual <https://mafft.cbrc.jp/alignment/software/manual/manual.html>, algorithm guide <https://mafft.cbrc.jp/alignment/software/algorithms/>, and trimAl `-backtrans` usage <https://trimal.readthedocs.io/en/latest/usage.html>.
