# Rejected Ascaphus candidate attempt

This immutable evidence directory records the first 50-candidate QC screen. It
differs from the promoted reference-review set at one tip: `XP_075436319.1`
(*Ascaphus truei*) occupied the second Anura slot before it was replaced by
`XP_053555561.1` (*Bombina bombina*).

The raw Pfam 37.1 and SMART 9.0 output contains two overlapping calls for one
BRCT interval at residues 1598–1677, but no second distinct C-terminal BRCT
repeat. The declared gate requires an N-terminal RING and two distinct
C-terminal BRCT intervals, so this candidate was rejected before planning.
The BLAST and InterProScan files contain all 50 proteins from that historical
attempt; they are retained to make the exclusion independently auditable and
must not be mixed with the promoted QC files in the parent directory.

`checksums.sha256` binds the exact initial manifest, candidate FASTA, pre-QC
table, shared human query, raw tool outputs, and software-version receipt. The
frozen provider archive used to reconstruct these inputs had SHA-256
`30f1ef86dbe2989981bed005c837e7538226cae24b800aef44a32a589bfe04c0`.
