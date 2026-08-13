# NCBI Taxonomy exact-name resolution

Use an official NCBI Taxonomy dump snapshot to turn organism names into current TaxIDs. Treat name resolution as a validation step, not a fuzzy search.

## Contents

- [Official snapshot](#official-snapshot)
- [Exact-match policy](#exact-match-policy)
- [Local resolver](#local-resolver)
- [TaxID status and lineage](#taxid-status-and-lineage)
- [Provenance and failures](#provenance-and-failures)

## Official snapshot

Prefer NCBI `new_taxdump`, which contains `names.dmp`, `nodes.dmp`, `merged.dmp`, `delnodes.dmp`, and lineage tables:

- current archive: <https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz>
- current official MD5: <https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz.md5>
- monthly archives: <https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/>
- field definitions: <https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/taxdump_readme.txt>

Use a dated monthly archive for a publication-grade frozen analysis. The current archive changes over time. Verify the provider MD5 after download, compute a local SHA-256, and extract with a trusted archive tool into a new directory. Never mix files from different snapshots.

The bundled resolver performs no download or extraction. Supply already-extracted `names.dmp` and `nodes.dmp` from one verified snapshot.

## Exact-match policy

Resolve an organism name only when all of these are true:

1. `names.dmp` `name_txt` equals the supplied name character for character;
2. `name class` equals `scientific name`;
3. the match maps to exactly one TaxID;
4. that TaxID exists in `nodes.dmp` from the same snapshot;
5. a supplied TaxID, when present, equals the resolved TaxID.

Equality is case-, whitespace-, punctuation-, and Unicode-sensitive. Do not case-fold, collapse spaces, remove authorship, remove `Candidatus`, accept a substring, or silently trim leading/trailing whitespace. A diagnostic match to `synonym`, `equivalent name`, `common name`, or `misspelling` is not an automatic TaxID assignment.

Zero exact scientific-name matches are unresolved. Multiple exact matches are ambiguous even when the text looks like a binomial; report every TaxID and NCBI `unique name`, then require explicit disambiguation. Never choose the first record.

## Local resolver

Resolve every candidate row and verify its existing `taxon_id`:

```text
python3 <skill-root>/scripts/ncbi_taxonomy.py \
  --names <snapshot>/names.dmp \
  --nodes <snapshot>/nodes.dmp \
  --snapshot <archive-date-or-label> \
  --source-url <exact-official-NCBI-archive-URL> \
  --retrieved-at <UTC-date-or-timestamp> \
  --input candidates.tsv \
  --id-column accession --name-column species --taxid-column taxon_id \
  --out taxonomy_resolution.tsv
```

For one organism, replace `--input ...` with `--name "Homo sapiens" --expected-taxid 9606 --record-id query`.

The output records the input and matched name, name class, TaxID, direct parent TaxID, rank, snapshot metadata, and SHA-256 of both dump files. The resolver refuses overwrite and returns non-zero on malformed dumps, missing names, ambiguity, missing nodes, or mismatched TaxIDs.

Enable the same check inside `gene_to_tree.py plan` with an optional request object:

```json
{
  "taxonomy": {
    "enabled": true,
    "source": "ncbi-taxdump",
    "match_mode": "exact-scientific-name",
    "names_dmp": "taxonomy/names.dmp",
    "nodes_dmp": "taxonomy/nodes.dmp",
    "snapshot": "new_taxdump_2026-08-01",
    "source_url": "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/new_taxdump_2026-08-01.zip",
    "retrieved_at": "2026-08-13T00:00:00Z"
  }
}
```

When enabled, every candidate `species`/`taxon_id` pair is validated before selection. Dump hashes become decision-bearing, `taxonomy_resolution.tsv` is emitted, and the manifest records both dump files as inputs. This validation status applies only to the name, TaxID, direct parent, and rank fields returned from the supplied dumps; it does not validate candidate `lineage`, `clade`, ingroup/outgroup role, or outgroup suitability.

## TaxID status and lineage

For a directly supplied TaxID, a host-side acquisition workflow should interpret the same snapshot in this order:

1. present in `nodes.dmp`: current;
2. present as `old_tax_id` in `merged.dmp`: follow the complete merge chain to a current node and record the chain;
3. present in `delnodes.dmp`: deleted without a supplied replacement; stop;
4. absent from all three: unknown in that snapshot; stop.

Never guess a replacement for a deleted TaxID. For full lineage, traverse `tax_id -> parent_tax_id` in `nodes.dmp` to TaxID 1, or cross-check `taxidlineage.dmp` from the same `new_taxdump` snapshot. Detect missing parents, cycles, and non-root self-parent nodes. NCBI lineage is a classification and does not replace an inferred phylogeny or justify an outgroup by itself.

The bundled resolver implements exact name assignment plus direct `nodes.dmp` presence, parent, and rank validation. It does not parse `merged.dmp`, `delnodes.dmp`, or a full lineage table. Resolve obsolete/deleted direct TaxIDs and build full lineage in the host workflow before handing records to the planner, and retain that separate evidence.

## Provenance and failures

The downloader/host workflow should record archive URL, snapshot date/label, retrieval time, provider MD5, and local archive SHA-256. The bundled resolver records the asserted archive URL/snapshot/retrieval time plus individual `names.dmp` and `nodes.dmp` SHA-256 values, raw input name, exact match mode, matched name class, current TaxID, rank, and resolver version. Because NCBI does not embed a shared snapshot identifier in both files, the resolver cannot independently prove that two user-supplied extracted files came from the same archive; the archive verifier must establish that provenance before invocation.

The host workflow must fail closed on provider checksum failure or mixed snapshots. The bundled resolver fails closed on malformed delimiters, non-ASCII-numeric TaxIDs, a non-NCBI archive URL, no exact scientific-name match, ambiguity, missing `nodes.dmp` record, or TaxID disagreement. Suggested trimmed, case-insensitive, or secondary-name matches may be shown only as diagnostics and must never silently enter the reference table.
