# Synthetic offline-review fixture

Every identifier and sequence in this directory is synthetic. The fixture is
small enough to inspect by eye and exercises one decision at a time:

- one query self-hit;
- two records from the same taxon, where the canonical record wins;
- a low-coverage fragment;
- a paralog rejected by `ortholog-tree` mode;
- two eligible outgroups competing for one outgroup slot; and
- a second input bundle containing exactly the same records in a different
  order.

These fixtures intentionally retain request schema 0.1 to test the compatibility
migration. The shuffled bundle must produce the same selected accessions,
semantic `run_id`, and `plan_hash`. Raw-file checksums may still differ and
should remain recorded in the manifest.
