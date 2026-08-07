# Contributing

Contributions should preserve deterministic offline planning, explicit biological decision gates, portable Agent Skills instructions, and honest capability claims.

Before opening a pull request:

1. Keep `SKILL.md` concise and vendor-neutral; place detailed policies in one-level references.
2. Store every executable command as an argv array and never use `shell=True`.
3. Add synthetic fixtures only; do not commit restricted database records or unpublished sequences.
4. Add regression tests for new schemas, artifacts, reason codes, command semantics, and privacy properties.
5. Run `python3 -m unittest discover -s tests -v` on Python 3.10 or newer.
6. Update README capability claims, the script version, schemas, and agent metadata together.

Scientific changes should explain how orthology, taxonomic sampling, domain architecture, outgroup choice, alignment uncertainty, support semantics, or gene-tree/species-tree interpretation is protected.
