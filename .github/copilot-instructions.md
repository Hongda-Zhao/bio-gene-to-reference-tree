# Repository instructions

- Treat `skills/bio-gene-to-reference-tree/` as the sole canonical, installable Skill.
- When a request concerns this gene-tree workflow, use the closest task adapter in `.github/skills/`, then follow its links into the canonical Skill.
- Keep `.github/skills/` thin: do not copy scientific policy, commands, schemas, scripts, references, or assets into the adapters.
- When task routing changes, update the adapters, README diagram, and metadata tests together.
