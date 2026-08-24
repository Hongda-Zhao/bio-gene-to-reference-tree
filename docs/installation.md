# Installation

[Back to Install on the project homepage](../README.md#install)

[skills CLI](#skills-cli) · [Invoke](#invoke) · [Claude Marketplace](#claude-code-marketplace) · [Manual locations](#manual-locations) · [Cursor Remote Rule](#cursor-remote-rule)

Install the complete `skills/bio-gene-to-reference-tree/` package, not `SKILL.md` alone. Its relative links require the bundled `references/`, `scripts/`, and `assets/` directories.

## skills CLI

Interactive installation for a supported agent:

```bash
npx skills add Hongda-Zhao/bio-gene-to-reference-tree \
  --skill bio-gene-to-reference-tree
```

Explicit global installations:

```bash
# Codex
npx skills add Hongda-Zhao/bio-gene-to-reference-tree \
  --skill bio-gene-to-reference-tree --agent codex --global

# Cursor
npx skills add Hongda-Zhao/bio-gene-to-reference-tree \
  --skill bio-gene-to-reference-tree --agent cursor --global

# Claude Code
npx skills add Hongda-Zhao/bio-gene-to-reference-tree \
  --skill bio-gene-to-reference-tree --agent claude-code --global
```

Omit `--global` for project scope. The third-party `skills` CLI reports anonymous installation telemetry by default; set `DISABLE_TELEMETRY=1` if you do not want the installation counted.

## Invoke

| Client | Command |
|---|---|
| Codex | `$bio-gene-to-reference-tree ...` |
| Cursor 2.4+ | `/bio-gene-to-reference-tree ...` |
| Claude Code | `/bio-gene-to-reference-tree ...` |

Cursor users can type `/` in Agent chat to discover the Skill. `Option+Enter` on macOS or `Alt+Enter` on Windows keeps it active as a Custom Mode for the session.

## Claude Code Marketplace

This repository includes a thin `.claude-plugin/marketplace.json` wrapper:

```text
/plugin marketplace add Hongda-Zhao/bio-gene-to-reference-tree
/plugin install bio-gene-to-reference-tree@hongda-zhao-bio-skills
```

Run `/reload-plugins` if the installation summary requests it. A marketplace installation also exposes the namespaced command:

```text
/bio-gene-to-reference-tree:bio-gene-to-reference-tree ...
```

See the official [Claude Code Skills](https://code.claude.com/docs/en/skills) and [plugin marketplace](https://code.claude.com/docs/en/plugin-marketplaces) documentation.

## Manual locations

Copy the complete Skill directory to one supported location:

| Client | Project-level destination | User-level destination |
|---|---|---|
| Codex | `.codex/skills/bio-gene-to-reference-tree/` | `~/.codex/skills/bio-gene-to-reference-tree/` |
| Cursor | `.cursor/skills/bio-gene-to-reference-tree/` or `.agents/skills/bio-gene-to-reference-tree/` | `~/.cursor/skills/bio-gene-to-reference-tree/` or `~/.agents/skills/bio-gene-to-reference-tree/` |
| Claude Code | `.claude/skills/bio-gene-to-reference-tree/` | `~/.claude/skills/bio-gene-to-reference-tree/` |

Use `.agents/skills/bio-gene-to-reference-tree/` only with Cursor or another client that explicitly supports that shared directory; a Claude Code project installation belongs in `.claude/skills/bio-gene-to-reference-tree/`.

The same canonical Skill is copied to every supported location; no forked Cursor- or Claude-specific prompt is required.

## Cursor Remote Rule

Cursor's [Skills documentation](https://cursor.com/docs/skills) describes a logged-in UI route:

1. Open **Customize → Rules → Add Rule → Remote Rule (Github)**.
2. Enter `https://github.com/Hongda-Zhao/bio-gene-to-reference-tree`.
3. Confirm the result under **Customize → Skills**.

If the UI does not import the nested Skill, use the CLI or copy the complete directory manually. Restart Cursor if a newly installed Skill does not appear. Claude Code watches an existing `.claude/skills/` directory live; restart it only if that top-level directory did not exist when the session began.
