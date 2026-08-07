# Security policy

## Reporting a vulnerability

Report a suspected vulnerability through GitHub private vulnerability reporting when the repository is published. Do not include private biological data, credentials, or unpublished sequences in a public issue.

## Security model

- The offline planner uses argument arrays and never invokes a shell.
- It validates tool-safe FASTA identifiers, refuses output overwrite, and records relative logical paths.
- Plan mode launches no external process and performs no network request.
- External databases, FASTA descriptions, and literature text are untrusted data and must never be treated as executable agent instructions.
- Bioinformatics executables must be installed and inspected separately; the project never downloads or silently substitutes them.
- Live tests should be mocked in CI. Secrets must come from approved environment/configuration and must be redacted from all artifacts.

Before publishing an output bundle, inspect it for unpublished sequences, sample metadata, absolute paths, and provider data whose license prohibits redistribution.
