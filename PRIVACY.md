# Privacy

This repository contains no telemetry or hosted service. The bundled `plan`
command performs no network requests and does not inspect credentials. The
separate third-party `skills` installer recommended in the README reports
anonymous installation telemetry unless it is run with
`DISABLE_TELEMETRY=1`; that installer behavior is not repository telemetry.

Users and host agents remain responsible for every external database, literature, structure, or iTOL request. Obtain explicit permission before transmitting an unpublished sequence, unpublished tree, or sensitive sample metadata. Treat permissions as endpoint-specific: approval for BLAST does not automatically authorize iTOL or another service.

Keep credentials in approved environment or client configuration, redact them from logs, and never write authorization headers, cookies, API keys, home-directory paths, or unpublished sequence content into manifests. Do not commit generated unpublished outputs automatically or paste private data into public issues.

Downloaded records and external services are governed by their own privacy, licensing, retention, and rate-limit policies.
