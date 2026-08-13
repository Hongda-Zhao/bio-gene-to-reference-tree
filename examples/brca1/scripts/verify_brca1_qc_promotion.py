#!/usr/bin/env python3
"""Rebuild BRCA1 QC summaries and audit their historical promotion byte for byte.

The original copy utility and copy time were not retained.  This verifier does
not invent them: it reruns the deterministic summarizer from the committed
inputs, compares both regenerated outputs with their promoted destinations,
and writes a post-hoc reconciliation receipt.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path


UTC_TIMESTAMP = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")


class PromotionAuditError(RuntimeError):
    """Raised when promoted outputs cannot be reproduced exactly."""


def fail(message: str) -> None:
    raise PromotionAuditError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require_regular_file(path: Path, label: str) -> None:
    if path.is_symlink() or not path.is_file() or path.stat().st_size == 0:
        fail(f"{label} must be a non-empty regular file: {path}")


def logical_path(path: Path) -> str:
    value = path.as_posix()
    if path.is_absolute() or value.startswith("../") or "/../" in value:
        fail(f"public receipt paths must be repository-relative: {path}")
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summarizer", type=Path, required=True)
    parser.add_argument("--candidate-fasta", type=Path, required=True)
    parser.add_argument("--candidate-table", type=Path, required=True)
    parser.add_argument("--blast-tsv", type=Path, required=True)
    parser.add_argument("--interpro-tsv", type=Path, required=True)
    parser.add_argument("--promoted-qc", type=Path, required=True)
    parser.add_argument("--promoted-candidates", type=Path, required=True)
    parser.add_argument("--verified-at-utc", required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not UTC_TIMESTAMP.fullmatch(args.verified_at_utc):
        fail("--verified-at-utc must use YYYY-MM-DDTHH:MM:SSZ")
    if args.output.exists() or args.output.is_symlink():
        fail(f"refusing to overwrite existing receipt: {args.output}")

    inputs = (
        (args.summarizer, "summarizer"),
        (args.candidate_fasta, "candidate FASTA"),
        (args.candidate_table, "pre-QC candidate table"),
        (args.blast_tsv, "BLAST TSV"),
        (args.interpro_tsv, "InterProScan TSV"),
        (args.promoted_qc, "promoted QC table"),
        (args.promoted_candidates, "promoted candidate table"),
    )
    for path, label in inputs:
        require_regular_file(path, label)

    with tempfile.TemporaryDirectory(prefix="brca1-qc-promotion-") as temporary:
        summary_dir = Path(temporary) / "summary"
        completed = subprocess.run(
            [
                sys.executable,
                str(args.summarizer),
                "--candidate-fasta",
                str(args.candidate_fasta),
                "--candidate-table",
                str(args.candidate_table),
                "--blast-tsv",
                str(args.blast_tsv),
                "--interpro-tsv",
                str(args.interpro_tsv),
                "--out-dir",
                str(summary_dir),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            fail(
                "deterministic QC regeneration failed: "
                + (completed.stderr.strip() or completed.stdout.strip())
            )
        regenerated = {
            "candidate_qc.tsv": summary_dir / "candidate_qc.tsv",
            "candidates.tsv": summary_dir / "candidates.tsv",
        }
        promoted = {
            "candidate_qc.tsv": args.promoted_qc,
            "candidates.tsv": args.promoted_candidates,
        }
        rows: list[dict[str, object]] = []
        for name in ("candidate_qc.tsv", "candidates.tsv"):
            source = regenerated[name]
            destination = promoted[name]
            require_regular_file(source, f"regenerated {name}")
            source_bytes = source.read_bytes()
            destination_bytes = destination.read_bytes()
            if source_bytes != destination_bytes:
                fail(f"promoted output is not byte-identical to regeneration: {destination}")
            digest = hashlib.sha256(source_bytes).hexdigest()
            rows.append(
                {
                    "historical_source_path": f"qc/summary/{name}",
                    "promoted_destination_path": logical_path(destination),
                    "sha256": digest,
                    "byte_size": len(source_bytes),
                    "verification": "byte-identical-deterministic-regeneration",
                }
            )

    receipt = {
        "schema_version": "1.0",
        "record_type": "post-hoc-qc-artifact-promotion-reconciliation",
        "status": "verified",
        "historical_copy_utility": "not-recorded",
        "historical_copy_time_utc": "not-recorded",
        "prospective_replay_claimed": False,
        "verified_at_utc": args.verified_at_utc,
        "regeneration_inputs": {
            logical_path(path): file_sha256(path)
            for path in (
                args.summarizer,
                args.candidate_fasta,
                args.candidate_table,
                args.blast_tsv,
                args.interpro_tsv,
            )
        },
        "promotions": rows,
        "interpretation": (
            "The promoted files reproduce byte for byte from the recorded summarizer "
            "inputs. The historical copy program and time were not retained and are not invented."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(receipt, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )


if __name__ == "__main__":
    try:
        main()
    except PromotionAuditError as error:
        print(f"BRCA1 QC promotion audit failed: {error}", file=sys.stderr)
        raise SystemExit(2) from error
