#!/usr/bin/env python3
"""Publish a fail-closed, redacted subset of native IQ-TREE outputs.

The raw calculation directory is treated as private.  This program validates
the complete directory before creating its output directory, copies only a
small explicit artifact whitelist, and records content hashes for every
published file.  IQ-TREE checkpoints, model archives, and wrapper checksums
are bound into a raw-output receipt but are intentionally not published.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import stat
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


PUBLIC_ALIGNMENT_PATH = "alignment/alignment.trimmed.balanced.faa"

# Longest suffixes must be tested first (``.best_scheme.nex`` before ``.nex``).
IQTREE_SUFFIXES = tuple(
    sorted(
        (
            ".log",
            ".iqtree",
            ".treefile",
            ".contree",
            ".ufboot",
            ".splits.nex",
            ".suptree",
            ".best_scheme.nex",
            ".bionj",
            ".mldist",
        ),
        key=len,
        reverse=True,
    )
)
CORE_SUFFIXES = (".log", ".iqtree", ".treefile", ".contree", ".ufboot")
TEXT_SUFFIXES = (
    ".log",
    ".iqtree",
    ".nex",
    ".treefile",
    ".contree",
    ".ufboot",
    ".txt",
    ".suptree",
    ".bionj",
    ".mldist",
)
OMITTED_CHECKSUM = "checksums.sha256"
SOFTWARE_VERSION = "software_version.txt"

SAFE_PUBLISHED_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\Z")
HOST_LINE = re.compile(
    r"^([ \t]*Host[ \t]*:[ \t]*)\S+",
    flags=re.MULTILINE,
)
# Match a Unix absolute path without mistaking the slashes in an HTTP URL or
# an IQ-TREE support label such as 99/100 for a path.
UNIX_ABSOLUTE_PATH = re.compile(
    r"(?<![A-Za-z0-9._:/-])/(?!/)[^\s\x00\"'<>]+"
)
WINDOWS_ABSOLUTE_PATH = re.compile(
    r"(?<![A-Za-z0-9._-])(?:[A-Za-z]:\\|\\\\[A-Za-z0-9])[^\s\x00\"'<>]*"
)
HOME_SHORTHAND_PATH = re.compile(r"(?<![A-Za-z0-9._-])~/[^\s\x00\"'<>]+")
FILE_URI = re.compile(r"\bfile://[^\s\x00\"'<>]+", flags=re.IGNORECASE)
LOGIN_ADDRESS = re.compile(
    r"(?<![A-Za-z0-9_.+-])[A-Za-z0-9_.+-]+@"
    r"(?:[A-Za-z0-9-]+\.)*[A-Za-z0-9-]+(?![A-Za-z0-9_.-])"
)
PBS_LABEL = re.compile(
    r"\b(?:PBS_JOBID|PBS_JOB_ID|Job[ _-]?ID)\b[ \t]*[:=][ \t]*\S+",
    flags=re.IGNORECASE,
)
PBS_DOTTED_ID = re.compile(
    r"(?<![A-Za-z0-9])\d{4,}\.[A-Za-z][A-Za-z0-9._-]*(?![A-Za-z0-9])"
)
PBS_LOG_NAME = re.compile(
    r"(?<![A-Za-z0-9._-])[A-Za-z][A-Za-z0-9._-]*\.[oe]\d{4,}(?!\d)"
)


class PublicationError(RuntimeError):
    """Raised when an IQ-TREE output directory is unsafe to publish."""


@dataclass(frozen=True)
class PublishedArtifact:
    name: str
    original: bytes
    public: bytes
    redactions: dict[str, int]

    @property
    def original_sha256(self) -> str:
        return hashlib.sha256(self.original).hexdigest()

    @property
    def public_sha256(self) -> str:
        return hashlib.sha256(self.public).hexdigest()


def fail(message: str) -> None:
    raise PublicationError(message)


def iqtree_suffix(name: str) -> str | None:
    for suffix in IQTREE_SUFFIXES:
        if name.endswith(suffix) and len(name) > len(suffix):
            return suffix
    return None


def is_text_artifact(name: str) -> bool:
    return name.endswith(TEXT_SUFFIXES)


def reject_private_residue(text: str, filename: str) -> None:
    checks = (
        (UNIX_ABSOLUTE_PATH, "absolute Unix path"),
        (WINDOWS_ABSOLUTE_PATH, "absolute Windows or UNC path"),
        (HOME_SHORTHAND_PATH, "home-relative path"),
        (FILE_URI, "file URI"),
        (LOGIN_ADDRESS, "login address"),
        (PBS_LABEL, "scheduler job identifier"),
        (PBS_DOTTED_ID, "scheduler job identifier"),
        (PBS_LOG_NAME, "scheduler log identifier"),
    )
    for pattern, label in checks:
        if pattern.search(text):
            fail(f"{filename}: remaining {label} after permitted redactions")


def sanitize_text(
    data: bytes, *, filename: str, private_alignment_path: str
) -> tuple[bytes, dict[str, int]]:
    try:
        text = data.decode("utf-8", errors="strict")
    except UnicodeDecodeError:
        fail(f"{filename}: expected UTF-8 text")

    redactions: dict[str, int] = {}
    path_count = text.count(private_alignment_path)
    if path_count:
        text = text.replace(private_alignment_path, PUBLIC_ALIGNMENT_PATH)
        redactions["private_alignment_path"] = path_count

    text, host_count = HOST_LINE.subn(
        lambda match: f"{match.group(1)}<batch-compute-node>", text
    )
    if host_count:
        redactions["host_line"] = host_count

    reject_private_residue(text, filename)
    return text.encode("utf-8"), redactions


def raw_manifest_sha256(records: list[tuple[str, int, str]]) -> str:
    """Hash filename/size/digest records with an unambiguous encoding."""

    digest = hashlib.sha256()
    for name, byte_size, file_sha256 in sorted(records):
        name_bytes = name.encode("utf-8")
        digest.update(len(name_bytes).to_bytes(8, byteorder="big"))
        digest.update(name_bytes)
        digest.update(byte_size.to_bytes(8, byteorder="big"))
        digest.update(bytes.fromhex(file_sha256))
    return digest.hexdigest()


def render_redactions_tsv(artifacts: list[PublishedArtifact]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(
        (
            "filename",
            "original_sha256",
            "public_sha256",
            "redaction_types",
            "redaction_count",
        )
    )
    for artifact in sorted(artifacts, key=lambda item: item.name):
        redaction_types = ";".join(
            f"{kind}:{count}" for kind, count in sorted(artifact.redactions.items())
        )
        writer.writerow(
            (
                artifact.name,
                artifact.original_sha256,
                artifact.public_sha256,
                redaction_types or "none",
                sum(artifact.redactions.values()),
            )
        )
    return buffer.getvalue().encode("utf-8")


def validate_paths(raw_dir: Path, public_dir: Path, private_alignment_path: str) -> None:
    if "\x00" in private_alignment_path or "\n" in private_alignment_path:
        fail("private alignment path contains a forbidden control character")
    if not private_alignment_path.startswith("/"):
        fail("private alignment path must be an absolute Unix path")
    if private_alignment_path.endswith("/"):
        fail("private alignment path must identify a file, not a directory")
    if raw_dir.is_symlink():
        fail("raw output directory must not be a symlink")
    if not raw_dir.is_dir():
        fail("raw output directory does not exist or is not a directory")
    if public_dir.exists() or public_dir.is_symlink():
        fail("public output directory already exists; overwrite is forbidden")
    if not public_dir.parent.is_dir():
        fail("parent of public output directory must already exist")
    if public_dir.parent.is_symlink():
        fail("parent of public output directory must not be a symlink")

    raw_resolved = raw_dir.resolve()
    parent_resolved = public_dir.parent.resolve()
    if raw_resolved == parent_resolved or raw_resolved in parent_resolved.parents:
        fail("public output directory must not be created inside the raw directory")


def collect_artifacts(
    raw_dir: Path, private_alignment_path: str
) -> tuple[list[PublishedArtifact], list[tuple[str, int, str]], dict[str, int]]:
    entries = sorted(os.scandir(raw_dir), key=lambda entry: entry.name)
    if not entries:
        fail("raw output directory is empty")

    copied: dict[str, tuple[str | None, bytes]] = {}
    raw_manifest: list[tuple[str, int, str]] = []
    omitted_iqtree: list[tuple[str, str]] = []
    omitted = {
        "iqtree_checkpoint": 0,
        "iqtree_model": 0,
        "wrapper_checksums": 0,
    }

    for entry in entries:
        if entry.is_symlink():
            fail(f"{entry.name}: symlinks are forbidden")
        mode = entry.stat(follow_symlinks=False).st_mode
        if not stat.S_ISREG(mode):
            fail(f"{entry.name}: only regular files are permitted")
        if not SAFE_PUBLISHED_NAME.fullmatch(entry.name):
            fail(f"{entry.name!r}: unsafe raw artifact filename")
        if PBS_DOTTED_ID.search(entry.name) or PBS_LOG_NAME.search(entry.name):
            fail("raw artifact filename contains a scheduler identifier")

        data = Path(entry.path).read_bytes()
        raw_manifest.append(
            (entry.name, len(data), hashlib.sha256(data).hexdigest())
        )

        if entry.name == OMITTED_CHECKSUM:
            omitted["wrapper_checksums"] += 1
            continue
        if entry.name.endswith(".ckp.gz"):
            omitted["iqtree_checkpoint"] += 1
            omitted_iqtree.append((entry.name, ".ckp.gz"))
            continue
        if entry.name.endswith(".model.gz"):
            omitted["iqtree_model"] += 1
            omitted_iqtree.append((entry.name, ".model.gz"))
            continue

        suffix = None if entry.name == SOFTWARE_VERSION else iqtree_suffix(entry.name)
        if entry.name != SOFTWARE_VERSION and suffix is None:
            fail(f"{entry.name}: unknown IQ-TREE output; explicit review is required")
        copied[entry.name] = (suffix, data)

    for suffix in CORE_SUFFIXES:
        matches = [name for name, (kind, _) in copied.items() if kind == suffix]
        if len(matches) != 1:
            fail(f"expected exactly one core {suffix} artifact, found {len(matches)}")
    if SOFTWARE_VERSION not in copied:
        fail(f"missing core artifact: {SOFTWARE_VERSION}")

    treefile_name = next(
        name for name, (suffix, _) in copied.items() if suffix == ".treefile"
    )
    prefix = treefile_name[: -len(".treefile")]
    for name, suffix in omitted_iqtree:
        if name != f"{prefix}{suffix}":
            fail(f"{name}: omitted artifact does not share the core IQ-TREE prefix")
    for name, (suffix, _) in copied.items():
        if name == SOFTWARE_VERSION:
            continue
        if name != f"{prefix}{suffix}":
            fail(f"{name}: artifact does not share the core IQ-TREE prefix")

    artifacts: list[PublishedArtifact] = []
    for name, (_, original) in sorted(copied.items()):
        if not original:
            fail(f"{name}: empty publication artifact")
        if is_text_artifact(name):
            public, redactions = sanitize_text(
                original,
                filename=name,
                private_alignment_path=private_alignment_path,
            )
        else:
            public, redactions = original, {}
        artifacts.append(PublishedArtifact(name, original, public, redactions))

    return artifacts, raw_manifest, omitted


def publish(raw_dir: Path, public_dir: Path, private_alignment_path: str) -> None:
    validate_paths(raw_dir, public_dir, private_alignment_path)
    artifacts, raw_manifest, omitted = collect_artifacts(
        raw_dir, private_alignment_path
    )

    redactions_tsv = render_redactions_tsv(artifacts)
    receipt = {
        "schema_version": 1,
        "raw_manifest_sha256": raw_manifest_sha256(raw_manifest),
        "raw_manifest_encoding": (
            "UTF-8 filename length (8-byte big-endian), filename bytes, "
            "byte size (8-byte big-endian), raw SHA-256 digest bytes; "
            "records sorted by filename"
        ),
        "raw_file_manifest": [
            {"filename": name, "byte_size": byte_size, "sha256": file_sha256}
            for name, byte_size, file_sha256 in sorted(raw_manifest)
        ],
        "raw_regular_file_count": len(raw_manifest),
        "published_file_count": len(artifacts),
        "redacted_file_count": sum(bool(item.redactions) for item in artifacts),
        "redacted_occurrence_count": sum(
            sum(item.redactions.values()) for item in artifacts
        ),
        "omitted_file_counts": omitted,
    }
    receipt_bytes = (json.dumps(receipt, indent=2, sort_keys=True) + "\n").encode(
        "utf-8"
    )

    # All validation and transformations happen before the output directory is
    # created.  Exclusive creation plus a previously nonexistent leaf prevents
    # accidental overwrite of a prior public result.
    public_dir.mkdir(mode=0o755)
    for artifact in artifacts:
        with (public_dir / artifact.name).open("xb") as handle:
            handle.write(artifact.public)
    with (public_dir / "redactions.tsv").open("xb") as handle:
        handle.write(redactions_tsv)
    with (public_dir / "raw-output-receipt.json").open("xb") as handle:
        handle.write(receipt_bytes)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", required=True, type=Path)
    parser.add_argument("--public-dir", required=True, type=Path)
    parser.add_argument(
        "--private-alignment-path",
        required=True,
        help="exact absolute alignment path to replace in IQ-TREE text outputs",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        publish(args.raw_dir, args.public_dir, args.private_alignment_path)
    except (OSError, PublicationError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print("Published a validated, redacted IQ-TREE artifact set.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
