#!/usr/bin/env python3
"""Provenance report: database checksums + tool/pipeline versions for a run."""
import argparse
import hashlib
import json
import os
from pathlib import Path


def md5sum(path, chunk=8 << 20):
    # ponytail: md5 is the real reproducibility guarantee; slow on multi-GB gz.
    #           If throughput ever matters, fall back to size+mtime here.
    h = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description="Write a run provenance report.")
    parser.add_argument("--db", nargs="*", default=[], help="Database files to checksum")
    parser.add_argument("--versions", help="Collected versions.yml (tool versions)")
    parser.add_argument("--workflow-json", help="JSON blob of workflow.* metadata")
    parser.add_argument("--output", required=True, help="Output provenance report")
    args = parser.parse_args()

    lines = ["# Provenance report", ""]

    if args.workflow_json:
        meta = json.loads(Path(args.workflow_json).read_text())
        lines.append("## Run")
        for key in sorted(meta):
            lines.append(f"{key}: {meta[key]}")
        lines.append("")

    lines.append("## Databases")
    for db in args.db:
        p = Path(db)
        lines.append(
            f"{p.name}\tsize_bytes={os.path.getsize(p)}\tmd5={md5sum(p)}"
        )
    lines.append("")

    if args.versions and Path(args.versions).exists():
        lines.append("## Tool versions")
        lines.append(Path(args.versions).read_text().rstrip())
        lines.append("")

    Path(args.output).write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
