#!/usr/bin/env python3
"""Verify every SHA-256 entry in MANIFEST_SHA256.txt."""

from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "MANIFEST_SHA256.txt"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    failures: list[str] = []
    checked = 0
    for raw in MANIFEST.read_text(encoding="utf-8").splitlines():
        if not raw or raw.startswith("#"):
            continue
        expected, relative = raw.split("  ", 1)
        path = ROOT / relative
        if not path.is_file():
            failures.append(f"MISSING  {relative}")
            continue
        actual = digest(path)
        checked += 1
        if actual != expected:
            failures.append(f"MISMATCH {relative}: {actual}")
    if failures:
        raise SystemExit("\n".join(failures))
    print(f"PASS: {checked} files match MANIFEST_SHA256.txt")


if __name__ == "__main__":
    main()
