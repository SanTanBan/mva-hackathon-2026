"""Small helpers shared by the pipeline scripts."""
import gzip
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def open_text(path):
    """Open a plain or gzip/bgzip text file for streaming."""
    path = Path(path)
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def assert_outside_repo(path):
    """Patient-derived output must never be written inside the git repository."""
    resolved = Path(path).resolve()
    if resolved == REPO or REPO in resolved.parents:
        raise SystemExit(f"Refusing to write patient-derived output inside the repository: {resolved}")
    return resolved


def add_chr(chrom):
    """'15' -> 'chr15', 'MT' -> 'chrM'; already-prefixed names are left alone."""
    if chrom.startswith("chr"):
        return chrom
    return "chrM" if chrom in ("M", "MT") else "chr" + chrom


def is_autosome(chrom):
    c = chrom[3:] if chrom.startswith("chr") else chrom
    return c.isdigit() and 1 <= int(c) <= 22
