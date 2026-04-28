"""Extract index entries for known PPI pairs from index_fixed.tsv."""

import re
from pathlib import Path

from loguru import logger
import typer

from pirna_ppi.config import INTERIM_DATA_DIR

app = typer.Typer()

_BASE_DIR = INTERIM_DATA_DIR / "piRNA_PPI_20260210" / "pairs_msa_generation_output"
_SEG_RE = re.compile(r"^(.+)-seg\d+$")


def _base_id(pid: str) -> str:
    m = _SEG_RE.match(pid)
    return m.group(1) if m else pid


@app.command()
def main(
    known_ppi_file: Path = typer.Option(_BASE_DIR / "known_ppi.tsv", help="TSV with known PPI pairs (p1, p2)"),
    index_file: Path = typer.Option(_BASE_DIR / "index_fixed.tsv", help="Full index TSV to filter"),
    output_file: Path = typer.Option(_BASE_DIR / "known_ppi_index_fixed.tsv", help="Output filtered index TSV"),
) -> None:
    # Load known PPI pairs (protein IDs without seg suffix)
    known_pairs: set[tuple[str, str]] = set()
    with known_ppi_file.open() as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                known_pairs.add((parts[0], parts[1]))
    logger.info("Loaded {} known PPI pairs from {}", len(known_pairs), known_ppi_file)

    # Filter index entries whose filename matches a known pair (seg-aware)
    matched = 0
    with index_file.open() as fin, output_file.open("w") as fout:
        for line in fin:
            path_str = line.split("\t")[0]
            stem = Path(path_str).stem  # e.g. KWMTBOMO02358-seg1_vs_KWMTBOMO05469-seg1
            parts = stem.split("_vs_")
            if len(parts) != 2:
                continue
            p1_base = _base_id(parts[0])
            p2_base = _base_id(parts[1])
            if (p1_base, p2_base) in known_pairs or (p2_base, p1_base) in known_pairs:
                fout.write(line)
                matched += 1

    logger.success("Matched {} index entries -> {}", matched, output_file)


if __name__ == "__main__":
    app()
