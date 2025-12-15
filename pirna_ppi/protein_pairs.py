from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import comb
from pathlib import Path
from typing import Iterator

import typer
from loguru import logger
from tqdm import tqdm

try:
    from Bio import SeqIO
except ModuleNotFoundError as e:
    raise ModuleNotFoundError(
        "Missing dependency 'biopython'. Install it with: pixi add biopython (or pip install biopython)"
    ) from e

from pirna_ppi.config import INTERIM_DATA_DIR, RAW_DATA_DIR

app = typer.Typer(add_completion=False)


@dataclass(frozen=True, slots=True)
class ProteinRecord:
    index: int
    protein_id: str
    sequence: str


def read_fasta_records(path: Path) -> Iterator[ProteinRecord]:
    with path.open("r", encoding="utf-8") as f:
        for index, record in enumerate(SeqIO.parse(f, "fasta")):
            yield ProteinRecord(
                index=index,
                protein_id=str(record.id),
                sequence=str(record.seq).replace(" ", "").upper(),
            )


def write_pairs_fasta(records: list[ProteinRecord], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    n = len(records)
    pair_count = comb(n, 2) if n >= 2 else 0
    with output_path.open("w", encoding="utf-8") as w:
        for a, b in tqdm(combinations(records, 2), total=pair_count):
            pair_id = f"{a.protein_id}_{b.protein_id}"
            pair_seq = f"{a.sequence}:{b.sequence}"
            w.write(f">{pair_id}\n{pair_seq}\n")


@app.command()
def main(
    input_path: Path = RAW_DATA_DIR / "piRNA_PPI_protein.fasta",
    output_path: Path = INTERIM_DATA_DIR / "protein_pairs.fasta",
):
    if not input_path.exists():
        raise typer.BadParameter(f"Input FASTA not found: {input_path}")

    logger.info(f"Reading FASTA: {input_path}")
    records = list(read_fasta_records(input_path))

    if not records:
        raise typer.BadParameter(f"No FASTA records found in: {input_path}")

    empty_seq = [r for r in records if not r.sequence]
    if empty_seq:
        raise typer.BadParameter(
            f"Found {len(empty_seq)} records with empty sequence in: {input_path}"
        )

    n = len(records)
    pair_count = comb(n, 2) if n >= 2 else 0
    logger.info(f"Loaded {n} proteins. Pairs to write: {pair_count}")

    if pair_count == 0:
        raise typer.BadParameter("Need at least 2 protein records to create pairs.")

    if pair_count > 50_000_000:
        logger.warning(
            "The number of pairs is very large; this may take a long time and produce a huge file."
        )

    logger.info(f"Writing pairs to: {output_path}")

    write_pairs_fasta(records=records, output_path=output_path)

    logger.success("Protein pair generation complete.")


if __name__ == "__main__":
    app()
