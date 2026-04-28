from __future__ import annotations

import csv
from dataclasses import dataclass
from itertools import combinations
from math import comb
from pathlib import Path

import typer
from loguru import logger
from tqdm import tqdm

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
except ModuleNotFoundError as e:
    raise ModuleNotFoundError(
        "Missing dependency 'biopython'. Install it with: pixi add biopython (or pip install biopython)"
    ) from e

from pirna_ppi.config import INTERIM_DATA_DIR

app = typer.Typer(add_completion=False)

# Input FASTA files from segmenter output
INPUT_FASTA_NAMES: list[str] = [
    "segmented.fasta",
    "unsegmented.fasta",
    "no_domain.fasta",
    "oversized_segments.fasta",
]

# Output file names
PAIRS_FASTA = "protein_pairs.fasta"
PASSED_FASTA = "passed_sequences.fasta"
FAILED_FASTA = "failed_sequences.fasta"
SUMMARY_TSV = "summary.tsv"
LOG_FILE = "protein_pairs.log"


@dataclass(frozen=True, slots=True)
class ProteinRecord:
    index: int
    protein_id: str
    sequence: str
    source_file: str

    @property
    def base_protein_id(self) -> str:
        """Extract base protein ID (strip '-segN' suffix)."""
        return self.protein_id.rsplit("-seg", 1)[0]

    @property
    def length(self) -> int:
        return len(self.sequence)


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_fasta_records(path: Path, start_index: int = 0) -> list[ProteinRecord]:
    """Read all records from a single FASTA file."""
    records: list[ProteinRecord] = []
    with path.open("r", encoding="utf-8") as f:
        for i, rec in enumerate(SeqIO.parse(f, "fasta")):
            records.append(
                ProteinRecord(
                    index=start_index + i,
                    protein_id=str(rec.id),
                    sequence=str(rec.seq).replace(" ", "").upper(),
                    source_file=path.name,
                )
            )
    return records


def load_all_records(input_dir: Path) -> list[ProteinRecord]:
    """Load records from all input FASTA files in *input_dir*."""
    all_records: list[ProteinRecord] = []
    for name in INPUT_FASTA_NAMES:
        fasta_path = input_dir / name
        if not fasta_path.exists():
            logger.warning(f"FASTA not found, skipped: {fasta_path}")
            continue
        recs = read_fasta_records(fasta_path, start_index=len(all_records))
        logger.info(f"  {name}: {len(recs)} records")
        all_records.extend(recs)
    return all_records


def filter_by_length(
    records: list[ProteinRecord], max_length: int
) -> tuple[list[ProteinRecord], list[ProteinRecord]]:
    """Split records into (passed, failed) based on *max_length*."""
    passed, failed = [], []
    for r in records:
        (passed if r.length <= max_length else failed).append(r)
    return passed, failed


def write_records_fasta(records: list[ProteinRecord], path: Path) -> None:
    """Write ProteinRecords as a FASTA file."""
    bio_records = [
        SeqRecord(Seq(r.sequence), id=r.protein_id, description=f"source={r.source_file}")
        for r in records
    ]
    with path.open("w", encoding="utf-8") as f:
        SeqIO.write(bio_records, f, "fasta")


def write_pairs_fasta(records: list[ProteinRecord], output_path: Path) -> int:
    """Write all valid pairs (skip same-protein segments) and return count."""
    n = len(records)
    total_comb = comb(n, 2) if n >= 2 else 0
    written = 0
    with output_path.open("w", encoding="utf-8") as w:
        for a, b in tqdm(combinations(records, 2), total=total_comb):
            if a.base_protein_id == b.base_protein_id:
                continue
            pair_id = f"{a.protein_id}_{b.protein_id}"
            pair_seq = f"{a.sequence}:{b.sequence}"
            w.write(f">{pair_id}\n{pair_seq}\n")
            written += 1
    return written


def write_summary(
    output_dir: Path,
    all_records: list[ProteinRecord],
    passed: list[ProteinRecord],
    failed: list[ProteinRecord],
    max_length: int,
    n_pairs: int,
) -> None:
    """Write a TSV summary table."""
    path = output_dir / SUMMARY_TSV

    def _stats(recs: list[ProteinRecord]) -> dict:
        if not recs:
            return {"n_entries": 0, "n_proteins": 0, "len_min": 0, "len_max": 0, "len_mean": 0.0}
        lengths = [r.length for r in recs]
        return {
            "n_entries": len(recs),
            "n_proteins": len({r.base_protein_id for r in recs}),
            "len_min": min(lengths),
            "len_max": max(lengths),
            "len_mean": round(sum(lengths) / len(lengths), 1),
        }

    rows = [
        {"category": "all_input", **_stats(all_records)},
        {"category": f"passed (len<={max_length})", **_stats(passed)},
        {"category": f"failed (len>{max_length})", **_stats(failed)},
        {"category": "pairs_written", "n_entries": n_pairs, "n_proteins": "", "len_min": "", "len_max": "", "len_mean": ""},
    ]
    fieldnames = ["category", "n_entries", "n_proteins", "len_min", "len_max", "len_mean"]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

@app.command()
def main(
    input_dir: Path = INTERIM_DATA_DIR / "piRNA_PPI_20260210" / "proteins_segment_output",
    output_dir: Path = INTERIM_DATA_DIR / "piRNA_PPI_20260210" / "protein_pairs_output",
    max_length: int = 600,
):
    output_dir.mkdir(parents=True, exist_ok=True)

    # Setup file logger
    log_path = output_dir / LOG_FILE
    log_id = logger.add(log_path, format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}", level="DEBUG")

    try:
        if not input_dir.is_dir():
            raise typer.BadParameter(f"Input directory not found: {input_dir}")

        logger.info(f"Input directory : {input_dir}")
        logger.info(f"Output directory: {output_dir}")
        logger.info(f"Max length      : {max_length}")

        # 1. Load all records
        logger.info("Loading FASTA records ...")
        all_records = load_all_records(input_dir)
        if not all_records:
            raise typer.BadParameter("No FASTA records found in input directory.")

        empty_seq = [r for r in all_records if not r.sequence]
        if empty_seq:
            raise typer.BadParameter(f"Found {len(empty_seq)} records with empty sequence.")

        n_proteins = len({r.base_protein_id for r in all_records})
        logger.info(f"Total: {len(all_records)} segments from {n_proteins} proteins")

        # 2. Filter by length
        passed, failed = filter_by_length(all_records, max_length)
        logger.info(f"Passed length filter (<={max_length}): {len(passed)} segments")
        logger.info(f"Failed length filter (>{max_length}) : {len(failed)} segments")

        # 3. Write passed / failed FASTA
        write_records_fasta(passed, output_dir / PASSED_FASTA)
        logger.info(f"Written: {PASSED_FASTA}")
        write_records_fasta(failed, output_dir / FAILED_FASTA)
        logger.info(f"Written: {FAILED_FASTA}")

        # 4. Generate pairs
        if len(passed) < 2:
            raise typer.BadParameter("Need at least 2 passed records to create pairs.")

        logger.info("Generating pairs ...")
        n_pairs = write_pairs_fasta(passed, output_dir / PAIRS_FASTA)
        logger.info(f"Written: {PAIRS_FASTA} ({n_pairs} pairs)")

        # 5. Summary
        write_summary(output_dir, all_records, passed, failed, max_length, n_pairs)
        logger.info(f"Written: {SUMMARY_TSV}")

        logger.success(f"Done — {n_pairs} pairs written, same-protein pairs skipped.")
    finally:
        logger.remove(log_id)


if __name__ == "__main__":
    app()
