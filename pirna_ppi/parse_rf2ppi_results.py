import re
from pathlib import Path

import polars as pl
from loguru import logger
from tqdm import tqdm
import typer

from pirna_ppi.config import INTERIM_DATA_DIR

app = typer.Typer()

# Pattern: ProteinA-segX_vs_ProteinB-segY.a3m
SEG_PATTERN = re.compile(r"^(.+)-seg(\d+)_vs_(.+)-seg(\d+)\.a3m$")


def parse_segment_name(msa_path: str) -> tuple[str, str, str, str]:
    """Extract protein names and segment IDs from MSA filename."""
    filename = Path(msa_path).name
    m = SEG_PATTERN.match(filename)
    if not m:
        raise ValueError(f"Cannot parse filename: {filename}")
    protein_a, seg_a, protein_b, seg_b = m.groups()
    return protein_a, seg_a, protein_b, seg_b


def parse_log_file(log_path: Path) -> pl.DataFrame:
    """Parse RF2-PPI log file into a Polars DataFrame with segment info."""
    df = pl.read_csv(log_path, separator="\t")

    records = []
    for row in df.iter_rows(named=True):
        msa_path = row["Input_MSA"]
        prob = row["Interaction_probability"]
        protein_a, seg_a, protein_b, seg_b = parse_segment_name(msa_path)
        records.append({
            "Protein_A": protein_a,
            "Protein_B": protein_b,
            "segment_a": int(seg_a),
            "segment_b": int(seg_b),
            "Interaction_probability": prob,
            "Input_MSA": msa_path,
        })

    return pl.DataFrame(records)


def build_protein_pair_summary(seg_df: pl.DataFrame) -> pl.DataFrame:
    """Aggregate segment-level results to protein pair level."""
    return (
        seg_df
        .group_by(["Protein_A", "Protein_B"])
        .agg(
            pl.col("Interaction_probability").max().alias("Max_Prob"),
            pl.len().alias("Num_segs"),
        )
        .sort(["Protein_A", "Protein_B"])
    )


def build_segment_detail(seg_df: pl.DataFrame) -> pl.DataFrame:
    """Build segment-level detail table."""
    return (
        seg_df
        .with_columns(
            (pl.col("Protein_A") + "-seg" + pl.col("segment_a").cast(pl.Utf8))
            .alias("segment1"),
            (pl.col("Protein_B") + "-seg" + pl.col("segment_b").cast(pl.Utf8))
            .alias("segment2"),
        )
        .select(["segment1", "segment2", "Interaction_probability"])
    )


@app.command()
def main(
    log_path: Path = typer.Option(
        INTERIM_DATA_DIR
        / "piRNA_PPI_20260210"
        / "segmented_paired_msa_index_fixed_550aa.tsv.log",
        help="RF2-PPI log file path",
    ),
    output_dir: Path = typer.Option(
        INTERIM_DATA_DIR / "piRNA_PPI_20260210" / "RF2-PPI_phrase_output",
        help="Output directory",
    ),
):
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Parse log ---
    logger.info(f"Parsing log file: {log_path}")
    seg_df = parse_log_file(log_path)
    logger.info(f"Parsed {len(seg_df)} segment pair records")

    # --- Protein pair summary ---
    pair_df = build_protein_pair_summary(seg_df)
    pair_output = output_dir / "protein_pair_summary.tsv"
    pair_df.write_csv(pair_output, separator="\t")
    logger.success(f"Protein pair summary ({len(pair_df)} pairs) -> {pair_output}")

    # --- Segment detail ---
    detail_df = build_segment_detail(seg_df)
    detail_output = output_dir / "segment_interaction_detail.tsv"
    detail_df.write_csv(detail_output, separator="\t")
    logger.success(f"Segment detail ({len(detail_df)} rows) -> {detail_output}")

    logger.success("All done.")


if __name__ == "__main__":
    app()
