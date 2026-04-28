import shutil
from pathlib import Path

import polars as pl
from loguru import logger
from tqdm import tqdm
import typer

from pirna_ppi.config import INTERIM_DATA_DIR

app = typer.Typer()

_DEFAULT_LOG = (
    INTERIM_DATA_DIR
    / "piRNA_PPI_20260210/pairs_msa_generation_output/index_fixed.tsv.log"
)
_DEFAULT_OUT = (
    INTERIM_DATA_DIR
    / "piRNA_PPI_20260210/pairs_msa_generation_output/paired_msa_filtered"
)


@app.command()
def main(
    log_path: Path = _DEFAULT_LOG,
    output_dir: Path = _DEFAULT_OUT,
    threshold: float = 0.3,
):
    df = pl.read_csv(log_path, separator="\t")
    filtered = df.filter(pl.col("Interaction_probability") > threshold)
    logger.info(
        f"Total entries: {len(df)}, filtered (>{threshold}): {len(filtered)}"
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    moved, skipped = 0, 0
    for src_str in tqdm(filtered["Input_MSA"], desc="Copying files"):
        src = Path(src_str)
        if not src.exists():
            logger.warning(f"File not found, skipped: {src}")
            skipped += 1
            continue
        shutil.copy2(str(src), output_dir / src.name)
        moved += 1

    logger.success(f"Done. Copied: {moved}, skipped: {skipped}")


if __name__ == "__main__":
    app()
