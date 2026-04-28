"""Repeat RF2-PPI prediction N times and analyze Max_Prob > threshold distribution."""

import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl
from loguru import logger
import typer

from pirna_ppi.config import INTERIM_DATA_DIR, PROCESSED_DATA_DIR
from pirna_ppi.parse_rf2ppi_results import parse_log_file, build_protein_pair_summary

app = typer.Typer()

_INDEX_FILE = (
    INTERIM_DATA_DIR
    / "piRNA_PPI_20260210"
    / "pairs_msa_generation_output"
    / "known_ppi_index_fixed.tsv"
)
_DEFAULT_OUTPUT_DIR = PROCESSED_DATA_DIR / "piRNA_PPI_20260210" / "known_ppi_rf2ppi_repeat"

_SIF = "/home/gyk/sifs/SE3nv-20230612.sif"
_RF2PPI_SCRIPT = "/home/RoseTTAFold2-PPI/src/predict_list_PPI.py"
_MODEL_FILE = "/home/RoseTTAFold2-PPI/src/models/RF2-PPI.pt"
_DATA_BIND = "/home/gyk/projects/pirna_ppi/data:/work/users"
_CODE_BIND = "/home/gyk/projects/pirna_ppi/RoseTTAFold2-PPI:/home/RoseTTAFold2-PPI"
_INDEX_IN_CONTAINER = "/work/users/interim/piRNA_PPI_20260210/pairs_msa_generation_output/known_ppi_index_fixed.tsv"


def _run_rf2ppi(index_file: Path) -> None:
    """Run RF2-PPI prediction via singularity."""
    cmd = [
        "singularity", "exec",
        "--bind", _DATA_BIND,
        "--bind", _CODE_BIND,
        "--nv", _SIF,
        "/bin/bash", "-c",
        f"cd /work/users && python {_RF2PPI_SCRIPT} -list_fn {_INDEX_IN_CONTAINER} -model_file {_MODEL_FILE}",
    ]
    logger.info("Running RF2-PPI: {}", " ".join(cmd))
    subprocess.run(cmd, check=True)


def _parse_results(log_path: Path, run_output_dir: Path) -> pl.DataFrame:
    """Parse RF2-PPI log and save per-run results."""
    run_output_dir.mkdir(parents=True, exist_ok=True)
    seg_df = parse_log_file(log_path)
    pair_df = build_protein_pair_summary(seg_df)
    pair_df.write_csv(run_output_dir / "protein_pair_summary.tsv", separator="\t")
    return pair_df


def _plot_distribution(counts: list[int], threshold: float, output_dir: Path) -> None:
    """Plot line chart of Max_Prob > threshold count across runs."""
    fig, ax = plt.subplots(figsize=(10, 5))
    runs = list(range(1, len(counts) + 1))
    ax.plot(runs, counts, marker="o", linewidth=1.5, markersize=4)
    ax.set_xlabel("Run")
    ax.set_ylabel(f"Pairs with Max_Prob > {threshold}")
    ax.set_title(f"RF2-PPI repeat runs: pairs with Max_Prob > {threshold}")
    ax.set_xticks(runs)
    ax.grid(True, linestyle="--", alpha=0.5)
    fig.tight_layout()
    out_path = output_dir / "max_prob_distribution.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    logger.success("Saved distribution plot -> {}", out_path)


@app.command()
def main(
    n_repeats: int = typer.Option(50, help="Number of repeat runs"),
    index_file: Path = typer.Option(_INDEX_FILE, help="Index TSV passed to RF2-PPI"),
    log_path: Path = typer.Option(
        _INDEX_FILE.parent / "known_ppi_index_fixed.tsv.log",
        help="RF2-PPI output log file (overwritten each run)",
    ),
    output_dir: Path = typer.Option(_DEFAULT_OUTPUT_DIR, help="Root output directory"),
    threshold: float = typer.Option(0.7, help="Max_Prob threshold for counting"),
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    counts: list[int] = []

    for run_idx in range(1, n_repeats + 1):
        logger.info("=== Run {}/{} ===", run_idx, n_repeats)
        run_dir = output_dir / f"run_{run_idx:03d}"

        # Step 1: Run RF2-PPI
        _run_rf2ppi(index_file)

        # Step 2: Parse results
        pair_df = _parse_results(log_path, run_dir)

        # Count pairs above threshold
        n_above = pair_df.filter(pl.col("Max_Prob") > threshold).height
        counts.append(n_above)
        logger.info("Run {}: {} pairs with Max_Prob > {}", run_idx, n_above, threshold)

    # Aggregate summary
    summary_df = pl.DataFrame({
        "run": list(range(1, n_repeats + 1)),
        f"n_above_{threshold}".replace(".", "_"): counts,
    })
    summary_path = output_dir / "repeat_summary.tsv"
    summary_df.write_csv(summary_path, separator="\t")
    logger.success("Summary table -> {}", summary_path)

    # Plot
    _plot_distribution(counts, threshold, output_dir)
    logger.success("All {} runs complete.", n_repeats)


if __name__ == "__main__":
    app()
