import json
from pathlib import Path

from loguru import logger
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import typer

app = typer.Typer()


def _per_pair_summary(
    af2: pl.DataFrame, afmm: pl.DataFrame, contact_threshold: float = 0.5
) -> pl.DataFrame:
    """Aggregate per-network model rows to per-pair stats and merge."""

    def agg(df: pl.DataFrame, prefix: str) -> pl.DataFrame:
        return df.group_by("pair").agg(
            pl.col("interaction_probability").max().alias(f"{prefix}_max"),
            pl.col("interaction_probability").mean().alias(f"{prefix}_mean"),
            (pl.col("interaction_probability") > contact_threshold)
            .sum()
            .alias(f"{prefix}_n_pass"),
            pl.col("mean_plddt").mean().alias(f"{prefix}_mean_plddt"),
            pl.col("ptm").max().alias(f"{prefix}_ptm_max"),
        )

    a = agg(af2, "af2")
    b = agg(afmm, "afmm")
    merged = a.join(b, on="pair", how="full", coalesce=True)
    merged = merged.with_columns(
        pl.max_horizontal("af2_max", "afmm_max").alias("max10"),
        (pl.col("af2_n_pass").fill_null(0) + pl.col("afmm_n_pass").fill_null(0))
        .alias("n_pass_total"),
    )
    return merged.sort("max10", descending=True)


def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return float("nan")
    rx = pl.Series(x).rank(method="average").to_numpy()
    ry = pl.Series(y).rank(method="average").to_numpy()
    return _pearson(rx, ry)


def _confusion_at(
    af2: np.ndarray, afmm: np.ndarray, threshold: float
) -> dict[str, int | float]:
    """Confusion matrix and Cohen's kappa for binary agreement at threshold."""
    a = af2 > threshold
    b = afmm > threshold
    tt = int(((a) & (b)).sum())
    ff = int(((~a) & (~b)).sum())
    tf = int(((a) & (~b)).sum())
    ft = int(((~a) & (b)).sum())
    n = tt + ff + tf + ft
    po = (tt + ff) / n if n else float("nan")
    pa = ((tt + tf) * (tt + ft) + (ff + tf) * (ff + ft)) / (n * n) if n else float("nan")
    kappa = (po - pa) / (1 - pa) if pa not in (1.0, float("nan")) else float("nan")
    return {
        "threshold": threshold,
        "both_positive": tt,
        "af2_only": tf,
        "afmm_only": ft,
        "both_negative": ff,
        "agreement": po,
        "cohen_kappa": kappa,
    }


def _plot_scatter(
    summary: pl.DataFrame, threshold: float, out_path: Path, pearson: float, spearman: float
) -> None:
    df = summary.drop_nulls(["af2_max", "afmm_max"])
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(df["af2_max"], df["afmm_max"], s=4, alpha=0.3, edgecolors="none")
    ax.plot([0, 1], [0, 1], color="grey", lw=0.8, ls="--")
    ax.axvline(threshold, color="red", lw=0.6, ls=":")
    ax.axhline(threshold, color="red", lw=0.6, ls=":")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("AF2 max interaction probability (5 models)")
    ax.set_ylabel("AFmm max interaction probability (5 models)")
    ax.set_title(f"AF2 vs AFmm (n={len(df)})\nPearson={pearson:.3f}  Spearman={spearman:.3f}")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _plot_max10_hist(summary: pl.DataFrame, threshold: float, out_path: Path) -> None:
    vals = summary["max10"].drop_nulls().to_numpy()
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(vals, bins=50, color="steelblue", edgecolor="white")
    ax.axvline(threshold, color="red", lw=1.0, ls="--", label=f"threshold={threshold}")
    ax.set_xlabel("max interaction probability across 10 models")
    ax.set_ylabel("# pairs")
    ax.set_title(f"Distribution of max10 (n={len(vals)})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _plot_agreement_heatmap(conf: dict, out_path: Path) -> None:
    mat = np.array(
        [[conf["both_negative"], conf["afmm_only"]], [conf["af2_only"], conf["both_positive"]]]
    )
    fig, ax = plt.subplots(figsize=(4.5, 4))
    im = ax.imshow(mat, cmap="Blues")
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["AFmm ≤ thr", "AFmm > thr"])
    ax.set_yticklabels(["AF2 ≤ thr", "AF2 > thr"])
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f"{mat[i, j]:,}", ha="center", va="center", color="black")
    ax.set_title(
        f"Agreement at thr={conf['threshold']}\n"
        f"agreement={conf['agreement']:.3f}, kappa={conf['cohen_kappa']:.3f}"
    )
    fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


@app.command()
def main(
    af2_tsv: Path = typer.Option(..., help="Per-model TSV for AF2"),
    afmm_tsv: Path = typer.Option(..., help="Per-model TSV for AFmm"),
    output_dir: Path = typer.Option(..., help="Directory to write outputs"),
    threshold: float = typer.Option(0.6, help="Cutoff for max10 filter (paper default 0.6)"),
):
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(exist_ok=True)

    logger.info(f"Polars thread pool size: {pl.thread_pool_size()}")
    logger.info(f"Loading {af2_tsv}")
    af2 = pl.read_csv(af2_tsv, separator="\t")
    logger.info(f"Loading {afmm_tsv}")
    afmm = pl.read_csv(afmm_tsv, separator="\t")

    # Per-pair summary (intersect on pairs present in both networks for stats)
    summary = _per_pair_summary(af2, afmm)
    summary_path = output_dir / "pair_summary.tsv"
    summary.write_csv(summary_path, separator="\t")
    logger.success(f"Pair summary ({len(summary)} pairs) -> {summary_path}")

    # Stats only on pairs present in both networks
    both = summary.drop_nulls(["af2_max", "afmm_max"])
    x = both["af2_max"].to_numpy()
    y = both["afmm_max"].to_numpy()
    pearson = _pearson(x, y)
    spearman = _spearman(x, y)
    conf_05 = _confusion_at(x, y, 0.5)
    conf_06 = _confusion_at(x, y, 0.6)

    stats = {
        "n_pairs_total": len(summary),
        "n_pairs_in_both": len(both),
        "pearson": pearson,
        "spearman": spearman,
        "agreement_at_0.5": conf_05,
        "agreement_at_0.6": conf_06,
        "filter_threshold": threshold,
        "n_kept": int((summary["max10"] > threshold).sum()),
        "n_excluded": int((summary["max10"] <= threshold).sum()),
    }
    stats_path = output_dir / "consistency_stats.json"
    stats_path.write_text(json.dumps(stats, indent=2))
    logger.success(f"Consistency stats -> {stats_path}")

    # Filtered subsets
    kept = summary.filter(pl.col("max10") > threshold)
    excluded = summary.filter(pl.col("max10") <= threshold)
    kept.write_csv(output_dir / "pairs_kept.tsv", separator="\t")
    excluded.write_csv(output_dir / "pairs_excluded.tsv", separator="\t")
    logger.success(f"Kept: {len(kept)}, Excluded: {len(excluded)}")

    # Plots
    _plot_scatter(summary, threshold, fig_dir / "scatter_af2_vs_afmm_max.png", pearson, spearman)
    _plot_max10_hist(summary, threshold, fig_dir / "hist_max10.png")
    _plot_agreement_heatmap(conf_05, fig_dir / "agreement_heatmap_0.5.png")
    _plot_agreement_heatmap(conf_06, fig_dir / "agreement_heatmap_0.6.png")
    logger.success(f"Figures -> {fig_dir}")


if __name__ == "__main__":
    app()
