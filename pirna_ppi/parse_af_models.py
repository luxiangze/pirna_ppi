import multiprocessing as mp
import os
from pathlib import Path
import pickle
import re

from loguru import logger
import numpy as np
import polars as pl
from tqdm import tqdm
import typer

from pirna_ppi.config import INTERIM_DATA_DIR, PROCESSED_DATA_DIR

app = typer.Typer()

# Match: <pair>_all_rank_<rank>_<arch>_model_<model>_seed_000.pickle
PICKLE_PATTERN = re.compile(
    r"^(?P<pair>.+)_all_rank_(?P<rank>\d{3})_"
    r"(?P<arch>alphafold2_(?:ptm|multimer_v3))"
    r"_model_(?P<model>\d)_seed_\d+\.pickle$"
)


def _softmax(x: np.ndarray, axis: int = -1) -> np.ndarray:
    x = x.astype(np.float32)
    x = x - x.max(axis=axis, keepdims=True)
    e = np.exp(x)
    return e / e.sum(axis=axis, keepdims=True)


def _contact_prob_under(distogram: dict, distance_cutoff: float) -> np.ndarray:
    """Compute per-residue-pair contact probability: sum of distogram bins under cutoff (Å)."""
    bin_edges = np.asarray(distogram["bin_edges"], dtype=np.float32)
    logits = np.asarray(distogram["logits"], dtype=np.float32)
    probs = _softmax(logits, axis=-1)
    # 64 bins, 63 edges. bin k upper edge = bin_edges[k] for k<63, else +inf.
    upper = np.concatenate([bin_edges, [np.inf]])
    mask = upper <= distance_cutoff
    return probs[..., mask].sum(axis=-1)


# Task tuple: (pkl_path_str, L1, L2_or_neg1, pair, network, rank, model_id)
# L2 == -1 means "derive L2 = total - L1 at parse time" (legacy single-chain a3m).
Task = tuple[str, int, int, str, str, int, int]


def _read_a3m_chain_split(a3m_path: Path) -> tuple[int, int] | None:
    """Parse `#L1,L2\t1,1` header from a ColabFold a3m. Return (L1, L2) or None."""
    try:
        with open(a3m_path, "r") as fh:
            first = fh.readline().strip()
    except OSError:
        return None
    if not first.startswith("#"):
        return None
    # Format: "#L1,L2\t1,1" (multi-chain) or "#T\t1" (single-chain).
    body = first[1:].split("\t", 1)[0]
    if "," not in body:
        return None
    try:
        l1, l2 = (int(x) for x in body.split(",", 1))
    except ValueError:
        return None
    return l1, l2


def _worker_init() -> None:
    """Limit BLAS/OpenMP threads inside each worker to avoid oversubscription."""
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[var] = "1"
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except Exception:  # noqa: BLE001
        pass


def _process_pickle(task: Task) -> dict | None:
    """Load one pickle and return per-model record.

    Chain B's true start may differ from L1 because ColabFold inserts a 9-G
    linker for `alphafold2_ptm` complexes. We use chain B = [total - L2, total).
    """
    pkl_path, l1, l2_known, pair, network, rank, model_id = task
    try:
        with open(pkl_path, "rb") as fh:
            d = pickle.load(fh)
        contact = _contact_prob_under(d["distogram"], distance_cutoff=12.0)
        total = contact.shape[0]
        l2 = l2_known if l2_known > 0 else total - l1
        if l1 <= 0 or l2 <= 0 or l1 + l2 > total:
            return {"__error__": f"Invalid L1={l1} L2={l2} total={total}: {pair}"}
        b_start = total - l2  # skip any linker between chain A end and chain B start
        inter = contact[:l1, b_start:]
        plddt = np.asarray(d["plddt"], dtype=np.float32)
        return {
            "pair": pair,
            "network": network,
            "rank": rank,
            "model_id": model_id,
            "L1": l1,
            "L2": l2,
            "linker": total - l1 - l2,
            "interaction_probability": float(inter.max()),
            "mean_plddt": float(plddt.mean()),
            "ptm": float(np.asarray(d.get("ptm", np.nan))),
            "ranking_confidence": float(np.asarray(d.get("ranking_confidence", np.nan))),
        }
    except Exception as exc:  # noqa: BLE001
        return {"__error__": f"{Path(pkl_path).name}: {exc}"}


def _load_pair_lengths(index_tsv: Path) -> dict[str, int]:
    """Map pair-name (basename without .a3m) -> L1 from index_fixed.tsv."""
    df = pl.read_csv(index_tsv, separator="\t", has_header=False, new_columns=["msa", "L1"])
    return {Path(p).stem: int(L1) for p, L1 in zip(df["msa"], df["L1"])}


def _collect_tasks(cf_dir: Path, l1_map: dict[str, int]) -> list[Task]:
    """Build per-pickle tasks.

    For each pair we first try to read L1/L2 from the ColabFold output
    `<pair>.a3m` (the `#L1,L2\t1,1` header is authoritative). If absent, we
    fall back to L1 from `index_fixed.tsv` and let the worker derive L2.
    """
    tasks: list[Task] = []
    missing_l1: set[str] = set()
    a3m_split_cache: dict[str, tuple[int, int] | None] = {}
    for pkl in cf_dir.glob("*_all_rank_*.pickle"):
        m = PICKLE_PATTERN.match(pkl.name)
        if not m:
            continue
        pair = m["pair"]
        # Per-pair a3m header (cached across the 5 ranks).
        if pair not in a3m_split_cache:
            a3m_split_cache[pair] = _read_a3m_chain_split(cf_dir / f"{pair}.a3m")
        split = a3m_split_cache[pair]
        if split is not None:
            l1, l2 = split
        elif pair in l1_map:
            l1, l2 = l1_map[pair], -1
        else:
            missing_l1.add(pair)
            continue
        arch = m["arch"]
        network = "AF2" if arch == "alphafold2_ptm" else "AFmm"
        tasks.append((str(pkl), l1, l2, pair, network, int(m["rank"]), int(m["model"])))
    if missing_l1:
        logger.warning(
            f"{len(missing_l1)} pairs without a3m chain-split header AND missing from index_fixed.tsv (skipped)"
        )
    return tasks


@app.command()
def main(
    cf_dir: Path = typer.Option(..., help="ColabFold output directory (AF2 or AFmm)"),
    output_tsv: Path = typer.Option(..., help="Output TSV path for per-model scores"),
    index_tsv: Path = typer.Option(
        INTERIM_DATA_DIR
        / "piRNA_PPI_20260210/pairs_msa_generation_output/index_fixed.tsv",
        help="index_fixed.tsv providing per-pair L1",
    ),
    n_workers: int = typer.Option(
        0, help="Parallel worker processes (0 = os.cpu_count())"
    ),
    chunksize: int = typer.Option(64, help="Tasks per IPC batch for imap_unordered"),
):
    if n_workers <= 0:
        n_workers = os.cpu_count() or 1
    logger.info(f"Loading L1 from {index_tsv}")
    l1_map = _load_pair_lengths(index_tsv)
    logger.info(f"Indexed {len(l1_map)} pairs")

    logger.info(f"Scanning pickles in {cf_dir}")
    tasks = _collect_tasks(cf_dir, l1_map)
    logger.info(f"Total tasks: {len(tasks)} | workers={n_workers} | chunksize={chunksize}")

    records: list[dict] = []
    n_fail = 0
    ctx = mp.get_context("fork")
    with ctx.Pool(processes=n_workers, initializer=_worker_init) as pool:
        for rec in tqdm(
            pool.imap_unordered(_process_pickle, tasks, chunksize=chunksize),
            total=len(tasks),
            desc="Parsing models",
        ):
            if rec is None:
                continue
            if "__error__" in rec:
                n_fail += 1
                if n_fail <= 20:
                    logger.warning(rec["__error__"])
                continue
            records.append(rec)
    if n_fail:
        logger.warning(f"{n_fail} tasks failed (first 20 logged)")

    df = pl.DataFrame(records).sort(["pair", "model_id"])
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(output_tsv, separator="\t")
    logger.success(f"Wrote {len(df)} rows -> {output_tsv}")


if __name__ == "__main__":
    app()


# Default output dirs (for reference in notebooks):
_DEFAULT_AF_OUT = PROCESSED_DATA_DIR / "piRNA_PPI_20260210" / "af_models"
