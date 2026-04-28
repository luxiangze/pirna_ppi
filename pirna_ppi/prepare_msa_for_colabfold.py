"""Prepare RF2-PPI paired MSA a3m files for ColabFold paired-MSA mode.

Patches each `<pair>.a3m` by prepending a `#L1,L2\\t1,1` chain-split header so
ColabFold knows the input is two chains and applies a 200 aa residue-index gap
(and, for `alphafold2_ptm`, also inserts its standard 9-G linker).

Body is left untouched (the RF2-PPI paired MSA already concatenates each homolog
as `[chain_A][chain_B]` of length L1+L2, which is exactly what ColabFold expects
for paired alignments).
"""

import multiprocessing as mp
import os
from pathlib import Path

from loguru import logger
import polars as pl
from tqdm import tqdm
import typer

from pirna_ppi.config import INTERIM_DATA_DIR

app = typer.Typer()


def _load_pair_lengths(index_tsv: Path) -> dict[str, int]:
    """Map pair-name (basename without .a3m) -> L1 from index_fixed.tsv."""
    df = pl.read_csv(index_tsv, separator="\t", has_header=False, new_columns=["msa", "L1"])
    return {Path(p).stem: int(L1) for p, L1 in zip(df["msa"], df["L1"])}


def _read_first_query_len(a3m_path: Path) -> int:
    """Return the length (in chars, gaps included) of the first sequence entry.

    a3m structure: lines starting with '>' are headers; subsequent non-'>' lines
    until the next '>' or EOF concatenate to one sequence.
    """
    seq_chars = 0
    seen_first_header = False
    with open(a3m_path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                if seen_first_header:
                    break  # second header => first sequence finished
                seen_first_header = True
                continue
            seq_chars += len(line.rstrip("\n"))
    return seq_chars


def _patch_one(args: tuple[str, str, int]) -> dict | None:
    """Patch a single a3m: prepend '#L1,L2\\t1,1' header. Returns audit record or error."""
    src, dst, l1 = args
    src_path = Path(src)
    dst_path = Path(dst)
    try:
        total = _read_first_query_len(src_path)
        if total <= l1:
            return {"__error__": f"{src_path.name}: query len {total} <= L1 {l1}"}
        l2 = total - l1
        dst_path.parent.mkdir(parents=True, exist_ok=True)
        # Stream copy with prepended header (avoids loading full file into memory).
        with open(src_path, "rb") as fin, open(dst_path, "wb") as fout:
            fout.write(f"#{l1},{l2}\t1,1\n".encode())
            # Skip an existing leading '#' header line if present (idempotent re-runs).
            first_chunk = fin.read(4096)
            if first_chunk.startswith(b"#"):
                nl = first_chunk.find(b"\n")
                if nl >= 0:
                    first_chunk = first_chunk[nl + 1 :]
            fout.write(first_chunk)
            while chunk := fin.read(1 << 20):
                fout.write(chunk)
        return {"pair": src_path.stem, "L1": l1, "L2": l2}
    except Exception as exc:  # noqa: BLE001
        return {"__error__": f"{src_path.name}: {exc}"}


def _worker_init() -> None:
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[var] = "1"


@app.command()
def main(
    input_dir: Path = typer.Option(..., help="Directory of RF2-PPI paired MSA a3m files"),
    output_dir: Path = typer.Option(..., help="Directory to write patched a3m files"),
    index_tsv: Path = typer.Option(
        INTERIM_DATA_DIR
        / "piRNA_PPI_20260210/pairs_msa_generation_output/index_fixed.tsv",
        help="index_fixed.tsv providing per-pair L1",
    ),
    n_workers: int = typer.Option(
        0, help="Parallel worker processes (0 = os.cpu_count())"
    ),
    chunksize: int = typer.Option(64, help="Tasks per IPC batch"),
):
    if n_workers <= 0:
        n_workers = os.cpu_count() or 1
    logger.info(f"Loading L1 from {index_tsv}")
    l1_map = _load_pair_lengths(index_tsv)
    logger.info(f"Indexed {len(l1_map)} pairs")

    a3m_files = sorted(input_dir.glob("*.a3m"))
    logger.info(f"Found {len(a3m_files)} a3m files in {input_dir}")

    tasks: list[tuple[str, str, int]] = []
    missing: list[str] = []
    for f in a3m_files:
        if f.stem not in l1_map:
            missing.append(f.stem)
            continue
        tasks.append((str(f), str(output_dir / f.name), l1_map[f.stem]))
    if missing:
        logger.warning(f"{len(missing)} a3m files missing in index_fixed.tsv (skipped)")

    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Patching {len(tasks)} files | workers={n_workers} | chunksize={chunksize}")

    n_ok, n_fail = 0, 0
    ctx = mp.get_context("fork")
    with ctx.Pool(processes=n_workers, initializer=_worker_init) as pool:
        for rec in tqdm(
            pool.imap_unordered(_patch_one, tasks, chunksize=chunksize),
            total=len(tasks),
            desc="Patching a3m",
        ):
            if rec is None:
                continue
            if "__error__" in rec:
                n_fail += 1
                if n_fail <= 20:
                    logger.warning(rec["__error__"])
            else:
                n_ok += 1
    logger.success(f"Patched {n_ok} / {len(tasks)} files (failed: {n_fail}) -> {output_dir}")


if __name__ == "__main__":
    app()
