"""Compare 4 ColabFold runs (baseline/fixed x AF2/AFmm) for the chain-split test.

Reads `outputs/{variant}_{tag}/` and prints a comparison table to verify whether
adding `#L1,L2\\t1,1` to the input a3m drives ColabFold to:
- expose the chain split in its output a3m header,
- (AFmm) produce a 2-chain PDB with a 200 aa residue-index gap,
- (both) push the inter-chain max contact off the sequence boundary.
"""

import json
from pathlib import Path
import pickle

import numpy as np
import polars as pl

from pirna_ppi.parse_af_models import _contact_prob_under, _read_a3m_chain_split

PAIR = "KWMTBOMO00189-seg1_vs_KWMTBOMO15589-seg2"
L1_KNOWN = 155
L2_KNOWN = 25
HERE = Path(__file__).parent
RUNS = [
    ("baseline", "AF2",  "alphafold2_ptm"),
    ("baseline", "AFmm", "alphafold2_multimer_v3"),
    ("fixed",    "AF2",  "alphafold2_ptm"),
    ("fixed",    "AFmm", "alphafold2_multimer_v3"),
]


def _find_rank1_pickle(out_dir: Path, arch: str) -> Path | None:
    cands = sorted(out_dir.glob(f"{PAIR}_all_rank_001_{arch}_model_*_seed_*.pickle"))
    return cands[0] if cands else None


def _find_rank1_pdb(out_dir: Path, arch: str) -> Path | None:
    # ColabFold may write either relaxed (rank_001_*) or unrelaxed pdbs.
    for pat in (f"{PAIR}_relaxed_rank_001_{arch}_*.pdb",
                f"{PAIR}_unrelaxed_rank_001_{arch}_*.pdb"):
        cands = sorted(out_dir.glob(pat))
        if cands:
            return cands[0]
    return None


def _pdb_chain_summary(pdb: Path) -> tuple[set[str], dict[str, list[int]]]:
    """Return (chain_ids, chain_id -> sorted unique residue numbers)."""
    chains: dict[str, set[int]] = {}
    with open(pdb) as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            cid = line[21]
            try:
                resnum = int(line[22:26])
            except ValueError:
                continue
            chains.setdefault(cid, set()).add(resnum)
    return set(chains), {k: sorted(v) for k, v in chains.items()}


def analyze_one(variant: str, tag: str, arch: str) -> dict:
    out_dir = HERE / "outputs" / f"{variant}_{tag}"
    rec: dict = {"variant": variant, "tag": tag}
    pkl = _find_rank1_pickle(out_dir, arch)
    a3m = out_dir / f"{PAIR}.a3m"
    rec["a3m_header"] = a3m.read_text().splitlines()[0] if a3m.exists() else "<missing>"
    split = _read_a3m_chain_split(a3m) if a3m.exists() else None
    rec["a3m_split"] = f"{split[0]},{split[1]}" if split else "(none)"
    if pkl is None:
        rec["status"] = "no_pickle"
        return rec
    d = pickle.load(open(pkl, "rb"))
    contact = _contact_prob_under(d["distogram"], 12.0)
    total = contact.shape[0]
    # Use authoritative split if available; else fall back to known L1+L2.
    l1, l2 = (split if split else (L1_KNOWN, L2_KNOWN))
    b_start = total - l2
    inter = contact[:l1, b_start:]
    rec["pickle_total"] = total
    rec["linker"] = total - l1 - l2
    rec["inter_max"] = float(inter.max())
    rec["inter_mean"] = float(inter.mean())
    flat = inter.ravel()
    k = max(1, flat.size // 100)
    rec["inter_top1pct_mean"] = float(np.sort(flat)[-k:].mean())
    rec["frac_gt_0.5"] = float((flat > 0.5).mean())
    border_val = float(contact[l1 - 1, b_start])
    rec["border_(L1-1, B0)"] = border_val
    i_max, j_max = np.unravel_index(inter.argmax(), inter.shape)
    rec["argmax_(i, j)"] = f"({int(i_max)},{int(j_max)})"
    rec["ptm"] = float(np.asarray(d.get("ptm", np.nan)))
    rec["ranking_confidence"] = float(np.asarray(d.get("ranking_confidence", np.nan)))
    pdb = _find_rank1_pdb(out_dir, arch)
    if pdb:
        chains, residues = _pdb_chain_summary(pdb)
        rec["pdb_chains"] = ",".join(sorted(chains))
        if "B" in chains and "A" in chains:
            gap = residues["B"][0] - residues["A"][-1]
            rec["A_last/B_first/gap"] = (
                f"{residues['A'][-1]}/{residues['B'][0]}/{gap}"
            )
        else:
            rec["A_last/B_first/gap"] = "(single chain)"
    else:
        rec["pdb_chains"] = "<no pdb>"
        rec["A_last/B_first/gap"] = "<no pdb>"
    return rec


def main():
    rows = [analyze_one(v, t, a) for v, t, a in RUNS]
    df = pl.DataFrame(rows)
    print("\n=== Chain-split test summary ===\n")
    with pl.Config(tbl_cols=-1, tbl_width_chars=200, fmt_str_lengths=80,
                   set_tbl_hide_dataframe_shape=True):
        print(df)
    out = HERE / "summary.json"
    out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
