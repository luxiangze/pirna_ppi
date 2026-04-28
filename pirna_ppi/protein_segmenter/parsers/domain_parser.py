"""Parser for DPAM domain files."""

from pathlib import Path

import polars as pl

from ..models import Domain


def parse_range_string(range_str: str) -> list[tuple[int, int]]:
    """
    Parse a range string like '1-215,266-275,331-340' into list of tuples.

    Args:
        range_str: Range string in format 'start-end,start-end,...'

    Returns:
        List of (start, end) tuples.
    """
    ranges = []
    for part in range_str.split(","):
        if "-" in part:
            start, end = part.split("-")
            ranges.append((int(start), int(end)))
        else:
            # Single residue
            pos = int(part)
            ranges.append((pos, pos))
    return ranges


def parse_dpam_domains(domains_path: Path) -> dict[str, list[Domain]]:
    """
    Parse DPAM domains file and return domains grouped by protein.

    Args:
        domains_path: Path to the DPAM domains file.

    Returns:
        Dictionary mapping protein IDs to lists of Domain objects.
    """
    df = pl.read_csv(domains_path, separator="\t", null_values=["na"])

    # Filter out low-quality domains
    df = df.filter(~pl.col("Judge").is_in(["low_confidence", "partial_domain"]))

    protein_domains: dict[str, list[Domain]] = {}

    for row in df.iter_rows(named=True):
        protein_id = row["Protein"]
        domain_name = row["Domain"]
        range_str = row["Range"]

        ranges = parse_range_string(range_str)
        domain = Domain(name=domain_name, ranges=ranges)

        if protein_id not in protein_domains:
            protein_domains[protein_id] = []
        protein_domains[protein_id].append(domain)

    return protein_domains


def parse_failed_list(failed_path: Path) -> set[str]:
    """
    Parse the DPAM failed list file.

    Args:
        failed_path: Path to the failed list file.

    Returns:
        Set of protein IDs that failed DPAM prediction.
    """
    failed_proteins = set()
    with open(failed_path) as f:
        for line in f:
            line = line.strip()
            if line:
                # Format: "KWMTBOMO13492\tPDB contains NaN"
                protein_id = line.split("\t")[0]
                if protein_id:
                    failed_proteins.add(protein_id)
    return failed_proteins
