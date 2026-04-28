"""Parser for InterProScan output files."""

from pathlib import Path

import polars as pl

from ..models import Domain


def parse_iprscan(
    iprscan_path: Path,
    priority: list[str] | None = None,
) -> dict[str, list[Domain]]:
    """
    Parse InterProScan TSV output and return domains grouped by protein.

    Args:
        iprscan_path: Path to the iprscan output TSV file.
        priority: List of analysis methods in priority order (e.g., ['Pfam', 'Gene3D']).
                  If None, defaults to ['Pfam', 'Gene3D'].

    Returns:
        Dictionary mapping protein IDs to lists of Domain objects.
    """
    if priority is None:
        priority = ["Pfam", "Gene3D"]

    # iprscan TSV columns (standard format):
    # 0: protein_id, 1: md5, 2: length, 3: analysis, 4: signature_accession,
    # 5: signature_description, 6: start, 7: end, 8: score, 9: status, ...
    column_names = [
        "protein_id",
        "md5",
        "length",
        "analysis",
        "signature_accession",
        "signature_description",
        "start",
        "end",
        "score",
        "status",
        "date",
        "interpro_accession",
        "interpro_description",
        "go_annotations",
        "pathways",
    ]

    df = pl.read_csv(
        iprscan_path,
        separator="\t",
        has_header=False,
        new_columns=column_names,
        truncate_ragged_lines=True,
        quote_char=None,
    )

    # Exclude non-domain analysis methods from fallback
    _non_domain_analyses = {"MobiDBLite", "Phobius", "SignalP_EUK", "SignalP_GRAM_NEGATIVE",
                            "SignalP_GRAM_POSITIVE", "TMHMM", "Coils"}

    # Group by protein and collect domains
    protein_domains: dict[str, list[Domain]] = {}

    for protein_id in df["protein_id"].unique().to_list():
        protein_df = df.filter(pl.col("protein_id") == protein_id)

        # Try to get domains from highest priority analysis first
        domains = _extract_domains_by_priority(protein_df, priority)

        # Fallback: use all domain-type analyses if priority ones are absent
        if not domains:
            fallback_df = protein_df.filter(~pl.col("analysis").is_in(_non_domain_analyses))
            if len(fallback_df) > 0:
                domains = _extract_domains_from_df(fallback_df)

        if domains:
            domains = _merge_domains(domains)
            protein_domains[protein_id] = domains

    return protein_domains


def _extract_domains_by_priority(
    protein_df: pl.DataFrame, priority: list[str]
) -> list[Domain]:
    """Extract domains from the highest-priority analysis available."""
    for analysis in priority:
        analysis_df = protein_df.filter(pl.col("analysis") == analysis)
        if len(analysis_df) > 0:
            return _extract_domains_from_df(analysis_df)
    return []


def _extract_domains_from_df(df: pl.DataFrame) -> list[Domain]:
    """Extract Domain objects from a filtered iprscan DataFrame."""
    domains = []
    for row in df.iter_rows(named=True):
        start = int(row["start"])
        end = int(row["end"])
        name = row["signature_accession"]
        domains.append(Domain(name=name, ranges=[(start, end)]))
    return domains


def _merge_domains(domains: list[Domain]) -> list[Domain]:
    """
    Merge overlapping domains with the same name.

    Args:
        domains: List of Domain objects.

    Returns:
        List of merged Domain objects.
    """
    if not domains:
        return []

    # Group by name
    name_to_ranges: dict[str, list[tuple[int, int]]] = {}
    for domain in domains:
        if domain.name not in name_to_ranges:
            name_to_ranges[domain.name] = []
        name_to_ranges[domain.name].extend(domain.ranges)

    # Merge overlapping ranges for each name
    merged_domains = []
    for name, ranges in name_to_ranges.items():
        merged_ranges = _merge_ranges(ranges)
        merged_domains.append(Domain(name=name, ranges=merged_ranges))

    # Sort by start position
    merged_domains.sort(key=lambda d: d.start)

    return merged_domains


def _merge_ranges(ranges: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """
    Merge overlapping or adjacent ranges.

    Args:
        ranges: List of (start, end) tuples.

    Returns:
        List of merged (start, end) tuples.
    """
    if not ranges:
        return []

    # Sort by start position
    sorted_ranges = sorted(ranges, key=lambda x: x[0])

    merged = [sorted_ranges[0]]
    for start, end in sorted_ranges[1:]:
        last_start, last_end = merged[-1]
        # Merge if overlapping or adjacent (within 1 residue)
        if start <= last_end + 1:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))

    return merged
