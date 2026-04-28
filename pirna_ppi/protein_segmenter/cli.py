"""Command-line interface for protein segmenter."""

from pathlib import Path

from loguru import logger
from tqdm import tqdm
import typer

from pirna_ppi.config import INTERIM_DATA_DIR, RAW_DATA_DIR

from .config import SegmenterConfig
from .core.segmenter import ProteinSegmenter
from .models import Segment
from .parsers.domain_parser import parse_dpam_domains
from .parsers.fasta_parser import parse_fasta
from .parsers.iprscan_parser import parse_iprscan

app = typer.Typer()


def write_segments_to_fasta(segments: list[Segment], output_path: Path) -> None:
    """
    Write segments to a single FASTA file.

    Args:
        segments: List of Segment objects.
        output_path: Path to output FASTA file.
    """
    with open(output_path, "w") as f:
        for seg in segments:
            f.write(seg.to_fasta() + "\n")


def _write_summary(summary_path: Path, fasta_groups: dict[str, list[Segment]]) -> None:
    """
    Write a summary TSV with statistics for each FASTA output file.

    Columns: file, n_entries, n_proteins, len_min, len_max, len_mean, len_median
    """
    with open(summary_path, "w") as f:
        f.write("file\tn_entries\tn_proteins\tlen_min\tlen_max\tlen_mean\tlen_median\n")
        for filename, segs in fasta_groups.items():
            if not segs:
                f.write(f"{filename}\t0\t0\t-\t-\t-\t-\n")
                continue
            lengths = sorted(len(s.sequence) for s in segs)
            n_entries = len(lengths)
            n_proteins = len({s.protein_id for s in segs})
            mid = n_entries // 2
            median = (lengths[mid] if n_entries % 2 else (lengths[mid - 1] + lengths[mid]) / 2)
            f.write(
                f"{filename}\t{n_entries}\t{n_proteins}\t{min(lengths)}\t{max(lengths)}"
                f"\t{sum(lengths) / n_entries:.1f}\t{median:.1f}\n"
            )


@app.command()
def main(
    fasta_path: Path = RAW_DATA_DIR / "piRNA_PPI_protein.fasta",
    dpam_domains_path: Path = INTERIM_DATA_DIR / "DPAM_input" / "piRNA_domains",
    iprscan_path: Path = INTERIM_DATA_DIR / "iprscan5-output.tsv",
    structure_dir: Path = INTERIM_DATA_DIR / "DPAM_input" / "piRNA",
    output_dir: Path = INTERIM_DATA_DIR / "segmented_proteins",
    max_segment_length: int = 750,
):
    config = SegmenterConfig(
        max_segment_length=max_segment_length,
        fasta_path=fasta_path,
        dpam_domains_path=dpam_domains_path,
        iprscan_path=iprscan_path,
        structure_dir=structure_dir,
        output_dir=output_dir,
    )
    run_segmentation(config)


def _explain_not_segmented(protein, config: "SegmenterConfig") -> str:
    """Return a human-readable reason why a protein was not segmented."""
    if protein.length <= config.max_segment_length:
        return (
            f"protein length ({protein.length}) <= threshold ({config.max_segment_length})"
        )
    if protein.has_structure:
        return (
            f"all {len(protein.domains)} domains merged into one element by contact/PAE metrics"
        )
    return (
        f"only {len(protein.domains)} domain(s) found after splitting; "
        "all domains grouped into a single segment"
    )


def run_segmentation(config: SegmenterConfig) -> None:
    """
    Run the protein segmentation pipeline.

    Args:
        config: Segmenter configuration.
    """
    config.output_dir.mkdir(parents=True, exist_ok=True)
    log_id = logger.add(config.output_dir / "segmenter.log", level="DEBUG")

    logger.info("Loading protein sequences...")
    proteins = parse_fasta(config.fasta_path)
    logger.info(f"Loaded {len(proteins)} proteins")

    logger.info("Loading DPAM domains...")
    dpam_domains = parse_dpam_domains(config.dpam_domains_path)
    logger.info(f"Loaded domains for {len(dpam_domains)} proteins")

    logger.info("Loading iprscan domains...")
    iprscan_domains = parse_iprscan(config.iprscan_path, config.iprscan_priority)
    logger.info(f"Loaded iprscan domains for {len(iprscan_domains)} proteins")

    # Assign domains to proteins: DPAM first, fallback to iprscan
    n_dpam, n_iprscan, n_none = 0, 0, 0
    for protein_id, protein in proteins.items():
        if protein_id in dpam_domains:
            protein.domains = dpam_domains[protein_id]
            pdb_path = config.structure_dir / f"{protein_id}.pdb"
            pae_path = config.structure_dir / f"{protein_id}.json"
            protein.has_structure = pdb_path.exists() and pae_path.exists()
            n_dpam += 1
            logger.debug(f"{protein_id}: using DPAM domains ({len(protein.domains)} domains)")
        elif protein_id in iprscan_domains:
            protein.domains = iprscan_domains[protein_id]
            protein.has_structure = False
            n_iprscan += 1
            logger.debug(f"{protein_id}: using iprscan domains ({len(protein.domains)} domains)")
        else:
            n_none += 1
            logger.debug(f"{protein_id}: no domains found")
    logger.info(f"Domain sources: {n_dpam} DPAM, {n_iprscan} iprscan, {n_none} none")

    # Fallback: if a protein has only 1 domain and length > threshold,
    # try iprscan domains which may provide more segments
    n_fallback = 0
    for protein_id, protein in proteins.items():
        if (
            len(protein.domains) == 1
            and protein.length > config.max_segment_length
            and protein_id in iprscan_domains
            and len(iprscan_domains[protein_id]) > 1
        ):
            old_source = "DPAM" if protein_id in dpam_domains else "iprscan"
            protein.domains = iprscan_domains[protein_id]
            protein.has_structure = False
            n_fallback += 1
            logger.info(
                f"{protein_id}: single-domain ({old_source}) + length {protein.length} > threshold "
                f"{config.max_segment_length}, fallback to iprscan "
                f"({len(protein.domains)} domains)"
            )
    if n_fallback:
        logger.info(f"Fallback to iprscan for {n_fallback} single-domain oversized proteins")

    # Run segmentation
    segmenter = ProteinSegmenter(config)
    all_results = []  # (protein_id, protein, segments)
    summary_rows = []

    for protein_id, protein in tqdm(proteins.items(), desc="Segmenting proteins"):
        try:
            segments = segmenter.segment_protein(protein)
            all_results.append((protein_id, protein, segments))
        except Exception as e:
            logger.error(f"{protein_id}: segmentation failed: {e}")
            all_results.append((protein_id, protein, None))

    # Prepare output directory
    config.output_dir.mkdir(parents=True, exist_ok=True)

    # Classify results and build summary
    segmented_segs = []   # qualified segments (length < threshold)
    unsegmented_segs = [] # whole proteins that don't need segmentation
    oversized_segs = []   # segments exceeding threshold after segmentation
    no_domain_segs = []   # whole proteins with no domain annotations

    for protein_id, protein, segments in all_results:
        if segments is None:
            summary_rows.append((
                protein_id, protein.length, len(protein.domains),
                0, "error", "", "none", "",
            ))
            continue

        n_seg = len(segments)
        seg_lengths = [len(s.sequence) for s in segments]
        ranges_str = ";".join(f"{s.start}-{s.end}" for s in segments)

        if n_seg > 1:
            status = "segmented"
            qualified = [s for s in segments if len(s.sequence) < config.max_segment_length]
            oversized = [s for s in segments if len(s.sequence) >= config.max_segment_length]

            if qualified:
                segmented_segs.extend(qualified)
            if oversized:
                oversized_segs.extend(oversized)
                logger.warning(
                    f"{protein_id}: {len(oversized)}/{n_seg} segments exceed threshold"
                )

            output_to = []
            if qualified:
                output_to.append("segmented")
            if oversized:
                output_to.append("oversized")
            output_str = ",".join(output_to) if output_to else "none"

            logger.info(
                f"{protein_id}: segmented into {n_seg} parts "
                f"(length={protein.length}, seg_lengths={seg_lengths})"
            )
        else:
            if not protein.domains:
                status = "no_domain"
                no_domain_segs.extend(segments)
                output_str = "no_domain"
                logger.info(
                    f"{protein_id}: not segmented → no domain annotations "
                    f"(length={protein.length}), output to no_domain.fasta"
                )
            else:
                status = "not_segmented"
                unsegmented_segs.extend(segments)
                output_str = "unsegmented"
                reason = _explain_not_segmented(protein, config)
                logger.info(
                    f"{protein_id}: not segmented → {reason} "
                    f"(length={protein.length}, n_domains={len(protein.domains)}, "
                    f"has_structure={protein.has_structure}), output to unsegmented.fasta"
                )

        summary_rows.append((
            protein_id, protein.length, len(protein.domains),
            n_seg, status, ranges_str, output_str,
            ";".join(str(l) for l in seg_lengths),
        ))

    # Write FASTA files
    write_segments_to_fasta(segmented_segs, config.output_dir / "segmented.fasta")
    write_segments_to_fasta(unsegmented_segs, config.output_dir / "unsegmented.fasta")
    write_segments_to_fasta(oversized_segs, config.output_dir / "oversized_segments.fasta")
    write_segments_to_fasta(no_domain_segs, config.output_dir / "no_domain.fasta")
    logger.info(
        f"Wrote {len(segmented_segs)} segmented, {len(unsegmented_segs)} unsegmented, "
        f"{len(oversized_segs)} oversized, {len(no_domain_segs)} no_domain segments"
    )

    # Write segment detail TSV
    detail_path = config.output_dir / "segment_detail.tsv"
    with open(detail_path, "w") as f:
        f.write("protein_id\tlength\tn_domains\tn_segments\tstatus"
                "\tsegment_ranges\toutput_to\tsegment_lengths\n")
        for row in summary_rows:
            f.write("\t".join(str(v) for v in row) + "\n")
    logger.info(f"Wrote segment detail to {detail_path}")

    # Write summary TSV (statistics for each fasta file)
    summary_path = config.output_dir / "summary.tsv"
    _write_summary(summary_path, {
        "segmented.fasta": segmented_segs,
        "unsegmented.fasta": unsegmented_segs,
        "oversized_segments.fasta": oversized_segs,
        "no_domain.fasta": no_domain_segs,
    })
    logger.info(f"Wrote summary to {summary_path}")

    n_error = sum(1 for r in summary_rows if r[4] == "error")
    logger.success(
        f"Done: {len(segmented_segs)} segmented, {len(unsegmented_segs)} unsegmented, "
        f"{len(oversized_segs)} oversized, {len(no_domain_segs)} no_domain, {n_error} errors"
    )

    logger.remove(log_id)


if __name__ == "__main__":
    app()
