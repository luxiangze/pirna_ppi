"""Configuration for protein segmenter."""

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class SegmenterConfig:
    """Configuration parameters for protein segmentation."""

    # Protein length threshold: segment only if protein length > this value
    max_segment_length: int = 750

    # Contact calculation parameters
    min_sequence_separation: int = 10  # Minimum sequence distance for contacts
    max_contact_distance: float = 8.0  # Maximum 3D distance (Angstrom) for contacts

    # Input paths
    fasta_path: Path = field(default_factory=lambda: Path("data/raw/piRNA_PPI_protein.fasta"))
    dpam_domains_path: Path = field(default_factory=lambda: Path("data/interim/DPAM_input/piRNA_domains"))
    iprscan_path: Path = field(default_factory=lambda: Path("data/interim/iprscan5-output.tsv"))
    structure_dir: Path = field(default_factory=lambda: Path("data/interim/DPAM_input/piRNA"))

    # Output directory
    output_dir: Path = field(default_factory=lambda: Path("data/interim/segmented_proteins"))

    # iprscan analysis priority
    iprscan_priority: list = field(default_factory=lambda: ["Pfam", "Gene3D"])
