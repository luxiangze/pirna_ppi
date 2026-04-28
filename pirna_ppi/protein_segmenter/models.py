"""Data models for protein segmenter."""

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class Domain:
    """Represents a protein domain with its range."""

    name: str
    ranges: list[tuple[int, int]]  # List of (start, end) tuples, 1-indexed

    @property
    def total_length(self) -> int:
        """Calculate total length of domain (sum of all ranges)."""
        return sum(end - start + 1 for start, end in self.ranges)

    @property
    def start(self) -> int:
        """Get the minimum start position."""
        return min(s for s, _ in self.ranges)

    @property
    def end(self) -> int:
        """Get the maximum end position."""
        return max(e for _, e in self.ranges)

    def get_residue_indices(self) -> set[int]:
        """Get all residue indices covered by this domain."""
        indices = set()
        for start, end in self.ranges:
            indices.update(range(start, end + 1))
        return indices


@dataclass
class Protein:
    """Represents a protein with its sequence and domains."""

    protein_id: str
    sequence: str
    domains: list[Domain] = field(default_factory=list)
    has_structure: bool = False

    @property
    def length(self) -> int:
        """Get protein sequence length."""
        return len(self.sequence)

    @property
    def cumulative_domain_length(self) -> int:
        """Calculate cumulative length of all domains."""
        return sum(d.total_length for d in self.domains)


@dataclass
class Segment:
    """Represents a protein segment containing one or more domains."""

    protein_id: str
    segment_id: int
    domains: list[Domain]
    start: int  # 1-indexed
    end: int  # 1-indexed
    sequence: str

    @property
    def length(self) -> int:
        """Get segment length."""
        return self.end - self.start + 1

    def to_fasta_header(self) -> str:
        """Generate FASTA header for this segment."""
        domain_names = ",".join(d.name for d in self.domains)
        return f">{self.protein_id}-seg{self.segment_id} range={self.start}-{self.end} domains={domain_names}"

    def to_fasta(self) -> str:
        """Generate FASTA format string."""
        return f"{self.to_fasta_header()}\n{self.sequence}"


@dataclass
class ElementPairMetrics:
    """Metrics for a pair of elements (domains)."""

    element1: Domain
    element2: Domain
    n_contact: int = 0  # Number of inter-element contacts
    l_min: int = 0  # Minimum length of the two elements
    avg_pae: float = 0.0  # Average PAE between elements

    @property
    def should_merge(self) -> bool:
        """Determine if elements should be merged based on the formula."""
        if self.l_min == 0:
            return False
        # N_contact / L_min + N_contact / 100 > Avg_PAE / 8
        left_side = self.n_contact / self.l_min + self.n_contact / 100
        right_side = self.avg_pae / 8
        return left_side > right_side
