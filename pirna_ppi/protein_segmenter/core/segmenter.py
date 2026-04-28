"""Core protein segmentation algorithm.

Algorithm:
    1. If protein length <= threshold, no segmentation needed.
    2. If structure files exist, use structure-based segmentation:
       - Compute pairwise domain metrics (Ncontact, lmin, avg_pae).
       - Merge condition: Ncontact/lmin + Ncontact/100 > avg_pae/8.
       - Domains meeting merge condition are merged into one element.
    3. Without structure, each domain becomes an independent element.
"""

import logging

import numpy as np

from ..config import SegmenterConfig
from ..models import Domain, ElementPairMetrics, Protein, Segment
from ..parsers.pae_parser import calculate_average_pae_between_elements, parse_pae
from ..parsers.pdb_parser import calculate_distance_matrix, parse_pdb_coordinates
from .contact_calculator import calculate_inter_element_contacts

logger = logging.getLogger(__name__)


class ProteinSegmenter:
    """Segments proteins based on domain structure and inter-domain interactions."""

    def __init__(self, config: SegmenterConfig):
        self.config = config

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def segment_protein(self, protein: Protein) -> list[Segment]:
        """Segment a protein into elements based on domain structure."""
        # No domains → whole protein as single segment
        if not protein.domains:
            logger.debug(
                "%s: no domains, keeping whole protein (length=%d)",
                protein.protein_id, protein.length,
            )
            return [self._create_whole_protein_segment(protein)]

        # Protein length <= threshold → no segmentation, trim to domain range
        if protein.length <= self.config.max_segment_length:
            logger.debug(
                "%s: length %d <= threshold %d, no segmentation needed",
                protein.protein_id, protein.length, self.config.max_segment_length,
            )
            return [self._create_segment_from_domains(protein, protein.domains, 1)]

        # Segment based on structure availability
        logger.debug(
            "%s: length %d > threshold %d, attempting segmentation "
            "(n_domains=%d, has_structure=%s)",
            protein.protein_id, protein.length, self.config.max_segment_length,
            len(protein.domains), protein.has_structure,
        )
        if protein.has_structure:
            return self._segment_with_structure(protein)
        return self._segment_without_structure(protein)

    # ------------------------------------------------------------------
    # Structure-based segmentation
    # ------------------------------------------------------------------

    def _segment_with_structure(self, protein: Protein) -> list[Segment]:
        """Merge domains by contact/PAE metrics; each merged group = one element."""
        pdb_path = self.config.structure_dir / f"{protein.protein_id}.pdb"
        pae_path = self.config.structure_dir / f"{protein.protein_id}.json"

        if not pdb_path.exists() or not pae_path.exists():
            logger.warning(f"{protein.protein_id}: structure files missing, fallback to domain-based")
            return self._segment_without_structure(protein)

        try:
            coordinates = parse_pdb_coordinates(pdb_path)
            distance_matrix, residue_numbers = calculate_distance_matrix(coordinates)
            pae_matrix = parse_pae(pae_path)
        except Exception as e:
            logger.warning(f"{protein.protein_id}: error loading structure: {e}, fallback to domain-based")
            return self._segment_without_structure(protein)

        domains = protein.domains
        n = len(domains)

        # Build merge adjacency matrix
        should_merge = np.zeros((n, n), dtype=bool)
        for i in range(n):
            for j in range(i + 1, n):
                metrics = self._calculate_pair_metrics(
                    domains[i], domains[j], distance_matrix, residue_numbers, pae_matrix,
                )
                decision = "MERGE" if metrics.should_merge else "KEEP"
                logger.debug(
                    "%s: domain pair (%s, %s) → %s "
                    "(n_contact=%d, l_min=%d, avg_pae=%.2f, "
                    "score=%.4f vs threshold=%.4f)",
                    protein.protein_id, domains[i].name, domains[j].name, decision,
                    metrics.n_contact, metrics.l_min, metrics.avg_pae,
                    (metrics.n_contact / metrics.l_min + metrics.n_contact / 100)
                    if metrics.l_min > 0 else 0,
                    metrics.avg_pae / 8,
                )
                if metrics.should_merge:
                    should_merge[i, j] = should_merge[j, i] = True

        # Group merged domains → each group becomes one element/segment
        groups = self._union_find_groups(domains, should_merge)
        n_groups = len(groups)
        logger.debug(
            "%s: structure-based grouping produced %d group(s) from %d domain(s)",
            protein.protein_id, n_groups, n,
        )
        return self._groups_to_segments(protein, groups)

    # ------------------------------------------------------------------
    # Domain-based segmentation (no structure)
    # ------------------------------------------------------------------

    def _segment_without_structure(self, protein: Protein) -> list[Segment]:
        """Each domain becomes an independent element/segment."""
        domains = sorted(protein.domains, key=lambda d: d.start)
        if not domains:
            return [self._create_whole_protein_segment(protein)]

        logger.debug(
            "%s: no structure, each of %d domain(s) becomes an independent segment",
            protein.protein_id, len(domains),
        )
        return self._groups_to_segments(protein, [[d] for d in domains])

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _calculate_pair_metrics(
        self,
        domain1: Domain,
        domain2: Domain,
        distance_matrix: np.ndarray,
        residue_numbers: list[int],
        pae_matrix: np.ndarray,
    ) -> ElementPairMetrics:
        """Calculate Ncontact, lmin, avg_pae for a domain pair."""
        indices1 = domain1.get_residue_indices()
        indices2 = domain2.get_residue_indices()

        n_contact = calculate_inter_element_contacts(
            distance_matrix, residue_numbers, indices1, indices2,
            self.config.min_sequence_separation, self.config.max_contact_distance,
        )
        l_min = min(domain1.total_length, domain2.total_length)
        avg_pae = calculate_average_pae_between_elements(pae_matrix, indices1, indices2)

        return ElementPairMetrics(
            element1=domain1, element2=domain2,
            n_contact=n_contact, l_min=l_min, avg_pae=avg_pae,
        )

    @staticmethod
    def _union_find_groups(domains: list[Domain], should_merge: np.ndarray) -> list[list[Domain]]:
        """Group domains using union-find based on merge adjacency matrix."""
        n = len(domains)
        parent = list(range(n))

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for i in range(n):
            for j in range(i + 1, n):
                if should_merge[i, j]:
                    pi, pj = find(i), find(j)
                    if pi != pj:
                        parent[pi] = pj

        groups: dict[int, list[int]] = {}
        for i in range(n):
            groups.setdefault(find(i), []).append(i)

        result = [sorted([domains[i] for i in idx], key=lambda d: d.start) for idx in groups.values()]
        result.sort(key=lambda g: g[0].start)
        return result

    def _groups_to_segments(self, protein: Protein, groups: list[list[Domain]]) -> list[Segment]:
        """Convert domain groups to Segment objects, one per group."""
        return [
            self._create_segment_from_domains(protein, group, seg_id)
            for seg_id, group in enumerate(groups, 1)
        ]

    def _create_segment_from_domains(
        self, protein: Protein, domains: list[Domain], segment_id: int,
    ) -> Segment:
        """Create a Segment spanning from min domain start to max domain end."""
        if not domains:
            return self._create_whole_protein_segment(protein)

        start = max(1, min(d.start for d in domains))
        end = min(protein.length, max(d.end for d in domains))
        sequence = protein.sequence[start - 1 : end]

        return Segment(
            protein_id=protein.protein_id,
            segment_id=segment_id,
            domains=sorted(domains, key=lambda d: d.start),
            start=start, end=end, sequence=sequence,
        )

    def _create_whole_protein_segment(self, protein: Protein) -> Segment:
        """Create a segment containing the whole protein."""
        return Segment(
            protein_id=protein.protein_id,
            segment_id=1,
            domains=protein.domains,
            start=1, end=protein.length, sequence=protein.sequence,
        )
