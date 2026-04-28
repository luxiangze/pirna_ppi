"""Calculate inter-element contacts from structure data."""

import numpy as np

from ..parsers.pdb_parser import calculate_distance_matrix


def calculate_inter_element_contacts(
    distance_matrix: np.ndarray,
    residue_numbers: list[int],
    element1_indices: set[int],
    element2_indices: set[int],
    min_sequence_separation: int = 10,
    max_contact_distance: float = 8.0,
) -> int:
    """
    Calculate number of inter-element contacts.

    A contact is defined as:
    - Sequence separation >= min_sequence_separation
    - 3D distance < max_contact_distance

    Args:
        distance_matrix: Pairwise distance matrix.
        residue_numbers: List of residue numbers corresponding to matrix indices.
        element1_indices: Set of residue indices for element 1 (1-indexed).
        element2_indices: Set of residue indices for element 2 (1-indexed).
        min_sequence_separation: Minimum sequence distance for contacts.
        max_contact_distance: Maximum 3D distance (Angstrom) for contacts.

    Returns:
        Number of inter-element contacts.
    """
    # Create mapping from residue number to matrix index
    res_to_idx = {res: idx for idx, res in enumerate(residue_numbers)}

    n_contacts = 0
    for res1 in element1_indices:
        for res2 in element2_indices:
            # Check sequence separation
            if abs(res1 - res2) < min_sequence_separation:
                continue

            # Check if both residues are in the distance matrix
            if res1 not in res_to_idx or res2 not in res_to_idx:
                continue

            idx1 = res_to_idx[res1]
            idx2 = res_to_idx[res2]

            # Check 3D distance
            if distance_matrix[idx1, idx2] < max_contact_distance:
                n_contacts += 1

    return n_contacts
