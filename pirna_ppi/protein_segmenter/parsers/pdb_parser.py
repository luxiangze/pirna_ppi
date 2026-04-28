"""Parser for PDB structure files."""

from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser


def parse_pdb_coordinates(pdb_path: Path) -> dict[int, np.ndarray]:
    """
    Parse PDB file and extract CA atom coordinates for each residue.

    Args:
        pdb_path: Path to the PDB file.

    Returns:
        Dictionary mapping residue number (1-indexed) to CA coordinates (x, y, z).
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    coordinates = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                # Skip hetero atoms
                if residue.id[0] != " ":
                    continue
                res_num = residue.id[1]
                if "CA" in residue:
                    ca_atom = residue["CA"]
                    coordinates[res_num] = np.array(ca_atom.get_coord())

    return coordinates


def calculate_distance_matrix(coordinates: dict[int, np.ndarray]) -> tuple[np.ndarray, list[int]]:
    """
    Calculate pairwise distance matrix from CA coordinates.

    Args:
        coordinates: Dictionary mapping residue number to CA coordinates.

    Returns:
        Tuple of (distance_matrix, residue_numbers).
        distance_matrix[i][j] is the distance between residue_numbers[i] and residue_numbers[j].
    """
    residue_numbers = sorted(coordinates.keys())
    n = len(residue_numbers)

    # Build coordinate matrix
    coord_matrix = np.array([coordinates[res_num] for res_num in residue_numbers])

    # Calculate pairwise distances
    diff = coord_matrix[:, np.newaxis, :] - coord_matrix[np.newaxis, :, :]
    distance_matrix = np.sqrt(np.sum(diff**2, axis=2))

    return distance_matrix, residue_numbers
