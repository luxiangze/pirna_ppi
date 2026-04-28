"""Parser for AlphaFold2 PAE (Predicted Aligned Error) JSON files."""

import json
from pathlib import Path

import numpy as np


def parse_pae(pae_path: Path) -> np.ndarray:
    """
    Parse AlphaFold2 PAE JSON file and return the PAE matrix.

    Args:
        pae_path: Path to the PAE JSON file.

    Returns:
        2D numpy array of PAE values. pae[i][j] is the expected position error
        of residue j when aligned on residue i.
    """
    with open(pae_path) as f:
        data = json.load(f)

    # AlphaFold2 format: list with single dict containing 'predicted_aligned_error'
    if isinstance(data, list) and len(data) > 0:
        pae_data = data[0].get("predicted_aligned_error")
    elif isinstance(data, dict):
        pae_data = data.get("predicted_aligned_error")
    else:
        raise ValueError(f"Unexpected PAE JSON format in {pae_path}")

    if pae_data is None:
        raise ValueError(f"No 'predicted_aligned_error' found in {pae_path}")

    return np.array(pae_data)


def calculate_average_pae_between_elements(
    pae_matrix: np.ndarray,
    element1_indices: set[int],
    element2_indices: set[int],
) -> float:
    """
    Calculate average PAE between two elements (domains).

    Args:
        pae_matrix: PAE matrix (N x N), 0-indexed.
        element1_indices: Set of residue indices for element 1 (1-indexed).
        element2_indices: Set of residue indices for element 2 (1-indexed).

    Returns:
        Average PAE value between the two elements.
    """
    # Convert to 0-indexed
    idx1 = [i - 1 for i in element1_indices]
    idx2 = [i - 1 for i in element2_indices]

    # Get PAE values between elements (both directions)
    pae_values = []
    for i in idx1:
        for j in idx2:
            if 0 <= i < pae_matrix.shape[0] and 0 <= j < pae_matrix.shape[1]:
                pae_values.append(pae_matrix[i, j])
                pae_values.append(pae_matrix[j, i])

    if not pae_values:
        return float("inf")

    return np.mean(pae_values)
