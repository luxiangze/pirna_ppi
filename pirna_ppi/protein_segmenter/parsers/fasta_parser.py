"""Parser for FASTA files using Biopython."""

from pathlib import Path

from Bio import SeqIO

from ..models import Protein


def parse_fasta(fasta_path: Path) -> dict[str, Protein]:
    """
    Parse a FASTA file and return a dictionary of Protein objects.

    Args:
        fasta_path: Path to the FASTA file.

    Returns:
        Dictionary mapping protein IDs to Protein objects.
    """
    proteins = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        protein_id = record.id
        sequence = str(record.seq)
        proteins[protein_id] = Protein(protein_id=protein_id, sequence=sequence)
    return proteins
