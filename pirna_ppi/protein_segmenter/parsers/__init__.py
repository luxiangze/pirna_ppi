"""Parsers for various input file formats."""

from .domain_parser import parse_dpam_domains
from .fasta_parser import parse_fasta
from .iprscan_parser import parse_iprscan
from .pae_parser import parse_pae
from .pdb_parser import parse_pdb_coordinates

__all__ = [
    "parse_dpam_domains",
    "parse_fasta",
    "parse_iprscan",
    "parse_pae",
    "parse_pdb_coordinates",
]
