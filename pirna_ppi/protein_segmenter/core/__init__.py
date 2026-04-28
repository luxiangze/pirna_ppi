"""Core algorithms for protein segmentation."""

from .contact_calculator import calculate_inter_element_contacts
from .segmenter import ProteinSegmenter

__all__ = ["calculate_inter_element_contacts", "ProteinSegmenter"]
