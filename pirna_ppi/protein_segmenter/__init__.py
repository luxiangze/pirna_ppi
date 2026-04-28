"""Protein segmenter module for splitting large proteins into segments."""

from .config import SegmenterConfig
from .models import Domain, Protein, Segment

__all__ = ["SegmenterConfig", "Domain", "Protein", "Segment"]
