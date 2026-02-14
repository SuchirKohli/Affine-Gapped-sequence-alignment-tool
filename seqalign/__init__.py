"""
SeqAlign - Nucleotide Sequence Alignment Tool

A Python package for pairwise nucleotide sequence alignment
using the Smith-Waterman algorithm with affine gap penalties.
"""

from .fasta_parser import parse_fasta
from .alignment import smith_waterman, AlignmentResult
from .display import format_alignment

__version__ = "1.0.0"
__all__ = ["parse_fasta", "smith_waterman", "AlignmentResult", "format_alignment"]
