"""
Command-line interface for the SeqAlign tool.

Takes user input for:
  1. Paths to two FASTA files
  2. Scoring parameters: match, mismatch, gap opening, gap extension
"""

import sys
import os

# Allow running from the project root
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from seqalign import parse_fasta, smith_waterman, format_alignment


def get_float(prompt: str, default: float) -> float:
    """Prompt the user for a float value with a default."""
    raw = input(f"{prompt} [{default}]: ").strip()
    if not raw:
        return default
    try:
        return float(raw)
    except ValueError:
        print(f"  Invalid number '{raw}', using default {default}")
        return default


def main():
    print("=" * 60)
    print("  SeqAlign — Nucleotide Sequence Alignment Tool")
    print("  Smith-Waterman | Affine Gap Penalties")
    print("=" * 60)
    print()

    # Input 1: FASTA files
    file1 = input("Path to FASTA file 1: ").strip()
    if not os.path.isfile(file1):
        print(f"Error: File '{file1}' not found.")
        sys.exit(1)

    file2 = input("Path to FASTA file 2: ").strip()
    if not os.path.isfile(file2):
        print(f"Error: File '{file2}' not found.")
        sys.exit(1)

    # Input 2: Scoring parameters
    print("\n--- Scoring Parameters ---")
    print("  (Press Enter to use the default value)\n")

    match    = get_float("  Match score",          default=5.0)
    mismatch = get_float("  Mismatch penalty",     default=-4.0)
    gap_open = get_float("  Gap opening penalty",  default=-12.0)
    gap_ext  = get_float("  Gap extension penalty", default=-2.0)

    # Parse sequences
    print("\nParsing FASTA files...")
    try:
        header1, seq1 = parse_fasta(file1)
        header2, seq2 = parse_fasta(file2)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")
        sys.exit(1)

    print(f"  Seq1: {header1} ({len(seq1)} bp)")
    print(f"  Seq2: {header2} ({len(seq2)} bp)")

    # Align
    print("\nAligning sequences...")
    print(f"  Parameters: match={match}, mismatch={mismatch}, "
          f"gap_open={gap_open}, gap_ext={gap_ext}")
    print()

    result = smith_waterman(
        seq1, seq2,
        match=match,
        mismatch=mismatch,
        gap_open=gap_open,
        gap_ext=gap_ext,
    )

    # Display
    print(format_alignment(result))


if __name__ == "__main__":
    main()
