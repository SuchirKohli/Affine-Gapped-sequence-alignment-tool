"""
Alignment display and formatting utilities.
"""

from .alignment import AlignmentResult


def format_alignment(result: AlignmentResult, line_width: int = 60) -> str:
    """
    Format an AlignmentResult into a human-readable string.

    The output shows position numbers, a midline showing matches (|), 
    mismatches (.), and gaps ( ), and summary statistics.

    Parameters
    ----------
    result : AlignmentResult
        The alignment result to format.
    line_width : int
        Number of alignment columns per line (default 60).

    Returns
    -------
    str
        Formatted alignment string.
    """
    lines = []

    # Header
    lines.append("=" * 70)
    lines.append("  LOCAL ALIGNMENT (Smith-Waterman, Affine Gap Penalties)")
    lines.append("=" * 70)
    lines.append("")
    lines.append(f"  Score: {result.score:.1f}")
    lines.append(
        f"  Identity: {result.identity:.1f}% "
        f"({_count_matches(result)}/{result.alignment_length})"
    )
    lines.append(f"  Gaps: {result.gaps}/{result.alignment_length}")
    lines.append(f"  Mismatches: {result.mismatches}/{result.alignment_length}")
    lines.append(f"  Alignment length: {result.alignment_length}")
    lines.append("")
    lines.append(
        f"  Seq1 region: {result.start_seq1} – {result.end_seq1}"
    )
    lines.append(
        f"  Seq2 region: {result.start_seq2} – {result.end_seq2}"
    )
    lines.append("")
    lines.append("-" * 70)
    lines.append("")

    # Build midline
    midline = []
    for a, b in zip(result.aligned_seq1, result.aligned_seq2):
        if a == b and a != "-":
            midline.append("|")
        elif a == "-" or b == "-":
            midline.append(" ")
        else:
            midline.append(".")
    midline_str = "".join(midline)

    # Print in blocks
    pos1 = result.start_seq1  # current position in seq1
    pos2 = result.start_seq2  # current position in seq2

    for start in range(0, result.alignment_length, line_width):
        end = min(start + line_width, result.alignment_length)
        chunk1 = result.aligned_seq1[start:end]
        chunk2 = result.aligned_seq2[start:end]
        chunk_mid = midline_str[start:end]

        # Count actual residues (non-gap) to compute end positions
        residues1 = sum(1 for c in chunk1 if c != "-")
        residues2 = sum(1 for c in chunk2 if c != "-")

        end_pos1 = pos1 + residues1 - 1
        end_pos2 = pos2 + residues2 - 1

        # Format position labels with padding
        label_width = max(len(str(end_pos1)), len(str(end_pos2)), 4)

        lines.append(
            f"  Seq1  {pos1:>{label_width}}  {chunk1}  {end_pos1}"
        )
        lines.append(
            f"        {' ' * label_width}  {chunk_mid}"
        )
        lines.append(
            f"  Seq2  {pos2:>{label_width}}  {chunk2}  {end_pos2}"
        )
        lines.append("")

        pos1 = end_pos1 + 1
        pos2 = end_pos2 + 1

    lines.append("=" * 70)
    return "\n".join(lines)


def _count_matches(result: AlignmentResult) -> int:
    """Count the number of matches in the alignment."""
    return sum(
        1
        for a, b in zip(result.aligned_seq1, result.aligned_seq2)
        if a == b and a != "-"
    )
