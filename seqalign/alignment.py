"""
Local alignment with affine gap penalties.

Implements the three-matrix approach:
  H[i][j] = main scoring matrix
  E[i][j] = best score ending with a gap in seq1 (insertion)
  F[i][j] = best score ending with a gap in seq2 (deletion)
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class AlignmentResult:
    """Container for alignment results."""

    aligned_seq1: str
    aligned_seq2: str
    score: float
    start_seq1: int
    end_seq1: int
    start_seq2: int
    end_seq2: int
    identity: float
    gaps: int
    mismatches: int
    alignment_length: int


def _score_func(a: str, b: str, match: float, mismatch: float) -> float:
    """Returns the match/mismatch score for two nucleotides."""
    return match if a == b else mismatch


def smith_waterman(
    seq1: str,
    seq2: str,
    match: float = 5.0,
    mismatch: float = -4.0,
    gap_open: float = -12.0,
    gap_ext: float = -2.0,
) -> AlignmentResult:
    """
    Perform local pairwise alignment using Smith-Waterman with affine gaps.

    Parameters
    ----------
    seq1 : str
        First nucleotide sequence.
    seq2 : str
        Second nucleotide sequence.
    match : float
        Score for a matching pair (should be > 0). Default: 5.
    mismatch : float
        Score for a mismatching pair (should be ≤ 0). Default: -4.
    gap_open : float
        Penalty for opening a new gap (should be ≤ 0). Default: -12.
    gap_ext : float
        Penalty for extending an existing gap (should be ≤ 0). Default: -2.

    Returns
    -------
    AlignmentResult
        The best local alignment found.
    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    n = len(seq1)
    m = len(seq2)

    if n == 0 or m == 0:
        raise ValueError("Both sequences must be non-empty.")

    NEG_INF = float("-inf")

    # Initialise matrices 
    H = np.zeros((n + 1, m + 1), dtype=np.float64)   # main scoring matrix
    E = np.full((n + 1, m + 1), NEG_INF, dtype=np.float64)  # gap in seq1
    F = np.full((n + 1, m + 1), NEG_INF, dtype=np.float64)  # gap in seq2

    # Traceback constants
    STOP = 0
    DIAG = 1  # match / mismatch
    UP   = 2  # gap in seq2 (consume seq1[i])
    LEFT = 3  # gap in seq1 (consume seq2[j])

    traceback = np.zeros((n + 1, m + 1), dtype=np.int8)

    # Fill matrices 
    max_score = 0.0
    max_i, max_j = 0, 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Score for aligning seq1[i-1] with seq2[j-1]
            diag_score = H[i - 1, j - 1] + _score_func(
                seq1[i - 1], seq2[j - 1], match, mismatch
            )

            # E: best score ending with gap in seq1 (insertion in seq2)
            E[i, j] = max(
                H[i, j - 1] + gap_open + gap_ext,
                E[i, j - 1] + gap_ext,
            )

            # F: best score ending with gap in seq2 (insertion in seq1)
            F[i, j] = max(
                H[i - 1, j] + gap_open + gap_ext,
                F[i - 1, j] + gap_ext,
            )

            # H: main matrix (local alignment, so minimum is 0)
            H[i, j] = max(0, diag_score, E[i, j], F[i, j])

            # Traceback pointer
            if H[i, j] == 0:
                traceback[i, j] = STOP
            elif H[i, j] == diag_score:
                traceback[i, j] = DIAG
            elif H[i, j] == F[i, j]:
                traceback[i, j] = UP
            else:
                traceback[i, j] = LEFT

            # Track global maximum
            if H[i, j] > max_score:
                max_score = H[i, j]
                max_i, max_j = i, j

    # Traceback 
    aligned1 = []
    aligned2 = []
    i, j = max_i, max_j

    while i > 0 and j > 0 and H[i, j] > 0:
        if traceback[i, j] == DIAG:
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif traceback[i, j] == UP:
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
        elif traceback[i, j] == LEFT:
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1
        else:  
            break

    aligned1.reverse()
    aligned2.reverse()

    aligned_seq1 = "".join(aligned1)
    aligned_seq2 = "".join(aligned2)

    # Computing statistics
    alignment_length = len(aligned_seq1)
    matches = sum(
        1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != "-"
    )
    mismatches = sum(
        1
        for a, b in zip(aligned_seq1, aligned_seq2)
        if a != b and a != "-" and b != "-"
    )
    gaps = aligned_seq1.count("-") + aligned_seq2.count("-")
    identity = (matches / alignment_length * 100) if alignment_length > 0 else 0.0

    # Start positions are 1-indexed (biological convention)
    start_seq1 = i + 1
    start_seq2 = j + 1
    end_seq1 = max_i
    end_seq2 = max_j

    return AlignmentResult(
        aligned_seq1=aligned_seq1,
        aligned_seq2=aligned_seq2,
        score=max_score,
        start_seq1=start_seq1,
        end_seq1=end_seq1,
        start_seq2=start_seq2,
        end_seq2=end_seq2,
        identity=identity,
        gaps=gaps,
        mismatches=mismatches,
        alignment_length=alignment_length,
    )
