"""
FASTA file parser for nucleotide sequences.
"""


def parse_fasta(filepath: str) -> tuple[str, str]:
    """
    Parse a FASTA file and return the header and sequence.

    Parameters
    ----------
    filepath : str
        Path to the FASTA file.

    Returns
    -------
    tuple[str, str]
        A tuple of (header, sequence) where header is the description line
        (without '>') and sequence is the concatenated nucleotide string
        in uppercase.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file is not in valid FASTA format or contains
        invalid nucleotide characters.
    """
    header = None
    sequence_lines = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    # Only read the first sequence in the file
                    break
                header = line[1:].strip()
            else:
                sequence_lines.append(line.upper())

    if header is None:
        raise ValueError(
            f"Invalid FASTA file '{filepath}': no header line starting with '>' found."
        )

    sequence = "".join(sequence_lines)

    if not sequence:
        raise ValueError(
            f"Invalid FASTA file '{filepath}': no sequence data found."
        )

    # Validate nucleotide characters
    valid_chars = set("ACGTURYKMSWBDHVN")  # IUPAC nucleotide codes
    invalid = set(sequence) - valid_chars
    if invalid:
        raise ValueError(
            f"Invalid nucleotide characters found in '{filepath}': {invalid}"
        )

    return header, sequence
