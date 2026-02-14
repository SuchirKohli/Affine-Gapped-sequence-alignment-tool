# Affine-Gapped Sequence Alignment Tool

A modular Python tool for **pairwise nucleotide sequence alignment** using the **Smith-Waterman algorithm** with **affine gap penalties**.

## Algorithm

The tool implements **local alignment** (finds the best matching sub-region between two sequences) using a three-matrix dynamic programming approach:

| Matrix | Purpose |
|--------|---------|
| **H** (main) | Overall best score at each position; minimum is 0 (local alignment) |
| **E** | Best score ending with a gap in Sequence 1 (insertion) |
| **F** | Best score ending with a gap in Sequence 2 (deletion) |

### Recurrence Relations

```
E(i,j) = max( H(i, j-1) + gap_open + gap_ext,  E(i, j-1) + gap_ext )
F(i,j) = max( H(i-1, j) + gap_open + gap_ext,  F(i-1, j) + gap_ext )
H(i,j) = max( 0,  H(i-1,j-1) + score(a,b),  E(i,j),  F(i,j) )
```

Where:
- `score(a,b)` = **match** if `a == b`, else **mismatch**
- The **gap cost** for a gap of length *k* = `gap_open + k × gap_ext`

The best alignment is found by tracing back from the maximum-scoring cell in **H** until a cell with score 0 is reached.

### Why Affine Gap Penalties?

A simple linear gap penalty treats every gap position equally. In biology, however, **one long insertion/deletion** is more likely than many small ones. Affine gap penalties model this by making the **first gap position expensive** (gap_open) and **each extension cheap** (gap_ext).

## Project Structure

```
├── run_alignment.py           # Interactive CLI (main entry point)
├── seqalign/                  # Core Python package
│   ├── __init__.py            # Package exports
│   ├── fasta_parser.py        # FASTA file parsing & validation
│   ├── alignment.py           # Smith-Waterman (affine gaps)
│   └── display.py             # Alignment formatting
├── examples/
│   ├── seq1.fasta             # Example DNA sequence 1
│   ├── seq2.fasta             # Example DNA sequence 2
│   └── alignment_demo.ipynb   # Jupyter notebook with demos
├── environment.yml            # Conda environment specification
├── README.md
└── .gitignore
```

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/SuchirKohli/Affine-Gapped-sequence-alignment-tool.git
cd Affine-Gapped-sequence-alignment-tool
```

### 2. Create the conda environment
```bash
conda env create -f environment.yml
conda activate seqalign
```

## Usage

### Interactive CLI (Recommended)

The tool prompts you for **all inputs** — FASTA file paths and scoring parameters:

```bash
python run_alignment.py
```

You will be asked for:
1. **Path to FASTA file 1**
2. **Path to FASTA file 2**
3. **Match score** (default: 5)
4. **Mismatch penalty** (default: -4)
5. **Gap opening penalty** (default: -12)
6. **Gap extension penalty** (default: -2)

Example session:
```
============================================================
  SeqAlign — Nucleotide Sequence Alignment Tool
  Smith-Waterman | Affine Gap Penalties
============================================================

Path to FASTA file 1: examples/seq1.fasta
Path to FASTA file 2: examples/seq2.fasta

--- Scoring Parameters ---
  (Press Enter to use the default value)

  Match score [5.0]: 5
  Mismatch penalty [-4.0]: -4
  Gap opening penalty [-12.0]: -12
  Gap extension penalty [-2.0]: 0

Parsing FASTA files...
  Seq1: Seq1_Human_BRCA1_exon (186 bp)
  Seq2: Seq2_Mouse_Brca1_exon (186 bp)

Aligning sequences...
  Parameters: match=5.0, mismatch=-4.0, gap_open=-12.0, gap_ext=0.0

(alignment output follows...)
```

### In Python

```python
from seqalign import parse_fasta, smith_waterman, format_alignment

# Parse FASTA files
header1, seq1 = parse_fasta("examples/seq1.fasta")
header2, seq2 = parse_fasta("examples/seq2.fasta")

# Perform alignment
result = smith_waterman(
    seq1, seq2,
    match=5,
    mismatch=-4,
    gap_open=-12,
    gap_ext=-2
)

# Display the result
print(format_alignment(result))
```

### In Jupyter Notebook

Open and run the demo notebook:
```bash
jupyter notebook examples/alignment_demo.ipynb
```

## Scoring Parameters

| Parameter | Description | Default | Notes |
|-----------|-------------|---------|-------|
| `match` | Score for identical nucleotides | `5` | Positive value |
| `mismatch` | Penalty for different nucleotides | `-4` | Negative (or zero) |
| `gap_open` | Penalty for starting a new gap | `-12` | Negative (or zero) |
| `gap_ext` | Penalty for extending an existing gap | `-2` | Negative (or zero) |

> **Note:** Gap penalties follow the same sign convention as [EBI LALIGN](https://www.ebi.ac.uk/jdispatcher/psa/lalign?stype=dna) — all penalties are **negative** values.

## Validation

The tool's output has been validated against the [EBI LALIGN online tool](https://www.ebi.ac.uk/jdispatcher/psa/lalign?stype=dna&gapext=0). See the demo notebook (`examples/alignment_demo.ipynb`) for comparison examples.

## Dependencies

- Python ≥ 3.11
- NumPy
- Jupyter / Notebook

All dependencies are managed via the provided `environment.yml`.
