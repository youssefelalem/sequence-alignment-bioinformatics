# 🧬 Sequence Alignment - Bioinformatics

Implementation of classic sequence alignment algorithms in Python for DNA sequence comparison.

## 📋 Description

This project implements the two fundamental algorithms for biological sequence alignment:

- **Needleman-Wunsch** (Global Alignment) - Aligns entire sequences end-to-end
- **Smith-Waterman** (Local Alignment) - Finds the best matching subsequences

## 🚀 Features

- ✅ Read sequences from FASTA files
- ✅ Global alignment with Needleman-Wunsch algorithm
- ✅ Local alignment with Smith-Waterman algorithm  
- ✅ Backtracking to extract optimal alignments
- ✅ Score matrix visualization
- ✅ Validated against EMBOSS Needle & Water tools

## 📁 Project Structure

```
sequence-alignment-bioinformatics/
├── alignment.py        # Main Python script with all algorithms
├── sequences.fasta     # Sample DNA sequences in FASTA format
└── README.md           # This file
```

## 🔧 Requirements

- Python 3.x
- NumPy

```bash
pip install numpy
```

## 💻 Usage

```bash
python alignment.py
```

### Example Output

```
Séquence 1: ATGCGTACGTTAGC
Séquence 2: ATGCCGTCGTTAGG

============================================================
         ALIGNEMENT GLOBAL (Needleman-Wunsch)
============================================================

Score final: 10

=== Alignement Global ===
Seq1: ATGCGTACGTTAGC
Seq2: ATGCCGTCGTTAGG

============================================================
          ALIGNEMENT LOCAL (Smith-Waterman)
============================================================

Score maximal: 10

=== Alignement Local ===
Seq1: ATGCGTACGTTAG
Seq2: ATGCCGTCGTTAG
```

## ⚙️ Parameters

Default scoring scheme:
| Parameter | Value |
|-----------|-------|
| Match     | +1    |
| Mismatch  | 0     |
| Gap       | -1    |

You can modify these parameters in the function calls:
```python
matrix = score_matrix_alignement_global(seq1, seq2, match=2, mismatch=-1, gap=-2)
```

## 📊 Algorithms

### Needleman-Wunsch (Global)
- Initializes borders with gap penalties
- Fills matrix using dynamic programming
- Backtracking from bottom-right corner

### Smith-Waterman (Local)
- Initializes borders with zeros
- Minimum score is 0 (never goes negative)
- Backtracking from maximum score position until reaching 0

## 🔗 Validation

Results validated against professional tools:
- [EMBOSS Needle](https://www.ebi.ac.uk/Tools/psa/emboss_needle/) (Global)
- [EMBOSS Water](https://www.ebi.ac.uk/Tools/psa/emboss_water/) (Local)

## 👤 Author

**EL ALEM YOUSSEF**

## 📄 License

This project is for educational purposes - Bioinformatics TP3.
