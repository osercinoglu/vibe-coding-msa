# Multiple Sequence Alignment (MSA) and Analysis Script

A self-contained Python script for performing progressive multiple sequence alignment (MSA) on protein sequences using only standard scientific libraries (numpy).

## Features

- **Progressive MSA**: Uses UPGMA-based guide tree construction with pairwise Needleman-Wunsch alignment
- **Pairwise Identity Matrix**: Computes sequence identity for all sequence pairs (0–1 scale)
- **Average Pairwise Identity**: Reports mean identity across all pairwise comparisons
- **Conservation Analysis**: Identifies fully conserved alignment columns (ignoring gaps)
- **BLOSUM62 Scoring**: Uses simplified BLOSUM62 matrix for biologically meaningful alignment
- **Output**: Writes aligned sequences to FASTA format for downstream analysis

## Requirements

- Python 3.11+
- `numpy` (standard scientific library)

No internet access or external tools (e.g., Clustal, MAFFT) required.

## Installation

No installation needed beyond Python 3.11 and numpy. If numpy is not available:

```bash
pip install numpy
```

## Usage

```bash
python msa_task.py --input input.fasta --output aligned.fasta
```

Or with short options:

```bash
python msa_task.py -i input.fasta -o aligned.fasta
```

### Example

```bash
python msa_task.py --input sequences.fasta --output aligned_output.fasta
```

## Output

The script prints analysis results to stderr and produces two outputs:

### Console Output (stderr)

```
Pairwise Sequence Identity Matrix:
                 seq1      seq2      seq3
seq1           1.0000    0.5234    0.4821
seq2           0.5234    1.0000    0.6123
seq3           0.4821    0.6123    1.0000

Average Pairwise Identity: 0.5526

Conservation Analysis:
  Fully conserved columns: 45 / 287
  Fraction conserved: 0.1568
```

### Output FASTA File (`aligned.fasta`)

Contains all sequences aligned with gaps (`-`) inserted for alignment:

```fasta
>seq1
MKGKDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWK-TMSAKEKGKFED
>seq2
MQDRVKRP-----MNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEA
>seq3
M--HIKKP------LNAFMLYMKEMRANVVAESTLKESAAINQILGRRWHALSREEQAKYYEL
```

## Algorithm Details

### Progressive Alignment (UPGMA-like)

1. **Distance Calculation**: Compute pairwise distances as `1 - identity` between all sequences
2. **Hierarchical Clustering**: Iteratively merge the closest sequence pairs using average linkage
3. **Guide Tree**: Build a binary tree representing the merge order
4. **Progressive Alignment**: Align sequences along the guide tree using the Needleman-Wunsch algorithm
5. **Refinement**: Optional iterative refinement to improve alignment quality

### Needleman-Wunsch Algorithm

- **Global alignment** of two sequences (or profiles)
- **Gap penalties**: 
  - Gap open: -10
  - Gap extension: -1
- **Scoring**: BLOSUM62-based matrix for amino acid substitutions
- **Traceback**: Dynamic programming matrix traceback for optimal alignment

### Pairwise Identity Calculation

```
Identity = (Number of matching positions) / (Total valid positions excluding gaps)
```

Gaps are excluded from the calculation to focus on aligned residues.

### Conservation Scoring

Columns are considered "fully conserved" if:
- All positions contain the same amino acid
- No gaps or ambiguous residues (`-` or `.`) are present

## Input Format

Standard FASTA format with sequences on one or multiple lines:

```fasta
>sequence_header_1
MKVLIVFASFAVSSAYDIENLSAETKK
>sequence_header_2
MAVLIVFTSFSVSNAVDVENLSAETSK
```

## Interpretation

### Pairwise Identity Matrix

- **Diagonal (1.0)**: Each sequence is identical to itself
- **Off-diagonal values**: Similarity between sequence pairs
  - Close to 1.0: High similarity (similar sequences)
  - Close to 0.0: Low similarity (divergent sequences)

### Average Pairwise Identity

- **High values (>0.7)**: Sequences are well-conserved as a family
- **Low values (<0.3)**: Sequences are quite divergent
- Used to assess overall sequence family conservation

### Fraction Conserved

- **High values (>0.5)**: Well-conserved alignment (possibly orthologous sequences)
- **Low values (<0.1)**: Highly divergent sequences with few conserved positions
- Useful for identifying functionally important residues

## Performance Characteristics

- **Time Complexity**: O(N² × L²) where N = number of sequences, L = sequence length
- **Space Complexity**: O(N × L) for storing alignments + O(L²) for DP matrix
- **Typical Performance**: 
  - 4 sequences × 85 residues: < 1 second
  - 10 sequences × 200 residues: < 5 seconds
  - 50 sequences × 1000 residues: ~1-2 minutes

## Limitations and Considerations

1. **No external tools**: All alignment done in pure Python (may be slower than optimized tools like Clustal/MAFFT)
2. **Progressive alignment**: Greedy strategy; does not guarantee globally optimal MSA (heuristic)
3. **UPGMA tree**: Assumes constant evolutionary rates (may not be accurate for distantly related sequences)
4. **BLOSUM62 approximation**: Simplified version; full matrix available in BioPython
5. **No profiles**: Treats all sequences equally in progressive alignment (not profile-to-sequence)

For production bioinformatics pipelines with many sequences or sequences with large gaps, consider:
- **Clustal Omega**: Fast, accurate, handles many sequences
- **MAFFT**: Very fast, good for large datasets
- **T-Coffee**: Accurate, good for structurally complex sequences
- **Muscle**: Good balance of speed and accuracy

## Example Workflow

```bash
# Align sequences
python msa_task.py --input my_proteins.fasta --output my_proteins_aligned.fasta

# Output analysis is printed to console
# - Check average pairwise identity to assess sequence conservation
# - Use conservation score to identify potential active site residues
# - Use aligned.fasta for phylogenetic analysis, domain identification, etc.
```

## Troubleshooting

**Issue**: Script terminates with "No sequences found"
- **Solution**: Ensure FASTA file is properly formatted with `>` headers

**Issue**: Very low pairwise identities (all near 0)
- **Cause**: Sequences may be from different families or too divergent
- **Solution**: Verify that sequences are truly homologous; re-check input

**Issue**: Memory error with many large sequences
- **Cause**: O(N² × L²) memory usage becomes prohibitive
- **Solution**: Align in smaller batches or use external MSA tools

## Implementation Notes

- **FASTA Parser**: Handles multi-line sequences and various header formats
- **Numpy for Speed**: Uses numpy arrays for dynamic programming tables
- **Gap Handling**: Gaps (`-`) and dots (`.`) are treated equivalently as missing data
- **Error Handling**: Validates input files and reports clear error messages

## Citation

If you use this script in published research, please cite:

- Needleman, S.B. and Wunsch, C.D. (1970). J. Mol. Biol., 48(3):443-453
- Henikoff, S. and Henikoff, J.G. (1992). PNAS, 89(22):10915-10919

## License

This script is provided as-is for educational and research purposes.
