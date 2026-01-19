# Quick Start Guide

## One-Minute Setup

1. **Ensure you have Python 3.11+ and numpy installed**:
   ```bash
   python --version  # Should be 3.11 or higher
   pip list | grep numpy  # Should show numpy
   ```

2. **Run the script**:
   ```bash
   python msa_task.py --input input.fasta --output aligned.fasta
   ```

3. **Check the results**:
   - Console output shows:
     * Pairwise sequence identity matrix
     * Average pairwise identity
     * Number and fraction of fully conserved columns
   - Output file `aligned.fasta` contains the aligned sequences

## What the Script Does

The script implements a **progressive multiple sequence alignment (MSA)** algorithm:

```
1. Read unaligned protein sequences from FASTA file
   ↓
2. Compute pairwise sequence distances
   ↓
3. Build a guide tree using hierarchical clustering (UPGMA)
   ↓
4. Progressively align sequences along the tree
   ↓
5. Refine the alignment (optional)
   ↓
6. Analyze conservation and identity
   ↓
7. Output aligned sequences to FASTA
```

## Understanding the Output

### Pairwise Identity Matrix

```
                 seq1      seq2      seq3      seq4
seq1           1.0000    0.2763    0.0000    0.0000
seq2           0.2763    1.0000    0.0000    0.0000
seq3           0.0000    0.0000    1.0000    0.0000
seq4           0.0000    0.0000    0.0000    1.0000
```

- **Diagonal = 1.0**: Each sequence is 100% identical to itself
- **Off-diagonal values**: Percentage of identical residues between sequence pairs
  - 0.50 = 50% identity
  - 0.27 = 27% identity
  - 0.00 = No identity (completely different)

### Average Pairwise Identity

```
Average Pairwise Identity: 0.0461
```

This is the mean identity across all pairwise comparisons. For your 4 sequences, the average similarity is 4.61%.

### Conservation Analysis

```
Fully conserved columns: 2 / 71
Fraction conserved: 0.0282
```

- **Fully conserved columns**: Alignment positions where ALL sequences have the same amino acid (ignoring gaps)
- **2 out of 71**: Only 2 positions out of 71 aligned columns are identical across all 4 sequences
- **Fraction 0.0282**: 2.82% of the alignment is fully conserved

## Using the Aligned Output

The output file `aligned.fasta` can be used for:

1. **Phylogenetic analysis**: Build evolutionary trees with FastTree, RAxML, etc.
2. **Domain identification**: Identify conserved domains (InterProScan, Pfam)
3. **Secondary structure prediction**: Predict 3D structure (AlphaFold, I-TASSER)
4. **Motif discovery**: Find conserved sequence patterns
5. **Evolution studies**: Compare divergence rates and selective pressures
6. **Quality assessment**: Identify orthologs vs. paralogs

## Example Commands

```bash
# Basic usage
python msa_task.py --input proteins.fasta --output proteins_aligned.fasta

# Short option syntax
python msa_task.py -i proteins.fasta -o proteins_aligned.fasta

# See help
python msa_task.py --help

# Check Python version
python --version

# Check if numpy is installed
python -c "import numpy; print(numpy.__version__)"
```

## Interpreting Results from Your Data

Your analysis output shows:

```
Average Pairwise Identity: 0.0461
Fraction conserved: 0.0282
```

**Interpretation**:
- These 4 sequences are **very distantly related** (only ~4.6% average identity)
- They are likely from **different protein families** or **highly divergent orthologs**
- Very **few functionally critical residues** are conserved (only 2.8%)
- The sequences may have undergone **significant evolutionary divergence**

If you expected higher identity, check:
1. Are all sequences truly from the same family?
2. Should you be using amino acid substitution matrices rather than strict identity?
3. Consider re-analyzing with E-value cutoffs if these are BLAST results

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| "No sequences found" | Invalid FASTA format | Check file starts with `>`, has sequences on next lines |
| Very low identity (all 0) | Sequences too different | Verify sequences are homologous |
| Script hangs/slow | Large dataset | Normal for >100 sequences × >5000 residues; use external tools instead |
| ImportError: numpy | numpy not installed | Run `pip install numpy` |
| Python version error | Python < 3.11 | Update to Python 3.11+ |

## Next Steps

1. **Examine the alignment**: Open `aligned.fasta` in a text editor or alignment viewer (Jalview, SeaView)
2. **Visualize conservation**: Use WebLogo to create sequence logos from aligned sequences
3. **Downstream analysis**: 
   - Build phylogenetic tree
   - Identify orthologs using conservation patterns
   - Run structure prediction on aligned sequences
4. **Parameter tuning**: Edit gap penalties in script if needed (`gap_open`, `gap_extend` in `needleman_wunsch`)

## Performance Notes

For your test data (4 sequences, ~85 residues each):
- Runtime: < 1 second
- Memory: < 10 MB

For larger datasets:
- 10 sequences × 500 residues: ~1-2 seconds
- 50 sequences × 1000 residues: ~30-60 seconds
- 100+ sequences or >5000 residues: Consider external tools (Clustal, MAFFT)

## Additional Resources

- **FASTA format**: https://www.ncbi.nlm.nih.gov/blast/fasta.shtml
- **Sequence alignment theory**: https://en.wikipedia.org/wiki/Sequence_alignment
- **BLOSUM matrices**: https://en.wikipedia.org/wiki/BLOSUM
- **Alignment viewers**: 
  - Jalview: http://www.jalview.org
  - SeaView: http://doua.prabi.fr/software/seaview
  - WebLogo: http://weblogo.threeplusone.com

## Questions?

Refer to the full `README.md` for:
- Algorithm details
- Limitations and considerations
- Implementation notes
- Proper citations for published work
