# MSA Script - Quick Reference Card

## One-Liner Usage
```bash
python msa_task.py --input input.fasta --output aligned.fasta
```

## Complete Command Structure
```bash
python msa_task.py -i INPUT_FILE -o OUTPUT_FILE
python msa_task.py --input INPUT_FILE --output OUTPUT_FILE
python msa_task.py -h  # Show help
```

## What Gets Output

**To Console (stderr):**
```
Pairwise Sequence Identity Matrix:
              seq1      seq2      seq3
seq1        1.0000    0.5230    0.4821
seq2        0.5230    1.0000    0.6123
seq3        0.4821    0.6123    1.0000

Average Pairwise Identity: 0.5393

Conservation Analysis:
  Fully conserved columns: 42 / 287
  Fraction conserved: 0.1464
```

**To File (output.fasta):**
```
>seq1
MKVLIVFASF-AVSSAYDIENLSAETKK
>seq2
MAVLIVFTSFSVSNAVDVENLSAETSK-
>seq3
MKVLIVFTSFSVQNAV-DVENLSAETSK
```

## Understanding the Output

| Metric | Interpretation |
|--------|-----------------|
| Identity 0.95 | 95% of positions match |
| Identity 0.50 | 50% similarity |
| Identity 0.05 | Very distant (5%) |
| **Avg Identity 0.7+** | Close homologs |
| **Avg Identity 0.3-0.7** | Divergent but related |
| **Avg Identity <0.3** | Very distant or unrelated |
| **Fraction Conserved >0.5** | Well-conserved alignment |
| **Fraction Conserved <0.1** | Highly variable sequences |

## File Descriptions

| File | Read/Write | Content |
|------|-----------|---------|
| input.fasta | READ | Your unaligned sequences |
| output.fasta | WRITE | Aligned sequences with gaps |

## Input Format (FASTA)
```
>Sequence_Header_1
MKVLIVFASFAVSSAYDIENLSAETKK
RELPEKKKM
>Sequence_Header_2
MAVLIVFTSFSVSNAVDVENLSAETSKRELPEKKKM
```

## Output Format (FASTA - Aligned)
```
>Sequence_Header_1
MKVLIVFASF-AVSSAYDIENLSAETKK
RELPEKKKM-
>Sequence_Header_2
MAVLIVFTSFSVSNAVDVENLSAETSKR
ELPEKKKM--
```

## Algorithms Used (In Order)

1. **Read FASTA** - Parse input file
2. **Compute Distances** - O(n² × L) for n sequences, length L
3. **Build UPGMA Tree** - O(n³) hierarchical clustering
4. **Progressive Alignment** - Needleman-Wunsch O(n × L²)
5. **Refinement** - Re-align sequences to profile
6. **Analysis** - Compute matrices and statistics
7. **Output** - Write FASTA and print metrics

## Needleman-Wunsch Algorithm
- **Type**: Global pairwise alignment
- **Scoring**: BLOSUM62 matrix + affine gaps
- **Gap Cost**: Open = -10, Extend = -1
- **Output**: Optimal alignment of two sequences

## Key Scores in BLOSUM62
```
Match (identity):              +4 to +11
Similar amino acids:           +1 to +2
Different amino acids:         -1
Gap opening:                   -10
Gap extension:                 -1
```

## Common Issues & Quick Fixes

| Problem | Check | Fix |
|---------|-------|-----|
| No sequences found | FASTA format | Ensure `>` headers before sequences |
| All zeros identity | Sequence relationship | Verify sequences are homologous |
| Script too slow | Dataset size | Reduce sequences or use external tools |
| NumPy error | Installation | `pip install numpy` |
| Python version | Python installed | Need Python 3.11+ |

## Performance Expectations

```
4 sequences × 85 residues  → <1 second, <50 MB
10 sequences × 500 residues → 1-2 seconds
50 sequences × 1000 residues → 30-60 seconds
100+ sequences → Consider Clustal/MAFFT instead
```

## Modifying the Script

Most common edits (in msa_task.py):

```python
# Change gap penalties (line ~75)
gap_open=-15    # More conservative (fewer gaps)
gap_extend=-0.5 # More liberal (longer gaps)

# Skip refinement for speed (line ~437)
# aligned_seqs = refine_alignment(aligned_seqs, iterations=1)

# Change number of refinement iterations
aligned_seqs = refine_alignment(aligned_seqs, iterations=3)
```

## Pairwise Identity Formula

```
Identity = Matches / Valid Positions

Valid Positions = Total positions excluding gaps
Matches = Number of identical residues at each position

Example:
Seq1: MKVL-IFAS
Seq2: MKVLFIFAT
      *** **** = 8 matches / 9 valid = 0.889 (88.9% identity)
```

## Conservation Formula

```
Conservation = Fully Conserved Columns / Total Valid Columns

Valid Columns = Columns with NO gaps in any sequence
Fully Conserved = Columns where ALL sequences have SAME amino acid

Example: 42 conserved / 287 total = 0.1464 (14.64% conserved)
```

## FASTA File Validation

A valid FASTA file has:
- ✓ Lines starting with `>` as headers
- ✓ Sequence on following lines (single or multiple)
- ✓ No spaces in sequences (replace with nothing)
- ✓ One or more sequences
- ✓ Same case letters (handles both upper and lower)

## Results Output Locations

```
Script console (stderr):
  → Pairwise identity matrix
  → Average pairwise identity
  → Conservation analysis

Output file (fasta):
  → Aligned sequences with gaps
  → Headers preserved from input
  → Ready for downstream analysis
```

## Using Output for Downstream Analysis

```bash
# Phylogenetic tree
FastTree -protein < aligned.fasta > tree.nwk

# Structure prediction  
AlphaFold2 -fasta aligned.fasta -output results/

# Domain identification
InterProScan -i aligned.fasta -f json

# Sequence logo
WebLogo -f aligned.fasta -o logo.png

# Multiple tools
clustalo -i aligned.fasta --outfmt=phylip > tree.phy
```

## Citation Info

In publications, cite:
- Needleman & Wunsch (1970) - Global alignment algorithm
- Henikoff & Henikoff (1992) - BLOSUM62 matrix

Example text: "Multiple sequence alignment was performed using progressive alignment with Needleman-Wunsch algorithm (Needleman & Wunsch, 1970) and BLOSUM62 scoring matrix (Henikoff & Henikoff, 1992)."

## Typical Workflow

```
1. Prepare input sequences
   input.fasta
   ↓
2. Run alignment
   python msa_task.py -i input.fasta -o aligned.fasta
   ↓
3. Check metrics in console output
   - Identity scores make sense?
   - Conservation level reasonable?
   ↓
4. Examine alignment in aligned.fasta
   - Look for gapped regions
   - Identify conserved blocks
   ↓
5. Use aligned.fasta for downstream analysis
   - Build phylogenetic trees
   - Predict structure
   - Identify domains
   - Run signature detection
```

## Debugging Checklist

- [ ] Python installed? `python --version` → Need 3.11+
- [ ] NumPy installed? `python -c "import numpy"`
- [ ] FASTA format valid? `head -3 input.fasta` → See `>` headers?
- [ ] File readable? Check permissions with `ls -l input.fasta`
- [ ] Sequences homologous? Check BLAST results or literature
- [ ] Output file writable? Check directory permissions
- [ ] Try example: `cp test_input.fasta mytest.fasta`

## Quick Stats Commands

```bash
# Count sequences in FASTA
grep -c "^>" input.fasta

# Count residues per sequence
grep -v "^>" input.fasta | wc -c

# Check for non-standard characters
grep -v "^>" input.fasta | grep -o . | sort | uniq -c

# Verify alignment (all same length)
awk '/^>/{if(seq)print length(seq); seq=""} !/^>/{seq=seq$0} END{print length(seq)}' aligned.fasta
```

## Key Formulas Quick Reference

| Formula | Code |
|---------|------|
| Pairwise Identity | `matches / valid_positions` |
| Average Identity | `sum(off_diagonal) / n_pairs` |
| Fraction Conserved | `conserved_cols / total_cols` |
| DP Score | `max(match_score, gap_score)` |
| Gap Cost | `gap_open + (length-1)*gap_extend` |

## Memory Usage Estimate

```
Input sequences: n sequences, L residues each
DP matrix: L × L ≈ L² numbers (8 bytes each)
Aligned output: n × L_aligned (1 byte per position)
Total: ~8*L² + n*L_aligned bytes

Example: 4 × 85 residues
DP: 8*85² ≈ 58 KB
Output: 4*85 ≈ 0.3 KB
Total: ~60 KB
```

---

**For more details, see the .md documentation files!**
- Quick help → QUICKSTART.md
- Full explanation → README.md
- Implementation → TECHNICAL.md
- Navigation → INDEX.md
