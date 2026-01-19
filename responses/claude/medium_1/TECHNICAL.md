# Technical Implementation Guide

## Code Architecture

The script is organized into logical sections:

```
msa_task.py
├── FASTA I/O
│   ├── read_fasta()
│   └── write_fasta()
├── Scoring Functions
│   ├── get_blosum62_score()
│   ├── pairwise_identity()
│   └── needleman_wunsch()
├── Progressive MSA
│   ├── guide_tree_upgma()
│   ├── progressive_align()
│   ├── refine_alignment()
│   └── build_consensus()
├── Analysis
│   ├── compute_pairwise_identity_matrix()
│   ├── compute_average_pairwise_identity()
│   └── compute_conservation()
└── Main CLI Interface
    └── main()
```

## Core Algorithms

### 1. Needleman-Wunsch Global Alignment

**Purpose**: Optimally align two sequences or sequence profiles globally.

**Key Components**:
- **DP Table** (m+1 × n+1): Stores optimal scores for all subproblems
- **Traceback Table**: Records direction of optimal path (diagonal/up/left)
- **Gap Penalties**: 
  - Open: -10 (cost to start a new gap)
  - Extend: -1 (cost to extend existing gap)

**Algorithm Steps**:

```python
# Initialize
DP[0,0] = 0
DP[i,0] = gap_open + (i-1) * gap_extend  # Gap in seq2
DP[0,j] = gap_open + (j-1) * gap_extend  # Gap in seq1

# Fill DP table
for i in range(1, m+1):
    for j in range(1, n+1):
        score_match = DP[i-1,j-1] + similarity(seq1[i-1], seq2[j-1])
        score_insert = DP[i,j-1] + gap_penalty
        score_delete = DP[i-1,j] + gap_penalty
        DP[i,j] = max(score_match, score_insert, score_delete)

# Traceback from DP[m,n] to DP[0,0]
```

**Time/Space Complexity**: O(m×n) both

**BLOSUM62 Integration**:
```python
score = get_blosum62_score(aa1, aa2)
# Returns:
#  + High positive for similar amino acids
#  + Lower for dissimilar residues
#  - Negative for gaps
```

### 2. UPGMA Guide Tree Construction

**Purpose**: Build a hierarchical clustering tree for progressive alignment.

**Algorithm**:

1. **Distance Matrix**: Compute all pairwise distances as `1 - identity`
2. **Iterative Clustering**:
   - Find closest pair of clusters using average linkage
   - Merge them into new cluster
   - Remove old clusters, add new one
   - Repeat until all sequences merged

**Average Linkage Formula**:
```
dist(C_i, C_j) = (1/n_i * n_j) * Σ dist(s_a, s_b)
                 for s_a in C_i, s_b in C_j
```

**Time Complexity**: O(n³) for full hierarchical clustering

**Output**: Tree structure `[(cluster_i, cluster_j, distance), ...]`

### 3. Progressive Alignment Strategy

**Purpose**: Combine individual pairwise alignments into full MSA.

**Process**:
```
Input: 4 sequences [S1, S2, S3, S4]
Guide Tree: [(0,1), (2,3), (2,3)]

Step 1: Align S1 ↔ S2 → Aligned pair A_12
Step 2: Align S3 ↔ S4 → Aligned pair A_34  
Step 3: Align A_12 ↔ A_34 → Final MSA
```

**Key Points**:
- Early alignments influence later ones (order matters!)
- Guide tree order is crucial (determined by UPGMA)
- Can produce suboptimal MSA (heuristic, not optimal)

### 4. Alignment Refinement

**Purpose**: Improve alignment quality through iterative re-alignment.

**Strategy**:
```
For each sequence in alignment:
  1. Remove sequence from alignment
  2. Build consensus from remaining sequences
  3. Realign sequence to consensus
  4. Update alignment
Repeat for N iterations (default: 1)
```

**Benefit**: Can correct errors from initial progressive alignment

## Data Structure Design

### Sequence Representation

**Raw sequences**: Strings
```python
seq = "MKVLIVFASFAVSSAYDIENLSAETKK"
```

**Aligned sequences**: Strings with gap characters
```python
aligned = "MKVLIVFASF-AVSSAYDIENLSAETKK"
#                  ↑ gap inserted
```

### Identity Matrix

```python
identity_matrix = np.ndarray((n_sequences, n_sequences))
# [i,j] = pairwise identity between seq_i and seq_j
# Values: 0.0 to 1.0
# Symmetric: matrix[i,j] = matrix[j,i]
```

## Scoring System: BLOSUM62

**Biological Background**: BLOSUM62 (BLOcks SUbstitution Matrix)
- Derived from alignments of evolutionarily divergent proteins
- 62 = sequences averaged over 62% identity
- Standard for protein sequence alignment

**Simplified Version in Script**:
```python
# Diagonal scores (identity)
('A', 'A'): 4,
('R', 'R'): 5,
('W', 'W'): 11,  # Tryptophan highest (rarest/most specific)

# Common substitutions
('A', 'S'): 1,   # Both small, polar
('R', 'K'): 2,   # Both basic
('D', 'E'): 2,   # Both acidic

# Gaps
('-', '*'): -4   # Gap penalty
```

**Full BLOSUM62 Matrix** (20×20):
```
     A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
A    4  -1  -2  -2   0  -1  -1   0  -2  -1  -1  -1  -1  -2  -1   1   0  -3  -2   0
R   -1   5   0  -2  -3   1   0  -2   0  -3  -2   2  -1  -3  -2  -1  -1  -3  -2  -3
...
```

## Pairwise Identity Calculation

**Formula**:
```
identity = Σ(seq1[i] == seq2[i]) / n_valid_positions

where:
  - seq1[i] == seq2[i]: Matching amino acids
  - n_valid_positions: Total positions excluding gaps
```

**Implementation**:
```python
def pairwise_identity(seq1, seq2):
    matches = 0
    valid = 0
    for aa1, aa2 in zip(seq1, seq2):
        if aa1 not in ['-', '.'] and aa2 not in ['-', '.']:
            valid += 1
            if aa1.upper() == aa2.upper():
                matches += 1
    return matches / valid if valid > 0 else 0.0
```

**Why ignore gaps**:
- Gaps are artifacts of alignment, not true evolutionary information
- Inclusion would unfairly penalize alignments with more gaps
- Ignoring gaps focuses on conserved residues

## Conservation Analysis

**Concept**: Identify alignment columns that are identical across all sequences.

**Calculation**:
```python
for column_pos in range(alignment_length):
    residues_in_column = set()
    has_gap = False
    
    for sequence in aligned_sequences:
        aa = sequence[column_pos]
        if aa in ['-', '.']:
            has_gap = True
            break
        residues_in_column.add(aa)
    
    if not has_gap and len(residues_in_column) == 1:
        count_conserved += 1
    if not has_gap:
        count_valid += 1

fraction = count_conserved / count_valid
```

**Biological Significance**:
- **Highly conserved columns** (>0.8 fraction): Likely functional sites
- **Moderately conserved** (0.3-0.8): Structurally important
- **Low conservation** (<0.3): Rapidly evolving regions

## Gap Penalties

**Affine Gap Model** (used here):
```
Gap cost = gap_open + (gap_length - 1) * gap_extend

Examples:
  Single gap: -10 + 0*(-1) = -10
  Two gaps in row: -10 + 1*(-1) = -11
  Three gaps: -10 + 2*(-1) = -12
```

**Why affine gaps**:
- Single gap opening cost discourages excessive gapping
- Lower extension cost allows long gaps (common in biology)
- More realistic than linear gap costs

**Tuning Parameters**:
```python
# Conservative (fewer gaps)
gap_open=-15, gap_extend=-2

# Liberal (more gaps allowed)
gap_open=-5, gap_extend=-0.5

# Current (balanced)
gap_open=-10, gap_extend=-1
```

## Complexity Analysis

### Time Complexity

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| Read FASTA | O(total_length) | Linear scan |
| Pairwise distances | O(n² × L) | n=sequences, L=length |
| UPGMA tree | O(n³) | Naive O(n²) per iteration × n iterations |
| Needleman-Wunsch × n | O(n × L²) | Single alignment O(L²), × n alignments |
| Conservation calc | O(n × L) | Single pass through alignment |
| **Total** | **O(n³ + n × L²)** | Dominated by tree building for large n |

### Space Complexity

| Data Structure | Space |
|---|---|
| Sequence storage | O(n × L) |
| DP table (NW) | O(L²) |
| Distance matrix | O(n²) |
| Identity matrix | O(n²) |
| Aligned sequences | O(n × L_aligned) |
| **Total** | **O(n × L + L² + n²)** |

**For typical data** (n=10, L=500):
- Space: ~50 KB
- Time: ~100 ms (pairwise bottleneck)

## Error Handling

```python
# File I/O errors
try:
    with open(filename) as f:
        # read/write
except IOError as e:
    print(f"Error: {e}", file=sys.stderr)
    sys.exit(1)

# Format validation
if not sequences:
    print("Error: No sequences found", file=sys.stderr)
    sys.exit(1)

# Safe indexing
for seq in aligned_seqs:
    if pos < len(seq):  # Boundary check
        aa = seq[pos]
```

## Performance Optimization Tips

### If Script is Slow

1. **Reduce sequences**: Group by species/family first
2. **Reduce length**: Use domain-only sequences (e.g., HMM profiles)
3. **Skip refinement**: Remove `refine_alignment()` call
4. **Lower gap penalty**: Fewer gaps = less DP table filling

### Vectorization Opportunities

Could speed up with:
- NumPy vectorized DP table filling
- Numba JIT compilation for tight loops
- Parallel pairwise distance computation

## Extending the Code

### Add Custom Scoring Matrix

```python
def get_custom_matrix_score(aa1, aa2):
    custom_matrix = {
        # your custom scores
    }
    return custom_matrix.get((aa1, aa2), -1)

# Use in needleman_wunsch:
score = get_custom_matrix_score(seq1[i-1], seq2[j-1])
```

### Add Local Alignment (Smith-Waterman)

```python
# Modify NW: Initialize all cells to 0 (not penalties)
dp[i, 0] = 0  # Instead of gap_open + (i-1)*gap_extend
dp[0, j] = 0  # Instead of gap_open + (j-1)*gap_extend

# Traceback: Start from maximum value (not DP[m,n])
max_score = np.max(dp)
i, j = np.unravel_index(np.argmax(dp), dp.shape)
```

### Output Conservation Score Per Position

```python
def get_conservation_scores(aligned_seqs):
    scores = []
    for pos in range(len(aligned_seqs[0])):
        residues = [seq[pos] for seq in aligned_seqs if pos < len(seq)]
        unique = len(set(residues))
        score = 1 - (unique - 1) / len(residues)  # 1=conserved, 0=variable
        scores.append(score)
    return scores
```

## Dependencies Breakdown

### numpy
- **Use**: Dynamic programming DP table storage and operations
- **Alternative**: Could use lists, but much slower
- **Specific functions**:
  - `np.full()`: Initialize DP table
  - `np.ndarray`: 2D array creation
  - `np.unravel_index()`: Find position of max value

### Standard library (no additional dependencies needed!)
- `argparse`: Command-line argument parsing
- `sys`: File handle and error reporting
- `typing`: Type hints for clarity
- `collections`: Not explicitly used, but could for efficiency

## Testing Recommendations

```bash
# Test 1: Basic functionality
python msa_task.py -i test.fasta -o output.fasta

# Test 2: Check output format
head output.fasta

# Test 3: Verify alignment lengths
python -c "
from msa_task import read_fasta
headers, seqs = zip(*read_fasta('output.fasta'))
print(f'Length: {[len(s) for s in seqs]}')
print(f'All equal: {len(set(len(s) for s in seqs)) == 1}')
"

# Test 4: Sanity checks
python -c "
from msa_task import read_fasta, pairwise_identity
headers, seqs = zip(*read_fasta('output.fasta'))
identity = pairwise_identity(seqs[0], seqs[0])
print(f'Self-identity: {identity}')  # Should be 1.0
"
```

## References

### Key Papers

1. **Needleman-Wunsch (1970)**
   - Global sequence alignment
   - Needleman, S.B. and Wunsch, C.D. "A general method applicable to the search for similarities in the amino acid sequence of two proteins"
   - J. Mol. Biol., 48(3):443-453

2. **BLOSUM Matrices (1992)**
   - Substitution scoring matrices
   - Henikoff, S. and Henikoff, J.G. "Amino acid substitution matrices from protein blocks"
   - PNAS, 89(22):10915-10919

3. **UPGMA (1965)**
   - Hierarchical clustering
   - Sokal, R.R. and Michener, C.D. "A statistical method for evaluating systematic relationships"
   - University of Kansas Science Bulletin, 38:1409-1438

### Python Resources

- NumPy Documentation: https://numpy.org/doc/
- Python argparse: https://docs.python.org/3/library/argparse.html
- Type Hints: https://docs.python.org/3/library/typing.html
