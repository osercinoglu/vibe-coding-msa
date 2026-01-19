from Bio import SeqIO
from itertools import combinations
import numpy as np

input_fasta = "input.fasta"
aligned_fasta = "aligned.fasta"

# 1. Read original sequences
records = list(SeqIO.parse(input_fasta, "fasta"))
print(f"Number of sequences: {len(records)}")

lengths = [len(r.seq) for r in records]
print(f"Min length: {min(lengths)}")
print(f"Max length: {max(lengths)}")
print(f"Mean length: {sum(lengths)/len(lengths):.2f}")

# 2. Read aligned sequences
aligned_records = list(SeqIO.parse(aligned_fasta, "fasta"))
if len(aligned_records) != len(records):
    raise ValueError("Number of sequences changed between input and alignment.")

alignment_length = len(aligned_records[0].seq)
print(f"Alignment length: {alignment_length}")

# 3. Compute pairwise sequence identity on the alignment
def pairwise_identity(seq1, seq2, ignore_gaps=True):
    """
    Simple percent identity: matches / positions (optionally excluding gap-only positions).
    seq1 and seq2 must be same length (aligned).
    """
    if len(seq1) != len(seq2):
        raise ValueError("Aligned sequences must have same length.")

    matches = 0
    total = 0
    for a, b in zip(seq1, seq2):
        if ignore_gaps and (a == "-" or b == "-"):
            continue
        if a == b:
            matches += 1
        if not ignore_gaps or (a != "-" and b != "-"):
            total += 1

    return 0.0 if total == 0 else matches / total

n = len(aligned_records)
id_matrix = np.zeros((n, n))

for i, j in combinations(range(n), 2):
    pid = pairwise_identity(str(aligned_records[i].seq), str(aligned_records[j].seq))
    id_matrix[i, j] = id_matrix[j, i] = pid

# 4. Simple summary on pairwise identities
upper_vals = id_matrix[np.triu_indices(n, k=1)]
if len(upper_vals) > 0:
    print(f"Pairwise identity: mean={upper_vals.mean():.3f}, "
          f"min={upper_vals.min():.3f}, max={upper_vals.max():.3f}")
else:
    print("Only one sequence; no pairwise identities to compute.")

# 5. Optional: print a small identity table
names = [r.id for r in aligned_records]
print("\nPairwise identity matrix:")
header = "\t" + "\t".join(names)
print(header)
for i, name in enumerate(names):
    row = [name] + [f"{id_matrix[i, j]:.2f}" for j in range(n)]
    print("\t".join(row))

# 6. Optional: per-column conservation
# Fraction of most frequent non-gap residue at each column
alignment_array = [list(str(r.seq)) for r in aligned_records]
column_conservation = []

for col_idx in range(alignment_length):
    col = [row[col_idx] for row in alignment_array]
    non_gaps = [c for c in col if c != "-"]
    if not non_gaps:
        column_conservation.append(0.0)
        continue
    counts = {}
    for aa in non_gaps:
        counts[aa] = counts.get(aa, 0) + 1
    max_freq = max(counts.values()) / len(non_gaps)
    column_conservation.append(max_freq)

print(f"\nMean column conservation (0–1): {np.mean(column_conservation):.3f}")

