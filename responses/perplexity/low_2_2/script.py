from Bio import SeqIO
from itertools import combinations
import numpy as np

input_fasta = "input.fasta"
records = list(SeqIO.parse(input_fasta, "fasta"))

print(f"Number of sequences: {len(records)}")
lengths = [len(r.seq) for r in records]
print(f"Min length: {min(lengths)}")
print(f"Max length: {max(lengths)}")
print(f"Mean length: {sum(lengths)/len(lengths):.2f}")

def ungapped_identity(seq1, seq2):
    # Truncate to shortest length and compute identity
    L = min(len(seq1), len(seq2))
    if L == 0:
        return 0.0
    matches = sum(a == b for a, b in zip(seq1[:L], seq2[:L]))
    return matches / L

n = len(records)
id_matrix = np.zeros((n, n))

for i, j in combinations(range(n), 2):
    pid = ungapped_identity(str(records[i].seq), str(records[j].seq))
    id_matrix[i, j] = id_matrix[j, i] = pid

names = [r.id for r in records]
print("\nPairwise identity matrix (ungapped, truncated to shortest):")
header = "\t" + "\t".join(names)
print(header)
for i, name in enumerate(names):
    row = [name] + [f"{id_matrix[i, j]:.2f}" for j in range(n)]
    print("\t".join(row))

upper_vals = id_matrix[np.triu_indices(n, k=1)]
if len(upper_vals) > 0:
    print(f"\nPairwise identity: mean={upper_vals.mean():.3f}, "
          f"min={upper_vals.min():.3f}, max={upper_vals.max():.3f}")

