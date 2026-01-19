import subprocess
from Bio import AlignIO, SeqIO
import numpy as np
from itertools import combinations

# -----------------------------
# Input / output files
# -----------------------------
input_fasta = "input.fasta"
aligned_fasta = "aligned.fasta"

# -----------------------------
# 1. Run MUSCLE to create MSA
# -----------------------------
cmd = [
    "muscle",
    "-align", input_fasta,
    "-output", aligned_fasta
]

print("Running MUSCLE...")
subprocess.run(cmd, check=True)
print(f"Alignment written to {aligned_fasta}")

# -----------------------------
# 2. Read the alignment
# -----------------------------
alignment = AlignIO.read(aligned_fasta, "fasta")
n_seqs = len(alignment)
aln_length = alignment.get_alignment_length()

print(f"\nNumber of sequences: {n_seqs}")
print(f"Alignment length: {aln_length}")

# -----------------------------
# 3. Pairwise percent identity
# -----------------------------
def percent_identity(seq1, seq2):
    """
    Compute percent identity between two aligned sequences.
    Gaps ('-') are ignored.
    """
    matches = 0
    compared = 0

    for a, b in zip(seq1, seq2):
        if a == "-" or b == "-":
            continue
        compared += 1
        if a == b:
            matches += 1

    if compared == 0:
        return 0.0

    return 100.0 * matches / compared

# Compute all pairwise identities
pairwise_ids = []

for rec1, rec2 in combinations(alignment, 2):
    pid = percent_identity(str(rec1.seq), str(rec2.seq))
    pairwise_ids.append(pid)

pairwise_ids = np.array(pairwise_ids)

# -----------------------------
# 4. Summary statistics
# -----------------------------
print("\nPairwise percent identity statistics:")
print(f"  Mean: {pairwise_ids.mean():.2f}%")
print(f"  Std:  {pairwise_ids.std():.2f}%")
print(f"  Min:  {pairwise_ids.min():.2f}%")
print(f"  Max:  {pairwise_ids.max():.2f}%")

# Optional: show full matrix
print("\nPairwise identities:")
for (rec1, rec2), pid in zip(combinations(alignment, 2), pairwise_ids):
    print(f"{rec1.id} vs {rec2.id}: {pid:.2f}%")

