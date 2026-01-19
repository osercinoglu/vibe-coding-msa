#!/usr/bin/env python

import subprocess
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
import itertools
import numpy as np

# -----------------------------
# Input / output files
# -----------------------------
input_fasta = "input.fasta"
aligned_fasta = "aligned.fasta"

# -----------------------------
# 1. Run MUSCLE
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
# 2. Read alignment
# -----------------------------
alignment = AlignIO.read(aligned_fasta, "fasta")
nseq = len(alignment)
aln_len = alignment.get_alignment_length()

print(f"\nAlignment contains {nseq} sequences")
print(f"Alignment length: {aln_len}")

# -----------------------------
# 3. Pairwise sequence identity
# -----------------------------
def pairwise_identity(seq1, seq2):
    matches = 0
    compared = 0
    for a, b in zip(seq1, seq2):
        if a == "-" or b == "-":
            continue
        compared += 1
        if a == b:
            matches += 1
    return matches / compared if compared > 0 else 0.0

ids = [rec.id for rec in alignment]
pid_matrix = np.zeros((nseq, nseq))

for i, j in itertools.combinations(range(nseq), 2):
    pid = pairwise_identity(alignment[i].seq, alignment[j].seq)
    pid_matrix[i, j] = pid
    pid_matrix[j, i] = pid

np.fill_diagonal(pid_matrix, 1.0)

print("\nPairwise identity matrix:")
print("\t" + "\t".join(ids))
for i, row in enumerate(pid_matrix):
    values = "\t".join(f"{v:.3f}" for v in row)
    print(f"{ids[i]}\t{values}")

# -----------------------------
# 4. Average pairwise identity
# -----------------------------
upper_triangle = [
    pid_matrix[i, j]
    for i in range(nseq)
    for j in range(i + 1, nseq)
]

avg_pid = sum(upper_triangle) / len(upper_triangle)
print(f"\nAverage pairwise identity: {avg_pid:.3f}")

# -----------------------------
# 5. Fully conserved columns
# -----------------------------
fully_conserved = 0
for col in range(aln_len):
    residues = [
        rec.seq[col]
        for rec in alignment
        if rec.seq[col] != "-"
    ]
    if len(residues) == 0:
        continue
    if len(set(residues)) == 1:
        fully_conserved += 1

fraction_conserved = fully_conserved / aln_len
print(f"Fully conserved columns: {fully_conserved}/{aln_len}")
print(f"Fraction fully conserved: {fraction_conserved:.3f}")

