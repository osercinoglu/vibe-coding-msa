#!/usr/bin/env python3
"""
msa_task.py

Simple multiple sequence alignment and identity statistics in pure Python
(+numpy).

Run as:
    python msa_task.py --input input.fasta --output aligned.fasta
"""

import argparse
import sys
from itertools import combinations
import numpy as np


# ------------------------- FASTA I/O ------------------------- #

def read_fasta(path):
    """Read a protein FASTA file, return list of (header, sequence)."""
    records = []
    header = None
    seq_lines = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.replace(" ", "").upper())
        if header is not None:
            records.append((header, "".join(seq_lines)))
    return records


def write_fasta(path, records):
    """Write list of (header, aligned_sequence) to FASTA."""
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            # Wrap to 60 chars for readability
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")


# ------------------------- Needleman–Wunsch ------------------------- #

def nw_global_align(seq1, seq2, match=1, mismatch=-1, gap=-1):
    """
    Global alignment (Needleman–Wunsch) of two sequences.
    Returns aligned_seq1, aligned_seq2.
    Very simple scoring: +match, mismatch, gap penalties.
    """
    n = len(seq1)
    m = len(seq2)

    # DP matrix
    score = np.zeros((n + 1, m + 1), dtype=int)
    # Traceback: 0 = diag, 1 = up, 2 = left
    traceback = np.zeros((n + 1, m + 1), dtype=np.int8)

    # Initialization
    for i in range(1, n + 1):
        score[i, 0] = score[i - 1, 0] + gap
        traceback[i, 0] = 1  # up
    for j in range(1, m + 1):
        score[0, j] = score[0, j - 1] + gap
        traceback[0, j] = 2  # left

    # Fill
    for i in range(1, n + 1):
        c1 = seq1[i - 1]
        for j in range(1, m + 1):
            c2 = seq2[j - 1]
            s = match if c1 == c2 else mismatch
            diag = score[i - 1, j - 1] + s
            up = score[i - 1, j] + gap
            left = score[i, j - 1] + gap
            best = diag
            tb = 0
            if up > best:
                best = up
                tb = 1
            if left > best:
                best = left
                tb = 2
            score[i, j] = best
            traceback[i, j] = tb

    # Traceback
    aligned1 = []
    aligned2 = []
    i, j = n, m
    while i > 0 or j > 0:
        tb = traceback[i, j]
        if i > 0 and j > 0 and tb == 0:
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or tb == 1):
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
        else:
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1

    return "".join(reversed(aligned1)), "".join(reversed(aligned2))


# ------------------------- Profile alignment ------------------------- #

def align_seq_to_profile(profile_seqs, new_seq, match=1, mismatch=-1, gap=-1):
    """
    Align a single sequence to an existing alignment profile.

    profile_seqs: list of aligned sequences (all same length).
    new_seq: unaligned sequence string.

    Strategy:
    - Define a consensus-like reference by majority vote (excluding gaps).
    - Align new_seq to this reference with NW.
    - Propagate gaps introduced in the reference alignment back into
      all profile sequences.

    Returns updated_profile_seqs (including new sequence as aligned).
    """
    # Build a simple consensus (excluding gaps)
    profile_length = len(profile_seqs[0])
    consensus_chars = []
    for col in range(profile_length):
        column = [s[col] for s in profile_seqs]
        residues = [c for c in column if c != "-"]
        if not residues:
            consensus_chars.append("-")
        else:
            # Majority vote
            counts = {}
            for r in residues:
                counts[r] = counts.get(r, 0) + 1
            consensus_chars.append(max(counts.items(), key=lambda x: x[1])[0])
    consensus = "".join(consensus_chars)

    # Align new_seq to consensus
    aln_cons, aln_new = nw_global_align(consensus, new_seq,
                                        match=match, mismatch=mismatch, gap=gap)

    # Propagate gaps (where aln_cons has '-') to all profile sequences
    new_profile = []
    for seq in profile_seqs:
        seq_idx = 0
        aligned_seq_chars = []
        for c in aln_cons:
            if c == "-":
                aligned_seq_chars.append("-")
            else:
                aligned_seq_chars.append(seq[seq_idx])
                seq_idx += 1
        new_profile.append("".join(aligned_seq_chars))

    # Add aligned new sequence
    new_profile.append(aln_new)
    return new_profile


def progressive_msa(seqs, match=1, mismatch=-1, gap=-1):
    """
    Very simple progressive MSA:
    - Take first sequence as seed.
    - Align each remaining sequence to current profile one by one.
    seqs: list of sequence strings (unaligned).
    Returns list of aligned sequences (same order as input).
    """
    if len(seqs) == 1:
        return seqs[:]

    # Start profile with the first two sequences
    aln1, aln2 = nw_global_align(seqs[0], seqs[1],
                                 match=match, mismatch=mismatch, gap=gap)
    profile = [aln1, aln2]

    for s in seqs[2:]:
        profile = align_seq_to_profile(profile, s,
                                       match=match, mismatch=mismatch, gap=gap)

    # profile currently holds aligned sequences in the order they were added:
    # [seq0, seq1, seq2, ...] already matching original seqs order.
    return profile


# ------------------------- Identity statistics ------------------------- #

def pairwise_identity(seq_a, seq_b):
    """
    Compute pairwise sequence identity (%) between two aligned sequences.
    Identity = matches / positions where neither is a gap.
    """
    assert len(seq_a) == len(seq_b)
    matches = 0
    comp = 0
    for a, b in zip(seq_a, seq_b):
        if a == "-" or b == "-":
            continue
        comp += 1
        if a == b:
            matches += 1
    if comp == 0:
        return 0.0
    return 100.0 * matches / comp


def compute_identity_matrix(aligned_seqs):
    """Return NxN identity matrix for aligned sequences."""
    n = len(aligned_seqs)
    mat = np.zeros((n, n), dtype=float)
    for i in range(n):
        mat[i, i] = 100.0
        for j in range(i + 1, n):
            pid = pairwise_identity(aligned_seqs[i], aligned_seqs[j])
            mat[i, j] = pid
            mat[j, i] = pid
    return mat


def average_pairwise_identity(identity_matrix):
    """Average over all unique i < j pairs."""
    n = identity_matrix.shape[0]
    if n < 2:
        return 100.0
    vals = []
    for i in range(n):
        for j in range(i + 1, n):
            vals.append(identity_matrix[i, j])
    if not vals:
        return 0.0
    return float(np.mean(vals))


def fraction_fully_conserved_columns(aligned_seqs):
    """
    Fraction of alignment columns that are fully conserved.
    A column is fully conserved if all sequences have the same residue
    and none is a gap.
    """
    n = len(aligned_seqs)
    if n == 0:
        return 0.0
    length = len(aligned_seqs[0])
    conserved = 0
    for col in range(length):
        column = [s[col] for s in aligned_seqs]
        if "-" in column:
            continue
        # all equal?
        first = column[0]
        if all(c == first for c in column[1:]):
            conserved += 1
    return conserved / length if length > 0 else 0.0


# ------------------------- Main CLI ------------------------- #

def main():
    parser = argparse.ArgumentParser(
        description="Simple MSA and identity statistics (pure Python + numpy)."
    )
    parser.add_argument("--input", required=True,
                        help="Input unaligned protein FASTA")
    parser.add_argument("--output", required=True,
                        help="Output aligned FASTA")
    args = parser.parse_args()

    records = read_fasta(args.input)
    if not records:
        print("No sequences found in input FASTA.", file=sys.stderr)
        sys.exit(1)

    headers = [h for h, _ in records]
    seqs = [s for _, s in records]

    # Multiple sequence alignment
    aligned_seqs = progressive_msa(seqs, match=1, mismatch=-1, gap=-1)

    # Sanity check lengths
    aln_len = len(aligned_seqs[0])
    if any(len(s) != aln_len for s in aligned_seqs):
        print("Internal error: aligned sequences have different lengths.",
              file=sys.stderr)
        sys.exit(1)

    # Write alignment
    write_fasta(args.output, list(zip(headers, aligned_seqs)))

    # Pairwise identity matrix
    id_mat = compute_identity_matrix(aligned_seqs)
    avg_pid = average_pairwise_identity(id_mat)
    frac_conserved = fraction_fully_conserved_columns(aligned_seqs)

    # Print results
    print("Pairwise sequence identity matrix (%)")
    # Header row
    print("\t" + "\t".join(headers))
    for i, h in enumerate(headers):
        row_vals = "\t".join(f"{id_mat[i, j]:.2f}" for j in range(len(headers)))
        print(f"{h}\t{row_vals}")

    print(f"\nAverage pairwise identity: {avg_pid:.2f}%")
    print(f"Fraction of fully conserved columns: {frac_conserved:.4f}")


if __name__ == "__main__":
    main()

