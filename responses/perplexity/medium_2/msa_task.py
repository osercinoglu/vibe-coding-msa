#!/usr/bin/env python3
import argparse
import sys
from itertools import combinations
import numpy as np


def read_fasta(path):
    """Read a FASTA file, return list of (header, sequence)."""
    records = []
    header = None
    seq_chunks = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks)))
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        records.append((header, "".join(seq_chunks)))
    return records


def write_fasta(records, path, width=60):
    """Write records = list of (header, sequence) to FASTA."""
    with open(path, "w") as f:
        for h, s in records:
            f.write(f">{h}\n")
            for i in range(0, len(s), width):
                f.write(s[i:i+width] + "\n")


def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1):
    """
    Global alignment (Needleman-Wunsch) with simple linear scoring.
    Returns aligned_seq1, aligned_seq2.
    """
    n = len(seq1)
    m = len(seq2)
    # Score matrix
    score = np.zeros((n + 1, m + 1), dtype=int)
    # Traceback matrix: 0 = diag, 1 = up, 2 = left
    trace = np.zeros((n + 1, m + 1), dtype=np.int8)

    # Initialize
    for i in range(1, n + 1):
        score[i, 0] = score[i - 1, 0] + gap
        trace[i, 0] = 1  # up
    for j in range(1, m + 1):
        score[0, j] = score[0, j - 1] + gap
        trace[0, j] = 2  # left

    # Fill
    for i in range(1, n + 1):
        ci = seq1[i - 1]
        for j in range(1, m + 1):
            cj = seq2[j - 1]
            diag = score[i - 1, j - 1] + (match if ci == cj else mismatch)
            up = score[i - 1, j] + gap
            left = score[i, j - 1] + gap
            best = max(diag, up, left)
            score[i, j] = best
            if best == diag:
                trace[i, j] = 0
            elif best == up:
                trace[i, j] = 1
            else:
                trace[i, j] = 2

    # Traceback
    i, j = n, m
    aln1 = []
    aln2 = []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and trace[i, j] == 0:
            aln1.append(seq1[i - 1])
            aln2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or trace[i, j] == 1):
            aln1.append(seq1[i - 1])
            aln2.append("-")
            i -= 1
        else:
            aln1.append("-")
            aln2.append(seq2[j - 1])
            j -= 1

    return "".join(reversed(aln1)), "".join(reversed(aln2))


def profile_to_string(profile_rows):
    """Convert a list of equal-length aligned sequences to a 'consensus-like' string
    used only for alignment guidance (not a real consensus)."""
    # Here we just take the most frequent non-gap character in each column;
    # if all gaps, use '-'.
    if not profile_rows:
        return ""
    length = len(profile_rows[0])
    cols = []
    for c in range(length):
        col_chars = [row[c] for row in profile_rows]
        non_gaps = [x for x in col_chars if x != "-"]
        if not non_gaps:
            cols.append("-")
        else:
            # simple majority
            unique, counts = np.unique(non_gaps, return_counts=True)
            cols.append(unique[np.argmax(counts)])
    return "".join(cols)


def align_profile_with_sequence(profile_rows, seq):
    """
    Align an existing profile (list of aligned sequences) with a new unaligned sequence.
    Returns a new list of profile rows (with gaps inserted) and aligned new sequence.
    """
    # Turn profile into a guide string
    guide = profile_to_string(profile_rows)
    guide_aln, seq_aln = needleman_wunsch(guide, seq)

    # Now insert gaps into each profile row according to guide_aln vs original guide
    new_profile = []
    for row in profile_rows:
        new_row = []
        idx = 0  # index into original row (same as into guide)
        for g_char in guide_aln:
            if g_char == "-":
                # gap introduced in guide -> gap in all existing rows
                new_row.append("-")
            else:
                new_row.append(row[idx])
                idx += 1
        new_profile.append("".join(new_row))

    return new_profile, seq_aln


def progressive_msa(seqs):
    """
    Very simple progressive MSA:
    - Start from the first sequence as the initial profile.
    - Iteratively align each next sequence to the profile.
    seqs: list of sequences (strings).
    Returns a list of aligned sequences in the same order.
    """
    if not seqs:
        return []
    # Initialize profile with first sequence
    profile = [seqs[0]]

    # Iteratively add sequences
    for s in seqs[1:]:
        profile, s_aln = align_profile_with_sequence(profile, s)
        profile.append(s_aln)

    return profile


def pairwise_identity(seq_a, seq_b):
    """
    Compute pairwise identity between two aligned or unaligned sequences.
    For unaligned input, sequences will be padded to the same length with gaps
    before calling this function; here we assume equal length.
    Identity is: matches / positions where neither is a gap.
    """
    if len(seq_a) != len(seq_b):
        raise ValueError("Sequences must be same length for identity computation")
    matches = 0
    comparable = 0
    for a, b in zip(seq_a, seq_b):
        if a == "-" or b == "-":
            continue
        comparable += 1
        if a == b:
            matches += 1
    if comparable == 0:
        return 0.0
    return matches / comparable


def compute_pairwise_identity_matrix(aligned_seqs):
    """
    Compute pairwise identity matrix from a list of aligned sequences (strings).
    Returns a numpy array (n x n) and list of indices (0..n-1).
    """
    n = len(aligned_seqs)
    mat = np.zeros((n, n), dtype=float)
    for i in range(n):
        mat[i, i] = 1.0
        for j in range(i + 1, n):
            pid = pairwise_identity(aligned_seqs[i], aligned_seqs[j])
            mat[i, j] = pid
            mat[j, i] = pid
    return mat


def fraction_fully_conserved_columns(aligned_seqs):
    """
    Fraction of columns where all non-gap residues are identical and at least one residue.
    """
    if not aligned_seqs:
        return 0.0
    n = len(aligned_seqs)
    L = len(aligned_seqs[0])
    conserved = 0
    for pos in range(L):
        column = [aligned_seqs[i][pos] for i in range(n)]
        non_gaps = [c for c in column if c != "-"]
        if not non_gaps:
            continue
        if all(c == non_gaps[0] for c in non_gaps):
            conserved += 1
    if L == 0:
        return 0.0
    return conserved / L


def main():
    parser = argparse.ArgumentParser(description="Simple MSA and identity statistics.")
    parser.add_argument("--input", required=True, help="Input unaligned protein FASTA")
    parser.add_argument("--output", required=True, help="Output aligned FASTA")
    args = parser.parse_args()

    records = read_fasta(args.input)
    if not records:
        print("No sequences found in input.", file=sys.stderr)
        sys.exit(1)

    headers = [h for h, s in records]
    seqs = [s for h, s in records]

    # Multiple sequence alignment
    aligned_seqs = progressive_msa(seqs)

    # Sanity: ensure all same length
    aln_len = len(aligned_seqs[0])
    if any(len(s) != aln_len for s in aligned_seqs):
        print("Internal error: aligned sequences have different lengths.", file=sys.stderr)
        sys.exit(1)

    # Write aligned FASTA
    write_fasta(list(zip(headers, aligned_seqs)), args.output)

    # Pairwise identity matrix
    pid_mat = compute_pairwise_identity_matrix(aligned_seqs)

    # Average pairwise identity (over unique pairs)
    n = len(aligned_seqs)
    if n > 1:
        vals = []
        for i, j in combinations(range(n), 2):
            vals.append(pid_mat[i, j])
        avg_pid = float(np.mean(vals))
    else:
        avg_pid = 1.0

    # Fraction of fully conserved columns
    frac_conserved = fraction_fully_conserved_columns(aligned_seqs)

    # Print results
    np.set_printoptions(precision=3, suppress=True)
    print("Pairwise sequence identity matrix (rows/cols in input order):")
    print(pid_mat)
    print()
    print(f"Average pairwise identity: {avg_pid:.4f}")
    print(f"Fraction of fully conserved alignment columns: {frac_conserved:.4f}")


if __name__ == "__main__":
    main()

