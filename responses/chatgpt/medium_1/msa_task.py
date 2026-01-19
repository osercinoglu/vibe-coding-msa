#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path
import numpy as np

def read_fasta(path):
    sequences = []
    header = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequences.append((header, "".join(seq)))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line)
        if header is not None:
            sequences.append((header, "".join(seq)))
    return sequences

def write_fasta(sequences, path):
    with open(path, "w") as f:
        for header, seq in sequences:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

def run_muscle(input_fasta, output_fasta):
    try:
        subprocess.run(
            ["muscle", "-align", str(input_fasta), "-output", str(output_fasta)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except FileNotFoundError:
        sys.exit(
            "ERROR: MUSCLE not found.\n"
            "Install it with:\n"
            "  conda install -c bioconda muscle\n"
        )
    except subprocess.CalledProcessError:
        sys.exit("ERROR: MUSCLE failed to run.")

def pairwise_identity(seq1, seq2):
    matches = 0
    length = 0
    for a, b in zip(seq1, seq2):
        if a == "-" or b == "-":
            continue
        length += 1
        if a == b:
            matches += 1
    return matches / length if length > 0 else 0.0

def compute_identity_matrix(alignment):
    n = len(alignment)
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            mat[i, j] = pairwise_identity(alignment[i][1], alignment[j][1])
    return mat

def fully_conserved_fraction(alignment):
    seqs = [seq for _, seq in alignment]
    n_seq = len(seqs)
    aln_len = len(seqs[0])

    conserved = 0
    valid_cols = 0

    for i in range(aln_len):
        column = [seq[i] for seq in seqs]
        residues = [r for r in column if r != "-"]
        if not residues:
            continue
        valid_cols += 1
        if len(set(residues)) == 1 and len(residues) == n_seq:
            conserved += 1

    return conserved / valid_cols if valid_cols > 0 else 0.0

def main():
    parser = argparse.ArgumentParser(description="MSA and basic identity statistics")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output aligned FASTA file")
    args = parser.parse_args()

    input_fasta = Path(args.input)
    output_fasta = Path(args.output)

    if not input_fasta.exists():
        sys.exit(f"ERROR: Input file not found: {input_fasta}")

    # Run MUSCLE
    run_muscle(input_fasta, output_fasta)

    # Read alignment
    alignment = read_fasta(output_fasta)

    if len(alignment) < 2:
        sys.exit("ERROR: Need at least two sequences for statistics.")

    # Compute statistics
    id_matrix = compute_identity_matrix(alignment)
    avg_identity = np.mean(id_matrix[np.triu_indices(len(alignment), k=1)])
    conserved_frac = fully_conserved_fraction(alignment)

    names = [h for h, _ in alignment]

    print("\nPairwise sequence identity matrix:")
    print(" " * 12 + " ".join(f"{n:>8}" for n in names))
    for i, name in enumerate(names):
        row = " ".join(f"{id_matrix[i, j]:8.3f}" for j in range(len(names)))
        print(f"{name:>12} {row}")

    print("\nAverage pairwise identity: {:.3f}".format(avg_identity))
    print("Fully conserved columns (ignoring gaps): {:.3f}".format(conserved_frac))
    print(f"\nAlignment written to: {output_fasta}")

if __name__ == "__main__":
    main()

