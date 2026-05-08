#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path
import numpy as np
import pandas as pd


def read_fasta(path):
    sequences = []
    names = []
    seq = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq:
                    sequences.append("".join(seq))
                    seq = []
                names.append(line[1:].split()[0])
            else:
                seq.append(line)
        if seq:
            sequences.append("".join(seq))

    if len(sequences) < 2:
        raise ValueError("FASTA must contain at least two sequences")

    return names, sequences


def run_muscle(input_fasta, output_fasta):
    cmd = [
        "muscle",
        "-align", str(input_fasta),
        "-output", str(output_fasta)
    ]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    except FileNotFoundError:
        sys.exit("ERROR: MUSCLE not found in PATH.")
    except subprocess.CalledProcessError:
        sys.exit("ERROR: MUSCLE failed.")


def read_alignment(path):
    names, seqs = read_fasta(path)
    length = {len(s) for s in seqs}
    if len(length) != 1:
        raise ValueError("Alignment sequences are not the same length")
    return names, np.array([list(s) for s in seqs])


def pairwise_identity(aln):
    n = aln.shape[0]
    pid = np.zeros((n, n))

    for i in range(n):
        for j in range(i, n):
            matches = 0
            compared = 0
            for a, b in zip(aln[i], aln[j]):
                if a == "-" or b == "-":
                    continue
                compared += 1
                if a == b:
                    matches += 1
            pid[i, j] = pid[j, i] = matches / compared if compared > 0 else 0.0
    return pid


def fraction_fully_conserved(aln):
    conserved = 0
    total = aln.shape[1]

    for col in aln.T:
        residues = [r for r in col if r != "-"]
        if len(residues) == 0:
            continue
        if len(set(residues)) == 1:
            conserved += 1

    return conserved / total


def main():
    parser = argparse.ArgumentParser(description="MSA + identity statistics")
    parser.add_argument("--input", required=True, help="Input FASTA")
    parser.add_argument("--output", required=True, help="Output aligned FASTA")
    args = parser.parse_args()

    input_fasta = Path(args.input)
    output_fasta = Path(args.output)

    # Sanity check
    read_fasta(input_fasta)

    # Alignment
    run_muscle(input_fasta, output_fasta)

    # Load alignment
    names, aln = read_alignment(output_fasta)

    # Pairwise identity
    pid = pairwise_identity(aln)
    pid_df = pd.DataFrame(pid, index=names, columns=names)

    avg_pid = pid[np.triu_indices_from(pid, k=1)].mean()
    frac_cons = fraction_fully_conserved(aln)

    # Output
    print("\nPairwise identity matrix:")
    print(pid_df.round(3))

    print(f"\nAverage pairwise identity: {avg_pid:.3f}")
    print(f"Fraction of fully conserved columns: {frac_cons:.3f}")


if __name__ == "__main__":
    main()

