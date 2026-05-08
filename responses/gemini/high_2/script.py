import subprocess
import os
import argparse
import sys
from io import StringIO
from typing import List, Tuple

import numpy as np
import pandas as pd

def parse_fasta(file_path: str) -> List[Tuple[str, str]]:
    """
    Parses a FASTA file into a list of (header, sequence) tuples.
    """
    sequences = []
    try:
        with open(file_path, 'r') as f:
            header, seq = None, []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if header:
                        sequences.append((header, "".join(seq)))
                    header = line[1:]
                    seq = []
                else:
                    seq.append(line)
            if header:
                sequences.append((header, "".join(seq)))
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        sys.exit(1)
    
    return sequences

def run_muscle_alignment(input_path: str, output_path: str):
    """
    Runs MUSCLE alignment via subprocess. 
    Assumes MUSCLE is installed and available in the PATH.
    """
    try:
        # Standard command for MUSCLE v3.8 or v5
        subprocess.run(
            ["muscle", "-in", input_path, "-out", output_path],
            check=True,
            capture_output=True,
            text=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error during alignment: {e.stderr}")
        sys.exit(1)
    except FileNotFoundError:
        print("Error: 'muscle' executable not found in PATH.")
        sys.exit(1)

def calculate_pairwise_identity(seq1: str, seq2: str) -> float:
    """
    Calculates the percent identity between two aligned sequences.
    Gaps are treated as non-matches.
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length for pairwise identity.")
    
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-' and b != '-')
    length = len(seq1)
    return (matches / length) * 100 if length > 0 else 0.0

def compute_metrics(aligned_seqs: List[Tuple[str, str]]) -> Tuple[pd.DataFrame, float, float]:
    """
    Computes identity matrix, average identity, and conservation fraction.
    """
    headers = [s[0] for s in aligned_seqs]
    seqs = [s[1] for s in aligned_seqs]
    num_seqs = len(seqs)
    aln_len = len(seqs[0])
    
    # 1. Pairwise Identity Matrix
    matrix = np.zeros((num_seqs, num_seqs))
    for i in range(num_seqs):
        for j in range(num_seqs):
            matrix[i, j] = calculate_pairwise_identity(seqs[i], seqs[j])
    
    df_matrix = pd.DataFrame(matrix, index=headers, columns=headers)
    
    # 2. Average Pairwise Identity (excluding diagonal)
    mask = np.ones(matrix.shape, dtype=bool)
    np.fill_diagonal(mask, 0)
    avg_identity = matrix[mask].mean()
    
    # 3. Fraction of fully conserved columns
    conserved_count = 0
    for col in range(aln_len):
        residues = set(s[col] for s in seqs)
        if len(residues) == 1 and '-' not in residues:
            conserved_count += 1
    conservation_fraction = conserved_count / aln_len if aln_len > 0 else 0.0
    
    return df_matrix, avg_identity, conservation_fraction

def main():
    parser = argparse.ArgumentParser(description="Align protein sequences and compute quality metrics.")
    parser.add_argument("--input", required=True, help="Input FASTA file path")
    parser.add_argument("--outdir", default=".", help="Directory to save output files")
    args = parser.parse_args()

    # Validation
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' does not exist.")
        sys.exit(1)

    raw_sequences = parse_fasta(args.input)
    if len(raw_sequences) < 3:
        print("Error: Input must contain at least 3 sequences for meaningful MSA.")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    aligned_fasta = os.path.join(args.outdir, "aligned.fasta")
    report_path = os.path.join(args.outdir, "alignment_report.csv")

    print(f"Aligning {len(raw_sequences)} sequences...")
    run_muscle_alignment(args.input, aligned_fasta)

    print("Computing metrics...")
    aligned_data = parse_fasta(aligned_fasta)
    identity_matrix, avg_id, cons_frac = compute_metrics(aligned_data)

    # Save outputs
    identity_matrix.to_csv(report_path)
    
    with open(os.path.join(args.outdir, "summary.txt"), "w") as f:
        f.write(f"Average Pairwise Identity: {avg_id:.2f}%\n")
        f.write(f"Fraction of fully conserved columns: {cons_frac:.4f}\n")

    print(f"Results saved to {args.outdir}")
    print(f"Average Identity: {avg_id:.2f}%")
    print(f"Conservation: {cons_frac:.4f}")

if __name__ == "__main__":
    main()
