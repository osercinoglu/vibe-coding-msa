import subprocess
import os
import argparse
import sys
import re
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
                if not line:
                    continue
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

def get_muscle_version() -> int:
    """
    Detects MUSCLE major version to determine correct CLI flags.
    """
    try:
        result = subprocess.run(["muscle", "-version"], capture_output=True, text=True)
        version_str = result.stdout or result.stderr
        match = re.search(r"v(\d+)", version_str)
        if match:
            return int(match.group(1))
        return 3  # Default fallback to v3
    except FileNotFoundError:
        print("Error: 'muscle' executable not found in PATH.")
        sys.exit(1)

def run_muscle_alignment(input_path: str, output_path: str):
    """
    Runs MUSCLE alignment using version-appropriate command-line flags.
    """
    version = get_muscle_version()
    
    if version >= 5:
        # MUSCLE v5 syntax
        cmd = ["muscle", "-align", input_path, "-output", output_path]
    else:
        # MUSCLE v3 syntax
        cmd = ["muscle", "-in", input_path, "-out", output_path]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error during alignment: {e.stderr}")
        sys.exit(1)

def calculate_pairwise_identity(seq1: str, seq2: str) -> float:
    """
    Calculates percent identity excluding double-gap positions.
    """
    if len(seq1) != len(seq2):
        raise ValueError("Aligned sequences must be of equal length.")
    
    matches = 0
    comparable_positions = 0
    
    for a, b in zip(seq1, seq2):
        if a == '-' and b == '-':
            continue 
        comparable_positions += 1
        if a == b:
            matches += 1
            
    return (matches / comparable_positions) * 100 if comparable_positions > 0 else 0.0

def compute_metrics(aligned_seqs: List[Tuple[str, str]]) -> Tuple[pd.DataFrame, float, float]:
    """
    Computes identity matrix, average identity, and conservation fraction.
    """
    headers = [s[0] for s in aligned_seqs]
    seqs = [s[1] for s in aligned_seqs]
    num_seqs = len(seqs)
    aln_len = len(seqs[0])
    
    # 1. Identity Matrix
    matrix = np.zeros((num_seqs, num_seqs))
    for i in range(num_seqs):
        for j in range(num_seqs):
            if i == j:
                matrix[i, j] = 100.0
            else:
                matrix[i, j] = calculate_pairwise_identity(seqs[i], seqs[j])
    
    df_matrix = pd.DataFrame(matrix, index=headers, columns=headers)
    
    # 2. Average Pairwise Identity
    if num_seqs > 1:
        mask = ~np.eye(num_seqs, dtype=bool)
        avg_identity = matrix[mask].mean()
    else:
        avg_identity = 100.0
    
    # 3. Conservation Fraction (Strictly fully conserved residues, no gaps)
    conserved_cols = 0
    for col in range(aln_len):
        residues = {s[col] for s in seqs}
        if len(residues) == 1 and '-' not in residues:
            conserved_cols += 1
            
    conservation_fraction = conserved_cols / aln_len if aln_len > 0 else 0.0
    
    return df_matrix, avg_identity, conservation_fraction

def main():
    parser = argparse.ArgumentParser(description="Protein Alignment and Metrics Utility")
    parser.add_argument("--input", required=True, help="Input FASTA file path")
    parser.add_argument("--outdir", default=".", help="Directory to save output files")
    args = parser.parse_args()

    # Pre-run Validation
    if not os.path.exists(args.input):
        print(f"Error: File '{args.input}' not found.")
        sys.exit(1)

    raw_data = parse_fasta(args.input)
    if len(raw_data) < 3:
        print("Error: At least 3 sequences are required for this analysis.")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    aligned_fasta = os.path.join(args.outdir, "aligned.fasta")
    matrix_csv = os.path.join(args.outdir, "identity_matrix.csv")
    summary_txt = os.path.join(args.outdir, "summary.txt")

    print(f"Running alignment on {len(raw_data)} sequences...")
    run_muscle_alignment(args.input, aligned_fasta)

    print("Computing metrics...")
    aligned_data = parse_fasta(aligned_fasta)
    id_matrix, avg_id, cons_frac = compute_metrics(aligned_data)

    # Export
    id_matrix.to_csv(matrix_csv)
    with open(summary_txt, "w") as f:
        f.write(f"Average Pairwise Identity: {avg_id:.2f}%\n")
        f.write(f"Fraction of fully conserved columns: {cons_frac:.4f}\n")

    print(f"\nSuccess. Results saved to: {args.outdir}")
    print(f"Avg Identity: {avg_id:.2f}% | Conservation: {cons_frac:.4f}")

if __name__ == "__main__":
    main()
