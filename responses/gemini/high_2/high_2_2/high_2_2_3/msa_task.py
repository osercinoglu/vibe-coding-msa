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
    Determines the major version of MUSCLE installed.
    """
    try:
        result = subprocess.run(["muscle", "-version"], capture_output=True, text=True)
        version_str = result.stdout or result.stderr
        match = re.search(r"v(\d+)", version_str)
        if match:
            return int(match.group(1))
        return 3  # Fallback to v3
    except FileNotFoundError:
        print("Error: 'muscle' executable not found in PATH.")
        sys.exit(1)

def run_muscle_alignment(input_path: str, output_path: str):
    """
    Runs MUSCLE alignment using version-appropriate flags.
    """
    version = get_muscle_version()
    
    if version >= 5:
        # MUSCLE v5+ syntax
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
        raise ValueError("Aligned sequences must have identical lengths.")
    
    matches = 0
    comparable_length = 0
    
    for a, b in zip(seq1, seq2):
        if a == '-' and b == '-':
            continue 
        comparable_length += 1
        if a == b:
            matches += 1
            
    return (matches / comparable_length) * 100 if comparable_length > 0 else 0.0

def compute_metrics(aligned_seqs: List[Tuple[str, str]]) -> Tuple[pd.DataFrame, float, float]:
    """
    Computes the identity matrix and conservation metrics.
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
    
    # 2. Average Identity (off-diagonal)
    if num_seqs > 1:
        mask = ~np.eye(num_seqs, dtype=bool)
        avg_identity = matrix[mask].mean()
    else:
        avg_identity = 100.0
    
    # 3. Conservation Fraction
    conserved_count = 0
    for col in range(aln_len):
        residues = set(s[col] for s in seqs)
        if len(residues) == 1 and '-' not in residues:
            conserved_count += 1
    conservation_fraction = conserved_count / aln_len if aln_len > 0 else 0.0
    
    return df_matrix, avg_identity, conservation_fraction

def main():
    parser = argparse.ArgumentParser(description="Protein MSA and Quality Metrics Tool")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--outdir", default=".", help="Output directory")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: {args.input} not found.")
        sys.exit(1)

    # Initial check on input sequences 
    raw_sequences = parse_fasta(args.input)
    if len(raw_sequences) < 3:
        print("Error: Minimum 3 sequences required.")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    aligned_fasta = os.path.join(args.outdir, "aligned.fasta")
    report_path = os.path.join(args.outdir, "alignment_report.csv")

    print(f"Starting alignment for {len(raw_sequences)} sequences...")
    run_muscle_alignment(args.input, aligned_fasta)

    print("Analyzing alignment quality...")
    aligned_data = parse_fasta(aligned_fasta)
    identity_matrix, avg_id, cons_frac = compute_metrics(aligned_data)

    # Save outputs
    identity_matrix.to_csv(report_path)
    
    summary_path = os.path.join(args.outdir, "summary.txt")
    with open(summary_path, "w") as f:
        f.write(f"Average Pairwise Identity: {avg_id:.2f}%\n")
        f.write(f"Fraction of fully conserved columns: {cons_frac:.4f}\n")

    print(f"Done. Reports available in {args.outdir}")

if __name__ == "__main__":
    main()
