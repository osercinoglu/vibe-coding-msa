import argparse
import subprocess
import os
import sys
from typing import List, Dict, Tuple
import numpy as np
import pandas as pd

def validate_input(input_path: str) -> List[str]:
    """
    Validates that the input file exists and contains at least 3 sequences.
    Returns a list of sequence headers.
    """
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    with open(input_path, 'r') as f:
        headers = [line for line in f if line.startswith('>')]
    
    if len(headers) < 3:
        raise ValueError("The input FASTA must contain at least 3 sequences for MSA.")
    
    return headers

def run_msa_alignment(input_fasta: str, output_fasta: str) -> None:
    """
    Executes Clustal Omega via subprocess to perform multiple sequence alignment.
    Assumes clustalo is installed and in the system PATH.
    """
    try:
        # Using Clustal Omega as the default engine
        cmd = ["clustalo", "-i", input_fasta, "-o", output_fasta, "--force", "--outfmt=fa"]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError:
        print("Error: 'clustalo' executable not found. Please ensure Clustal Omega is installed.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Error during alignment: {e.stderr}")
        sys.exit(1)

def parse_fasta(fasta_path: str) -> Dict[str, str]:
    """
    Parses a FASTA file into a dictionary of {header: sequence}.
    """
    sequences = {}
    current_header = None
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_header = line[1:]
                sequences[current_header] = ""
            else:
                sequences[current_header] += line
    return sequences

def calculate_pairwise_identity(seq1: str, seq2: str) -> float:
    """
    Calculates the identity percentage between two aligned sequences, 
    excluding positions where both have gaps.
    """
    matches = 0
    valid_length = 0
    for a, b in zip(seq1, seq2):
        if a == '-' and b == '-':
            continue
        valid_length += 1
        if a == b:
            matches += 1
    
    return (matches / valid_length) * 100 if valid_length > 0 else 0.0

def compute_metrics(aligned_seqs: Dict[str, str]) -> Tuple[pd.DataFrame, float, float]:
    """
    Computes identity matrix, average identity, and conservation fraction.
    """
    headers = list(aligned_seqs.keys())
    n = len(headers)
    matrix = np.zeros((n, n))
    
    # 1. Pairwise Identity Matrix
    for i in range(n):
        for j in range(n):
            matrix[i, j] = calculate_pairwise_identity(
                aligned_seqs[headers[i]], 
                aligned_seqs[headers[j]]
            )
    
    df_matrix = pd.DataFrame(matrix, index=headers, columns=headers)
    
    # 2. Average Pairwise Identity (excluding self-comparison)
    upper_triangle = matrix[np.triu_indices(n, k=1)]
    avg_identity = np.mean(upper_triangle)
    
    # 3. Fraction of fully conserved columns
    seq_list = list(aligned_seqs.values())
    alignment_length = len(seq_list[0])
    conserved_cols = 0
    
    for col_idx in range(alignment_length):
        column = [seq[col_idx] for seq in seq_list]
        # Check if all characters in column are identical and not gaps
        if len(set(column)) == 1 and column[0] != '-':
            conserved_cols += 1
            
    conservation_frac = conserved_cols / alignment_length
    
    return df_matrix, avg_identity, conservation_frac

def main():
    parser = argparse.ArgumentParser(description="MSA Pipeline and Quality Metrics Tool")
    parser.add_argument("--input", required=True, help="Path to input protein FASTA file")
    parser.add_argument("--outdir", default=".", help="Directory to save output files")
    
    args = parser.parse_args()
    
    try:
        # Setup output paths
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
            
        aligned_fasta = os.path.join(args.outdir, "aligned.fasta")
        report_path = os.path.join(args.outdir, "alignment_report.csv")
        metrics_summary = os.path.join(args.outdir, "summary.txt")

        # 1. Validation
        print(f"Validating {args.input}...")
        validate_input(args.input)

        # 2. Alignment
        print("Running Multiple Sequence Alignment...")
        run_msa_alignment(args.input, aligned_fasta)
        print(f"Alignment saved to {aligned_fasta}")

        # 3. Metric Computation
        aligned_seqs = parse_fasta(aligned_fasta)
        matrix, avg_id, cons_frac = compute_metrics(aligned_seqs)

        # 4. Save Results
        matrix.to_csv(report_path)
        
        with open(metrics_summary, 'w') as f:
            f.write("MSA Quality Metrics Report\n")
            f.write("==========================\n")
            f.write(f"Average Pairwise Identity: {avg_id:.2f}%\n")
            f.write(f"Fraction of fully conserved columns: {cons_frac:.4f}\n")

        print(f"Metrics report saved to {report_path} and {metrics_summary}")
        print("\nSummary:")
        print(f" - Avg Identity: {avg_id:.2f}%")
        print(f" - Conservation: {cons_frac:.2%}")

    except Exception as e:
        print(f"Critical Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
