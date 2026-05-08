#!/usr/bin/env python3
"""
msa_task.py

Align a set of homologous protein sequences and compute alignment quality metrics.

Outputs:
- aligned.fasta : multiple sequence alignment
- pairwise_identity.csv : pairwise sequence identity matrix
- summary.txt : average pairwise identity and fraction of fully conserved columns
"""

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd


def read_fasta(path: Path) -> Dict[str, str]:
    """
    Read a FASTA file into a dictionary of {header: sequence}.
    """
    sequences = {}
    header = None
    seq_chunks = []

    try:
        with path.open() as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        sequences[header] = "".join(seq_chunks)
                    header = line[1:].strip()
                    seq_chunks = []
                else:
                    seq_chunks.append(line)
            if header is not None:
                sequences[header] = "".join(seq_chunks)
    except OSError as e:
        raise RuntimeError(f"Failed to read FASTA file: {e}")

    return sequences


def validate_sequences(seqs: Dict[str, str], min_n: int = 3) -> None:
    """
    Validate that we have at least min_n non-empty sequences.
    """
    if len(seqs) < min_n:
        raise ValueError(f"Need at least {min_n} sequences, found {len(seqs)}.")

    for name, seq in seqs.items():
        if not seq:
            raise ValueError(f"Sequence '{name}' is empty.")


def run_muscle(input_fasta: Path, output_fasta: Path) -> None:
    """
    Run MUSCLE (v5+) to produce a multiple sequence alignment.

    MUSCLE 5 syntax differs from MUSCLE 3:
      muscle -align input.fasta -output aligned.fasta
    """
    cmd = [
        "muscle",
        "-align",
        str(input_fasta),
        "-output",
        str(output_fasta),
    ]

    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except FileNotFoundError:
        raise RuntimeError("MUSCLE executable not found in PATH.")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"MUSCLE failed:\n{e.stderr.strip()}")


def read_aligned_fasta(path: Path) -> Dict[str, str]:
    """
    Read an aligned FASTA file. All sequences must have equal length.
    """
    seqs = read_fasta(path)
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        raise ValueError("Aligned sequences do not all have the same length.")
    return seqs


def pairwise_identity(seq_a: str, seq_b: str) -> float:
    """
    Compute pairwise sequence identity between two aligned sequences.
    Gaps ('-') are ignored.
    """
    matches = 0
    comparable = 0

    for a, b in zip(seq_a, seq_b):
        if a == "-" or b == "-":
            continue
        comparable += 1
        if a == b:
            matches += 1

    if comparable == 0:
        return np.nan

    return matches / comparable


def compute_pairwise_identity_matrix(seqs: Dict[str, str]) -> pd.DataFrame:
    """
    Compute a pairwise identity matrix for aligned sequences.
    """
    names = list(seqs.keys())
    n = len(names)
    matrix = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(n):
            matrix[i, j] = pairwise_identity(seqs[names[i]], seqs[names[j]])

    return pd.DataFrame(matrix, index=names, columns=names)


def fraction_fully_conserved_columns(seqs: Dict[str, str]) -> float:
    """
    Fraction of alignment columns that are fully conserved,
    ignoring gaps.
    """
    aligned = list(seqs.values())
    n_cols = len(aligned[0])

    conserved = 0
    considered = 0

    for col in range(n_cols):
        residues = {seq[col] for seq in aligned if seq[col] != "-"}
        if not residues:
            continue
        considered += 1
        if len(residues) == 1:
            conserved += 1

    if considered == 0:
        return np.nan

    return conserved / considered


def write_summary(outdir: Path, avg_identity: float, frac_conserved: float) -> None:
    """
    Write a text summary of alignment statistics.
    """
    summary_path = outdir / "summary.txt"
    with summary_path.open("w") as fh:
        fh.write(f"Average pairwise identity: {avg_identity:.4f}\n")
        fh.write(f"Fraction of fully conserved columns: {frac_conserved:.4f}\n")


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Align protein sequences and compute alignment quality metrics."
    )
    parser.add_argument("--input", required=True, help="Input unaligned FASTA file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    input_fasta = Path(args.input)
    outdir = Path(args.outdir)

    if not input_fasta.exists():
        sys.exit(f"ERROR: Input file does not exist: {input_fasta}")

    try:
        outdir.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        sys.exit(f"ERROR: Could not create output directory: {e}")

    try:
        seqs = read_fasta(input_fasta)
        validate_sequences(seqs, min_n=3)
    except Exception as e:
        sys.exit(f"ERROR while reading input FASTA: {e}")

    aligned_fasta = outdir / "aligned.fasta"

    try:
        run_muscle(input_fasta, aligned_fasta)
    except Exception as e:
        sys.exit(f"ERROR during alignment: {e}")

    try:
        aligned_seqs = read_aligned_fasta(aligned_fasta)
    except Exception as e:
        sys.exit(f"ERROR while reading aligned FASTA: {e}")

    pid_df = compute_pairwise_identity_matrix(aligned_seqs)

    upper = pid_df.values[np.triu_indices_from(pid_df.values, k=1)]
    avg_identity = np.nanmean(upper)

    frac_conserved = fraction_fully_conserved_columns(aligned_seqs)

    pid_df.to_csv(outdir / "pairwise_identity.csv")
    write_summary(outdir, avg_identity, frac_conserved)

    print("MSA and statistics completed successfully.")


if __name__ == "__main__":
    main()

