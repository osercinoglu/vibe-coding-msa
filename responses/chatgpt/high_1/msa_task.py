#!/usr/bin/env python3
"""
msa_task.py

Align a set of homologous protein sequences using an external MSA tool
(MUSCLE or Clustal Omega) and compute basic alignment quality metrics.

Outputs:
- aligned.fasta : multiple sequence alignment
- pairwise_identity.csv : pairwise identity matrix
- summary.txt : summary statistics

Requirements:
- Python 3.11
- numpy, pandas
- External MSA tool available on PATH (muscle or clustalo)
"""

from __future__ import annotations

import argparse
import os
import sys
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


def parse_fasta(path: Path) -> Dict[str, str]:
    """
    Parse a FASTA file into an ordered dictionary of {header: sequence}.

    Parameters
    ----------
    path : Path
        Path to FASTA file.

    Returns
    -------
    Dict[str, str]
        Mapping from sequence ID to sequence string.

    Raises
    ------
    ValueError
        If the file is malformed or empty.
    """
    sequences: Dict[str, str] = {}
    current_id = None
    current_seq: List[str] = []

    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                if current_id is None:
                    raise ValueError("FASTA format error: sequence before header")
                current_seq.append(line)

        if current_id is not None:
            sequences[current_id] = "".join(current_seq)

    if not sequences:
        raise ValueError("No sequences found in FASTA file")

    return sequences


def validate_input_fasta(path: Path) -> Dict[str, str]:
    """
    Validate that the input FASTA exists and contains >= 3 sequences.

    Parameters
    ----------
    path : Path

    Returns
    -------
    Dict[str, str]
        Parsed sequences.

    Raises
    ------
    SystemExit
        On validation failure.
    """
    if not path.exists():
        sys.exit(f"ERROR: Input FASTA does not exist: {path}")

    try:
        sequences = parse_fasta(path)
    except ValueError as e:
        sys.exit(f"ERROR parsing FASTA: {e}")

    if len(sequences) < 3:
        sys.exit("ERROR: At least 3 sequences are required for MSA")

    return sequences


def run_msa(
    input_fasta: Path,
    output_fasta: Path,
    tool: str = "muscle",
) -> None:
    """
    Run an external MSA tool to generate a multiple sequence alignment.

    Parameters
    ----------
    input_fasta : Path
        Input FASTA file.
    output_fasta : Path
        Output aligned FASTA file.
    tool : str
        MSA tool to use ("muscle" or "clustalo").

    Raises
    ------
    SystemExit
        If the external tool fails.
    """
    if tool == "muscle":
        cmd = [
            "muscle",
            "-align",
            str(input_fasta),
            "-output",
            str(output_fasta),
        ]
    elif tool == "clustalo":
        cmd = [
            "clustalo",
            "-i",
            str(input_fasta),
            "-o",
            str(output_fasta),
            "--force",
        ]
    else:
        sys.exit(f"Unsupported MSA tool: {tool}")

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError:
        sys.exit(
            f"ERROR: MSA tool '{tool}' not found on PATH. "
            "Please install it or adjust the code."
        )
    except subprocess.CalledProcessError as e:
        sys.exit(
            f"ERROR: MSA tool failed with exit code {e.returncode}\n"
            f"STDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}"
        )


def compute_pairwise_identity(
    aligned_seqs: Dict[str, str],
) -> pd.DataFrame:
    """
    Compute pairwise sequence identity from an aligned FASTA.

    Identity is defined as:
        (# identical residues excluding gap-gap columns) /
        (# positions where neither sequence has a gap)

    Parameters
    ----------
    aligned_seqs : Dict[str, str]
        Aligned sequences.

    Returns
    -------
    pd.DataFrame
        Square matrix of pairwise identities (0-1).
    """
    ids = list(aligned_seqs.keys())
    seqs = [aligned_seqs[i] for i in ids]
    n = len(seqs)
    length = len(seqs[0])

    for s in seqs:
        if len(s) != length:
            raise ValueError("Aligned sequences are not all the same length")

    mat = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(n):
            matches = 0
            compared = 0
            for a, b in zip(seqs[i], seqs[j]):
                if a == "-" or b == "-":
                    continue
                compared += 1
                if a == b:
                    matches += 1
            mat[i, j] = matches / compared if compared > 0 else np.nan

    return pd.DataFrame(mat, index=ids, columns=ids)


def fraction_fully_conserved_columns(
    aligned_seqs: Dict[str, str],
) -> float:
    """
    Compute the fraction of alignment columns that are fully conserved.

    A column is considered fully conserved if:
    - No gaps
    - All residues are identical

    Parameters
    ----------
    aligned_seqs : Dict[str, str]

    Returns
    -------
    float
        Fraction of fully conserved columns.
    """
    seqs = list(aligned_seqs.values())
    length = len(seqs[0])
    n_seqs = len(seqs)

    conserved = 0

    for pos in range(length):
        column = [seq[pos] for seq in seqs]
        if "-" in column:
            continue
        if len(set(column)) == 1:
            conserved += 1

    return conserved / length if length > 0 else 0.0


def write_summary(
    outdir: Path,
    pairwise_df: pd.DataFrame,
    frac_conserved: float,
) -> None:
    """
    Write summary statistics to a text file.

    Parameters
    ----------
    outdir : Path
    pairwise_df : pd.DataFrame
    frac_conserved : float
    """
    avg_identity = np.nanmean(pairwise_df.values)

    summary_path = outdir / "summary.txt"
    with summary_path.open("w") as fh:
        fh.write("MSA Summary Statistics\n")
        fh.write("======================\n\n")
        fh.write(f"Number of sequences: {pairwise_df.shape[0]}\n")
        fh.write(f"Average pairwise identity: {avg_identity:.4f}\n")
        fh.write(f"Fraction fully conserved columns: {frac_conserved:.4f}\n")


def main() -> None:
    """
    Main entry point.
    """
    parser = argparse.ArgumentParser(
        description="Run MSA and compute alignment quality metrics"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input unaligned protein FASTA file",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory",
    )
    parser.add_argument(
        "--tool",
        default="muscle",
        choices=["muscle", "clustalo"],
        help="MSA tool to use (default: muscle)",
    )

    args = parser.parse_args()

    input_fasta = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Validate input
    validate_input_fasta(input_fasta)

    aligned_fasta = outdir / "aligned.fasta"

    # Run MSA
    run_msa(input_fasta, aligned_fasta, tool=args.tool)

    # Load alignment
    try:
        aligned_seqs = parse_fasta(aligned_fasta)
    except ValueError as e:
        sys.exit(f"ERROR parsing aligned FASTA: {e}")

    # Compute metrics
    pairwise_df = compute_pairwise_identity(aligned_seqs)
    frac_conserved = fraction_fully_conserved_columns(aligned_seqs)

    # Write outputs
    pairwise_csv = outdir / "pairwise_identity.csv"
    pairwise_df.to_csv(pairwise_csv)

    write_summary(outdir, pairwise_df, frac_conserved)


if __name__ == "__main__":
    main()

