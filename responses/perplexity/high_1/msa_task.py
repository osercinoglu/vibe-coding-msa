#!/usr/bin/env python3
"""
msa_task.py

Align a set of homologous protein sequences from a FASTA file and compute
alignment quality metrics:

- Pairwise identity matrix.
- Average pairwise identity.
- Fraction of fully conserved columns.

Outputs:
- aligned.fasta: multiple sequence alignment in FASTA format.
- metrics.csv: CSV report with pairwise identities and summary metrics.

Usage:
    python msa_task.py --input input.fasta --outdir results/
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


@dataclass
class SequenceRecord:
    """Container for a single sequence record."""

    id: str
    description: str
    sequence: str


def parse_fasta(path: str) -> List[SequenceRecord]:
    """
    Parse a FASTA file into a list of SequenceRecord objects.

    Parameters
    ----------
    path : str
        Path to a FASTA file.

    Returns
    -------
    List[SequenceRecord]
        Parsed sequence records.

    Raises
    ------
    ValueError
        If the file contains no valid sequences.
    """
    records: List[SequenceRecord] = []
    current_id: str | None = None
    current_desc: str | None = None
    current_seq_parts: List[str] = []

    try:
        with open(path, "r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    # Save previous record
                    if current_id is not None:
                        sequence = "".join(current_seq_parts).upper()
                        if sequence:
                            records.append(
                                SequenceRecord(
                                    id=current_id,
                                    description=current_desc or current_id,
                                    sequence=sequence,
                                )
                            )
                    # Start new record
                    header = line[1:].strip()
                    if not header:
                        continue
                    parts = header.split(None, 1)
                    current_id = parts[0]
                    current_desc = header
                    current_seq_parts = []
                else:
                    current_seq_parts.append(line.strip())
        # Save last record
        if current_id is not None:
            sequence = "".join(current_seq_parts).upper()
            if sequence:
                records.append(
                    SequenceRecord(
                        id=current_id,
                        description=current_desc or current_id,
                        sequence=sequence,
                    )
                )
    except OSError as exc:
        raise ValueError(f"Failed to read FASTA file '{path}': {exc}") from exc

    if not records:
        raise ValueError(f"No valid sequences found in FASTA file '{path}'.")

    return records


def validate_sequences(records: List[SequenceRecord]) -> None:
    """
    Validate the sequence list for downstream analysis.

    Parameters
    ----------
    records : List[SequenceRecord]
        List of sequence records.

    Raises
    ------
    ValueError
        If there are fewer than 3 sequences or sequences are empty.
    """
    if len(records) < 3:
        raise ValueError(
            f"Expected at least 3 sequences, found {len(records)} in the input FASTA."
        )

    empty_ids = [rec.id for rec in records if not rec.sequence]
    if empty_ids:
        raise ValueError(
            "The following records have empty sequences: "
            + ", ".join(empty_ids)
        )


def run_muscle_msa(
    input_fasta: str, output_fasta: str, muscle_exe: str = "muscle"
) -> None:
    """
    Run MUSCLE to generate a multiple sequence alignment.

    Parameters
    ----------
    input_fasta : str
        Path to the input unaligned FASTA file.
    output_fasta : str
        Path where the aligned FASTA will be written.
    muscle_exe : str, optional
        MUSCLE executable name or full path, by default "muscle".

    Raises
    ------
    RuntimeError
        If MUSCLE fails or is not available.
    """
    cmd = [
        muscle_exe,
        "-align",
        input_fasta,
        "-output",
        output_fasta,
    ]

    try:
        completed = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            "Failed to run MUSCLE: executable not found. "
            "Ensure MUSCLE is installed and in your PATH."
        ) from exc
    except OSError as exc:
        raise RuntimeError(f"Failed to run MUSCLE: {exc}") from exc

    if completed.returncode != 0:
        raise RuntimeError(
            "MUSCLE alignment failed with exit code "
            f"{completed.returncode}.\nSTDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )

    if not os.path.isfile(output_fasta):
        raise RuntimeError(
            f"MUSCLE completed without errors but output file '{output_fasta}' "
            "was not created."
        )


def write_fasta(records: List[SequenceRecord], path: str) -> None:
    """
    Write a list of SequenceRecord objects to a FASTA file.

    Parameters
    ----------
    records : List[SequenceRecord]
        Sequence records to write.
    path : str
        Output FASTA path.
    """
    try:
        with open(path, "w", encoding="utf-8") as handle:
            for rec in records:
                handle.write(f">{rec.description}\n")
                seq = rec.sequence
                for i in range(0, len(seq), 60):
                    handle.write(seq[i : i + 60] + "\n")
    except OSError as exc:
        raise RuntimeError(f"Failed to write FASTA file '{path}': {exc}") from exc


def load_aligned_fasta(path: str) -> List[SequenceRecord]:
    """
    Load an aligned FASTA file where all sequences have equal length.

    Parameters
    ----------
    path : str
        Path to aligned FASTA.

    Returns
    -------
    List[SequenceRecord]
        Loaded sequence records.

    Raises
    ------
    ValueError
        If sequences have inconsistent alignment lengths.
    """
    records = parse_fasta(path)
    lengths = {len(rec.sequence) for rec in records}
    if len(lengths) != 1:
        raise ValueError(
            "Aligned FASTA contains sequences with inconsistent lengths: "
            + ", ".join(str(l) for l in sorted(lengths))
        )
    return records


def compute_pairwise_identity_matrix(
    records: List[SequenceRecord],
) -> Tuple[pd.DataFrame, float]:
    """
    Compute the pairwise identity matrix for an MSA.

    For each sequence pair (i, j), pairwise identity is defined as:
        matches / aligned_positions
    where:
        matches = number of positions where residues are identical.
        aligned_positions = number of positions where neither residue is a gap
        (i.e. ignoring positions with '-' in either sequence).

    Parameters
    ----------
    records : List[SequenceRecord]
        Aligned sequences, all of equal length.

    Returns
    -------
    Tuple[pd.DataFrame, float]
        - A square DataFrame (n x n) with pairwise identities (0-1).
        - Average pairwise identity over all unique sequence pairs (i < j).
    """
    n = len(records)
    ids = [rec.id for rec in records]
    seqs = np.array([list(rec.sequence) for rec in records])
    length = seqs.shape[1]

    matrix = np.zeros((n, n), dtype=float)

    # Precompute for efficiency
    for i in range(n):
        for j in range(i, n):
            if i == j:
                matrix[i, j] = 1.0
                continue
            s1 = seqs[i]
            s2 = seqs[j]
            non_gap_mask = (s1 != "-") & (s2 != "-")
            aligned_positions = np.count_nonzero(non_gap_mask)
            if aligned_positions == 0:
                identity = 0.0
            else:
                matches = np.count_nonzero((s1 == s2) & non_gap_mask)
                identity = matches / aligned_positions
            matrix[i, j] = identity
            matrix[j, i] = identity

    df = pd.DataFrame(matrix, index=ids, columns=ids)

    # Compute average over upper triangle (excluding diagonal)
    if n > 1:
        triu_indices = np.triu_indices(n, k=1)
        avg_identity = float(matrix[triu_indices].mean())
    else:
        avg_identity = 1.0

    return df, avg_identity


def compute_fraction_fully_conserved_columns(
    records: List[SequenceRecord],
) -> float:
    """
    Compute the fraction of fully conserved columns in an MSA.

    A column is considered fully conserved if:
    - At least one sequence has a non-gap character at that position.
    - All non-gap characters in that column are identical.

    Parameters
    ----------
    records : List[SequenceRecord]
        Aligned sequences, all of equal length.

    Returns
    -------
    float
        Fraction of fully conserved columns (0-1).
    """
    if not records:
        return 0.0

    seqs = np.array([list(rec.sequence) for rec in records])
    length = seqs.shape[1]

    conserved_count = 0
    for col_idx in range(length):
        col = seqs[:, col_idx]
        non_gap = col[col != "-"]
        if non_gap.size == 0:
            # Column with only gaps is not counted as conserved
            continue
        if np.all(non_gap == non_gap[0]):
            conserved_count += 1

    # Fraction over all columns (including non-conserved and gap-only columns)
    return conserved_count / float(length)


def ensure_outdir(path: str) -> None:
    """
    Create the output directory if it does not exist.

    Parameters
    ----------
    path : str
        Directory path.

    Raises
    ------
    RuntimeError
        If the directory cannot be created.
    """
    try:
        os.makedirs(path, exist_ok=True)
    except OSError as exc:
        raise RuntimeError(f"Failed to create output directory '{path}': {exc}") from exc


def write_metrics_report(
    outdir: str,
    pairwise_df: pd.DataFrame,
    avg_identity: float,
    conserved_fraction: float,
    filename: str = "metrics.csv",
) -> str:
    """
    Write a CSV report with pairwise identity matrix and summary metrics.

    The CSV will contain:
    - A block with the pairwise identity matrix.
    - A blank line.
    - A block with summary metrics (average pairwise identity, fraction conserved).

    Parameters
    ----------
    outdir : str
        Output directory.
    pairwise_df : pd.DataFrame
        Pairwise identity matrix.
    avg_identity : float
        Average pairwise identity.
    conserved_fraction : float
        Fraction of fully conserved columns.
    filename : str, optional
        Name of the CSV file, by default "metrics.csv".

    Returns
    -------
    str
        Path to the written CSV file.
    """
    out_path = os.path.join(outdir, filename)

    try:
        with open(out_path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)

            # Pairwise identity matrix section
            writer.writerow(["Pairwise_identity_matrix"])
            writer.writerow([""] + list(pairwise_df.columns))
            for idx, row in pairwise_df.iterrows():
                writer.writerow([idx] + [f"{val:.6f}" for val in row.values])

            # Blank line
            writer.writerow([])

            # Summary metrics section
            writer.writerow(["Summary_metrics"])
            writer.writerow(["Metric", "Value"])
            writer.writerow(["Average_pairwise_identity", f"{avg_identity:.6f}"])
            writer.writerow(
                ["Fraction_fully_conserved_columns", f"{conserved_fraction:.6f}"]
            )
    except OSError as exc:
        raise RuntimeError(f"Failed to write metrics report '{out_path}': {exc}") from exc

    return out_path


def parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments.

    Parameters
    ----------
    argv : List[str] | None
        Optional list of arguments (for testing). If None, uses sys.argv.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Align a set of homologous protein sequences from a FASTA file "
            "and compute alignment quality metrics."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to the input unaligned protein FASTA file.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for aligned FASTA and metrics report.",
    )

    return parser.parse_args(argv)


def main(argv: List[str] | None = None) -> None:
    """
    Main entry point: orchestrate MSA and metrics computation.

    Parameters
    ----------
    argv : List[str] | None
        Optional list of arguments (for testing). If None, uses sys.argv.
    """
    args = parse_args(argv)

    input_path = args.input
    outdir = args.outdir

    # Validate input path
    if not os.path.isfile(input_path):
        print(f"ERROR: Input FASTA file '{input_path}' does not exist or is not a file.", file=sys.stderr)
        sys.exit(1)

    # Ensure output directory exists
    try:
        ensure_outdir(outdir)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    # Read and validate input sequences
    try:
        records = parse_fasta(input_path)
        validate_sequences(records)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    # Write a temporary copy of the input FASTA (ensures clean formatting)
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_input = os.path.join(tmpdir, "input.fasta")
        try:
            write_fasta(records, tmp_input)
        except RuntimeError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            sys.exit(1)

        aligned_fasta_path = os.path.join(outdir, "aligned.fasta")

        # Run MUSCLE to generate MSA
        try:
            run_muscle_msa(tmp_input, aligned_fasta_path)
        except RuntimeError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            sys.exit(1)

    # Load aligned sequences
    try:
        aligned_records = load_aligned_fasta(aligned_fasta_path)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    # Compute metrics
    pairwise_df, avg_identity = compute_pairwise_identity_matrix(aligned_records)
    conserved_fraction = compute_fraction_fully_conserved_columns(aligned_records)

    # Write metrics report
    try:
        metrics_path = write_metrics_report(
            outdir,
            pairwise_df,
            avg_identity,
            conserved_fraction,
        )
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"Alignment written to: {aligned_fasta_path}")
    print(f"Metrics report written to: {metrics_path}")


if __name__ == "__main__":
    main()

