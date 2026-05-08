#!/usr/bin/env python3
"""
msa_task.py

Align a set of homologous protein sequences from an unaligned FASTA file and
compute alignment quality metrics.

Outputs:
- Multiple sequence alignment in aligned.fasta
- CSV report with:
    * Pairwise identity matrix
    * Average pairwise identity
    * Fraction of fully conserved columns
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import tempfile
import textwrap
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


@dataclass
class SequenceRecord:
    """Container for a single sequence record."""
    id: str
    sequence: str


def parse_fasta(path: Path) -> List[SequenceRecord]:
    """
    Parse a FASTA file into a list of SequenceRecord.

    Parameters
    ----------
    path : Path
        Path to the FASTA file.

    Returns
    -------
    List[SequenceRecord]
        Parsed sequence records.

    Raises
    ------
    ValueError
        If FASTA is malformed or contains no sequences.
    """
    records: List[SequenceRecord] = []
    current_id: str | None = None
    current_seq: List[str] = []

    try:
        with path.open("r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_id is not None:
                        seq = "".join(current_seq).replace(" ", "").upper()
                        if not seq:
                            raise ValueError(
                                f"Empty sequence for record '{current_id}'."
                            )
                        records.append(SequenceRecord(current_id, seq))
                    current_id = line[1:].strip()
                    if not current_id:
                        raise ValueError("Found FASTA header with empty ID.")
                    current_seq = []
                else:
                    if current_id is None:
                        raise ValueError(
                            "FASTA file does not start with a header line ('>')."
                        )
                    current_seq.append(line.strip())
        # Add last record
        if current_id is not None:
            seq = "".join(current_seq).replace(" ", "").upper()
            if not seq:
                raise ValueError(f"Empty sequence for record '{current_id}'.")
            records.append(SequenceRecord(current_id, seq))
    except OSError as exc:
        raise ValueError(f"Failed to read FASTA file '{path}': {exc}") from exc

    if not records:
        raise ValueError(f"No sequences found in FASTA file '{path}'.")
    return records


def validate_sequences(records: List[SequenceRecord]) -> None:
    """
    Validate that the sequence list contains at least 3 sequences.

    Parameters
    ----------
    records : List[SequenceRecord]
        Parsed sequence records.

    Raises
    ------
    ValueError
        If there are fewer than 3 sequences.
    """
    if len(records) < 3:
        raise ValueError(
            f"Expected at least 3 sequences, found {len(records)}."
        )


def run_muscle_msa(
    records: List[SequenceRecord],
    aligned_fasta_path: Path,
) -> None:
    """
    Run MUSCLE to generate a multiple sequence alignment.

    This function uses subprocess to call MUSCLE, assuming it is installed
    and available on the PATH.

    Parameters
    ----------
    records : List[SequenceRecord]
        Input unaligned sequences.
    aligned_fasta_path : Path
        Path where the aligned FASTA will be written.

    Raises
    ------
    RuntimeError
        If MUSCLE fails or returns a non-zero exit code.
    """
    # Write input sequences to a temporary FASTA file for MUSCLE
    with tempfile.TemporaryDirectory() as tmpdir:
        input_tmp = Path(tmpdir) / "input.fasta"
        with input_tmp.open("w", encoding="utf-8") as fh:
            for rec in records:
                fh.write(f">{rec.id}\n")
                # Wrap sequence to 60 chars per line for readability (optional)
                for i in range(0, len(rec.sequence), 60):
                    fh.write(rec.sequence[i : i + 60] + "\n")

        cmd = [
            "muscle",
            "-in",
            str(input_tmp),
            "-out",
            str(aligned_fasta_path),
        ]

        try:
            completed = subprocess.run(
                cmd,
                check=False,
                capture_output=True,
                text=True,
            )
        except FileNotFoundError as exc:
            raise RuntimeError(
                "MUSCLE executable not found. Please ensure MUSCLE is installed "
                "and available on your PATH."
            ) from exc

        if completed.returncode != 0:
            stderr_preview = textwrap.shorten(
                completed.stderr.strip(), width=300, placeholder="..."
            )
            raise RuntimeError(
                f"MUSCLE alignment failed with exit code {completed.returncode}."
                f"\nCommand: {' '.join(cmd)}"
                f"\nError output (truncated):\n{stderr_preview}"
            )


def parse_aligned_fasta(path: Path) -> List[SequenceRecord]:
    """
    Parse an aligned FASTA file (MSA) into a list of SequenceRecord.

    Parameters
    ----------
    path : Path
        Path to the aligned FASTA file.

    Returns
    -------
    List[SequenceRecord]
        Parsed aligned sequence records.

    Raises
    ------
    ValueError
        If the file cannot be parsed or contains no sequences.
    """
    return parse_fasta(path)


def sequences_to_matrix(
    records: List[SequenceRecord],
) -> Tuple[List[str], np.ndarray]:
    """
    Convert aligned sequences to a character matrix.

    Parameters
    ----------
    records : List[SequenceRecord]
        Aligned sequence records.

    Returns
    -------
    Tuple[List[str], np.ndarray]
        A tuple of sequence IDs and a 2D numpy character array of shape
        (n_sequences, alignment_length).

    Raises
    ------
    ValueError
        If sequences have inconsistent alignment lengths.
    """
    ids = [r.id for r in records]
    sequences = [r.sequence for r in records]
    lengths = {len(seq) for seq in sequences}
    if len(lengths) != 1:
        raise ValueError("Aligned sequences do not have the same length.")
    aln_len = lengths.pop()
    mat = np.array([list(seq) for seq in sequences], dtype="U1")
    assert mat.shape == (len(records), aln_len)
    return ids, mat


def compute_pairwise_identity_matrix(
    mat: np.ndarray,
    ids: List[str],
) -> pd.DataFrame:
    """
    Compute the pairwise sequence identity matrix from an alignment.

    Identity is computed as:
        matches / positions_compared
    where positions_compared excludes columns where either sequence has a gap.

    Parameters
    ----------
    mat : np.ndarray
        Character matrix of shape (n_sequences, alignment_length).
    ids : List[str]
        Sequence IDs.

    Returns
    -------
    pd.DataFrame
        DataFrame of shape (n_sequences, n_sequences) with identity values in [0, 1].
    """
    n_seq, aln_len = mat.shape
    identity = np.zeros((n_seq, n_seq), dtype=float)

    for i in range(n_seq):
        for j in range(i, n_seq):
            if i == j:
                identity[i, j] = 1.0
                continue
            seq_i = mat[i]
            seq_j = mat[j]
            non_gap_mask = (seq_i != "-") & (seq_j != "-")
            positions_compared = np.count_nonzero(non_gap_mask)
            if positions_compared == 0:
                pid = 0.0
            else:
                matches = np.count_nonzero(
                    (seq_i == seq_j) & non_gap_mask
                )
                pid = matches / positions_compared
            identity[i, j] = pid
            identity[j, i] = pid

    df = pd.DataFrame(identity, index=ids, columns=ids)
    return df


def compute_average_pairwise_identity(identity_df: pd.DataFrame) -> float:
    """
    Compute the average pairwise identity over all unique sequence pairs.

    Parameters
    ----------
    identity_df : pd.DataFrame
        Symmetric identity matrix with 1.0 on the diagonal.

    Returns
    -------
    float
        Average pairwise identity in [0, 1].
    """
    n = identity_df.shape[0]
    if n <= 1:
        return float("nan")
    # Extract upper triangle without diagonal
    vals = identity_df.values
    iu = np.triu_indices(n, k=1)
    upper_vals = vals[iu]
    if upper_vals.size == 0:
        return float("nan")
    return float(np.mean(upper_vals))


def compute_fraction_fully_conserved_columns(mat: np.ndarray) -> float:
    """
    Compute the fraction of fully conserved alignment columns.

    A column is considered fully conserved if:
    - All non-gap characters in the column are identical.
    - And there is at least one non-gap character.

    Parameters
    ----------
    mat : np.ndarray
        Character matrix of shape (n_sequences, alignment_length).

    Returns
    -------
    float
        Fraction of fully conserved columns in [0, 1].
    """
    n_seq, aln_len = mat.shape
    if aln_len == 0:
        return float("nan")

    conserved_count = 0
    for col_idx in range(aln_len):
        col = mat[:, col_idx]
        non_gap = col[col != "-"]
        if non_gap.size == 0:
            # Column of all gaps, ignore
            continue
        unique_chars = np.unique(non_gap)
        if unique_chars.size == 1:
            conserved_count += 1

    return conserved_count / aln_len


def write_identity_matrix_csv(
    identity_df: pd.DataFrame,
    path: Path,
) -> None:
    """
    Write the pairwise identity matrix to CSV.

    Parameters
    ----------
    identity_df : pd.DataFrame
        Pairwise identity matrix.
    path : Path
        Output CSV file path.
    """
    identity_df.to_csv(path, float_format="%.6f")


def append_summary_metrics(
    csv_path: Path,
    avg_identity: float,
    frac_conserved: float,
) -> None:
    """
    Append summary metrics to an existing CSV file as a small block below the matrix.

    The CSV is appended with a blank line, followed by two rows:
    - Metric,Value
    - Average_pairwise_identity,<value>
    - Fraction_fully_conserved_columns,<value>

    Parameters
    ----------
    csv_path : Path
        Path to the CSV file containing the identity matrix.
    avg_identity : float
        Average pairwise identity.
    frac_conserved : float
        Fraction of fully conserved columns.
    """
    with csv_path.open("a", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([])
        writer.writerow(["Metric", "Value"])
        writer.writerow(["Average_pairwise_identity", f"{avg_identity:.6f}"])
        writer.writerow(["Fraction_fully_conserved_columns", f"{frac_conserved:.6f}"])


def ensure_outdir(path: Path) -> None:
    """
    Ensure that the output directory exists, creating it if necessary.

    Parameters
    ----------
    path : Path
        Directory path.

    Raises
    ------
    RuntimeError
        If the directory cannot be created.
    """
    try:
        path.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise RuntimeError(f"Failed to create output directory '{path}': {exc}") from exc


def parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    """
    Parse command-line arguments.

    Parameters
    ----------
    argv : List[str] or None
        Argument list or None to use sys.argv[1:].

    Returns
    -------
    argparse.Namespace
        Parsed arguments with attributes 'input' and 'outdir'.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Align homologous protein sequences and compute alignment quality metrics."
        )
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Path to unaligned protein FASTA file.",
    )
    parser.add_argument(
        "--outdir",
        "-o",
        required=True,
        help="Output directory for aligned FASTA and report CSV.",
    )
    return parser.parse_args(argv)


def main(argv: List[str] | None = None) -> None:
    """
    Main entry point for the command-line tool.

    Parameters
    ----------
    argv : List[str] or None
        Argument list or None to use sys.argv[1:].
    """
    args = parse_args(argv)

    input_path = Path(args.input)
    outdir = Path(args.outdir)

    # Validate input file
    if not input_path.exists():
        print(
            f"ERROR: Input FASTA file '{input_path}' does not exist.",
            file=sys.stderr,
        )
        sys.exit(1)
    if not input_path.is_file():
        print(
            f"ERROR: Input path '{input_path}' is not a file.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Ensure output directory
    try:
        ensure_outdir(outdir)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    aligned_fasta_path = outdir / "aligned.fasta"
    report_csv_path = outdir / "alignment_report.csv"

    # Parse input FASTA and validate
    try:
        records = parse_fasta(input_path)
        validate_sequences(records)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    # Run MUSCLE to compute MSA
    try:
        run_muscle_msa(records, aligned_fasta_path)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        # Clean up potentially incomplete output
        if aligned_fasta_path.exists():
            try:
                aligned_fasta_path.unlink()
            except OSError:
                pass
        sys.exit(1)

    # Parse aligned FASTA
    try:
        aligned_records = parse_aligned_fasta(aligned_fasta_path)
    except ValueError as exc:
        print(
            f"ERROR: Failed to parse aligned FASTA '{aligned_fasta_path}': {exc}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Convert to matrix and compute metrics
    try:
        ids, mat = sequences_to_matrix(aligned_records)
        identity_df = compute_pairwise_identity_matrix(mat, ids)
        avg_identity = compute_average_pairwise_identity(identity_df)
        frac_conserved = compute_fraction_fully_conserved_columns(mat)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    # Write outputs
    try:
        write_identity_matrix_csv(identity_df, report_csv_path)
        append_summary_metrics(report_csv_path, avg_identity, frac_conserved)
    except OSError as exc:
        print(
            f"ERROR: Failed to write output report '{report_csv_path}': {exc}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Final user-facing message
    print(f"Alignment written to: {aligned_fasta_path}")
    print(f"Report written to:    {report_csv_path}")


if __name__ == "__main__":
    main()

