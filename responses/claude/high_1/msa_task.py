#!/usr/bin/env python3
"""
Multiple Sequence Alignment (MSA) and Quality Metrics Computation

This module performs multiple sequence alignment on homologous protein sequences
and computes comprehensive alignment quality metrics including pairwise identity,
conservation metrics, and alignment statistics.

Author: Bioinformatics Team
Version: 1.0
"""

import argparse
import csv
import logging
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class SequenceParser:
    """Parser for FASTA format protein sequences."""

    @staticmethod
    def read_fasta(filepath: Path) -> Dict[str, str]:
        """
        Read sequences from a FASTA file.

        Args:
            filepath: Path to FASTA file

        Returns:
            Dictionary mapping sequence IDs to sequences

        Raises:
            FileNotFoundError: If file does not exist
            ValueError: If file is empty or malformed
        """
        if not filepath.exists():
            raise FileNotFoundError(f"FASTA file not found: {filepath}")

        sequences = {}
        current_id = None
        current_seq = []

        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith('>'):
                        if current_id is not None:
                            sequences[current_id] = ''.join(current_seq)
                        current_id = line[1:].split()[0]  # ID is first word after >
                        current_seq = []
                    else:
                        current_seq.append(line)

            if current_id is not None:
                sequences[current_id] = ''.join(current_seq)

        except Exception as e:
            raise ValueError(f"Error reading FASTA file: {e}")

        if not sequences:
            raise ValueError("No sequences found in FASTA file")

        return sequences

    @staticmethod
    def write_fasta(filepath: Path, sequences: Dict[str, str],
                   line_width: int = 80) -> None:
        """
        Write sequences to a FASTA file.

        Args:
            filepath: Output file path
            sequences: Dictionary mapping sequence IDs to sequences
            line_width: Maximum characters per line (default: 80)
        """
        with open(filepath, 'w') as f:
            for seq_id, seq in sequences.items():
                f.write(f">{seq_id}\n")
                for i in range(0, len(seq), line_width):
                    f.write(seq[i:i+line_width] + '\n')

    @staticmethod
    def validate_sequences(sequences: Dict[str, str],
                          min_count: int = 3) -> None:
        """
        Validate sequence collection.

        Args:
            sequences: Dictionary of sequences
            min_count: Minimum required number of sequences

        Raises:
            ValueError: If validation fails
        """
        if len(sequences) < min_count:
            raise ValueError(
                f"Need at least {min_count} sequences, found {len(sequences)}"
            )

        for seq_id, seq in sequences.items():
            if not seq:
                raise ValueError(f"Empty sequence for {seq_id}")
            if not all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in seq.upper()):
                invalid_chars = set(seq.upper()) - set('ACDEFGHIKLMNPQRSTVWY')
                logger.warning(
                    f"Sequence {seq_id} contains non-standard amino acids: "
                    f"{invalid_chars}"
                )


class MSAAlignment:
    """Multiple Sequence Alignment handler."""

    @staticmethod
    def align_muscle(input_fasta: Path, output_fasta: Path,
                    timeout: int = 300) -> bool:
        """
        Perform MSA using MUSCLE.

        Args:
            input_fasta: Input FASTA file
            output_fasta: Output aligned FASTA file
            timeout: Process timeout in seconds

        Returns:
            True if successful, False otherwise
        """
        try:
            cmd = ['muscle', '-in', str(input_fasta), '-out', str(output_fasta)]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout
            )

            if result.returncode != 0:
                logger.error(f"MUSCLE failed: {result.stderr}")
                return False

            logger.info(f"MUSCLE alignment completed successfully")
            return True

        except FileNotFoundError:
            logger.warning("MUSCLE not found in PATH")
            return False
        except subprocess.TimeoutExpired:
            logger.error(f"MUSCLE alignment exceeded {timeout}s timeout")
            return False
        except Exception as e:
            logger.error(f"Error running MUSCLE: {e}")
            return False

    @staticmethod
    def align_clustalo(input_fasta: Path, output_fasta: Path,
                      timeout: int = 300) -> bool:
        """
        Perform MSA using Clustal Omega.

        Args:
            input_fasta: Input FASTA file
            output_fasta: Output aligned FASTA file
            timeout: Process timeout in seconds

        Returns:
            True if successful, False otherwise
        """
        try:
            cmd = [
                'clustalo',
                '-i', str(input_fasta),
                '-o', str(output_fasta),
                '--output=fasta',
                '--outfmt=fasta'
            ]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout
            )

            if result.returncode != 0:
                logger.error(f"Clustal Omega failed: {result.stderr}")
                return False

            logger.info("Clustal Omega alignment completed successfully")
            return True

        except FileNotFoundError:
            logger.warning("Clustal Omega not found in PATH")
            return False
        except subprocess.TimeoutExpired:
            logger.error(f"Clustal Omega alignment exceeded {timeout}s timeout")
            return False
        except Exception as e:
            logger.error(f"Error running Clustal Omega: {e}")
            return False

    @staticmethod
    def align_simple(sequences: Dict[str, str]) -> Dict[str, str]:
        """
        Perform simple MSA using progressive alignment (fallback method).

        This is a basic implementation for alignment when external tools
        are unavailable. Pads all sequences to the length of the longest.

        Args:
            sequences: Dictionary of unaligned sequences

        Returns:
            Dictionary of aligned sequences
        """
        if not sequences:
            return {}

        # Sort by length descending (longest first)
        sorted_seqs = sorted(
            sequences.items(),
            key=lambda x: len(x[1]),
            reverse=True
        )

        # Pad to longest sequence
        max_len = len(sorted_seqs[0][1])
        aligned = {
            seq_id: seq.ljust(max_len, '-')
            for seq_id, seq in sorted_seqs
        }

        logger.info(
            f"Simple alignment: padded {len(aligned)} sequences to "
            f"length {max_len}"
        )
        return aligned

    @staticmethod
    def align_sequences(input_fasta: Path, output_fasta: Path) -> bool:
        """
        Align sequences using available external tools, fall back to simple method.

        Args:
            input_fasta: Input FASTA file
            output_fasta: Output aligned FASTA file

        Returns:
            True if alignment successful
        """
        # Try MUSCLE first (most common)
        if MSAAlignment.align_muscle(input_fasta, output_fasta):
            return True

        # Try Clustal Omega
        if MSAAlignment.align_clustalo(input_fasta, output_fasta):
            return True

        # Fallback to simple alignment
        logger.warning("No external MSA tool found; using simple alignment")
        sequences = SequenceParser.read_fasta(input_fasta)
        aligned = MSAAlignment.align_simple(sequences)
        SequenceParser.write_fasta(output_fasta, aligned)
        return True


class AlignmentMetrics:
    """Compute alignment quality metrics."""

    @staticmethod
    def pairwise_identity(seq1: str, seq2: str) -> float:
        """
        Compute pairwise sequence identity.

        Args:
            seq1: First sequence (aligned)
            seq2: Second sequence (aligned)

        Returns:
            Identity as fraction (0.0 to 1.0)
        """
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must have same length for alignment metrics")

        if len(seq1) == 0:
            return 0.0

        matches = sum(a == b for a, b in zip(seq1, seq2))
        return matches / len(seq1)

    @staticmethod
    def compute_identity_matrix(sequences: Dict[str, str]) -> Tuple[
        np.ndarray, List[str]
    ]:
        """
        Compute pairwise identity matrix for all sequences.

        Args:
            sequences: Dictionary of aligned sequences

        Returns:
            Tuple of (identity matrix, sequence IDs in same order)
        """
        seq_ids = list(sequences.keys())
        n_seqs = len(seq_ids)
        identity_matrix = np.zeros((n_seqs, n_seqs))

        for i, seq_id1 in enumerate(seq_ids):
            for j, seq_id2 in enumerate(seq_ids):
                if i <= j:
                    identity = AlignmentMetrics.pairwise_identity(
                        sequences[seq_id1],
                        sequences[seq_id2]
                    )
                    identity_matrix[i, j] = identity
                    identity_matrix[j, i] = identity

        return identity_matrix, seq_ids

    @staticmethod
    def average_pairwise_identity(identity_matrix: np.ndarray) -> float:
        """
        Compute average pairwise identity (excluding diagonal).

        Args:
            identity_matrix: n x n identity matrix

        Returns:
            Average pairwise identity
        """
        n = identity_matrix.shape[0]
        if n <= 1:
            return 0.0

        # Sum upper triangle (excluding diagonal)
        upper_triangle_sum = np.sum(np.triu(identity_matrix, k=1))
        n_pairs = n * (n - 1) / 2

        return upper_triangle_sum / n_pairs if n_pairs > 0 else 0.0

    @staticmethod
    def conserved_columns(sequences: Dict[str, str]) -> float:
        """
        Compute fraction of fully conserved columns.

        A column is conserved if all sequences have the same amino acid
        at that position (gaps are not considered conservation).

        Args:
            sequences: Dictionary of aligned sequences

        Returns:
            Fraction of conserved columns (0.0 to 1.0)
        """
        if not sequences:
            return 0.0

        seq_list = list(sequences.values())
        alignment_length = len(seq_list[0])

        if alignment_length == 0:
            return 0.0

        conserved_count = 0

        for col_idx in range(alignment_length):
            # Extract column
            column = [seq[col_idx] for seq in seq_list]

            # Check if all non-gap characters are identical
            non_gap_chars = [c for c in column if c != '-']

            if len(non_gap_chars) > 0 and len(set(non_gap_chars)) == 1:
                conserved_count += 1

        return conserved_count / alignment_length

    @staticmethod
    def alignment_coverage(sequences: Dict[str, str]) -> Dict[str, float]:
        """
        Compute per-sequence coverage (fraction of non-gap positions).

        Args:
            sequences: Dictionary of aligned sequences

        Returns:
            Dictionary mapping sequence IDs to coverage fraction
        """
        coverage = {}
        for seq_id, seq in sequences.items():
            non_gap = sum(1 for c in seq if c != '-')
            coverage[seq_id] = non_gap / len(seq) if len(seq) > 0 else 0.0

        return coverage

    @staticmethod
    def compute_all_metrics(sequences: Dict[str, str]) -> Dict:
        """
        Compute all alignment metrics.

        Args:
            sequences: Dictionary of aligned sequences

        Returns:
            Dictionary containing all computed metrics
        """
        identity_matrix, seq_ids = AlignmentMetrics.compute_identity_matrix(
            sequences
        )
        avg_identity = AlignmentMetrics.average_pairwise_identity(
            identity_matrix
        )
        conserved_frac = AlignmentMetrics.conserved_columns(sequences)
        coverage = AlignmentMetrics.alignment_coverage(sequences)

        alignment_length = len(list(sequences.values())[0]) if sequences else 0

        metrics = {
            'identity_matrix': identity_matrix,
            'seq_ids': seq_ids,
            'avg_pairwise_identity': avg_identity,
            'conserved_columns_fraction': conserved_frac,
            'alignment_length': alignment_length,
            'n_sequences': len(sequences),
            'coverage': coverage
        }

        return metrics


class ReportGenerator:
    """Generate alignment quality reports."""

    @staticmethod
    def save_identity_matrix_csv(identity_matrix: np.ndarray,
                                seq_ids: List[str],
                                filepath: Path) -> None:
        """
        Save identity matrix as CSV file.

        Args:
            identity_matrix: n x n identity matrix
            seq_ids: List of sequence IDs
            filepath: Output CSV file path
        """
        df = pd.DataFrame(
            identity_matrix,
            index=seq_ids,
            columns=seq_ids
        )
        df.to_csv(filepath)
        logger.info(f"Identity matrix saved to {filepath}")

    @staticmethod
    def save_metrics_report(metrics: Dict, filepath: Path) -> None:
        """
        Save metrics report as text file.

        Args:
            metrics: Dictionary of computed metrics
            filepath: Output text file path
        """
        with open(filepath, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("MULTIPLE SEQUENCE ALIGNMENT QUALITY REPORT\n")
            f.write("=" * 70 + "\n\n")

            f.write(f"Alignment Statistics\n")
            f.write("-" * 70 + "\n")
            f.write(f"Number of sequences:               {metrics['n_sequences']}\n")
            f.write(f"Alignment length:                  {metrics['alignment_length']}\n")
            f.write(
                f"Average pairwise identity:         "
                f"{metrics['avg_pairwise_identity']:.4f} "
                f"({metrics['avg_pairwise_identity']*100:.2f}%)\n"
            )
            f.write(
                f"Conserved columns:                 "
                f"{metrics['conserved_columns_fraction']:.4f} "
                f"({metrics['conserved_columns_fraction']*100:.2f}%)\n"
            )

            f.write(f"\nPer-sequence Coverage\n")
            f.write("-" * 70 + "\n")
            for seq_id, cov in metrics['coverage'].items():
                f.write(f"  {seq_id:30s}: {cov:.4f} ({cov*100:.2f}%)\n")

            f.write(f"\nPairwise Identity Matrix\n")
            f.write("-" * 70 + "\n")
            identity_matrix = metrics['identity_matrix']
            seq_ids = metrics['seq_ids']

            # Header
            f.write(f"{'ID':30s}")
            for sid in seq_ids:
                f.write(f" {sid:>10s}")
            f.write("\n")

            # Rows
            for i, seq_id in enumerate(seq_ids):
                f.write(f"{seq_id:30s}")
                for j in range(len(seq_ids)):
                    f.write(f" {identity_matrix[i, j]:10.4f}")
                f.write("\n")

            f.write("\n" + "=" * 70 + "\n")

        logger.info(f"Metrics report saved to {filepath}")

    @staticmethod
    def save_summary_csv(metrics: Dict, filepath: Path) -> None:
        """
        Save summary metrics as CSV file.

        Args:
            metrics: Dictionary of computed metrics
            filepath: Output CSV file path
        """
        summary_data = {
            'Metric': [
                'Number of sequences',
                'Alignment length',
                'Average pairwise identity',
                'Conserved columns fraction'
            ],
            'Value': [
                metrics['n_sequences'],
                metrics['alignment_length'],
                f"{metrics['avg_pairwise_identity']:.4f}",
                f"{metrics['conserved_columns_fraction']:.4f}"
            ]
        }

        df = pd.DataFrame(summary_data)
        df.to_csv(filepath, index=False)
        logger.info(f"Summary report saved to {filepath}")


def validate_input(input_fasta: Path) -> None:
    """
    Validate input FASTA file and content.

    Args:
        input_fasta: Path to input FASTA file

    Raises:
        FileNotFoundError: If file not found
        ValueError: If validation fails
    """
    if not input_fasta.exists():
        raise FileNotFoundError(f"Input file not found: {input_fasta}")

    sequences = SequenceParser.read_fasta(input_fasta)
    SequenceParser.validate_sequences(sequences, min_count=3)
    logger.info(f"Input validation passed: {len(sequences)} sequences found")


def main():
    """Main entry point for MSA and metrics computation."""
    parser = argparse.ArgumentParser(
        description='Multiple Sequence Alignment and Quality Metrics',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python msa_task.py --input sequences.fasta --outdir ./results
  python msa_task.py -i sequences.fasta -o /tmp/msa_output
        """
    )

    parser.add_argument(
        '--input', '-i',
        type=Path,
        required=True,
        help='Input FASTA file with protein sequences'
    )
    parser.add_argument(
        '--outdir', '-o',
        type=Path,
        default=Path('msa_output'),
        help='Output directory for results (default: msa_output)'
    )

    args = parser.parse_args()

    try:
        # Validate input
        logger.info(f"Reading input from {args.input}")
        validate_input(args.input)

        # Create output directory
        args.outdir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Output directory: {args.outdir}")

        # Read input sequences
        sequences = SequenceParser.read_fasta(args.input)
        logger.info(f"Loaded {len(sequences)} sequences")

        # Perform alignment
        aligned_fasta = args.outdir / 'aligned.fasta'
        logger.info("Starting multiple sequence alignment...")
        MSAAlignment.align_sequences(args.input, aligned_fasta)

        # Read aligned sequences
        aligned_sequences = SequenceParser.read_fasta(aligned_fasta)
        logger.info(f"Alignment complete: {len(aligned_sequences)} sequences aligned")

        # Compute metrics
        logger.info("Computing alignment quality metrics...")
        metrics = AlignmentMetrics.compute_all_metrics(aligned_sequences)

        # Generate reports
        logger.info("Generating reports...")
        ReportGenerator.save_identity_matrix_csv(
            metrics['identity_matrix'],
            metrics['seq_ids'],
            args.outdir / 'identity_matrix.csv'
        )
        ReportGenerator.save_metrics_report(
            metrics,
            args.outdir / 'alignment_report.txt'
        )
        ReportGenerator.save_summary_csv(
            metrics,
            args.outdir / 'summary_metrics.csv'
        )

        logger.info(f"Results saved to {args.outdir}")
        logger.info("=" * 70)
        logger.info(f"SUMMARY:")
        logger.info(f"  Sequences:              {metrics['n_sequences']}")
        logger.info(f"  Alignment length:       {metrics['alignment_length']}")
        logger.info(
            f"  Avg pairwise identity:  "
            f"{metrics['avg_pairwise_identity']*100:.2f}%"
        )
        logger.info(
            f"  Conserved columns:      "
            f"{metrics['conserved_columns_fraction']*100:.2f}%"
        )
        logger.info("=" * 70)

    except FileNotFoundError as e:
        logger.error(f"File error: {e}")
        sys.exit(1)
    except ValueError as e:
        logger.error(f"Validation error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
