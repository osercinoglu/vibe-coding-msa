#!/usr/bin/env python3
"""
Multiple Sequence Alignment (MSA) and Analysis Script
Uses a progressive alignment strategy with pairwise alignment of sequences.
"""

import argparse
import sys
from typing import List, Tuple, Dict
import numpy as np


# ============================================================================
# FASTA I/O
# ============================================================================

def read_fasta(filename: str) -> List[Tuple[str, str]]:
    """Read a FASTA file and return list of (header, sequence) tuples."""
    sequences = []
    current_header = None
    current_seq = []
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if current_header is not None:
                        sequences.append((current_header, ''.join(current_seq)))
                    current_header = line[1:]  # Remove '>'
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Don't forget the last sequence
            if current_header is not None:
                sequences.append((current_header, ''.join(current_seq)))
    
    except IOError as e:
        print(f"Error reading file {filename}: {e}", file=sys.stderr)
        sys.exit(1)
    
    if not sequences:
        print("Error: No sequences found in FASTA file.", file=sys.stderr)
        sys.exit(1)
    
    return sequences


def write_fasta(filename: str, sequences: List[Tuple[str, str]]) -> None:
    """Write sequences to a FASTA file."""
    try:
        with open(filename, 'w') as f:
            for header, seq in sequences:
                f.write(f'>{header}\n')
                # Write sequence in blocks of 60 characters
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + '\n')
    except IOError as e:
        print(f"Error writing file {filename}: {e}", file=sys.stderr)
        sys.exit(1)


# ============================================================================
# Scoring matrices and functions
# ============================================================================

def get_blosum62_score(aa1: str, aa2: str) -> int:
    """
    Return BLOSUM62 score for two amino acids.
    Returns 0 for gaps and gap characters.
    """
    # Simplified BLOSUM62 matrix (diagonal and common substitutions)
    blosum62 = {
        ('A', 'A'): 4, ('R', 'R'): 5, ('N', 'N'): 6, ('D', 'D'): 6,
        ('C', 'C'): 9, ('Q', 'Q'): 5, ('E', 'E'): 5, ('G', 'G'): 7,
        ('H', 'H'): 8, ('I', 'I'): 4, ('L', 'L'): 4, ('K', 'K'): 5,
        ('M', 'M'): 5, ('F', 'F'): 6, ('P', 'P'): 7, ('S', 'S'): 4,
        ('T', 'T'): 5, ('W', 'W'): 11, ('Y', 'Y'): 7, ('V', 'V'): 4,
        # Common matches
        ('A', 'S'): 1, ('A', 'T'): 0, ('A', 'V'): 0,
        ('R', 'K'): 2, ('N', 'D'): 1, ('N', 'S'): 1,
        ('D', 'E'): 2, ('Q', 'E'): 2, ('Q', 'K'): 1,
        ('L', 'M'): 2, ('L', 'I'): 2, ('L', 'V'): 1,
        ('F', 'Y'): 3, ('F', 'W'): 1,
    }
    
    # Handle gaps
    if aa1 in ['-', '.'] or aa2 in ['-', '.']:
        return -4
    
    # Normalize to upper case
    aa1, aa2 = aa1.upper(), aa2.upper()
    
    # Look up score (symmetric)
    if (aa1, aa2) in blosum62:
        return blosum62[(aa1, aa2)]
    elif (aa2, aa1) in blosum62:
        return blosum62[(aa2, aa1)]
    else:
        return -1  # Default mismatch


def pairwise_identity(seq1: str, seq2: str) -> float:
    """
    Calculate pairwise sequence identity (0 to 1).
    Ignores gaps in the calculation.
    """
    if len(seq1) != len(seq2):
        return 0.0
    
    matches = 0
    valid_positions = 0
    
    for aa1, aa2 in zip(seq1, seq2):
        # Skip if either is a gap
        if aa1 in ['-', '.'] or aa2 in ['-', '.']:
            continue
        valid_positions += 1
        if aa1.upper() == aa2.upper():
            matches += 1
    
    if valid_positions == 0:
        return 0.0
    
    return matches / valid_positions


def needleman_wunsch(seq1: str, seq2: str, gap_open: int = -10, 
                     gap_extend: int = -1) -> Tuple[str, str, int]:
    """
    Needleman-Wunsch global alignment of two sequences.
    Returns (aligned_seq1, aligned_seq2, score).
    """
    m, n = len(seq1), len(seq2)
    
    # Initialize DP table and traceback table
    dp = np.full((m + 1, n + 1), -np.inf, dtype=float)
    traceback = np.zeros((m + 1, n + 1), dtype=int)  # 0:diag, 1:up, 2:left
    
    # Fill first row and column
    dp[0, 0] = 0
    for i in range(1, m + 1):
        dp[i, 0] = gap_open + (i - 1) * gap_extend
        traceback[i, 0] = 1  # from up
    for j in range(1, n + 1):
        dp[0, j] = gap_open + (j - 1) * gap_extend
        traceback[0, j] = 2  # from left
    
    # Fill DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score_match = dp[i - 1, j - 1] + get_blosum62_score(seq1[i - 1], seq2[j - 1])
            score_delete = dp[i - 1, j] + (gap_extend if traceback[i - 1, j] == 1 else gap_open)
            score_insert = dp[i, j - 1] + (gap_extend if traceback[i, j - 1] == 2 else gap_open)
            
            dp[i, j] = max(score_match, score_delete, score_insert)
            
            if dp[i, j] == score_match:
                traceback[i, j] = 0
            elif dp[i, j] == score_delete:
                traceback[i, j] = 1
            else:
                traceback[i, j] = 2
    
    # Traceback to get alignment
    aligned1, aligned2 = [], []
    i, j = m, n
    
    while i > 0 or j > 0:
        if traceback[i, j] == 0:  # diagonal
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif traceback[i, j] == 1:  # up
            aligned1.append(seq1[i - 1])
            aligned2.append('-')
            i -= 1
        else:  # left
            aligned1.append('-')
            aligned2.append(seq2[j - 1])
            j -= 1
    
    aligned1.reverse()
    aligned2.reverse()
    
    return ''.join(aligned1), ''.join(aligned2), int(dp[m, n])


# ============================================================================
# Progressive MSA
# ============================================================================

def guide_tree_upgma(sequences: List[str]) -> List[Tuple[int, int, float]]:
    """
    Build a simplified guide tree using UPGMA-like clustering.
    Returns list of (cluster1_idx, cluster2_idx, distance).
    """
    n = len(sequences)
    
    # Initialize clusters
    clusters = [[i] for i in range(n)]
    distances = {}
    
    # Compute pairwise distances
    for i in range(n):
        for j in range(i + 1, n):
            dist = 1.0 - pairwise_identity(sequences[i], sequences[j])
            distances[(i, j)] = dist
    
    tree = []
    
    # Hierarchical clustering
    while len(clusters) > 1:
        # Find closest pair
        min_dist = float('inf')
        merge_i, merge_j = -1, -1
        
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                # Average linkage
                avg_dist = 0
                pair_count = 0
                for seq_i in clusters[i]:
                    for seq_j in clusters[j]:
                        if seq_i < seq_j:
                            pair_count += 1
                            avg_dist += distances.get((seq_i, seq_j), 1.0)
                        else:
                            pair_count += 1
                            avg_dist += distances.get((seq_j, seq_i), 1.0)
                
                if pair_count > 0:
                    avg_dist /= pair_count
                
                if avg_dist < min_dist:
                    min_dist = avg_dist
                    merge_i, merge_j = i, j
        
        if merge_i == -1:
            break
        
        # Merge clusters
        new_cluster = clusters[merge_i] + clusters[merge_j]
        tree.append((merge_i, merge_j, min_dist))
        
        # Remove old clusters (largest index first to maintain indices)
        if merge_i > merge_j:
            del clusters[merge_i]
            del clusters[merge_j]
        else:
            del clusters[merge_j]
            del clusters[merge_i]
        
        clusters.append(new_cluster)
    
    return tree


def progressive_align(sequences: List[str]) -> List[str]:
    """
    Perform progressive multiple sequence alignment.
    Returns list of aligned sequences.
    """
    n = len(sequences)
    
    if n == 1:
        return sequences
    
    # Build guide tree
    tree = guide_tree_upgma(sequences)
    
    # Initialize aligned sequences
    aligned = [seq for seq in sequences]
    
    # Progressive alignment along tree
    next_cluster_idx = n
    
    for cluster_i, cluster_j, _ in tree:
        # Align representatives from clusters
        seq_i = aligned[cluster_i]
        seq_j = aligned[cluster_j]
        
        aligned_i, aligned_j, _ = needleman_wunsch(seq_i, seq_j)
        
        # Store new aligned cluster
        aligned.append(aligned_i)
        aligned[cluster_i] = aligned_i
        aligned[cluster_j] = aligned_j
        
        next_cluster_idx += 1
    
    return aligned[:n]


def refine_alignment(aligned_seqs: List[str], iterations: int = 2) -> List[str]:
    """
    Refine alignment by realigning each sequence to the profile.
    """
    for iteration in range(iterations):
        refined = []
        
        for idx, query in enumerate(aligned_seqs):
            # Build profile without this sequence
            profile_seqs = aligned_seqs[:idx] + aligned_seqs[idx + 1:]
            
            if not profile_seqs:
                refined.append(query)
                continue
            
            # Align query to consensus of profile
            consensus = build_consensus(profile_seqs)
            aligned_query, _, _ = needleman_wunsch(query.replace('-', ''), consensus)
            refined.append(aligned_query)
        
        aligned_seqs = refined
    
    return aligned_seqs


def build_consensus(aligned_seqs: List[str], threshold: float = 0.5) -> str:
    """
    Build consensus sequence from aligned sequences.
    Uses the most frequent residue at each position (ignoring gaps).
    """
    if not aligned_seqs or len(aligned_seqs[0]) == 0:
        return ""
    
    consensus = []
    aln_len = len(aligned_seqs[0])
    
    for pos in range(aln_len):
        residues = {}
        total = 0
        
        for seq in aligned_seqs:
            if pos < len(seq):
                aa = seq[pos].upper()
                if aa != '-' and aa != '.':
                    residues[aa] = residues.get(aa, 0) + 1
                    total += 1
        
        if residues:
            most_common = max(residues.items(), key=lambda x: x[1])[0]
            consensus.append(most_common)
        else:
            consensus.append('X')
    
    return ''.join(consensus)


# ============================================================================
# Analysis
# ============================================================================

def compute_pairwise_identity_matrix(aligned_seqs: List[str]) -> np.ndarray:
    """Compute pairwise identity matrix for aligned sequences."""
    n = len(aligned_seqs)
    matrix = np.zeros((n, n), dtype=float)
    
    for i in range(n):
        for j in range(n):
            matrix[i, j] = pairwise_identity(aligned_seqs[i], aligned_seqs[j])
    
    return matrix


def compute_average_pairwise_identity(matrix: np.ndarray) -> float:
    """Compute average pairwise identity from matrix (excluding diagonal)."""
    n = len(matrix)
    
    if n <= 1:
        return 1.0
    
    total = 0.0
    count = 0
    
    for i in range(n):
        for j in range(i + 1, n):
            total += matrix[i, j]
            count += 1
    
    return total / count if count > 0 else 0.0


def compute_conservation(aligned_seqs: List[str]) -> Tuple[float, int, int]:
    """
    Compute fraction of fully conserved columns.
    Returns (fraction_conserved, num_conserved, total_columns).
    Ignores columns with any gaps.
    """
    if not aligned_seqs or len(aligned_seqs[0]) == 0:
        return 0.0, 0, 0
    
    aln_len = len(aligned_seqs[0])
    num_conserved = 0
    total_columns = 0
    
    for pos in range(aln_len):
        residues = set()
        has_gap = False
        
        for seq in aligned_seqs:
            if pos < len(seq):
                aa = seq[pos].upper()
                if aa in ['-', '.']:
                    has_gap = True
                    break
                residues.add(aa)
        
        if not has_gap and residues:
            total_columns += 1
            if len(residues) == 1:
                num_conserved += 1
    
    fraction = num_conserved / total_columns if total_columns > 0 else 0.0
    
    return fraction, num_conserved, total_columns


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Multiple Sequence Alignment and Analysis"
    )
    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Input FASTA file'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output aligned FASTA file'
    )
    
    args = parser.parse_args()
    
    # Read input
    print("Reading FASTA file...", file=sys.stderr)
    fasta_sequences = read_fasta(args.input)
    headers = [h for h, _ in fasta_sequences]
    sequences = [s for _, s in fasta_sequences]
    
    print(f"Read {len(sequences)} sequences.", file=sys.stderr)
    
    # Perform MSA
    print("Performing progressive multiple sequence alignment...", file=sys.stderr)
    aligned_seqs = progressive_align(sequences)
    
    print("Refining alignment...", file=sys.stderr)
    aligned_seqs = refine_alignment(aligned_seqs, iterations=1)
    
    # Compute analyses
    print("\n" + "=" * 70, file=sys.stderr)
    print("ALIGNMENT ANALYSIS RESULTS", file=sys.stderr)
    print("=" * 70 + "\n", file=sys.stderr)
    
    # Pairwise identity matrix
    print("Pairwise Sequence Identity Matrix:", file=sys.stderr)
    print("-" * 70, file=sys.stderr)
    identity_matrix = compute_pairwise_identity_matrix(aligned_seqs)
    
    # Print matrix header
    header_format = "{:<12}" + "{:>10}" * len(headers)
    print(header_format.format("", *[h[:9] for h in headers]), file=sys.stderr)
    
    for i, header in enumerate(headers):
        row_format = "{:<12}" + "{:>10.4f}" * len(headers)
        print(row_format.format(header[:11], *identity_matrix[i]), file=sys.stderr)
    
    # Average pairwise identity
    avg_identity = compute_average_pairwise_identity(identity_matrix)
    print(f"\nAverage Pairwise Identity: {avg_identity:.4f}", file=sys.stderr)
    
    # Conservation analysis
    conservation_frac, num_conserved, total_cols = compute_conservation(aligned_seqs)
    print(f"\nConservation Analysis:", file=sys.stderr)
    print(f"  Fully conserved columns: {num_conserved} / {total_cols}", file=sys.stderr)
    print(f"  Fraction conserved: {conservation_frac:.4f}", file=sys.stderr)
    
    print("\n" + "=" * 70 + "\n", file=sys.stderr)
    
    # Write aligned sequences
    print(f"Writing aligned sequences to {args.output}...", file=sys.stderr)
    aligned_output = list(zip(headers, aligned_seqs))
    write_fasta(args.output, aligned_output)
    
    print("Done!", file=sys.stderr)


if __name__ == '__main__':
    main()
