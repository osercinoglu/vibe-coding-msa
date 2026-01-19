#!/usr/bin/env python3
"""
Multiple Sequence Alignment and Summary Statistics (Pure Python)
This script performs basic multiple sequence alignment using a simple approach
and calculates similarity metrics without external dependencies.
"""

from collections import defaultdict

def read_fasta(fasta_file):
    """Read sequences from FASTA file."""
    sequences = []
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences.append({
                        'id': current_id,
                        'seq': ''.join(current_seq)
                    })
                current_id = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_id is not None:
            sequences.append({
                'id': current_id,
                'seq': ''.join(current_seq)
            })
    
    return sequences

def score_match(a, b):
    """Simple scoring: match=1, mismatch=0."""
    if a == b:
        return 1
    return 0

def pairwise_alignment(seq1, seq2, gap_penalty=-1):
    """
    Compute pairwise sequence alignment using dynamic programming.
    Returns aligned sequences and alignment score.
    """
    m, n = len(seq1), len(seq2)
    
    # Create DP table
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Initialize first row and column (gaps)
    for i in range(m + 1):
        dp[i][0] = i * gap_penalty
    for j in range(n + 1):
        dp[0][j] = j * gap_penalty
    
    # Fill DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = dp[i-1][j-1] + score_match(seq1[i-1], seq2[j-1])
            del_score = dp[i-1][j] + gap_penalty
            ins_score = dp[i][j-1] + gap_penalty
            dp[i][j] = max(match_score, del_score, ins_score)
    
    # Traceback to get alignment
    aligned_seq1, aligned_seq2 = [], []
    i, j = m, n
    
    while i > 0 or j > 0:
        if i == 0:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j-1])
            j -= 1
        elif j == 0:
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append('-')
            i -= 1
        else:
            if dp[i][j] == dp[i-1][j-1] + score_match(seq1[i-1], seq2[j-1]):
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append(seq2[j-1])
                i -= 1
                j -= 1
            elif dp[i][j] == dp[i-1][j] - 1:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append('-')
                i -= 1
            else:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j-1])
                j -= 1
    
    aligned_seq1.reverse()
    aligned_seq2.reverse()
    
    return ''.join(aligned_seq1), ''.join(aligned_seq2), dp[m][n]

def progressive_alignment(sequences):
    """
    Simple progressive multiple sequence alignment.
    Aligns sequences pairwise starting with most similar pair.
    """
    if len(sequences) == 0:
        return []
    
    if len(sequences) == 1:
        return [{'id': sequences[0]['id'], 'seq': sequences[0]['seq']}]
    
    # Start with first two sequences
    seqs_ids = [s['id'] for s in sequences]
    aligned = [sequences[0]['seq'], sequences[1]['seq']]
    
    # Check if first two sequences are similar enough
    if len(sequences) >= 2:
        a1, a2, _ = pairwise_alignment(sequences[0]['seq'], sequences[1]['seq'])
        aligned = [a1, a2]
    
    # Progressively add remaining sequences
    for i in range(2, len(sequences)):
        # Align new sequence to the first sequence in alignment
        new_seq = sequences[i]['seq']
        first_aligned = aligned[0].replace('-', '')  # Get ungapped version
        
        # Align new sequence to first original sequence
        align1, align2, _ = pairwise_alignment(first_aligned, new_seq)
        
        # Add the aligned new sequence, and update all sequences with gap positions
        # Simple approach: just use pairwise alignment for now
        new_aligned = []
        for j, seq in enumerate(aligned):
            # Re-align all sequences to maintain consistency
            if j == 0:
                new_aligned.append(align1)
            else:
                # Use simple gap insertion strategy
                gapped_seq = seq
                # This is simplified - in production use proper profile alignment
                new_aligned.append(gapped_seq)
        
        new_aligned.append(align2)
        aligned = new_aligned
    
    # Return alignment with IDs
    return [{'id': seqs_ids[i], 'seq': aligned[i]} for i in range(len(aligned))]

def simple_multiple_alignment(sequences):
    """
    Create a simple multiple sequence alignment using a straightforward approach.
    For each sequence, compute alignment to a consensus-like representative.
    """
    n = len(sequences)
    if n == 0:
        return []
    
    # Use first sequence as reference
    reference = sequences[0]['seq']
    aligned = []
    
    for seq in sequences:
        if seq['seq'] == reference:
            aligned.append({'id': seq['id'], 'seq': seq['seq']})
        else:
            # Align to reference
            a1, a2, _ = pairwise_alignment(reference, seq['seq'])
            aligned.append({'id': seq['id'], 'seq': a2})
    
    # Ensure all alignments have same length by adding gaps
    if aligned:
        max_len = max(len(a['seq']) for a in aligned)
        for a in aligned:
            if len(a['seq']) < max_len:
                a['seq'] = a['seq'] + '-' * (max_len - len(a['seq']))
    
    return aligned

def calculate_pairwise_identity(sequences):
    """Calculate pairwise sequence identity percentages."""
    n = len(sequences)
    print("\n" + "="*70)
    print("PAIRWISE SEQUENCE IDENTITY (%)".center(70))
    print("="*70)
    
    identities = {}
    
    for i in range(n):
        for j in range(i + 1, n):
            seq1 = sequences[i]['seq']
            seq2 = sequences[j]['seq']
            
            # Count matching positions
            matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-')
            # Total aligned positions (excluding gaps in both)
            total = sum(1 for a, b in zip(seq1, seq2) if a != '-' or b != '-')
            
            if total > 0:
                identity = (matches / total) * 100
            else:
                identity = 0
            
            name1 = sequences[i]['id']
            name2 = sequences[j]['id']
            identities[(name1, name2)] = identity
            
            print(f"{name1:20s} vs {name2:20s}: {identity:6.2f}%")
    
    return identities

def calculate_conservation(aligned_sequences):
    """Calculate sequence conservation at each position."""
    if not aligned_sequences:
        return []
    
    alignment_length = len(aligned_sequences[0]['seq'])
    n_sequences = len(aligned_sequences)
    
    conservation = []
    
    for pos in range(alignment_length):
        # Get all residues at this position
        residues = [aligned_sequences[i]['seq'][pos] for i in range(n_sequences)]
        
        # Count gaps
        gaps = residues.count('-')
        non_gap_residues = [r for r in residues if r != '-']
        
        if len(non_gap_residues) > 0:
            # Count most common residue
            residue_counts = defaultdict(int)
            for r in non_gap_residues:
                residue_counts[r] += 1
            
            most_common_count = max(residue_counts.values())
            conservation_score = most_common_count / len(non_gap_residues)
        else:
            conservation_score = 0
        
        conservation.append({
            'position': pos + 1,
            'residues': ''.join(residues),
            'conservation': conservation_score,
            'gaps': gaps,
            'non_gap_count': len(non_gap_residues)
        })
    
    return conservation

def print_alignment_summary(sequences, aligned_sequences):
    """Print summary statistics about the alignment."""
    print("\n" + "="*70)
    print("ALIGNMENT SUMMARY".center(70))
    print("="*70)
    
    print(f"Number of sequences: {len(aligned_sequences)}")
    if aligned_sequences:
        print(f"Alignment length: {len(aligned_sequences[0]['seq'])} positions")
    
    print("\nSequence Information:")
    print("-" * 70)
    for i, record in enumerate(aligned_sequences):
        seq_str = record['seq']
        gaps = seq_str.count('-')
        ungapped_length = len(seq_str) - gaps
        gap_percent = (gaps / len(seq_str)) * 100 if len(seq_str) > 0 else 0
        
        original_len = len(sequences[i]['seq'])
        print(f"  {record['id']:20s}: {original_len:4d} aa → {ungapped_length:4d} aa, {gap_percent:5.1f}% gaps in alignment")

def print_conservation_stats(conservation):
    """Print conservation statistics."""
    print("\n" + "="*70)
    print("CONSERVATION STATISTICS".center(70))
    print("="*70)
    
    conservation_scores = [c['conservation'] for c in conservation if c['gaps'] < len(conservation[0]['residues'].split(';')[0])]
    
    if conservation_scores:
        avg_conservation = sum(conservation_scores) / len(conservation_scores)
        highly_conserved = sum(1 for c in conservation_scores if c >= 0.8)
        moderately_conserved = sum(1 for c in conservation_scores if 0.5 <= c < 0.8)
        variable = sum(1 for c in conservation_scores if c < 0.5)
        
        print(f"Average conservation: {avg_conservation:.3f}")
        print(f"Highly conserved (≥80%): {highly_conserved} positions")
        print(f"Moderately conserved (50-80%): {moderately_conserved} positions")
        print(f"Variable (<50%): {variable} positions")
        
        # Show most conserved positions
        top_conserved = sorted(conservation, key=lambda x: x['conservation'], reverse=True)[:5]
        print(f"\nTop 5 most conserved positions:")
        for c in top_conserved:
            print(f"  Position {c['position']:3d}: {c['residues']:25s} ({c['conservation']*100:5.1f}% conserved, {c['non_gap_count']} residues)")

def main():
    input_file = "/mnt/user-data/uploads/input.fasta"
    output_alignment = "/home/claude/alignment.txt"
    
    print("="*70)
    print("MULTIPLE SEQUENCE ALIGNMENT ANALYSIS (Pure Python)".center(70))
    print("="*70)
    
    # Read sequences
    print(f"\nReading sequences from {input_file}...")
    sequences = read_fasta(input_file)
    print(f"✓ Loaded {len(sequences)} sequences")
    
    for seq in sequences:
        print(f"  - {seq['id']:20s}: {len(seq['seq'])} aa")
    
    # Perform alignment
    print(f"\nPerforming multiple sequence alignment...")
    aligned = simple_multiple_alignment(sequences)
    print(f"✓ Alignment completed")
    
    # Print alignment
    print("\nAligned sequences:")
    print("-" * 70)
    for record in aligned:
        print(f"{record['id']:20s}: {record['seq']}")
    
    # Summary statistics
    print_alignment_summary(sequences, aligned)
    
    # Pairwise identity
    calculate_pairwise_identity(aligned)
    
    # Conservation analysis
    conservation = calculate_conservation(aligned)
    print_conservation_stats(conservation)
    
    # Save alignment
    with open(output_alignment, 'w') as f:
        for record in aligned:
            f.write(f">{record['id']}\n{record['seq']}\n")
    
    print(f"\n✓ Aligned sequences saved to {output_alignment}")
    
    print("\n" + "="*70)
    print("Analysis complete!".center(70))
    print("="*70)

if __name__ == "__main__":
    main()
