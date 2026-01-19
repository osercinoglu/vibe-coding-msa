import numpy as np
import pandas as pd
import argparse
import collections

def read_fasta(file_path):
    sequences = collections.OrderedDict()
    with open(file_path, 'r') as f:
        label = None
        seq = []
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if label: sequences[label] = "".join(seq)
                label = line[1:]
                seq = []
            else:
                seq.append(line)
        if label: sequences[label] = "".join(seq)
    return sequences

def needleman_wunsch(seq1, seq2, match=2, mismatch=-1, gap=-2):
    """Pairwise global alignment for profile building."""
    n, m = len(seq1), len(seq2)
    score = np.zeros((n + 1, m + 1))
    
    for i in range(n + 1): score[i, 0] = i * gap
    for j in range(m + 1): score[0, j] = j * gap
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            score[i, j] = max(score[i-1, j-1] + s, 
                              score[i-1, j] + gap, 
                              score[i, j-1] + gap)
    
    # Traceback
    align1, align2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score[i, j] == score[i-1, j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            align1.append(seq1[i-1]); align2.append(seq2[j-1])
            i -= 1; j -= 1
        elif i > 0 and score[i, j] == score[i-1, j] + gap:
            align1.append(seq1[i-1]); align2.append("-")
            i -= 1
        else:
            align2.append(seq2[j-1]); align1.append("-")
            j -= 1
    return "".join(align1[::-1]), "".join(align2[::-1])

def align_to_profile(profile_seqs, new_seq):
    """Simple progressive step: aligns a new sequence to the first sequence of the current MSA."""
    # In a simple implementation, we align the new sequence to the 'consensus' or first member
    ref_seq = profile_seqs[0].replace("-", "")
    p_aligned, s_aligned = needleman_wunsch(profile_seqs[0], new_seq)
    
    # Adjust all previous sequences in profile to match the new gaps introduced in the reference
    new_profile = []
    for seq in profile_seqs:
        updated_seq = ""
        p_idx = 0
        for char in p_aligned:
            if char == "-" and (p_idx >= len(seq) or seq[p_idx] != "-"):
                updated_seq += "-"
            else:
                updated_seq += seq[p_idx]
                p_idx += 1
        new_profile.append(updated_seq)
    
    new_profile.append(s_aligned)
    return new_profile

def calculate_metrics(aligned_seqs):
    names = list(aligned_seqs.keys())
    seqs = list(aligned_seqs.values())
    n = len(seqs)
    L = len(seqs[0])
    
    # Identity Matrix
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            matches = sum(1 for a, b in zip(seqs[i], seqs[j]) if a == b and a != "-")
            # Identity = matches / alignment length
            matrix[i, j] = matches / L
            
    # Conservation
    conserved_cols = 0
    for col in zip(*seqs):
        if "-" not in col and len(set(col)) == 1:
            conserved_cols += 1
    fraction_conserved = conserved_cols / L
    
    return names, matrix, fraction_conserved

def main():
    parser = argparse.ArgumentParser(description="Simple MSA Tool")
    parser.add_argument("--input", required=True, help="Input FASTA")
    parser.add_argument("--output", required=True, help="Output Aligned FASTA")
    args = parser.parse_args()

    # 1. Read
    raw_data = read_fasta(args.input)
    ids = list(raw_data.keys())
    seq_list = list(raw_data.values())

    # 2. Align (Progressive)
    current_alignment = [seq_list[0]]
    for i in range(1, len(seq_list)):
        current_alignment = align_to_profile(current_alignment, seq_list[i])
    
    msa_dict = dict(zip(ids, current_alignment))

    # 3. Save
    with open(args.output, 'w') as f:
        for header, seq in msa_dict.items():
            f.write(f">{header}\n{seq}\n")

    # 4. Metrics
    names, matrix, conservation = calculate_metrics(msa_dict)
    df_matrix = pd.DataFrame(matrix, index=names, columns=names)

    print("\n--- Pairwise Sequence Identity Matrix ---")
    print(df_matrix.round(4))
    print(f"\nAverage Pairwise Identity: {np.mean(matrix[np.triu_indices(len(names), k=1)]):.4f}")
    print(f"Fraction of fully conserved columns: {conservation:.4f}")
    print(f"\nAlignment saved to: {args.output}")

if __name__ == "__main__":
    main()
