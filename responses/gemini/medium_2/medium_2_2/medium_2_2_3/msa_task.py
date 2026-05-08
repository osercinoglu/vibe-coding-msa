import argparse
import numpy as np
import pandas as pd
from itertools import combinations

def read_fasta(filepath):
    """Reads a FASTA file and returns a list of (header, sequence) tuples."""
    sequences = []
    with open(filepath, 'r') as f:
        header, seq = None, []
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if header:
                    sequences.append((header, "".join(seq)))
                header, seq = line[1:], []
            else:
                seq.append(line)
        if header:
            sequences.append((header, "".join(seq)))
    return sequences

def needleman_wunsch(seq1, seq2, match=2, mismatch=-1, gap=-2):
    """Global alignment of two sequences returning (aligned_seq1, aligned_seq2)."""
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
            align1.append("-"); align2.append(seq2[j-1])
            j -= 1
    return "".join(align1[::-1]), "".join(align2[::-1])

def star_alignment(sequences):
    """Aligns all sequences to the first sequence using a Star Alignment strategy."""
    ref_name, ref_seq = sequences[0]
    msa = {ref_name: ref_seq}
    
    for i in range(1, len(sequences)):
        target_name, target_seq = sequences[i]
        # Align new sequence to the growing reference profile
        ref_aligned, target_aligned = needleman_wunsch(msa[ref_name], target_seq)
        
        new_msa = {}
        for name, existing_seq in msa.items():
            # Update existing sequences to include any new gaps added to reference
            updated_seq = ""
            ref_idx = 0
            for char in ref_aligned:
                if char == "-" and (ref_idx >= len(existing_seq) or msa[ref_name][ref_idx] != "-"):
                    updated_seq += "-"
                else:
                    if ref_idx < len(existing_seq):
                        updated_seq += existing_seq[ref_idx]
                        ref_idx += 1
            new_msa[name] = updated_seq
        
        new_msa[target_name] = target_aligned
        msa = new_msa

    # Final length normalization
    max_len = max(len(s) for s in msa.values())
    for name in msa:
        msa[name] = msa[name].ljust(max_len, "-")
    return msa

def calculate_metrics(aligned_dict):
    names = list(aligned_dict.keys())
    seqs = list(aligned_dict.values())
    num_seqs = len(seqs)
    
    # Pairwise Identity Matrix
    matrix = np.eye(num_seqs)
    for i, j in combinations(range(num_seqs), 2):
        s1, s2 = seqs[i], seqs[j]
        matches = sum(1 for a, b in zip(s1, s2) if a == b and a != "-")
        # Identify based on alignment length
        total_len = sum(1 for a, b in zip(s1, s2) if a != "-" or b != "-")
        identity = matches / total_len if total_len > 0 else 0
        matrix[i, j] = matrix[j, i] = identity
    
    df_matrix = pd.DataFrame(matrix, index=names, columns=names)
    avg_identity = matrix[np.triu_indices(num_seqs, k=1)].mean()
    
    # Conservation
    alignment_array = np.array([list(s) for s in seqs])
    total_cols = alignment_array.shape[1]
    conserved_cols = 0
    for col in range(total_cols):
        column_data = alignment_array[:, col]
        if "-" not in column_data and len(set(column_data)) == 1:
            conserved_cols += 1
            
    fraction_conserved = conserved_cols / total_cols if total_cols > 0 else 0
    return df_matrix, avg_identity, fraction_conserved

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # [cite_start]Read the provided input file [cite: 1]
    raw_data = read_fasta(args.input)
    if not raw_data: return

    aligned_dict = star_alignment(raw_data)
    matrix, avg_id, cons = calculate_metrics(aligned_dict)

    with open(args.output, 'w') as f:
        for name, seq in aligned_dict.items():
            f.write(f">{name}\n{seq}\n")

    print("\n### Pairwise Sequence Identity Matrix ###")
    print(matrix.round(3))
    print(f"\nAverage Pairwise Identity: {avg_id:.4f}")
    print(f"Fraction of Fully Conserved Columns: {cons:.4f}")

if __name__ == "__main__":
    main()
