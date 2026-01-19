from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import numpy as np

# 1. Read the sequences to verify content
sequences = list(SeqIO.parse("input.fasta", "fasta"))
print(f"Loaded {len(sequences)} sequences.\n")

# 2. Perform Multiple Sequence Alignment
# Note: This requires 'clustalw' or 'clustalo' installed on your OS.
# If you don't have it, you can use online tools or the 'Muscle' wrapper.
def perform_alignment(input_file):
    # This is a placeholder for the command line call
    # Example: clustalw2 -infile=input.fasta
    print("Aligning sequences...")
    # For this example, we assume the alignment is saved as 'aligned.fasta'
    # In a real scenario, you'd call a tool here.
    pass

# 3. Calculate Similarity Statistics
def calculate_identity_matrix(alignment):
    num_seqs = len(alignment)
    matrix = np.zeros((num_seqs, num_seqs))
    
    for i in range(num_seqs):
        for j in range(num_seqs):
            seq1 = alignment[i].seq
            seq2 = alignment[j].seq
            matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != "-")
            # Simple identity: matches / length of alignment
            matrix[i, j] = (matches / len(seq1)) * 100
            
    return matrix

# Assuming you have an aligned file named 'aligned.fasta'
# For the sake of this code, I'll use your input file to demonstrate the logic
try:
    # In practice, you would load 'aligned.fasta' here
    alignment = AlignIO.read("input.fasta", "fasta") 
    
    print("--- Summary Statistics ---")
    print(f"Alignment Length: {alignment.get_alignment_length()}")
    
    matrix = calculate_identity_matrix(alignment)
    print("\nPairwise Identity Matrix (%):")
    print(matrix)
except Exception as e:
    print(f"Note: To get a real MSA, run a tool like Clustal Omega on your file first.")
