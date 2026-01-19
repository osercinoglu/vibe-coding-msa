import subprocess
from Bio import AlignIO
from Bio.SeqIO import parse
import io

# 1. Run the alignment using MAFFT via subprocess
# This assumes 'mafft' is installed on your system
def run_mafft_alignment(input_file):
    print(f"Running MAFFT on {input_file}...")
    try:
        # Executes: mafft input.fasta
        result = subprocess.run(
            ['mafft', '--auto', input_file], 
            capture_output=True, 
            text=True, 
            check=True
        )
        # Parse the output string directly into an alignment object
        alignment = AlignIO.read(io.StringIO(result.stdout), "fasta")
        return alignment
    except FileNotFoundError:
        print("Error: 'mafft' command not found. Please install mafft (e.g., 'conda install mafft').")
        return None

# 2. Calculate a simple Similarity Score (Percent Identity)
def get_summary_stats(alignment):
    num_seqs = len(alignment)
    aln_len = alignment.get_alignment_length()
    
    print(f"\n--- Alignment Summary ---")
    print(f"Number of sequences: {num_seqs}")
    print(f"Alignment Length:    {aln_len}")
    
    # Calculate average pairwise identity
    identities = []
    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            seq1 = alignment[i].seq
            seq2 = alignment[j].seq
            matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != "-")
            identities.append(matches / aln_len)
            
    avg_identity = (sum(identities) / len(identities)) * 100
    print(f"Average Pairwise Identity: {avg_identity:.2f}%")

# Main Execution
alignment = run_mafft_alignment("input.fasta")
if alignment:
    get_summary_stats(alignment)
    # Save the alignment for future use
    AlignIO.write(alignment, "aligned_output.fasta", "fasta")
    print("\nAlignment saved to 'aligned_output.fasta'")
