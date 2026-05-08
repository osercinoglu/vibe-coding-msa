from Bio import AlignIO, SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator
import pandas as pd

# 1. Load the sequences
file_name = "input.fasta"
sequences = list(SeqIO.parse(file_name, "fasta"))
print(f"Loaded {len(sequences)} sequences.\n")

# 2. Perform the Alignment
# Note: This assumes 'clustalo' is installed and in your PATH. 
# For a quick Python-only alternative, you can use online tools, 
# but for local scripts, calling a wrapper is standard.
import subprocess

input_file = "input.fasta"
output_file = "aligned.fasta"

# Calling Clustal Omega via subprocess
subprocess.run(["clustalo", "-i", input_file, "-o", output_file, "--force"])

# 3. Read the Alignment
alignment = AlignIO.read(output_file, "fasta")
print("Multiple Sequence Alignment produced successfully.\n")

# 4. Calculate Similarity Statistics
# We use a DistanceCalculator to find the 'distance', then convert to identity
calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment)

# Convert distance to percentage identity (1 - distance) * 100
ids = [seq.id for seq in alignment]
matrix_data = []
for row in distance_matrix:
    matrix_data.append([(1 - d) * 100 for d in row])

df_identity = pd.DataFrame(matrix_data, index=ids, columns=ids)

print("--- Percentage Identity Matrix ---")
print(df_identity.round(2))

# 5. Simple Summary Stats
print(f"\nAlignment Length: {alignment.get_alignment_length()}")
