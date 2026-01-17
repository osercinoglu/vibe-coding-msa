## Promps and context on different levels to perform multiple sequence alignments (MSA) by LLMs.

### Low stringency

I have a small set of protein sequences in a FASTA file called input.fasta.
I want to get a reasonable multiple sequence alignment of these sequences and some basic summary numbers.

Please help me do this using Python in a typical scientific environment.
I’d like to:

Read the sequences from input.fasta.

Produce a multiple sequence alignment of these sequences.

See some simple summary statistics, like how similar the sequences are to each other.

You can write code or describe the steps, whatever you think is best. Just make sure it’s something I can actually run on my computer.

### Medium stringency

You are helping with practical bioinformatics programming.

Task: I have an unaligned protein FASTA file, input.fasta, containing several homologous sequences from the same family.
I want a Python script that:

Reads input.fasta.

Produces a multiple sequence alignment of all sequences.

Computes and prints:

The pairwise sequence identity matrix.

The average pairwise identity.

The fraction of alignment columns that are fully conserved (all same residue, ignoring gaps).

Saves the alignment as aligned.fasta.

Constraints:

Python 3.11 environment.

Only standard scientific libraries are guaranteed: numpy, pandas, matplotlib, seaborn.

Do not assume internet access or any non‑standard packages unless you show how to install them.

Please output a single, self‑contained Python script that I can run as:
python msa_task.py --input input.fasta --output aligned.fasta.

### High stringency
'''
ou are acting as a senior bioinformatics software engineer.

Objective: Write production‑style Python code to align a set of homologous protein sequences and compute alignment quality metrics.

Input: unaligned protein FASTA file input.fasta.
Output:

Multiple sequence alignment saved as aligned.fasta.

Text/CSV report with:

Pairwise identity matrix.

Average pairwise identity.

Fraction of fully conserved columns.

Environment:

Python 3.11.

Available libraries: numpy, pandas, matplotlib, seaborn (and standard library).

If you call an external MSA tool (like MUSCLE or Clustal Omega), assume it is installed and callable via subprocess.

Requirements:

Organize the code into small functions with docstrings.

Provide main() and if __name__ == "__main__": main().

Use argparse for --input and --outdir.

Validate that input.fasta exists and contains at least 3 sequences.

Handle common errors gracefully with informative messages.

Return only the content of msa_task.py.
'''
## Context


