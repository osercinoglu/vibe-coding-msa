from pathlib import Path
from itertools import combinations
from statistics import mean

from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo
from Bio.Align.Applications import MuscleCommandline  # or ClustalOmegaCommandline

# -------- settings --------
INPUT_FASTA = "input.fasta"
ALN_FASTA = "alignment.fasta"
MUSCLE_EXE = "muscle"  # change if needed, e.g. "muscle5" or full path
# --------------------------


def read_sequences(fasta_path):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    print(f"Read {len(records)} sequences from {fasta_path}")
    for r in records:
        print(f"  {r.id}: length {len(r.seq)}")
    return records


def run_muscle(input_fasta, output_fasta, muscle_exe="muscle"):
    """
    Run MUSCLE to generate a multiple sequence alignment.
    """
    muscle_cline = MuscleCommandline(
        cmd=muscle_exe, input=input_fasta, out=output_fasta
    )
    stdout, stderr = muscle_cline()
    print(f"Alignment written to {output_fasta}")


def load_alignment(aln_path):
    alignment = AlignIO.read(aln_path, "fasta")
    print(f"Alignment has {alignment.get_alignment_length()} columns")
    return alignment


def pairwise_identity(seq1, seq2, ignore_gaps=True):
    """
    Simple pairwise identity on aligned sequences.
    seq1, seq2: strings of equal length (including '-').
    """
    assert len(seq1) == len(seq2)
    matches = 0
    valid_positions = 0

    for a, b in zip(seq1, seq2):
        if ignore_gaps and ("-" in (a, b)):
            continue
        valid_positions += 1
        if a == b:
            matches += 1

    if valid_positions == 0:
        return 0.0
    return matches / valid_positions


def summarize_alignment(alignment):
    """
    Print basic stats:
    - consensus sequence
    - pairwise % identities
    - mean identity
    - per-column conservation (fraction of most common residue)
    """
    # Consensus (simple majority rule)
    summary = AlignInfo.SummaryInfo(alignment)
    consensus = summary.dumb_consensus(threshold=0.5, ambiguous="X")
    print("\nConsensus sequence:")
    print(str(consensus))

    # Pairwise identities
    print("\nPairwise sequence identities:")
    ids = []
    for rec1, rec2 in combinations(alignment, 2):
        pid = pairwise_identity(str(rec1.seq), str(rec2.seq)) * 100
        ids.append(pid)
        print(f"  {rec1.id} vs {rec2.id}: {pid:.1f}%")

    if ids:
        print(f"\nMean pairwise identity: {mean(ids):.1f}%")

    # Simple conservation per column (fraction of most common residue)
    aln_len = alignment.get_alignment_length()
    col_conservation = []
    for i in range(aln_len):
        column = alignment[:, i]
        # remove gaps
        residues = [c for c in column if c != "-"]
        if not residues:
            col_conservation.append(0.0)
            continue
        # frequency of the most common residue
        max_freq = max(residues.count(r) for r in set(residues))
        col_conservation.append(max_freq / len(residues))

    mean_cons = mean(col_conservation)
    print(f"\nMean column conservation: {mean_cons:.3f}")
    print("First 20 column conservation values:")
    print(", ".join(f"{c:.2f}" for c in col_conservation[:20]))


def main():
    # 1) Read sequences
    if not Path(INPUT_FASTA).exists():
        raise FileNotFoundError(f"{INPUT_FASTA} not found")
    records = read_sequences(INPUT_FASTA)

    # Optionally: write a cleaned FASTA (not needed here since input is fine)
    SeqIO.write(records, INPUT_FASTA, "fasta")

    # 2) Run MUSCLE to produce MSA
    run_muscle(INPUT_FASTA, ALN_FASTA, muscle_exe=MUSCLE_EXE)

    # 3) Load alignment
    alignment = load_alignment(ALN_FASTA)

    # 4) Summary statistics
    summarize_alignment(alignment)


if __name__ == "__main__":
    main()

