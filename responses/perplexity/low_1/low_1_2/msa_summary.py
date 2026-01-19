from pathlib import Path
from itertools import combinations
from statistics import mean
import subprocess

from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo

# -------- settings --------
INPUT_FASTA = "input.fasta"
ALN_FASTA = "alignment.fasta"
MUSCLE_EXE = "muscle"  # or full path to your muscle binary
# --------------------------


def read_sequences(fasta_path):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    print(f"Read {len(records)} sequences from {fasta_path}")
    for r in records:
        print(f"  {r.id}: length {len(r.seq)}")
    return records


def run_muscle_v5(input_fasta, output_fasta, muscle_exe="muscle"):
    """
    Run MUSCLE v5 from the command line using subprocess.
    For a small protein set, the 'align' command is appropriate.
    """
    cmd = [
        muscle_exe,
        "-align", input_fasta,
        "-output", output_fasta,
    ]
    print("Running:", " ".join(cmd))
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
        raise RuntimeError(
            f"MUSCLE failed with return code {result.returncode}"
        )
    print(f"Alignment written to {output_fasta}")


def load_alignment(aln_path):
    alignment = AlignIO.read(aln_path, "fasta")
    print(f"Alignment has {alignment.get_alignment_length()} columns")
    return alignment


def pairwise_identity(seq1, seq2, ignore_gaps=True):
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
    summary = AlignInfo.SummaryInfo(alignment)
    consensus = summary.dumb_consensus(threshold=0.5, ambiguous="X")
    print("\nConsensus sequence:")
    print(str(consensus))

    print("\nPairwise sequence identities:")
    ids = []
    for rec1, rec2 in combinations(alignment, 2):
        pid = pairwise_identity(str(rec1.seq), str(rec2.seq)) * 100
        ids.append(pid)
        print(f"  {rec1.id} vs {rec2.id}: {pid:.1f}%")

    if ids:
        print(f"\nMean pairwise identity: {mean(ids):.1f}%")

    aln_len = alignment.get_alignment_length()
    col_conservation = []
    for i in range(aln_len):
        column = alignment[:, i]
        residues = [c for c in column if c != "-"]
        if not residues:
            col_conservation.append(0.0)
            continue
        max_freq = max(residues.count(r) for r in set(residues))
        col_conservation.append(max_freq / len(residues))

    mean_cons = mean(col_conservation)
    print(f"\nMean column conservation: {mean_cons:.3f}")
    print("First 20 column conservation values:")
    print(", ".join(f"{c:.2f}" for c in col_conservation[:20]))


def main():
    if not Path(INPUT_FASTA).exists():
        raise FileNotFoundError(f"{INPUT_FASTA} not found")
    read_sequences(INPUT_FASTA)

    # Run MUSCLE v5
    run_muscle_v5(INPUT_FASTA, ALN_FASTA, muscle_exe=MUSCLE_EXE)

    # Load alignment and summarize
    alignment = load_alignment(ALN_FASTA)
    summarize_alignment(alignment)


if __name__ == "__main__":
    main()

