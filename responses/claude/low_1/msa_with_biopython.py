#!/usr/bin/env python3
"""
Multiple Sequence Alignment with Biopython
This is the recommended version to use on your own computer.
It provides better alignment quality and more features.

Installation:
    pip install biopython

For alignment tools, you'll also need:
    - Clustal Omega (recommended): http://www.clustal.org/omega/
    - Or ClustalW: http://www.clustal.org/
    
On macOS with Homebrew:
    brew install clustal-omega
    
On Ubuntu/Debian:
    sudo apt-get install clustal-omega
"""

from Bio import AlignIO, SeqIO
import os
import sys

def read_fasta(fasta_file):
    """Read sequences from FASTA file using Biopython."""
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(record)
    return sequences

def perform_alignment(fasta_file, output_alignment="alignment.aln"):
    """
    Perform multiple sequence alignment using Clustal Omega or ClustalW.
    """
    try:
        from Bio.Align.Applications import ClustalOmegaCommandLine
        cline = ClustalOmegaCommandLine(
            infile=fasta_file,
            outfile=output_alignment,
            verbose=True,
            auto=True
        )
        print(f"Running Clustal Omega...")
        print(f"Command: {cline}")
        cline()
        alignment = AlignIO.read(output_alignment, "clustal")
        print(f"✓ Alignment completed using Clustal Omega")
        return alignment
        
    except Exception as e:
        print(f"Clustal Omega not available ({e})")
        try:
            from Bio.Align.Applications import ClustalwCommandLine
            cline = ClustalwCommandLine(
                infile=fasta_file,
                outfile=output_alignment
            )
            print(f"Running ClustalW...")
            print(f"Command: {cline}")
            cline()
            alignment = AlignIO.read(output_alignment, "clustal")
            print(f"✓ Alignment completed using ClustalW")
            return alignment
            
        except Exception as e:
            print(f"ERROR: Neither Clustal Omega nor ClustalW available")
            print(f"\nPlease install one of these tools:")
            print(f"  - Clustal Omega: http://www.clustal.org/omega/")
            print(f"  - ClustalW: http://www.clustal.org/")
            print(f"\nOn macOS with Homebrew:")
            print(f"  brew install clustal-omega")
            print(f"\nOn Ubuntu/Debian:")
            print(f"  sudo apt-get install clustal-omega")
            raise

def calculate_pairwise_identity(alignment):
    """Calculate pairwise sequence identity percentages."""
    n_sequences = len(alignment)
    print("\n" + "="*70)
    print("PAIRWISE SEQUENCE IDENTITY (%)".center(70))
    print("="*70)
    
    identities = {}
    
    for i in range(n_sequences):
        for j in range(i + 1, n_sequences):
            seq1 = str(alignment[i].seq)
            seq2 = str(alignment[j].seq)
            
            # Count matching positions
            matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-')
            # Total aligned positions (excluding gaps)
            total = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-')
            
            if total > 0:
                identity = (matches / total) * 100
            else:
                identity = 0
                
            name1 = alignment[i].id
            name2 = alignment[j].id
            identities[(name1, name2)] = identity
            
            print(f"{name1:20s} vs {name2:20s}: {identity:6.2f}%")
    
    return identities

def calculate_conservation(alignment):
    """Calculate sequence conservation at each position."""
    n_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()
    
    conservation = []
    
    for pos in range(alignment_length):
        # Get all residues at this position
        residues = [str(alignment[i].seq[pos]) for i in range(n_sequences)]
        
        # Count gaps
        gaps = residues.count('-')
        non_gap_residues = [r for r in residues if r != '-']
        
        if len(non_gap_residues) > 0:
            # Count most common residue
            from collections import Counter
            most_common_count = Counter(non_gap_residues).most_common(1)[0][1]
            conservation_score = most_common_count / len(non_gap_residues)
        else:
            conservation_score = 0
            
        conservation.append({
            'position': pos + 1,
            'residues': ''.join(residues),
            'conservation': conservation_score,
            'gaps': gaps
        })
    
    return conservation

def print_alignment_summary(sequences, alignment):
    """Print summary statistics about the alignment."""
    print("\n" + "="*70)
    print("ALIGNMENT SUMMARY".center(70))
    print("="*70)
    
    print(f"Number of sequences: {len(alignment)}")
    print(f"Alignment length: {alignment.get_alignment_length()} positions")
    
    print("\nSequence Information:")
    print("-" * 70)
    for i, record in enumerate(alignment):
        seq_str = str(record.seq)
        gaps = seq_str.count('-')
        ungapped_length = len(seq_str) - gaps
        gap_percent = (gaps / len(seq_str)) * 100 if len(seq_str) > 0 else 0
        
        print(f"  {record.id:20s}: {ungapped_length:4d} aa, {gap_percent:5.1f}% gaps")

def print_conservation_stats(conservation):
    """Print conservation statistics."""
    print("\n" + "="*70)
    print("CONSERVATION STATISTICS".center(70))
    print("="*70)
    
    conservation_scores = [c['conservation'] for c in conservation]
    
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
            print(f"  Position {c['position']:3d}: {c['residues']:20s} ({c['conservation']*100:5.1f}% conserved)")

def main(fasta_file="input.fasta"):
    """Main analysis pipeline."""
    
    # Check if file exists
    if not os.path.exists(fasta_file):
        print(f"ERROR: File '{fasta_file}' not found")
        sys.exit(1)
    
    output_alignment = "alignment.aln"
    output_fasta = "aligned_sequences.fasta"
    
    print("="*70)
    print("MULTIPLE SEQUENCE ALIGNMENT ANALYSIS (with Biopython)".center(70))
    print("="*70)
    
    # Read sequences
    print(f"\nReading sequences from {fasta_file}...")
    sequences = read_fasta(fasta_file)
    print(f"✓ Loaded {len(sequences)} sequences")
    for seq in sequences:
        print(f"  - {seq.id:20s}: {len(seq.seq)} aa")
    
    # Perform alignment
    print(f"\nPerforming multiple sequence alignment...")
    alignment = perform_alignment(fasta_file, output_alignment)
    
    # Print alignment
    print("\nAligned sequences:")
    print("-" * 70)
    for record in alignment:
        print(f"{record.id:20s}: {str(record.seq)}")
    
    # Summary statistics
    print_alignment_summary(sequences, alignment)
    
    # Pairwise identity
    calculate_pairwise_identity(alignment)
    
    # Conservation analysis
    conservation = calculate_conservation(alignment)
    print_conservation_stats(conservation)
    
    # Save results
    AlignIO.write(alignment, output_fasta, "fasta")
    print(f"\n✓ Aligned sequences saved to {output_fasta}")
    print(f"✓ Alignment file saved to {output_alignment}")
    
    print("\n" + "="*70)
    print("Analysis complete!".center(70))
    print("="*70)

if __name__ == "__main__":
    # Get filename from command line or use default
    fasta_file = sys.argv[1] if len(sys.argv) > 1 else "input.fasta"
    main(fasta_file)
