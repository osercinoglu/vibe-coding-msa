# Multiple Sequence Alignment (MSA) Analysis Scripts

This folder contains two Python scripts for analyzing protein sequences with multiple sequence alignment.

## Quick Summary of Your Results

Your 4 protein sequences have been analyzed:
- **1aab_**: 83 amino acids
- **1j46_A**: 85 amino acids  
- **1k99_A**: 91 amino acids
- **2lef_A**: 86 amino acids

### Key Findings:

**Sequence Similarity:**
- Most similar pair: **1aab_ vs 1j46_A** at 22.22% identity
- Most other pairs: ~5-9% identity
- Overall: These are fairly diverse sequences (low sequence conservation)

**Alignment Conservation:**
- Average conservation: 41.7% (moderate)
- 2 highly conserved positions (≥80%)
- 46 moderately conserved positions (50-80%)
- 43 variable positions (<50%)

---

## Script Descriptions

### 1. `simple_msa.py` (Recommended for getting started)

**No external dependencies required** - runs with just Python 3.

This script implements a pure-Python multiple sequence alignment algorithm and works on any computer without needing to install additional tools.

**Usage:**
```bash
python3 simple_msa.py
```

**Pros:**
- No installation required beyond Python
- Works everywhere
- Good for quick analysis

**Cons:**
- Less sophisticated alignment algorithm than specialized tools
- Slower for large datasets

---

### 2. `msa_with_biopython.py` (Recommended for better quality)

**Requires Biopython and an external alignment tool** - produces higher quality alignments.

This is the professional-grade version that uses established bioinformatics algorithms.

**Installation:**

First, install Biopython:
```bash
pip install biopython
```

Then install an alignment tool. Choose one:

**Option A: Clustal Omega (Recommended - faster and better)**
```bash
# macOS with Homebrew
brew install clustal-omega

# Ubuntu/Debian
sudo apt-get install clustal-omega

# Or download from: http://www.clustal.org/omega/
```

**Option B: ClustalW (older but reliable)**
```bash
# macOS with Homebrew
brew install clustalw

# Ubuntu/Debian
sudo apt-get install clustalw

# Or download from: http://www.clustal.org/
```

**Usage:**
```bash
python3 msa_with_biopython.py input.fasta
```

**Pros:**
- Uses industry-standard alignment algorithms
- Better quality alignments
- Faster execution
- More robust

**Cons:**
- Requires additional software installation

---

## What the Scripts Do

Both scripts:

1. **Read your FASTA file** - Parse protein sequences in standard FASTA format
2. **Perform alignment** - Align all sequences to each other
3. **Calculate statistics:**
   - Pairwise identity percentages (how similar each pair is)
   - Alignment length and gap information
   - Conservation scores (which positions are most similar across sequences)
   - Highly conserved regions

4. **Generate output files:**
   - `alignment.txt` or `alignment.aln` - The aligned sequences in text format
   - `aligned_sequences.fasta` - The aligned sequences in FASTA format

---

## Understanding the Output

### Alignment View
Shows your sequences aligned with gaps (`-`) inserted to maximize similarity:
```
1aab_   : GKGDPKKPRGKMSSYAFFVQTS---REEHKKKHPDAS...
1j46_A  : -MQDRVK-R-PM--NAFIVWSRDQRRKMALENPRMR...
1k99_A  : MKKLKKHPDFPKKPLTPYFRFFMEKRAKYAKLHPEM...
2lef_A  : MHIKKPLNAFMLYMKEMRANVVAESTLKESAAINQI...
```

### Pairwise Identity
Shows percentage of matching amino acids between each pair:
- 22% means 22 out of 100 aligned positions match
- Higher is more similar

### Conservation
- **Highly conserved (≥80%)**: Position matches in 4/4 or 3/4 sequences
- **Moderately conserved (50-80%)**: Position matches in 2/4 sequences
- **Variable (<50%)**: Position differs in most sequences

---

## Tips for Your Own Use

### For your input.fasta:
Place your FASTA file in the same directory as the script, then run:
```bash
python3 simple_msa.py        # or
python3 msa_with_biopython.py input.fasta
```

### Creating FASTA files:
FASTA format is simple:
```
>sequence_name_1
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQV
>sequence_name_2
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQV
```

Just:
1. Start each sequence with `>` and a name (no spaces in the name)
2. Put the sequence on the next line(s)
3. Repeat for each sequence

### Interpreting low identity:
- Low sequence identity (<30%) is normal for:
  - Sequences from different organisms
  - Proteins with different functions
  - Very distant evolutionary relatives
- High conservation in some positions suggests functional importance

---

## Troubleshooting

**"ModuleNotFoundError: No module named 'Bio'"**
- Install Biopython: `pip install biopython`

**"Command 'clustal-omega' not found"** or similar
- Install Clustal Omega (see installation section above)
- Or use `simple_msa.py` instead (no external tools needed)

**Script runs but gives poor results**
- This can happen with very divergent sequences
- Try the Biopython version for better algorithms
- Check that your FASTA file is properly formatted

---

## For More Information

- Biopython: https://biopython.org/
- Clustal Omega: http://www.clustal.org/omega/
- Multiple Sequence Alignment basics: https://en.wikipedia.org/wiki/Multiple_sequence_alignment
- FASTA format: https://en.wikipedia.org/wiki/FASTA_format
