# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

A benchmarking experiment evaluating how well different LLMs (ChatGPT, Gemini Flash 2.5, Perplexity, Claude Haiku 4.5) perform at a standard bioinformatics task: multiple sequence alignment (MSA) of protein sequences. Input sequences come from the BAliBASE MSA benchmark (protein 3D structure-derived reference alignments). Each model was tested across three prompt stringency levels and two independent takes (runs).

Benchmark results are tracked in a Google Sheet (see `responses/chatgpt/README.md` and `responses/perplexity/README.md` for model-specific links).

## Repository Structure

- `input.fasta` — shared FASTA input for all runs (4 HMG-box domain protein sequences)
- `prompts.md` — the three prompt templates (low/medium/high stringency) plus optional benchmark context
- `responses/<model>/<stringency>_<take>/` — generated code and outputs per model per prompt level per take

### Directory naming convention

- `<stringency>_<N>` (e.g. `low_1`, `high_2`) — N is the take number (independent benchmark run)
- When a script failed and the LLM was asked to fix it, the corrected version is **nested** inside: `high_2/high_2_2/`, then `high_2/high_2_2/high_2_2_3/`, etc. Each nesting level = one additional prompting turn.
- Each directory typically contains: the generated Python script, `aligned.fasta`, and optionally a `results/` subdirectory with CSV/TXT metric reports
- `*_response.txt` files store the LLM's accompanying plain-text explanation

## Prompt Levels

**Low** — informal: "give me an alignment and basic stats, whatever you think is best"

**Medium** — structured: pairwise identity matrix + average + fraction of fully conserved columns; `--input`/`--output` CLI args; Python 3.11; only numpy/pandas/matplotlib/seaborn guaranteed

**High** — production-style: functions with docstrings, `main()`, `argparse` with `--input`/`--outdir`, input validation (≥3 sequences), graceful error handling, CSV/text reports in `results/`

## Running the Generated Scripts

```bash
# Medium stringency
python msa_task.py --input input.fasta --output aligned.fasta

# High stringency
python msa_task.py --input input.fasta --outdir results/
```

Low stringency scripts vary — check each script individually.

## Benchmark Findings

Results tracked per: Tool, Stringency, Take, Turn (prompting iteration), Outcome (Success/Error), Notes, Dependencies.

### Summary by model

**Claude Haiku 4.5** — Only take 1 exists. Ran the analysis directly within claude.ai rather than generating a standalone script; results are in `claude_response.txt`. No external tool dependency issues.

**Perplexity** — Generally reliable. Low/Medium: success within 1–2 turns. High take 2: required 3 turns due to MUSCLE command-line calling errors.

**ChatGPT** — Generally reliable. Low take 1 required a manual code replacement in the generated script. High take 2: required 2 turns (MUSCLE error on turn 1).

**Gemini Flash 2.5** — Most failures overall. High take 2 needed 5+ turns cycling through MUSCLE calling errors and syntax errors before succeeding. Medium take 2 needed 3 turns (NumPy error, then syntax error).

### Recurring failure modes

- **Biopython MUSCLE CLI** — Newer biopython versions removed the CLI wrapper for external tools; scripts calling `MuscleCommandline` fail. Affected: Perplexity low_1, Gemini low_1.
- **External tool not found** — Scripts assume `muscle`, `clustalo`, or `mafft` is on `PATH`; runs fail if the tool isn't installed. Gemini low_1 switched from muscle → mafft (also not found) before succeeding on turn 3.
- **MUSCLE command-line argument errors** — MUSCLE v5 changed its CLI interface; scripts written for v3/v4 argument style fail. Common in high stringency take 2 for Perplexity, ChatGPT, Gemini.
- **NumPy/syntax errors** — Appeared in Gemini medium/high take 2.

### Common dependencies across scripts

`biopython`, `muscle` (v3/v4 or v5 — interface differs), `clustalo`, `mafft`, `numpy`, `pandas`

## Adding New Benchmark Runs

1. Create `responses/<model>/<stringency>_<N>/`
2. Store the generated script, `aligned.fasta`, and any `results/` outputs
3. If the script required follow-up prompting, nest: `<dir>/<dir>_<N+1>/`
4. Optionally add `*_response.txt` with the LLM's explanation
5. Update the benchmark Google Sheet with Tool, Stringency, Take, Turn, Outcome, Notes, Dependencies
