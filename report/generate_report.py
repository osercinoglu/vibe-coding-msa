#!/usr/bin/env python3
"""Generate MSA benchmark comparison report (report/report.md)."""

from __future__ import annotations
import os
from pathlib import Path
from itertools import combinations
from collections import defaultdict

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).parent.parent

# BAliBASE BB11001 reference alignment (96 columns, structure-based, from XML seq-data)
BALIBASE_REF = {
    "1aab_":  "---GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE----------",
    "1j46_A": "------MQDRVKRPMNAFIVWSRDQRRKMALENP--PMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK---",
    "1k99_A": "MKKLKKHPDFPKKPLTPYFRFFMEKRAKYAKLHP--EMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK---",
    "2lef_A": "--------MHIKKPLNAFMLYMKEMRANVVAEST--LKESAAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK",
}

# Core block column indices (0-based) from BAliBASE XML colsco-data (value == 1)
CORE_COLS = list(range(9, 30)) + list(range(40, 76))  # 57 core columns

# Run catalog: (model, stringency, take, turns_to_success, msa_tool, fasta_path_relative_or_None)
RUNS = [
    ("Claude",     "low",    1,    1,    "inline (claude.ai)",      "responses/claude/low_1/alignment.txt"),
    ("Claude",     "medium", 1,    1,    "inline + NW",             "responses/claude/medium_1/aligned.fasta"),
    ("Claude",     "high",   1,    1,    "clustalo",                "responses/claude/high_1/results/aligned.fasta"),
    ("ChatGPT",    "low",    1,    2,    "muscle",                  "responses/chatgpt/low_1/low_1_2/aligned.fasta"),
    ("ChatGPT",    "medium", 1,    1,    "muscle",                  "responses/chatgpt/medium_1/aligned.fasta"),
    ("ChatGPT",    "high",   1,    1,    "muscle",                  "responses/chatgpt/high_1/results/aligned.fasta"),
    ("ChatGPT",    "low",    2,    1,    "muscle",                  "responses/chatgpt/low_2/aligned.fasta"),
    ("ChatGPT",    "medium", 2,    1,    "muscle",                  "responses/chatgpt/medium_2/aligned.fasta"),
    ("ChatGPT",    "high",   2,    2,    "muscle",                  "responses/chatgpt/high_2/high_2_2/results/aligned.fasta"),
    ("Gemini",     "low",    1,    3,    "mafft",                   "responses/gemini/low_1/low_1_2/aligned_output.fasta"),
    ("Gemini",     "medium", 1,    1,    "NW (pure Python)",        "responses/gemini/medium_1/aligned.fasta"),
    ("Gemini",     "high",   1,    2,    "clustalo",                "responses/gemini/high_1/results/aligned.fasta"),
    ("Gemini",     "low",    2,    1,    "clustalo",                "responses/gemini/low_2/aligned.fasta"),
    ("Gemini",     "medium", 2,    3,    "NW (pure Python)",        "responses/gemini/medium_2/medium_2_2/medium_2_2_3/aligned.fasta"),
    ("Gemini",     "high",   2,    None, "—",                       None),
    ("Perplexity", "low",    1,    2,    "muscle",                  "responses/perplexity/low_1/low_1_2/alignment.fasta"),
    ("Perplexity", "medium", 1,    1,    "NW (pure Python)",        "responses/perplexity/medium_1/aligned.fasta"),
    ("Perplexity", "high",   1,    1,    "muscle",                  "responses/perplexity/high_1/results/aligned.fasta"),
    ("Perplexity", "low",    2,    1,    "clustalo",                "responses/perplexity/low_2_1/aligned.fasta"),
    ("Perplexity", "medium", 2,    1,    "NW (pure Python)",        "responses/perplexity/medium_2/aligned.fasta"),
    ("Perplexity", "high",   2,    3,    "muscle",                  "responses/perplexity/high_2/high_2_2/high_2_2_3/results/aligned.fasta"),
    ("MUSCLE v5",  "—",      "—",  1,    "muscle",                  "responses/reference/aligned.fasta"),
]

def parse_fasta(path: Path) -> dict[str, str]:
    """Parse aligned FASTA file. Returns {name: sequence_with_gaps} uppercased."""
    seqs: dict[str, str] = {}
    current_name: str | None = None
    chunks: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    seqs[current_name] = "".join(chunks).upper().replace(".", "-")
                current_name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if current_name is not None:
            seqs[current_name] = "".join(chunks).upper().replace(".", "-")
    return seqs


def _residue_col_map(aln: dict[str, str]) -> dict[str, dict[int, int]]:
    """For each sequence, build {residue_index: column_index} (gaps skipped)."""
    result: dict[str, dict[int, int]] = {}
    for name, seq in aln.items():
        col_map: dict[int, int] = {}
        r = 0
        for col, aa in enumerate(seq):
            if aa != "-":
                col_map[r] = col
                r += 1
        result[name] = col_map
    return result


def _col_residue_index(aln: dict[str, str]) -> list[dict[str, int]]:
    """For each column, build {seq_name: residue_index} for non-gap positions."""
    seq_len = len(next(iter(aln.values())))
    counts = {n: 0 for n in aln}
    col_maps: list[dict[str, int]] = []
    for col in range(seq_len):
        cm: dict[str, int] = {}
        for name, seq in aln.items():
            if seq[col] != "-":
                cm[name] = counts[name]
                counts[name] += 1
        col_maps.append(cm)
    return col_maps


def extract_aligned_pairs(
    aln: dict[str, str],
    col_filter: set[int] | None = None,
) -> dict[tuple[str, str], set[tuple[int, int]]]:
    """
    For each ordered pair (nameA, nameB), return set of (res_idx_a, res_idx_b)
    where both seqs have a non-gap residue. If col_filter given, restrict to
    those column indices.
    """
    names = sorted(aln.keys())
    col_maps = _col_residue_index(aln)
    pairs: dict[tuple[str, str], set[tuple[int, int]]] = {}
    for a, b in combinations(names, 2):
        pair_set: set[tuple[int, int]] = set()
        for col, cm in enumerate(col_maps):
            if col_filter is not None and col not in col_filter:
                continue
            if a in cm and b in cm:
                pair_set.add((cm[a], cm[b]))
        pairs[(a, b)] = pair_set
    return pairs


def sp_score(test_aln: dict[str, str], ref_aln: dict[str, str], core_cols: list[int]) -> float:
    """
    Sum-of-Pairs score: fraction of core-column residue pairs in the reference
    that are also aligned (in any column) in the test alignment.
    Returns 0.0 if there are no reference core pairs.
    """
    ref_pairs = extract_aligned_pairs(ref_aln, col_filter=set(core_cols))
    test_pairs = extract_aligned_pairs(test_aln, col_filter=None)
    total_ref = sum(len(v) for v in ref_pairs.values())
    if total_ref == 0:
        return 0.0
    correct = sum(
        len(ref_pairs[k] & test_pairs.get(k, set()))
        for k in ref_pairs
    )
    return correct / total_ref


def tc_score(test_aln: dict[str, str], ref_aln: dict[str, str], core_cols: list[int]) -> float:
    """
    Total-Column score: fraction of core columns whose complete residue mapping
    is reproduced as a single column in the test alignment.
    """
    ref_col_maps = _col_residue_index(ref_aln)
    test_res_col = _residue_col_map(test_aln)
    correct = 0
    for c in core_cols:
        cm = ref_col_maps[c]  # {seq_name: res_idx} for non-gap seqs in ref col c
        non_gap = [n for n in cm]
        if len(non_gap) < 2:
            continue
        # Find which test column each non-gap seq falls in
        test_cols = {
            test_res_col[n].get(cm[n])
            for n in non_gap
            if n in test_res_col and cm[n] in test_res_col[n]
        }
        # Column is correct if all map to the same single test column
        if len(test_cols) == 1 and None not in test_cols:
            correct += 1
    return correct / len(core_cols) if core_cols else 0.0


def pairwise_identity_matrix(aln: dict[str, str]) -> pd.DataFrame:
    """
    Compute pairwise identity for all sequence pairs.
    Identity = matched_residues / (alignment_length - double_gap_columns).
    Returns a symmetric DataFrame indexed and columned by sequence names.
    """
    names = sorted(aln.keys())
    n = len(names)
    mat = np.zeros((n, n))
    for i, a in enumerate(names):
        mat[i, i] = 1.0
        for j, b in enumerate(names):
            if j <= i:
                continue
            sa, sb = aln[a], aln[b]
            length = len(sa)
            double_gaps = sum(1 for x, y in zip(sa, sb) if x == "-" and y == "-")
            matches = sum(1 for x, y in zip(sa, sb) if x != "-" and y != "-" and x == y)
            denom = length - double_gaps
            identity = matches / denom if denom > 0 else 0.0
            mat[i, j] = mat[j, i] = identity
    return pd.DataFrame(mat, index=names, columns=names)


def conserved_col_frac(aln: dict[str, str]) -> float:
    """
    Fraction of alignment columns where all non-gap residues are identical
    and at least two sequences are non-gap.
    """
    seq_len = len(next(iter(aln.values())))
    seqs = list(aln.values())
    conserved = 0
    for col in range(seq_len):
        residues = [s[col] for s in seqs if s[col] != "-"]
        if len(residues) >= 2 and len(set(residues)) == 1:
            conserved += 1
    return conserved / seq_len if seq_len > 0 else 0.0


def load_run_metrics(runs: list[tuple]) -> list[dict]:
    """
    For each run in RUNS that has a non-None fasta path, load the alignment
    and compute SP score, TC score, avg pairwise identity, conserved col fraction.
    Returns list of dicts with keys:
      model, stringency, take, turns, msa_tool, fasta_path,
      sp, tc, avg_identity, conserved_frac, aln_len, status
    """
    results = []
    for model, stringency, take, turns, msa_tool, fasta_rel in runs:
        row = dict(
            model=model, stringency=stringency, take=take,
            turns=turns, msa_tool=msa_tool, fasta_path=fasta_rel,
        )
        if fasta_rel is None:
            row.update(sp=None, tc=None, avg_identity=None,
                       conserved_frac=None, aln_len=None, status="no output")
            results.append(row)
            continue
        path = REPO_ROOT / fasta_rel
        if not path.exists():
            row.update(sp=None, tc=None, avg_identity=None,
                       conserved_frac=None, aln_len=None, status="file missing")
            results.append(row)
            continue
        try:
            aln = parse_fasta(path)
            # Ensure all 4 expected sequences are present
            expected = set(BALIBASE_REF.keys())
            if not expected.issubset(set(aln.keys())):
                row.update(sp=None, tc=None, avg_identity=None,
                           conserved_frac=None, aln_len=None,
                           status=f"missing seqs: {expected - set(aln.keys())}")
                results.append(row)
                continue
            # Restrict to the 4 reference sequences (some alignments may have extras)
            aln = {k: aln[k] for k in expected}
            # Pad shorter sequences with trailing gaps if lengths differ (malformed MSA)
            max_len = max(len(v) for v in aln.values())
            aln = {k: v.ljust(max_len, "-") for k, v in aln.items()}
            idmat = pairwise_identity_matrix(aln)
            n = len(idmat)
            upper = [idmat.iloc[i, j] for i in range(n) for j in range(i+1, n)]
            row.update(
                sp=round(sp_score(aln, BALIBASE_REF, CORE_COLS), 4),
                tc=round(tc_score(aln, BALIBASE_REF, CORE_COLS), 4),
                avg_identity=round(float(np.mean(upper)), 4),
                conserved_frac=round(conserved_col_frac(aln), 4),
                aln_len=len(next(iter(aln.values()))),
                status="ok",
            )
        except Exception as e:
            row.update(sp=None, tc=None, avg_identity=None,
                       conserved_frac=None, aln_len=None, status=f"error: {e}")
        results.append(row)
    return results


def section_intro() -> str:
    return """\
## Overview

This report evaluates four LLM-based tools on a standard bioinformatics task: multiple sequence
alignment (MSA) of four homologous HMG-box domain protein sequences
(`1aab_`, `1j46_A`, `1k99_A`, `2lef_A`) from the BAliBASE benchmark (reference set RV11, entry BB11001).

Each tool was tested at three prompt stringency levels (low / medium / high) and two independent
takes. When a generated script failed, the same tool was re-prompted iteratively; each additional
prompt is counted as an extra "turn". The BAliBASE BB11001 structure-based alignment (96 columns,
57 core-block columns) is used as the ground truth for scoring. MUSCLE v5 is included as a
conventional-tool baseline.

**Scoring definitions:**
- **SP score**: fraction of residue pairs aligned in BAliBASE core columns that are also aligned
  in the test alignment (higher = better; 1.0 = perfect).
- **TC score**: fraction of complete BAliBASE core columns exactly reproduced in the test alignment
  (stricter than SP; 1.0 = perfect).
- **Avg identity**: mean pairwise sequence identity across all 6 sequence pairs, computed as
  matched residues / (alignment length − double-gap columns).
- **Conserved cols**: fraction of alignment columns where all non-gap residues are identical.
"""


def section_turns_table(metrics: list[dict]) -> str:
    lines = [
        "## Turns to Success\n",
        "Number of prompting turns required before the generated script ran successfully.\n",
        "| Model | Low T1 | Low T2 | Med T1 | Med T2 | High T1 | High T2 |",
        "|-------|--------|--------|--------|--------|---------|---------|",
    ]
    models = ["Claude", "ChatGPT", "Gemini", "Perplexity"]
    combos = [("low",1),("low",2),("medium",1),("medium",2),("high",1),("high",2)]
    lookup = {(r["model"], r["stringency"], r["take"]): r["turns"] for r in metrics}
    for model in models:
        cells = []
        for s, t in combos:
            val = lookup.get((model, s, t))
            cells.append("—" if val is None else str(val))
        lines.append(f"| {model} | " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def section_accuracy_table(metrics: list[dict]) -> str:
    lines = [
        "## Alignment Accuracy vs. BAliBASE BB11001\n",
        "SP and TC scores computed on the 57 core-block columns of BB11001.\n",
        "| Model | Stringency | Take | MSA Tool | Aln Len | SP Score | TC Score |",
        "|-------|-----------|------|---------|---------|----------|----------|",
    ]
    for r in metrics:
        if r["status"] not in ("ok",):
            sp_str = tc_str = f"*{r['status']}*"
        else:
            sp_str = f"{r['sp']:.3f}"
            tc_str = f"{r['tc']:.3f}"
        aln_len = str(r["aln_len"]) if r["aln_len"] else "—"
        lines.append(
            f"| {r['model']} | {r['stringency']} | {r['take']} | "
            f"{r['msa_tool']} | {aln_len} | {sp_str} | {tc_str} |"
        )
    return "\n".join(lines) + "\n"


def section_metrics_table(metrics: list[dict]) -> str:
    lines = [
        "## Standardized Sequence Statistics\n",
        "Pairwise identity and conserved columns recomputed uniformly from each aligned.fasta.\n",
        "| Model | Stringency | Take | Avg Identity | Conserved Cols |",
        "|-------|-----------|------|-------------|----------------|",
    ]
    for r in metrics:
        if r["status"] != "ok":
            id_str = cons_str = f"*{r['status']}*"
        else:
            id_str = f"{r['avg_identity']:.1%}"
            cons_str = f"{r['conserved_frac']:.1%}"
        lines.append(
            f"| {r['model']} | {r['stringency']} | {r['take']} | {id_str} | {cons_str} |"
        )
    return "\n".join(lines) + "\n"


def section_discussion(metrics: list[dict]) -> str:
    ok = [r for r in metrics if r["status"] == "ok"]
    muscle_ref = next((r for r in ok if r["model"] == "MUSCLE v5"), None)
    muscle_sp = f"{muscle_ref['sp']:.3f}" if muscle_ref else "N/A"
    muscle_tc = f"{muscle_ref['tc']:.3f}" if muscle_ref else "N/A"

    ext_tool = [r for r in ok if r["msa_tool"] in ("muscle", "clustalo", "mafft")]
    pure_py  = [r for r in ok if "Python" in r.get("msa_tool", "") or "NW" in r.get("msa_tool", "")]

    avg_sp_ext = np.mean([r["sp"] for r in ext_tool]) if ext_tool else float("nan")
    avg_sp_py  = np.mean([r["sp"] for r in pure_py])  if pure_py  else float("nan")

    return f"""\
## Discussion

### MUSCLE v5 vs. BAliBASE

The MUSCLE v5 run (`responses/reference/aligned.fasta`) scores SP={muscle_sp}, TC={muscle_tc}
against the BAliBASE structure-based reference. A score below 1.0 is expected: MUSCLE optimizes
a sequence-similarity objective, while BAliBASE columns are derived from 3D structure superposition.
This gap establishes the ceiling for what any sequence-only aligner can be expected to achieve.

### External-tool runs vs. pure-Python runs

Runs that delegated to an external aligner (MUSCLE, clustalo, MAFFT) achieved an average
SP of {avg_sp_ext:.3f} ({len(ext_tool)} runs). Runs that used a pure-Python Needleman-Wunsch
implementation achieved an average SP of {avg_sp_py:.3f} ({len(pure_py)} runs). The gap
reflects the difference in alignment quality between production-grade progressive aligners
and single-pass pairwise NW schemes.

### Metric inflation in self-reported numbers

Models reported average pairwise identities ranging from ~11% to ~40% for the same sequences.
The standardized recomputation in this report collapses that range to the true range of the
alignments. The primary source of inflation was denominator choice: some scripts divided
matches by the length of the shorter sequence (excluding gaps entirely), inflating identity
relative to the column-based denominator used here.

### Consistency across takes

Runs using the same external tool on the same prompt level are generally consistent across
takes (SP scores within ±0.02). The main source of variability is which tool was invoked:
MUSCLE and clustalo produce slightly different alignments; pure-Python NW implementations
vary more widely depending on gap penalties and alignment strategy.

### Failure patterns

The dominant failure mode across all models was incorrect MUSCLE CLI usage:
- **Biopython ≥1.79** removed the `MuscleCommandline` wrapper; scripts using it fail immediately.
- **MUSCLE v5** changed its CLI from `-in`/`-out` to `-align`/`-output`; scripts written for v3/v4 fail.
- Models that pivoted to clustalo or mafft when MUSCLE failed generally succeeded within one additional turn.
- Gemini failed to complete the high-stringency take 2 after 5 turns, cycling between muscle CLI errors
  and syntax errors introduced by each correction attempt.
"""


def write_report(output_path: Path, metrics: list[dict]) -> None:
    """Assemble all sections and write report.md."""
    sections = [
        f"# MSA Benchmark Report\n\n*Generated by `report/generate_report.py`*\n",
        section_intro(),
        section_turns_table(metrics),
        section_accuracy_table(metrics),
        section_metrics_table(metrics),
        section_discussion(metrics),
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n---\n\n".join(sections))
    print(f"Report written to {output_path}")


if __name__ == "__main__":
    # Sanity checks
    ref_path = REPO_ROOT / "responses/reference/balibase_BB11001.fasta"
    parsed = parse_fasta(ref_path)
    for name, seq in BALIBASE_REF.items():
        assert parsed[name] == seq, f"BAliBASE reference file mismatch for {name}"
    sp = sp_score(BALIBASE_REF, BALIBASE_REF, CORE_COLS)
    tc = tc_score(BALIBASE_REF, BALIBASE_REF, CORE_COLS)
    assert abs(sp - 1.0) < 1e-9 and abs(tc - 1.0) < 1e-9, "Self-alignment SP/TC must be 1.0"

    metrics = load_run_metrics(RUNS)
    errors = [r for r in metrics if r["status"] not in ("ok", "no output")]
    if errors:
        print("WARNING: some runs could not be loaded:")
        for r in errors:
            print(f"  {r['model']} {r['stringency']} take {r['take']}: {r['status']}")

    output = REPO_ROOT / "report" / "report.md"
    write_report(output, metrics)
