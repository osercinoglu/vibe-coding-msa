# Vibe Coding in Scientific Computing: A Bioinformatics Case Study on LLM Reliability and Prompt Sensitivity

**Authors:** Onur Serçinoğlu, Department of Bioengineering, Faculty of Engineering, Gebze Technical University

**Keywords:** vibe coding, large language models, bioinformatics, multiple sequence alignment, benchmark, prompt engineering, API hallucination, reproducibility

## Abstract

Vibe coding — the practice of delegating code authorship to large language models (LLMs) without fully comprehending the generated implementation — is spreading from general software development into scientific computing, where output correctness cannot be inferred from the absence of runtime errors. We present a controlled benchmark evaluating four LLM tools (Claude Haiku 4.5, GPT 4.0, Gemini Flash 2.5, and Perplexity) on the canonical bioinformatics task of multiple sequence alignment (MSA), scored against structure-derived BAliBASE reference alignments using Sum-of-Pairs (SP) and Total Column (TC) metrics. Across 21 model conditions spanning three prompt stringency levels and two independent takes, SP scores ranged from 0.088 to 1.000 for the same input sequences — a more than ten-fold variation for an identical task. The dominant failure mode was deprecated and incorrect external-tool CLI usage, particularly the MUSCLE version-5 interface change and the removal of the Biopython `MuscleCommandline` wrapper, reflecting a structural mismatch between LLM training data cutoffs and current library versions. A counter-intuitive "medium stringency paradox" emerged: adding specificity to prompts — specifically, restricting the available library set — systematically degraded alignment quality by steering models toward toy pure-Python Needleman-Wunsch implementations (average SP = 0.302) rather than production-grade aligners (average SP = 0.988). Models also self-reported alignment statistics inflated by two to four times relative to standardized recomputation, due to inconsistent denominator choices in identity calculations. Critically, all of these failures are invisible without a formal domain benchmark — generated scripts ran without errors yet produced scientifically invalid output. We conclude with concrete recommendations for prompt design and structured validation protocols when deploying LLM-generated code in scientific workflows.

---

## Türkçe Başlık

Bilimsel Hesaplamada Hissel Kodlama: BDM Güvenilirliği ve İstem Hassasiyeti Üzerine Bir Biyoenformatik Örnek Çalışması

## Öz

Hissel kodlama — üretilen uygulamayı tam olarak kavramadan kod yazımını büyük dil modellerine (BDM) bırakma pratiği — çalışma zamanı hatalarının yokluğundan çıktı doğruluğunun çıkarılamadığı bilimsel hesaplama alanına doğru genel yazılım geliştirmeden yayılmaktadır. Bu çalışmada, dört BDM aracını (Claude Haiku 4.5, GPT 4.0, Gemini Flash 2.5 ve Perplexity) klasik bir biyoenformatik analizi olan çoklu dizi hizalaması (ÇDH) üzerinde, Çift Toplamı (SP) ve Toplam Sütun (TC) metrikleri kullanılarak yapıya dayalı BAliBASE referans hizalamalarına göre puanlayan kontrollü bir kıyaslama sunulmuştur. Üç istem katılık düzeyi ve iki bağımsız denemeyi kapsayan 21 model koşulunda, aynı girdi dizileri için SP skorları 0,088 ile 1,000 arasında değişmiş; özdeş bir görev için on kattan fazla bir değişkenlik gözlemlenmiştir. Baskın hata modu; BDM eğitim verisi kesim tarihleri ile güncel kütüphane sürümleri arasındaki yapısal uyumsuzluğu yansıtan, özellikle MUSCLE sürüm 5 arayüz değişikliği ve Biopython'un `MuscleCommandline` sarmalayıcısının kaldırılmasıyla somutlaşan, kullanımdan kalkmış ve hatalı harici araç CLI kullanımı olmuştur. Bunun yanı sıra beklenmedik bir "orta katılık paradoksu" gözlemlenmiştir: kullanılabilir kütüphane kümesinin kısıtlanması yoluyla istemlere özgüllük eklenmesi, modelleri üretim hizaleyicilerden (ortalama SP = 0,988) uzaklaştırarak saf Python Needleman-Wunsch gerçekleştirmelerine (ortalama SP = 0,302) yönlendirmiş ve hizalama kalitesini sistematik olarak düşürmüştür. Modeller ayrıca kimlik hesaplamalarında tutarsız payda seçimi nedeniyle standartlaştırılmış yeniden hesaplamaya kıyasla iki ila dört kat şişirilmiş hizalama istatistikleri raporlamıştır. Bu hataların tamamı, resmi bir alan kıyaslaması olmaksızın tespit edilmesi güçtür — üretilen betikler hata vermeden çalışmış, ancak bilimsel açıdan geçersiz çıktılar üretmiştir. Sonuç olarak, bilimsel iş akışlarında BDM tarafından üretilen kodun kullanımında istem tasarımı ve yapılandırılmış doğrulama protokollerine yönelik somut öneriler sunulmaktadır.

---

## 1. Introduction

The term "vibe coding," introduced by Andrei Karpathy in early 2025 [1], describes a mode of software development in which a programmer intentionally delegates code writing to an LLM and accepts the output without fully reading or understanding it — relying instead on whether the program "runs" as a correctness signal. Recent grey literature and qualitative studies document rapid adoption of this practice [2, 3]. For application developers building web interfaces or data pipelines, the cost of undetected errors is bounded by testing infrastructure and user feedback loops. For domain scientists, the calculus is different: correctness is not optional, and the absence of a crash does not imply the absence of a scientific error.

Bioinformatics presents an unusually well-instrumented stress test for vibe coding. LLMs are increasingly applied to bioinformatics code generation [8, 9], and multiple sequence alignment (MSA) is among the field's most canonical computational tasks. MSA is well-defined, widely implemented in mature tools, and — crucially — the BAliBASE benchmark [16, 17] provides a structure-based objective ground truth derived from three-dimensional protein structure superpositions. Unlike web or CRUD application code, where correctness is largely subjective or defined by unit tests, MSA correctness is externally measurable. This makes silent failures — code that executes without error but produces scientifically wrong output — directly detectable.

Three compounding risks motivate this investigation. First, non-determinism: even with deterministic generation settings, repeated LLM queries for the same task yield variable outputs [5], with non-determinism rates reaching 75% on some code generation benchmarks. Second, prompt sensitivity: minor phrasing changes cause accuracy swings exceeding 45% [6, 7], meaning the quality of generated code depends heavily on how the request is worded. Third, API hallucination and deprecation: comprehensive taxonomies of LLM code hallucinations [12] explicitly categorize deprecated API calls as a failure subcategory, and empirical studies find a 62% API misuse rate in GPT-4-generated code snippets [13]. The "Is Vibe Coding Safe?" study [4] addresses the security dimension of agent-generated code; we address scientific correctness — a complementary and underexplored dimension of the same problem.

Our contributions are: (1) a reproducible 24-condition benchmark of 4 LLM tools across 3 prompt stringency levels and 2 independent takes, scored against BAliBASE entry BB11001 using standardized SP and TC metrics; (2) a taxonomy of the observed failure modes instantiated concretely in bioinformatics tooling; (3) empirical demonstration of a "medium stringency paradox" in which adding precision to prompts systematically worsened scientific output; and (4) evidence of systematic metric inflation in self-reported alignment statistics, showing that plausible-sounding numbers can mask severe correctness deficits.

## 2. Related Work

The concept of vibe coding [1] describes intentional, rather than accidental, delegation of code authorship to LLMs, with the programmer acting as an orchestrator rather than a reviewer. Adoption surveys and qualitative interviews document that this practice is spreading rapidly across experience levels and domains [2, 3]. The trust dynamics are paradoxical: one recent critical review finds that 96% of developers report distrusting AI-generated code, yet adoption continues to accelerate [15]. Non-determinism compounds the trust problem: Ouyang et al. [5] empirically document non-determinism rates of 75.76%, 51.00%, and 47.56% across three widely-used code generation benchmark datasets, and Paleyes et al. [6] demonstrate accuracy swings exceeding 45% from minor prompt phrasing changes. Zi et al. [7] provide complementary evidence that prompt specificity affects code generation quality in ways that are not monotone — more detailed prompts do not reliably produce better code.

API hallucination is a structurally distinct failure mode from logical bugs or style issues. Zhang et al. [12] provide a comprehensive taxonomy of LLM code hallucinations that explicitly categorizes deprecated API calls as a recognized subcategory ("API Knowledge Conflict"). Zhong and Wang [13] find that 62% of GPT-4-generated code snippets contain API misuse when evaluated against current library documentation. Zhuo et al. [14] propose mitigation strategies including retrieval-augmented generation with up-to-date API documentation. Our study instantiates these abstract categories concretely: the Biopython `MuscleCommandline` wrapper, removed in version 1.79, and the MUSCLE v3/v4 argument flags (`-in`/`-out`), replaced by `-align`/`-output` in MUSCLE v5, appear as recurring failure patterns across multiple models and prompting turns.

The application of LLMs to bioinformatics code generation has received growing attention. BioCoder [8] is the most systematic benchmark to date, evaluating models on 2522 bioinformatics functions; GPT-4 achieves approximately 60% pass@k, indicating substantial but incomplete capability. Luo et al. [9] review one year of ChatGPT usage in bioinformatics and document persistent inconsistencies in code correctness. Rahman and Wong [10] critically evaluate ChatGPT's limitations for computational biology programming tasks, identifying brittleness under domain-specific constraints. Yin et al. [11] provide a broader evaluation of LLMs across bioinformatics research tasks. Our study differs from all of these by using formal SP and TC scoring against a structure-derived reference alignment — rather than pass/fail unit tests or human evaluation — allowing quantitative, reproducible comparison of alignment correctness across models and conditions.

## 3. Materials and Methods

### 3.1 Task and Sequences

The benchmark task was protein multiple sequence alignment of four HMG-box domain sequences (`1aab_`, `1j46_A`, `1k99_A`, `2lef_A`) drawn from the BAliBASE reference set RV11, entry BB11001 [16, 17]. The reference alignment spans 96 columns, of which 57 are designated core-block columns derived from three-dimensional protein structure superpositions. These core-block columns define the ground truth against which all model-generated alignments were scored, independently of any algorithmic objective function. BAliBASE is a standard in the MSA benchmarking literature precisely because structure-based reference alignments are not biased toward any particular sequence-alignment algorithm.

### 3.2 LLM Tools and Prompt Levels

Four tools were evaluated: Claude Haiku 4.5, ChatGPT, Gemini Flash 2.5, and Perplexity, each accessed via free-tier interfaces in January 2026. Three prompt stringency levels were applied (full text available in the public repository):

- **Low**: Informal, outcome-focused; no constraints on approach or tooling. The prompt requested "a reasonable multiple sequence alignment and some basic summary numbers," explicitly leaving approach open ("you can write code or describe the steps, whatever you think is best").
- **Medium**: Structured; requires specific metrics (pairwise identity matrix, average identity, and fraction of fully conserved columns), a `--input`/`--output` CLI interface, and Python 3.11. Critically, the prompt restricts dependencies to numpy, pandas, matplotlib, and seaborn — effectively prohibiting subprocess calls to external aligners without explicitly stating this prohibition.
- **High**: Production-style; requires functions with docstrings, a `main()` entry point, argparse with `--input`/`--outdir`, input validation for at least three sequences, graceful error handling, and CSV/text reports in a `results/` subdirectory. The key difference from medium: this prompt explicitly states that external MSA tools such as MUSCLE or Clustal Omega may be called via subprocess.

Two independent takes were run per (model, stringency) pair. When a generated script failed at runtime, the model was re-prompted with the error message iteratively; each additional exchange counts as one additional turn. Claude Haiku 4.5 was evaluated on take 1 only, and it ran the analysis directly within the claude.ai interface rather than generating a standalone script. The total evaluation set comprised 21 model conditions plus a MUSCLE v5 standalone baseline.

### 3.3 Evaluation Metrics

All metrics were computed uniformly from each `aligned.fasta` output by `report/generate_report.py` (available in the repository).

- **SP score** (Sum-of-Pairs): the fraction of residue pairs that are co-aligned in BAliBASE core-block columns and are reproduced as co-aligned in the test alignment [16]. The score ranges from 0 to 1, where 1.0 indicates exact recovery of all reference co-alignments.
- **TC score** (Total Column): the fraction of core-block columns that are reproduced exactly in the test alignment, with all residues in the correct positions and no spurious gap insertions. TC is strictly harder to satisfy than SP.
- **Turns to success**: the number of prompting exchanges required before an `aligned.fasta` file was successfully produced.
- **Standardized pairwise identity**: matches divided by (alignment length minus double-gap columns). This denominator was applied uniformly to avoid the inconsistencies in self-reported metrics, which used varying denominators across models.

### 3.4 Baseline

MUSCLE v5 [18] was run directly on the input sequences using its native command-line interface (`-align`/`-output` flags). It achieves SP=TC=1.000, exactly recovering the BAliBASE structure-based reference alignment at all 57 core-block columns. This confirms that the reference is achievable by a standard sequence aligner on this benchmark entry and establishes an empirical ceiling against which model-generated alignments are compared.

## 4. Results

### 4.1 Alignment Accuracy

SP scores span the full range from 0.088 (Gemini medium take 1, pure-Python Needleman-Wunsch) to 1.000 (ChatGPT all conditions, Perplexity low take 1, Perplexity high take 1, Perplexity high take 2). TC scores for Claude low and medium conditions, and for Gemini medium conditions, reached 0.000 — meaning no single core-block column was exactly reproduced. Table 1 presents the complete results across all 21 model conditions plus the baseline.

ChatGPT is the only model to achieve SP=TC=1.000 across all six tested conditions. Notably, it consistently delegated alignment to MUSCLE regardless of prompt stringency level, including under the medium prompt where other models interpreted the library constraint as a prohibition on external tools. Four runs that used Clustal Omega share a distinctive signature of SP=0.956/TC=0.912, appearing identically in Claude high take 1, Gemini high take 1, Gemini low take 2, and Perplexity low take 2 — reflecting Clustal Omega's slight, consistent divergence from the BAliBASE reference on this particular dataset.

Claude's two lowest-accuracy conditions (low take 1: SP=0.111; medium take 1: SP=0.152) were each produced in a single prompting turn — rapid convergence to an incorrect result. This illustrates that turn count is not a proxy for output quality in either direction. Gemini high take 2 is the only condition in the benchmark that produced no output at all, despite six prompting turns.

**Table 1: Alignment accuracy vs. BAliBASE BB11001 core-block columns.**

| Model | Stringency | Take | MSA Tool | Aln Len | SP Score | TC Score |
|-------|-----------|------|---------|---------|----------|----------|
| Claude | low | 1 | inline (claude.ai) | 91 | 0.111 | 0.000 |
| Claude | medium | 1 | inline + NW | 96 | 0.152 | 0.000 |
| Claude | high | 1 | clustalo | 96 | 0.956 | 0.912 |
| ChatGPT | low | 1 | muscle | 96 | 1.000 | 1.000 |
| ChatGPT | medium | 1 | muscle | 96 | 1.000 | 1.000 |
| ChatGPT | high | 1 | muscle | 96 | 1.000 | 1.000 |
| ChatGPT | low | 2 | muscle | 96 | 1.000 | 1.000 |
| ChatGPT | medium | 2 | muscle | 96 | 1.000 | 1.000 |
| ChatGPT | high | 2 | muscle | 96 | 1.000 | 1.000 |
| Gemini | low | 1 | mafft | 99 | 0.991 | 0.983 |
| Gemini | medium | 1 | NW (pure Python) | 97 | 0.088 | 0.000 |
| Gemini | high | 1 | clustalo | 96 | 0.956 | 0.912 |
| Gemini | low | 2 | clustalo | 96 | 0.956 | 0.912 |
| Gemini | medium | 2 | NW (pure Python) | 97 | 0.564 | 0.386 |
| Gemini | high | 2 | — | — | *no output* | *no output* |
| Perplexity | low | 1 | muscle | 96 | 1.000 | 1.000 |
| Perplexity | medium | 1 | NW (pure Python) | 115 | 0.412 | 0.228 |
| Perplexity | high | 1 | muscle | 96 | 1.000 | 1.000 |
| Perplexity | low | 2 | clustalo | 96 | 0.956 | 0.912 |
| Perplexity | medium | 2 | NW (pure Python) | 118 | 0.295 | 0.123 |
| Perplexity | high | 2 | muscle | 96 | 1.000 | 1.000 |
| MUSCLE v5 | — | — | muscle | 96 | 1.000 | 1.000 |

### 4.2 The Medium Stringency Paradox

A counter-intuitive pattern emerges from the medium stringency results: across all models, the medium prompt consistently produced worse alignments than either the low or high prompt. The five conditions that used pure-Python Needleman-Wunsch (NW) implementations — all arising from the medium prompt — achieved an average SP of 0.302. The fifteen conditions using external tools (MUSCLE, Clustal Omega, MAFFT) averaged SP=0.988.

The cause is the medium prompt's library constraint: by specifying that only numpy, pandas, matplotlib, and seaborn are guaranteed to be available, without mentioning subprocess or external tool calls, models interpreted this as a prohibition on calling production aligners. They instead implemented NW from scratch in pure Python. This interpretation is linguistically defensible — "only standard scientific libraries are guaranteed" does read as a closed list — but it had severe consequences for scientific correctness.

These pure-Python implementations compound inherent NW algorithmic limitations — including linear gap penalties that scatter indels, no guide tree, and no iterative refinement [20] — with implementation-level bugs. The worst-performing run (Gemini medium take 1, SP=0.088) aligned each new sequence against only the first sequence in the growing MSA rather than against a proper sequence profile, discarding all positional information accumulated from earlier alignment steps. The result was an alignment structurally inconsistent with the reference at nearly every core-block column.

The high-stringency prompt, by contrast, explicitly states that external MSA tools "may be called via subprocess" — and all high-stringency runs that completed achieved SP ≥ 0.956. The medium prompt, by adding specificity without specifying escape routes, actively degraded scientific output relative to the informal low prompt.

### 4.3 Turns to Success and Failure Patterns

**Table 2: Turns required before the generated script ran successfully. Dash (—) denotes either no take 2 data (Claude) or failure to produce output (Gemini high take 2).**

| Model | Low T1 | Low T2 | Med T1 | Med T2 | High T1 | High T2 |
|-------|--------|--------|--------|--------|---------|---------|
| Claude | 1 | — | 1 | — | 1 | — |
| ChatGPT | 2 | 1 | 1 | 1 | 1 | 2 |
| Gemini | 3 | 1 | 1 | 3 | 2 | — |
| Perplexity | 2 | 1 | 1 | 1 | 1 | 3 |

The dominant recurring failure mode was incorrect MUSCLE CLI usage, taking two distinct forms: (a) invocation of the `MuscleCommandline` wrapper from Biopython, which was removed in Biopython ≥1.79, causing an `AttributeError` or `ImportError` at runtime; and (b) use of MUSCLE v3/v4 argument flags (`-in`/`-out`) that were replaced by `-align`/`-output` in MUSCLE v5 [18], causing a command-line parsing error. These errors affected ChatGPT low take 1 (one additional turn to fix), Gemini low take 1 (which additionally attempted MAFFT via a missing PATH entry before pivoting to a working solution on turn 3), and Perplexity high take 2 (resolved on turn 3). Models that pivoted to Clustal Omega when MUSCLE failed generally succeeded within one additional turn, as Clustal Omega's CLI interface has remained stable.

Gemini high take 2 represents the most severe failure case in the benchmark: across six prompting turns, the model cycled between MUSCLE v5 CLI errors and Python syntax errors introduced by each successive correction attempt. This illustrates how iterative re-prompting can become trapped in a local failure loop, where each correction introduces a new class of error, and the model lacks sufficient context to escape the cycle.

### 4.4 Metric Inflation

Models self-reported average pairwise identities ranging from approximately 11% to 40% for the same input sequences. Standardized recomputation collapses this range to 9.8–22.7%. The primary source of inflation is denominator choice: several scripts divided matched residues by the length of the shorter sequence or by the number of non-gap positions in one sequence, rather than by the alignment length minus double-gap columns — the standard definition used in structural bioinformatics.

The practical consequence for vibe coding is significant. Perplexity medium take 1 self-reports a conserved-columns fraction of 19.1% — a plausible-sounding figure for a set of moderately conserved homologs — yet achieves SP=0.412 against the structural reference, meaning fewer than half of the reference co-alignments are recovered. Without a formal benchmark, this inflation is undetectable: the script runs, the output file is populated, the reported numbers fall within a reasonable range, and there is no signal to prompt verification.

## 5. Discussion

### 5.1 The "It Runs" Illusion

Of 21 non-baseline conditions, 13 produced a syntactically valid `aligned.fasta` without runtime error. Of those 13 apparent successes, five relied on pure-Python Needleman-Wunsch implementations with an average SP of 0.302 — incorrect by any bioinformatics standard, yet indistinguishable in output format from a correct alignment. A scientist who accepted these results at face value — trusting that a script that produced output without crashing must be working correctly — would report conclusions based on alignments that failed to recover the majority of structure-validated residue co-alignments.

This is a concrete instantiation of the trust dynamic identified in studies of vibe coding adoption [3]: runtime success serves as an implicit proxy for correctness, inflating confidence in results that have not been independently validated. Bioinformatics is unusual in possessing formal reference benchmarks like BAliBASE [16, 17] that can expose this class of silent failure. Most scientific computing domains — genomic data processing, statistical modeling, image analysis pipelines — lack equivalent ground truth. The implication is that the invisible-failure problem documented here is likely far more prevalent in scientific vibe coding practice than is currently appreciated, precisely because the infrastructure to detect it is absent in most domains.

### 5.2 API Deprecation as a Structural Risk

The MUSCLE CLI and Biopython failures observed in this benchmark are not random bugs. They are a predictable structural consequence of the mismatch between LLMs' training data cutoffs and the current state of actively-maintained scientific software libraries. Biopython removed the `MuscleCommandline` wrapper in version 1.79; MUSCLE changed its entire command-line argument interface between version 4 and version 5 [18]. LLMs trained on historical code corpora confidently generate usage patterns that were correct at training time but no longer function against current library releases.

This pattern is the concrete bioinformatics instantiation of the "API Knowledge Conflict" hallucination subcategory identified in the taxonomy of Zhang et al. [12] and the 62% API misuse rate documented empirically by Zhong and Wang [13]. High-stringency prompting provides partial mitigation: when MUSCLE fails under a high-stringency prompt, models with explicit subprocess permission pivot to Clustal Omega, which currently maintains a stable interface. But this resilience is opportunistic — it depends on having a working alternative that happens to be available — rather than structurally robust. Mitigation strategies proposed by Zhuo et al. [14], including retrieval-augmented generation with current API documentation, address the root cause more directly, but require deliberate infrastructure investment that is absent from typical scientific workflows. Pinned dependency environments and CI-style validation against known reference outputs are the practical near-term mitigations for production scientific use.

### 5.3 Prompt Stringency is Not Monotone

The medium stringency paradox is this study's most practically actionable finding. Adding explicit structure and precision to a prompt — specifying required output format, metric names, CLI interface, and library environment — can simultaneously remove the model's access to a correct solution if the new constraints inadvertently prohibit the tools needed to produce one.

This finding resonates with and extends the empirical results on prompt sensitivity [6, 7]: the relationship between prompt specificity and output quality is not monotone, and depends critically on which constraints are imposed and which are left implicit. For scientific computing tasks specifically, the recommendation follows directly from the observed failure mode: prompts should explicitly enumerate which external tools are permitted and state their expected version interfaces. Leaving library and tooling assumptions implicit — specifying some constraints (allowed libraries) without others (allowed subprocess calls) — is more dangerous than either leaving everything implicit (low prompt) or being explicit about everything (high prompt). The intermediate state creates interpretation ambiguities that models resolve through plausible-sounding but scientifically incorrect implementations.

The corollary for prompt design practice is concrete: a medium-stringency scientific prompt should include a clause analogous to the high-stringency prompt's explicit subprocess permission, e.g., "If a compiled external tool would produce better results, you may call it via subprocess." The absence of such a clause in a constrained prompt should be treated as a known risk factor for pure-Python fallback implementations.

### 5.4 Limitations

This study has several limitations that qualify its conclusions. The benchmark covers a single BAliBASE entry (BB11001) with four well-studied HMG-box sequences; results may not generalize to larger, more divergent protein families or to other classes of bioinformatics tasks. Only two independent takes per condition were collected, which is insufficient for robust statistical characterization of non-determinism rates (cf. Ouyang et al. [5], which uses substantially larger samples). Claude Haiku 4.5 was tested on take 1 only. The decision to stop re-prompting after a threshold number of turns introduces a subjective stopping criterion that affects the turns-to-success measure for Gemini high take 2. Finally, all models were accessed via free-tier interfaces during a specific two-week window; behavior may differ under paid-tier settings, with temperature parameters, or following model updates deployed after the evaluation period.

## 6. Conclusion

We presented a controlled benchmark of four LLM tools on protein multiple sequence alignment, scored uniformly against the BAliBASE structure-based reference. SP scores for identical input sequences ranged from 0.088 to 1.000 across 21 model conditions. Prompt stringency exhibited a non-monotone relationship with output quality: medium-stringency prompts systematically produced the worst results by steering models toward pure-Python Needleman-Wunsch implementations with an average SP of 0.302, while the informal low and explicit high prompts both yielded substantially better outcomes. API deprecation — specifically the MUSCLE v5 CLI interface change and the removal of the Biopython `MuscleCommandline` wrapper — was the dominant runtime failure mode, recurring across multiple models and prompting turns in a pattern consistent with the API Knowledge Conflict hallucination category. Self-reported alignment statistics were systematically inflated relative to standardized recomputation, producing plausible-sounding numbers that masked severe correctness deficits. None of these failure modes are visible from script execution alone.

For scientists adopting LLM-assisted coding in research workflows, we recommend: (1) include explicit version constraints and permitted external tools in prompts, naming expected CLI interfaces; (2) explicitly authorize subprocess calls to domain tools in any prompt that specifies library constraints; (3) benchmark generated code against known reference datasets before use in analysis — for bioinformatics, BAliBASE entries with known SP/TC scores provide a ready-made validation layer; and (4) prefer high-stringency prompts that explicitly invite appropriate tooling over medium-stringency prompts that inadvertently prohibit it. All experimental code, generated scripts, alignment outputs, and evaluation reports are available in the public repository. Future work should expand this evaluation across larger BAliBASE families with greater sequence divergence, additional LLM models and API tiers, and other standard bioinformatics tasks such as variant calling, differential expression analysis, and protein structure prediction pipelines.

## References

[1] A. Karpathy. "Vibe coding." *X (formerly Twitter)*, 2025. https://x.com/karpathy/status/1886192184808149289

[2] A. Fawzy, A. Tahir, and K. Blincoe. "Vibe Coding in Practice: Motivations, Challenges, and a Future Outlook — a Grey Literature Review." *arXiv*, 2025. https://arxiv.org/abs/2510.00328

[3] V. Pimenova, S. Fakhoury, C. Bird, M.-A. Storey, and M. Endres. "Good Vibrations? A Qualitative Study of Co-Creation, Communication, Flow, and Trust in Vibe Coding." *arXiv*, 2025. https://arxiv.org/abs/2509.12491

[4] S. Zhao et al. "Is Vibe Coding Safe? Benchmarking Vulnerability of Agent-Generated Code in Real-World Tasks." *arXiv*, 2025. https://arxiv.org/abs/2512.03262

[5] Ouyang et al. "An Empirical Study on the Non-Determinism of ChatGPT in Code Generation." *ACM Transactions on Software Engineering and Methodology*, 2024. https://doi.org/10.1145/3697010

[6] A. Paleyes, R. Sendyka, D. Robinson, C. Cabrera, and N.D. Lawrence. "Code Roulette: How Prompt Variability Affects LLM Code Generation." *arXiv*, 2025. https://arxiv.org/abs/2506.10204

[7] Y. Zi, H. Menon, and A. Guha. "More Than a Score: Probing the Impact of Prompt Specificity on LLM Code Generation." *arXiv*, 2025. https://arxiv.org/abs/2508.03678

[8] H. Dong et al. "BioCoder: A Benchmark for Bioinformatics Code Generation with Contextual Pragmatic Knowledge." *Bioinformatics*, 2024. https://doi.org/10.1093/bioinformatics/btae209

[9] Luo et al. "One Year of ChatGPT in Bioinformatics: Strengths, Inconsistencies, and Open Challenges." *PMC*, 2024. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11446534/

[10] C.R. Rahman and L. Wong. "How much can ChatGPT really help Computational Biologists in Programming?" *arXiv*, 2023. https://arxiv.org/abs/2309.09126

[11] H. Yin et al. "An Evaluation of Large Language Models in Bioinformatics Research." *arXiv*, 2024. https://arxiv.org/abs/2402.13714

[12] Z. Zhang, Y. Wang, C. Wang, J. Chen, and Z. Zheng. "LLM Hallucinations in Practical Code Generation: Phenomena, Mechanism, and Mitigation." *arXiv*, 2024. https://arxiv.org/abs/2409.20550

[13] L. Zhong and Z. Wang. "Can ChatGPT replace StackOverflow? A Study on Robustness and Reliability of Large Language Model Code Generation." *arXiv*, 2023. https://arxiv.org/abs/2308.10335

[14] T.Y. Zhuo et al. "Identifying and Mitigating API Misuse in Large Language Models." *arXiv*, 2025. https://arxiv.org/abs/2503.22821

[15] S. Baltes, T. Speith, B. Chiteri, S. Mohsenimofidi, S. Chakraborty, and D. Buschek. "On the Need to Rethink Trust in AI Assistants for Software Development: A Critical Review." *arXiv*, 2025. https://arxiv.org/abs/2504.12461

[16] J.D. Thompson, F. Plewniak, and O. Poch. "BAliBASE: A Benchmark Alignment Database for the Evaluation of Multiple Alignment Programs." *Bioinformatics*, 1999. https://doi.org/10.1093/bioinformatics/15.1.87

[17] J.D. Thompson, P. Koehl, R. Ripp, and O. Poch. "BAliBASE 3.0: Latest Developments in the Benchmark Alignment Database for the Evaluation and Comparison of Multiple Sequence Alignment Programs." *Proteins*, 2005. https://doi.org/10.1002/prot.20527

[18] R.C. Edgar. "Muscle5: High-Accuracy Alignment Ensembles Enable Unbiased Assessments of Sequence Homology and Phylogeny." *Nature Methods*, 2022. https://doi.org/10.1038/s41592-022-01570-8

[19] K. Katoh, K. Misawa, K. Kuma, and T. Miyata. "MAFFT: A Novel Method for Rapid Multiple Sequence Alignment Based on Fast Fourier Transform." *Nucleic Acids Research*, 2002. https://doi.org/10.1093/nar/gkf436

[20] S.B. Needleman and C.D. Wunsch. "A General Method Applicable to the Search for Similarities in the Amino Acid Sequence of Two Proteins." *Journal of Molecular Biology*, 1970. https://doi.org/10.1016/0022-2836(70)90057-4
