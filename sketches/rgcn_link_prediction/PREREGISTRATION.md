# Preregistration: SBGN-PD typed relations as inductive bias, a preregistered comparison

**Author.** Noor Muhammad ([github.com/NoorMuhammad1](https://github.com/NoorMuhammad1))
**Date registered.** 2026-07-05
**Status.** Registered, not yet executed
**Companion artifact.** [sketches/rgcn_link_prediction/](.) in this repo

---

## 1. Motivation

The [R-GCN link-prediction sketch](rgcn_sketch.py) in this directory demonstrates that the SBGN-PD BioCypher adapter's output feeds a relation-aware graph neural network end-to-end without special preprocessing. On the default glycolysis fixtures the training pipeline runs to convergence, but the resulting metrics are meaningless because the graph is 14 nodes with 12 arcs. The sketch is the pipeline. This preregistration extends the pipeline to a production-scale corpus and commits, before any experiment is run, to a protocol that would produce a defensible answer to the question the cover letter to Prof. Rousu poses. Pilot executions of the same pipeline at reduced synthetic scale are reported separately at [`PILOT_RESULTS.md`](PILOT_RESULTS.md) (v1) and [`PILOT_V2_RESULTS.md`](PILOT_V2_RESULTS.md) (v2).

That question in one sentence. How much of the link-prediction signal on curated pathway graphs is carried by the SBGN-PD relation typing (six typed arc classes) versus a single-relation collapse where all arcs share one type?

## 2. Research question and hypothesis

**Research question (RQ).** On a Reactome-scale SBGN-PD pathway graph emitted by the BioCypher adapter, does a relation-aware encoder that preserves the six SBGN-PD arc types outperform an otherwise-identical encoder that collapses all arcs to a single relation, on held-out link prediction?

**Primary hypothesis (H1).** The typed encoder T achieves higher raw MRR than the flat encoder F on the held-out test split, with per-seed difference detectable at alpha = 0.05 by a paired test across ten random seeds.

**Null hypothesis (H0).** T and F achieve indistinguishable raw MRR on the held-out test split under the same test.

**Falsification.** H1 is falsified if the paired difference on the pre-committed test statistic fails to reach significance, or if F outperforms T (see Section 10 for interpretation).

## 3. Data

**Source.** Reactome SBGN-ML pathway files. The exact Reactome release version, download date, and file list will be recorded in `data/manifest.json` before any experimental run and committed to this branch. No changes to the file list after that point.

**Pipeline.** The SBGN-PD BioCypher adapter at the tagged version `v0.2.0` (or later, tag pinned before run start) with default entity-matching threshold 0.7. All flags recorded in `pipeline_config.json`.

**Expected scale.** Order of 10^4 to 10^5 nodes and 10^5 to 10^6 arcs. If the emitted graph is smaller than 10^3 arcs, the preregistration is suspended and the deviation logged in Section 12.

**Reproducibility.** Adapter tag, Reactome dump timestamp, entity-matching threshold, and any parser flags are committed to this preregistration branch before any experimental run. Any deviation triggers Section 12.

## 4. Models

Two encoder architectures, identical decoder, matched hyperparameters, matched parameter counts where possible.

**Typed encoder (T).**
- Two-layer R-GCN with basis decomposition (Schlichtkrull et al. 2018)
- `num_relations = 6` (SBGN-PD arc classes preserved as distinct relations)
- `num_bases = 3` (basis compression rationale in Section 4.1)
- Hidden dimension 128
- Message normalization by per-target in-degree per relation
- Non-linearity ReLU
- Learned self-loop weight matrix

**Flat encoder (F).**
- Two-layer R-GCN with `num_relations = 1` (all SBGN-PD arcs collapsed to a single `interacts_with` relation while preserving edge direction to keep parity with T's basis-decomposed handling)
- `num_bases = 1`
- Hidden dimension 128
- All other hyperparameters identical to T

**Decoder (shared).**
- DistMult (Yang et al. 2015). Score `f(s, r, o) = <h_s, r, h_o>`
- Six relation vectors for T. One relation vector for F
- Same node embedding table (learned, dimension 128) shared between the encoder input and DistMult heads

**Loss.** Binary cross-entropy on positive triples versus 5 uniformly corrupted negatives per positive, corruption at the head-or-tail position with equal probability.

### 4.1 Parameter-count fairness

T has more parameters than F by construction because T carries three basis matrices and six relation vectors, while F carries one basis matrix and one relation vector. To rule out the possibility that any advantage of T comes from capacity rather than from relational typing, the comparison is run under two matched conditions.

**Condition A (matched-per-relation).** T uses hidden dim 128, num_bases 3. F uses hidden dim 128, num_bases 1. This is the primary comparison.

**Condition B (parameter-matched).** F's hidden dimension is increased until F's total parameter count matches T's within 5%. The resulting F-large is then compared against T under Condition B.

Both conditions are reported. H1 is considered supported only if T outperforms F under Condition A AND F-large under Condition B.

## 5. Training protocol

- Optimizer Adam, `lr = 1e-3`, `weight_decay = 1e-5`
- Full-graph batching (single forward per epoch). If the graph exceeds GPU memory, switch to neighbor-sampled mini-batches with `fanout = [20, 10]` per layer, and record the switch as a deviation
- Epochs 500 with early stopping on validation MRR (patience 30, minimum improvement 0.005)
- Seeds `0, 1, 2, 3, 4, 5, 6, 7, 8, 9`

## 6. Data splits

- Arc-level 80 / 10 / 10 split into train / validation / test
- Split is produced with fixed shuffle seed `S = 42` (pre-committed here) once before any training run, then frozen. Model-training seeds `0` through `9` are independent of `S`
- Stratified by relation type so each split contains the same proportion of each of the six SBGN-PD arc classes. Minimum floor. Each split must contain at least 30 arcs of each of the six relation classes. If any class falls below this floor after stratification, that class is excluded from the primary MRR computation and flagged in `deviations.md`
- Because the split is fixed across model seeds, the ten-seed variance reflects initialization and negative-sampling variance only, not test-set sampling variance. This is a deliberate choice to increase paired-test power. Generalisation to a different test split is not claimed
- All node IDs are present in every split (transductive setting)
- Split manifest committed to the preregistration branch as `data/split_manifest.json`

## 7. Evaluation

**Primary metric.** Raw MRR on the test split. For each true triple `(s, r, o)`, rank the true object `o` against every node in the graph under the DistMult score. Reciprocal rank averaged across the test set.

**Secondary metrics.** Hits@1, Hits@3, Hits@10 on the same rankings.

**Sensitivity check.** Filtered ranking (excluding other known-true triples with the same `(s, r, ?)` head from the ranking) reported alongside raw ranking. Filtered ranking is a sensitivity check only. The primary decision statistic is the raw MRR.

**Randomness control.** Ten seeds. All ten reported. No cherry-picking.

## 8. Analysis plan

- For each of T, F, and F-large, compute mean MRR ± 1.96 × standard error across the ten seeds
- Paired difference test on per-seed MRR difference `d_i = MRR_T(seed_i) - MRR_F(seed_i)`
  - Primary test. Exact one-sided Wilcoxon signed-rank on the ten paired differences, alternative `median(d_i) > 0`. Exact null distribution rather than asymptotic (n = 10 is too small for the normal approximation). Zero differences handled by the Pratt method (retained with sign assigned by rank averaging). The complementary alternative `median(d_i) < 0` is reported as a diagnostic for Case IV interpretation but does not enter the primary decision
  - Sensitivity test. Paired t-test if the differences are approximately normal by Shapiro-Wilk at alpha = 0.10
- Effect size. Cohen's `d_z = mean(d_i) / sd(d_i)` with Bessel-corrected sample SD (n minus one denominator) on the ten paired differences. `d_av = mean(d_i) / sqrt((var(MRR_T) + var(MRR_F)) / 2)` is reported as a sensitivity number. Bootstrap BCa confidence interval on `d_z` with 10,000 resamples
- Pre-committed alpha for H1 versus H0 decision. 0.05
- Sample-size justification. n = 10 seeds. At alpha = 0.05 exact one-sided Wilcoxon, power is approximately 0.75 to detect Cohen's `d_z` of 1.0 and approximately 0.60 to detect `d_z` of 0.8. If the true effect is smaller than `d_z` approximately 0.5, the study is underpowered and a non-rejection under Case III should be read as "no evidence of medium-or-larger effect" rather than "no effect." Compute-budget permitting, the study may be extended to n = 20 seeds. The decision to extend is made once, before any test-set MRR is inspected, based only on the observed variance across the first ten seeds' validation MRRs
- No multiple-comparison correction across primary and sensitivity tests. The primary test is Wilcoxon and is the sole decision statistic

## 9. Deliverables

- Trained model checkpoints for all runs at the best-validation-MRR epoch, committed as GitHub release assets
- `results.json` with per-seed metrics for T, F, F-large under both conditions
- `results.pdf` with the analysis figures (training-loss curves per seed, MRR bar chart with error bars, effect-size annotation)
- One-paragraph plain-language summary of the outcome, whether the primary test rejected H0

## 10. Interpretation of outcomes

Every possible outcome has a pre-committed reading.

**Case I. T > F significantly under both conditions.** Relation typing carries link-prediction signal beyond topology alone, and beyond additional capacity. The typed-arcs-carry-signal claim survives the parameter-matched control at corpus scale.

**Case II. T > F under Condition A but not F-large under Condition B.** Any advantage of T comes from capacity, not from relational typing. This would be a genuine null result on the framing and would be reported as such. The pipeline is still validated as a substrate for further work.

**Case III. T not significantly different from F under either condition.** Relation typing carries no measurable signal on this substrate. Either the arc types are redundant with node types on this data, or the R-GCN parameterization does not exploit them. Report and discuss. A null on this preregistration is a substantive biological finding, not a failure, and will be reported as such.

**Case IV. F significantly greater than T under Condition A or F-large significantly greater than T under Condition B.** Unusual. Would suggest overparameterization, basis-decomposition instability, or a substrate with no relational signal at all. Report with a follow-up ablation on `num_bases` for T. Distinguish Case IVa (F significantly greater than T under Condition A) from Case IVb (F-large significantly greater than T under Condition B while F ≤ T under Condition A). IVa suggests typing actively harms. IVb suggests capacity beats typing on this substrate.

**Case V. T greater than F-large under Condition B but T not greater than F under Condition A.** Unlikely but possible if F-large's added capacity harms it (overfitting to a substrate where typing gives no signal). Report as inconclusive on H1 with a follow-up ablation on F-large's hidden dim.

## 11. Compute budget

- Per-run estimates depend on realised graph scale. Small-end (10^5 arcs, full-graph): approximately 1 hour training plus 5 minutes evaluation per run. Thirty runs (three model variants times ten seeds) fits in about 32.5 hours. Large-end (10^6 arcs, neighbor-sampled): approximately 3 to 4 hours per run, thirty runs at 90 to 120 hours
- Target compute allocation. 20 to 40 GPU-hours on a single V100 (or L4, or A100 20 percent share) for the small-end case
- Contingency. If, after the first three T-under-Condition-A runs, the observed per-run time times thirty exceeds 40 hours, either reduce to n = 5 seeds per model (with the statistical-power renegotiation logged in `deviations.md`) or request additional compute. The decision is made once, before any test-set MRR is inspected
- Hard stop at 40 GPU-hours for the primary allocation. If the budget is exhausted with n less than 10 seeds per model, results are reported as an interim analysis and the primary test is not run until n = 10 is reached, possibly on a follow-up allocation. Interim analyses do not constitute a test of H1

## 12. Deviations policy

- Every deviation from this protocol is documented in `deviations.md` on the preregistration branch, with the reason and timestamp, before any downstream analysis reads the affected result
- No hyperparameter tuning against the test set. Hyperparameters are fixed in this document. The only tunable knob is the validation-based early-stopping decision
- No cherry-picking of seeds. All ten reported
- If the emitted graph is smaller than 10^3 arcs, the preregistration is suspended and the reason logged. Rerun on a larger corpus
- If a seed diverges or hits a numerical issue, the seed is reported as a failure and included in the mean-difference computation as a zero, not dropped
- **Data source.** If the pinned Reactome release is unavailable at run time, the next-earliest available release is substituted and logged. If the substitution is more than six months from the pinned release, the preregistration is re-frozen with the new release before runs begin
- **Adapter version.** Locked to a single git SHA committed before runs start. Section 3's "`v0.2.0` or later" clause is superseded once the SHA is written into `pipeline_config.json`
- **F-large width.** F-large's hidden dimension is chosen from {144, 160, 176, 192, 208, 224, 256} as the smallest value whose parameter count is within minus 5 percent to plus 5 percent of T's. The chosen value is committed to `pipeline_config.json` before any run
- **Number of relations.** Section 4 commits to `num_relations = 6` and stratified splits over six classes. If the adapter emits a different number of relation labels on the pinned Reactome release, the preregistration is re-frozen with the new label count and matching `num_bases` before runs begin, and the deviation is logged

## 13. Timeline

- **Preregistration frozen.** 2026-07-05 (this commit)
- **Data pull and pipeline verification.** Target August 2026
- **Training runs.** Target September to October 2026
- **Analysis and short report.** Target October 2026
- **Public archive of results.** Depends on outcome and downstream publication plans

## 14. Notes and out-of-scope items

- Alternative decoders (ComplEx, TransE, RESCAL). Out of scope. DistMult is fixed for this preregistration
- Alternative encoders (block-diagonal R-GCN, HetGNN, HGT, GAT). Out of scope
- Alternative pathway ontologies (BioPAX Level 3, PSI-MITAB). Out of scope
- LINCS L1000 to Reactome projection for a downstream drug-response task. Out of scope. Noted in the Rousu cover letter as an aspirational downstream target and is not part of this preregistration
- Filtered ranking as the primary decision statistic. Out of scope. Filtered is a sensitivity check only, and raw MRR remains the primary
- Kernel-methods baselines (Weisfeiler-Lehman kernel, subtree kernel, graph-kernel-based link prediction). Out of scope for this preregistration, though they are the natural next comparison and would be preregistered separately

## 15. Contact

Questions or requests for early access to the trained checkpoints. Email Noor Muhammad at [nmuhammad0900@gmail.com](mailto:nmuhammad0900@gmail.com), or open an issue on this repository.

---

*This preregistration was committed to the public repository before any of the training runs it describes were executed. Any future results that reference this protocol are checkable against this commit.*
