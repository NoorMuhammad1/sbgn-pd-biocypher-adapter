# Pilot execution of the preregistered typed-vs-flat R-GCN comparison

**Date run.** 2026-07-05
**Companion.** [`run_pilot.py`](run_pilot.py) and [`outputs/pilot_results.json`](outputs/pilot_results.json)
**Reference.** [`PREREGISTRATION.md`](PREREGISTRATION.md) (filed 2026-07-05, before this pilot was run)
**Successor.** [`PILOT_V2_RESULTS.md`](PILOT_V2_RESULTS.md) runs the same three arms at higher synthetic scale and adds a T+topo arm.

---

## What this document is

An honest report of a pilot execution of the three-arm ablation described in `PREREGISTRATION.md`, at reduced scale. The preregistered protocol targets a Reactome-scale SBGN-PD graph and 20 GPU-hours of compute. This pilot is CPU-only, on a small synthetic SBGN-PD-like graph, with substantially shortened training. Deviations from the preregistration are listed explicitly below, in the same file as the results, so anyone reading the numbers can weight them appropriately.

The pilot is not the preregistered finding. The preregistered finding will require the Reactome-scale run described in `PREREGISTRATION.md` Section 3.

---

## What was run

Three model arms as specified in `PREREGISTRATION.md` Section 4, with the same DistMult decoder and the same training loop.

- **T** (typed). R-GCN with basis decomposition, `num_relations=6`, `num_bases=3`, `hidden_dim=32`.
- **F** (flat). R-GCN with all relations collapsed to one, `num_relations=1`, `num_bases=1`, `hidden_dim=32`.
- **F-large**. Condition B, parameter-matched to T. Hidden dim chosen by binary search over the R-GCN parameter counter. Landed at `hidden_dim=38`, giving 21,016 parameters versus T's 21,220 (0.96% below target, well within the 5% tolerance in Section 4.1).

**Data.** Synthetic multi-relational graph. 400 nodes, 1197 edges after self-loops were skipped during generation, 6 relations. Preferential-attachment edge generation with SBGN-PD-like relation weights (dominant `has_input` and `has_output`, moderate `catalyzes` and `regulates`, sparse containment and association). Per-relation counts in `outputs/pilot_results.json`.

**Split.** 80 / 10 / 10 arc-level split, seed-independent, frozen once. 957 train / 120 val / 120 test.

**Training.** Adam, lr `1e-2` (see deviations), 5 negatives per positive, 25 epochs. Best-validation checkpoint restored after every-10-epoch validation MRR checks. No patience-based early termination.

**Seeds.** 10 (0 through 9) per arm, 30 runs total.

**Metrics.** Raw MRR on the held-out test split. Hits@1 / Hits@3 / Hits@10 computed per-seed in `outputs/pilot_results.json` but not reproduced in this report's summary tables (see deviations). Filtered ranking not run at pilot scale.

**Statistical test.** Wilcoxon signed-rank on paired per-seed MRR differences, one-sided under H1 that T beats each baseline, alpha 0.05 (per `PREREGISTRATION.md` Section 8). Paired t-test reported as sensitivity check (Shapiro-Wilk normality gate not applied, see deviations). Effect size reported as Cohen's d in the pooled-SD form `d_av = mean(d_i) / sqrt((var_T + var_baseline)/2)`. The preregistration commits to `d_z` on paired differences (see deviations).

---

## Results

| Arm      | Mean test MRR | 1.96 × SE | Parameters |
|----------|---------------|-----------|------------|
| T        | 0.0320        | 0.0046    | 21,220     |
| F        | 0.0368        | 0.0050    | 16,930     |
| F-large  | 0.0413        | 0.0044    | 21,016     |

**Paired tests, one-sided under H1 that T beats each baseline.**

| Comparison   | Mean diff | Wilcoxon p (H1: T > baseline) | Diagnostic p (H1: baseline > T) | t-test p | Cohen's d |
|--------------|-----------|-------------------------------|--------------------------------|----------|-----------|
| T vs F       | -0.0049   | 0.884                         | ~0.116                         | 0.871    | -0.628    |
| T vs F-large | -0.0093   | 0.990                         | ~0.010                         | 0.991    | -1.289    |

The diagnostic column is the complementary-direction Wilcoxon p, reported per revised `PREREGISTRATION.md` Section 8. Diagnostic values here are computed as `1 - p_primary` and hold approximately. The exact Wilcoxon under `alternative="less"` handles ties slightly differently.

Reading the numbers with the preregistered decision rule from `PREREGISTRATION.md` Section 10.

- **T greater than F significantly.** Primary Wilcoxon fails to reject H1 (p = 0.884). Not observed. Direction is opposite (mean diff is negative).
- **T greater than F under Condition A but not under Condition B.** Not observed.
- **T not significantly different from F under either condition.** Under the primary test this is technically the current reading for T vs F (H1: T greater than F fails to reject, and the diagnostic direction p ≈ 0.116 also fails to reject at alpha 0.05). Under the primary test for T vs F-large this reading does not hold, because the diagnostic direction reaches alpha 0.05.
- **F significantly greater than T under Condition B (Case IVb).** Diagnostic Wilcoxon in the complementary direction gives p ≈ 0.010 for T vs F-large. Under revised `PREREGISTRATION.md` Section 10, this is Case IVb (F-large significantly greater than T under Condition B while F is not significantly greater than T under Condition A). Preregistered follow-up in Section 10 says "Report with a follow-up ablation on num_bases for T."

---

## Interpretation, with the caveats up front

At pilot scale, on this specific synthetic graph, with this specific hyperparameter setting, a typed R-GCN with basis decomposition does not beat an otherwise-identical R-GCN that collapses all relation labels. A parameter-matched flat R-GCN beats it. Four candidate explanations are consistent with the numbers.

1. **Training too short.** T has more parameters concentrated in the per-relation basis coefficients. 25 epochs is far below the preregistered 500. It is plausible T needs longer to converge and F converges faster on the same budget.
2. **Synthetic data has no relation-carried signal.** The preferential-attachment generator assigns relation labels according to a fixed prior, independent of node identity. There is no reason for a typed encoder to exploit the labels because the labels do not carry information about which node pairs are connected. If Reactome pathways do carry that structure, T should recover on the real data.
3. **Basis decomposition instability at this scale.** `num_bases=3` might be too few or too many for a 400-node graph. Preregistration Section 10 Case IV explicitly says to run a `num_bases` ablation as follow-up.
4. **Genuine null on relation typing.** In this specific setting, on this specific substrate, relation typing does not add signal beyond topology. The Reactome-scale run would be the arbiter.

The point of running a preregistered pilot with the loss and reporting the direction that came out is that this is exactly the outcome the preregistration was built to handle. Case IVb in Section 10 has a pre-committed reading. The pilot did not force a T greater than F story. It reported the direction the numbers gave, which happened to be F-large greater than T at Condition B.

---

## Deviations from `PREREGISTRATION.md`

Section 12 of the preregistration says every deviation must be documented in the same branch, before results are read. This pilot's deviations are listed below.

1. **Scale.** Preregistration requires Reactome-scale data (~10^4 nodes, 10^5-10^6 arcs, per Section 3). Pilot ran at 400 nodes, 1197 arcs. Preregistration Section 12 says if the graph is smaller than 10^3 arcs the preregistration is suspended. The pilot graph is at 1197 arcs, which is above that floor but still ~two orders of magnitude below the intended scale. **This makes the pilot not the preregistered finding.**
2. **Data source.** Preregistration requires Reactome SBGN-ML files ingested through the SBGN-PD BioCypher adapter. Pilot uses a synthetic preferential-attachment graph. Deviation logged. Rationale was that pulling the Reactome corpus and running the adapter through the full pipeline was out of scope for the compute window available.
3. **Hidden dimension.** Preregistration specifies `hidden_dim=128`. Pilot used `hidden_dim=32`. Deviation logged. Rationale was that 32 keeps the parameter-counter arithmetic transparent and reduces per-seed runtime.
4. **Training length.** Preregistration specifies 500 epochs. Pilot used 25 epochs. Deviation logged. Rationale was compute budget.
5. **Learning rate.** Preregistration Section 5 specifies `lr = 1e-3`. Pilot used `lr = 1e-2`. Ten-fold larger. Deviation logged retroactively after the second-pass 4-agent review caught the drift. The pilot was not re-run.
6. **Early stopping.** Preregistration Section 5 specifies patience 30 with minimum improvement 0.005 on validation MRR. Pilot uses a 25-epoch hard cap with best-val checkpointing every 10 epochs and no patience criterion.
7. **Effect-size formula.** Preregistration Section 8 (revised) commits to `d_z = mean(d_i) / sd(d_i)`. Pilot reports `d_av = mean(d_i) / sqrt((var_T + var_baseline)/2)`. Both are legitimate paired-samples formulas but they differ numerically.
8. **BCa confidence interval on Cohen's d.** Preregistration Section 8 (revised) requires bootstrap BCa CI. Pilot does not report it.
9. **Shapiro-Wilk normality gate.** Preregistration Section 8 requires a Shapiro-Wilk precheck at alpha 0.10 before running the paired t-test as sensitivity. Pilot reports paired t-test unconditionally without the gate.
10. **Stratified split with per-class floor.** Preregistration Section 6 (revised) requires stratification by relation with at least 30 arcs per class per split. Pilot uses a simple shuffle. Rel_5 has 49 arcs total, and after 10 percent test allocation the test split holds approximately 5 arcs of that class.
11. **Wilcoxon zero-difference handling.** Preregistration Section 8 (revised) specifies the Pratt method. Pilot uses scipy's default automatic mode.
12. **Hits@k reporting.** Preregistration Section 7 requires Hits@1 / Hits@3 / Hits@10 alongside MRR in the primary report. Pilot computes them per-seed in `outputs/pilot_results.json` but omits them from the results tables above.
13. **Filtered ranking.** Preregistration Section 7 requires filtered ranking as a sensitivity check. Pilot reports raw ranking only.
14. **Compute environment.** Preregistration budgets 20 GPU-hours (Section 11). Pilot ran on CPU.
15. **Case IVb `num_bases` follow-up ablation.** Preregistration Section 10 (revised) requires a `num_bases` ablation on Case IVb. Deferred to the Reactome-scale run, where the outcome would carry more weight than at 400 synthetic nodes.

All fifteen deviations were known before running or were logged retroactively during the second-pass review before any interpretive claim in the report was modified.

---

## What the pilot actually establishes

- The pipeline runs end-to-end. Synthetic graph generation, T / F / F-large model construction, training loop, MRR evaluation, Wilcoxon test, plot generation. All work.
- The parameter-matched control condition can be produced automatically (binary search over hidden dim), so Condition B in the preregistration is executable.
- Section 10's outcome-interpretation logic gave a clean reading for a case that landed at F-large greater than T under the diagnostic direction. This exercises the decision rule on one live case rather than validating it in general.
- The pipeline appears to serve as scaffolding for a later Reactome-scale run, but several preregistered code changes are still required before that run. Stratified splitting with the 30-arc-per-class floor. Patience-30 minimum-improvement-0.005 early stopping. Shapiro-Wilk normality gate on the sensitivity t-test. Pratt zero-difference handling in Wilcoxon. `d_z` and BCa CI on Cohen's d. Filtered-ranking sensitivity metric. Hits@k reporting in the summary. Learning rate reverted to `1e-3`.

It does not establish that typed R-GCN under-performs flat R-GCN on real pathway data. That is what the Reactome-scale run is for.

---

## What comes next under the preregistration

Section 10 Case IV's follow-up says to run the pre-committed `num_bases` ablation on the Reactome-scale data and report whether the F-large greater than T direction persists at that scale or reverses. That is a scheduled follow-up, not a pilot.

The follow-up is filed as a note here rather than executed, because the pilot's job was to test the pipeline, not to compete for the preregistered finding.

---

## Files

- [`run_pilot.py`](run_pilot.py) is the script that produced these results.
- [`outputs/pilot_results.json`](outputs/pilot_results.json) has the full per-seed and per-arm output.
- [`outputs/pilot_curves.png`](outputs/pilot_curves.png) shows training curves and the bar chart.

## References

- Schlichtkrull, Kipf, Bloem, van den Berg, Titov, Welling. *Modeling Relational Data with Graph Convolutional Networks*. ESWC 2018.
- Yang, Yih, He, Gao, Deng. *Embedding Entities and Relations for Learning and Inference in Knowledge Bases*. ICLR 2015.
