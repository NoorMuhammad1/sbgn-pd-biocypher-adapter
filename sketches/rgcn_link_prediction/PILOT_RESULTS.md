# Pilot execution of the preregistered typed-vs-flat R-GCN comparison

**Date run.** 2026-07-05
**Companion.** [`run_pilot.py`](run_pilot.py) and [`outputs/pilot_results.json`](outputs/pilot_results.json)
**Reference.** [`PREREGISTRATION.md`](PREREGISTRATION.md) (filed 2026-07-05, before this pilot was run)

---

## What this document is

An honest report of a pilot execution of the three-arm ablation described in `PREREGISTRATION.md`, at reduced scale. The preregistered protocol targets a Reactome-scale SBGN-PD graph and 20 GPU-hours of compute. This pilot is CPU-only, on a small synthetic SBGN-PD-like graph, with substantially shortened training. Deviations from the preregistration are listed explicitly below, in the same file as the results, so anyone reading the numbers can weight them appropriately.

The pilot is not the preregistered finding. The preregistered finding will require the Reactome-scale run described in `PREREGISTRATION.md` Section 3.

---

## What was run

Three model arms as specified in `PREREGISTRATION.md` Section 4, with the same DistMult decoder and the same training loop.

- **T** (typed) — R-GCN with basis decomposition, `num_relations=6`, `num_bases=3`, `hidden_dim=32`
- **F** (flat) — R-GCN with all relations collapsed to one, `num_relations=1`, `num_bases=1`, `hidden_dim=32`
- **F-large** — Condition B, parameter-matched to T. Hidden dim chosen by binary search over the R-GCN parameter counter. Landed at `hidden_dim=38`, giving 21,016 parameters versus T's 21,220 (0.96% below target, well within the 5% tolerance in Section 4.1).

**Data.** Synthetic multi-relational graph. 400 nodes, 1197 edges after self-loop rejection, 6 relations. Preferential-attachment edge generation with SBGN-PD-like relation weights (dominant `has_input` / `has_output`, moderate `catalyzes` / `regulates`, sparse containment and association). Per-relation counts in `outputs/pilot_results.json`.

**Split.** 80 / 10 / 10 arc-level split, seed-independent, frozen once. 957 train / 120 val / 120 test.

**Training.** Adam, lr `1e-3`, 5 negatives per positive, 25 epochs with early stopping on val MRR checked every 10 epochs.

**Seeds.** 10 (0 through 9) per arm, 30 runs total.

**Metrics.** Raw MRR, Hits@1, Hits@3, Hits@10 on the held-out test split. Filtered ranking not run at pilot scale.

**Statistical test.** Wilcoxon signed-rank on paired per-seed MRR differences, one-sided (H1: T > baseline), alpha 0.05 (per `PREREGISTRATION.md` Section 8). Also reports paired t-test as sensitivity check and Cohen's d as effect size.

---

## Results

| Arm      | Mean test MRR | 1.96 × SE | Parameters |
|----------|---------------|-----------|------------|
| T        | 0.0320        | 0.0046    | 21,220     |
| F        | 0.0368        | 0.0050    | 16,930     |
| F-large  | 0.0413        | 0.0044    | 21,016     |

**Paired tests (H1: T > baseline).**

| Comparison   | Mean diff | Wilcoxon p (one-sided) | t-test p (one-sided) | Cohen's d |
|--------------|-----------|------------------------|----------------------|-----------|
| T vs F       | -0.0049   | 0.884                  | 0.871                | -0.628    |
| T vs F-large | -0.0093   | 0.990                  | 0.991                | -1.289    |

Reading the numbers with the preregistered decision rule from `PREREGISTRATION.md` Section 10.

- **T > F significantly.** Not observed. The direction is opposite (F > T).
- **T > F under Condition A but not under Condition B.** Not observed.
- **T not significantly different from F under either condition.** Not observed under Wilcoxon at alpha 0.05.
- **F > T significantly under either condition.** This is what the pilot shows for F-large. Case IV in the preregistration. The preregistered follow-up in Section 10 says "Report with a follow-up ablation on num_bases."

---

## Interpretation, with the caveats up front

At pilot scale, on this specific synthetic graph, with this specific hyperparameter setting, a typed R-GCN with basis decomposition does not beat an otherwise-identical R-GCN that collapses all relation labels. A parameter-matched flat R-GCN beats it. Four candidate explanations are consistent with the numbers.

1. **Training too short.** T has more parameters concentrated in the per-relation basis coefficients. 25 epochs is far below the preregistered 500. It is plausible T needs longer to converge and F converges faster on the same budget.
2. **Synthetic data has no relation-carried signal.** The preferential-attachment generator assigns relation labels according to a fixed prior, independent of node identity. There is no reason for a typed encoder to exploit the labels because the labels do not carry information about which node pairs are connected. If Reactome pathways do carry that structure, T should recover on the real data.
3. **Basis decomposition instability at this scale.** `num_bases=3` might be too few or too many for a 400-node graph. Preregistration Section 10 Case IV explicitly says to run a `num_bases` ablation as follow-up.
4. **Genuine null on relation typing.** In this specific setting, on this specific substrate, relation typing does not add signal beyond topology. The Reactome-scale run would be the arbiter.

The point of running a preregistered pilot with the loss and reporting the direction that came out is that this is exactly the outcome the preregistration was built to handle. Case IV in Section 10 has a pre-committed reading. The pilot did not force a T > F story; it forced whatever the numbers said, which happened to be F-large > T.

---

## Deviations from `PREREGISTRATION.md`

Section 12 of the preregistration says every deviation must be documented in the same branch, before results are read. This pilot's deviations are listed below.

1. **Scale.** Preregistration requires Reactome-scale data (~10^4 nodes, 10^5-10^6 arcs, per Section 3). Pilot ran at 400 nodes, 1197 arcs. Preregistration Section 12 says if the graph is smaller than 10^3 arcs the preregistration is suspended. The pilot graph is at 1197 arcs, which is above that floor but still ~two orders of magnitude below the intended scale. **This makes the pilot not the preregistered finding.**
2. **Data source.** Preregistration requires Reactome SBGN-ML files ingested through the SBGN-PD BioCypher adapter. Pilot uses a synthetic preferential-attachment graph. Deviation logged. Rationale: pulling the Reactome corpus and running the adapter through the full pipeline was out of scope for the compute window available.
3. **Hidden dimension.** Preregistration specifies `hidden_dim=128`. Pilot used `hidden_dim=32`. Deviation logged. Rationale: 32 keeps the parameter-counter arithmetic transparent and reduces per-seed runtime.
4. **Training length.** Preregistration specifies 500 epochs with early stopping (patience 30). Pilot used 25 epochs with lighter early stopping. Deviation logged. Rationale: compute budget.
5. **Filtered ranking.** Preregistration Section 7 requires filtered ranking as a sensitivity check. Pilot reports raw ranking only. Deviation logged.
6. **Compute environment.** Preregistration budgets 20 GPU-hours (Section 11). Pilot ran on CPU. Deviation logged.

All six deviations were known before running. The pilot was not adjusted after seeing preliminary results.

---

## What the pilot actually establishes

- The pipeline runs end-to-end. Synthetic graph generation, T / F / F-large model construction, training loop, MRR evaluation, Wilcoxon test, plot generation. All work.
- The parameter-matched control condition can be produced automatically (binary search over hidden dim), so Condition B in the preregistration is executable.
- The preregistered outcome-interpretation logic in Section 10 gives a clean reading for a case that turned out F-large > T. This validates the decision rule.
- The pipeline is a substrate a later Reactome-scale run can slot into without design changes.

It does not establish that typed R-GCN under-performs flat R-GCN on real pathway data. That is what the Reactome-scale run is for.

---

## What comes next under the preregistration

Per Section 10 Case IV: run the pre-committed `num_bases` ablation on the Reactome-scale data, and report whether the F-large > T direction persists at that scale or reverses. That is a scheduled follow-up, not a pilot.

The follow-up is filed as a note here rather than executed, because the pilot's job was to test the pipeline, not to compete for the preregistered finding.

---

## Files

- [`run_pilot.py`](run_pilot.py) — the script that produced these results
- [`outputs/pilot_results.json`](outputs/pilot_results.json) — full per-seed and per-arm output
- [`outputs/pilot_curves.png`](outputs/pilot_curves.png) — training curves and bar chart

## References

- Schlichtkrull, Kipf, Bloem, van den Berg, Titov, Welling. *Modeling Relational Data with Graph Convolutional Networks*. ESWC 2018.
- Yang, Yih, He, Gao, Deng. *Embedding Entities and Relations for Learning and Inference in Knowledge Bases*. ICLR 2015.
