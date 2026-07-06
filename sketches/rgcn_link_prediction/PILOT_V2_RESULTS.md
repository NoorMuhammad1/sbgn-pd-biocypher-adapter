# Extended pilot v2. T+topo arm added, still synthetic, still not the preregistered finding.

**Date run.** 2026-07-05
**Companion.** [`run_pilot_v2.py`](run_pilot_v2.py) and [`outputs/pilot_v2_results.json`](outputs/pilot_v2_results.json)
**Reference.** [`PREREGISTRATION.md`](PREREGISTRATION.md) and [`PILOT_RESULTS.md`](PILOT_RESULTS.md) (v1)
**Design rationale.** [`../reification_vs_missingness/NOTE.md`](../reification_vs_missingness/NOTE.md) and [`../topological_followup/DESIGN_NOTE.md`](../topological_followup/DESIGN_NOTE.md)

---

## What is new versus v1

Pilot v1 compared T (typed R-GCN with the six SBGN-PD relations preserved), F (flat R-GCN with all relations collapsed to one), and F-large (F with hidden dim increased until parameter count matches T within 5 percent) at 400 nodes / 1197 arcs / 25 epochs on a synthetic preferential-attachment graph. F-large won at Wilcoxon p equals 0.99, one-sided.

Pilot v2 does two things.

1. **Adds a T+topo arm.** A per-node per-relation degree signature (2 x num_relations = 12-dim vector per node, log-scaled in / out degrees per relation) is concatenated to the R-GCN output embedding before DistMult scoring. The DistMult relation embedding is widened to match. The signature is not learned. It is computed from the training triples only. This is a defensible cousin of the per-relation graph-magnitude fingerprint proposed in `../topological_followup/DESIGN_NOTE.md`, chosen because true Leinster graph magnitude is O(V^3) per relation and out of scope for the compute budget of this pilot.
2. **Runs at higher synthetic scale.** 800 nodes / 2400 arcs / 10 seeds / 20 epochs. Two-fold larger graph than v1, same number of seeds, same protocol.

Reactome-scale execution remains scheduled follow-up. The Reactome ContentService SBGN exporter endpoint was Cloudflare-gated at the time this pilot was run.

---

## Results

Graph after preferential-attachment generation. 800 nodes, 2400 arcs, six relations with SBGN-PD-like weights. Arc-level 80/10/10 split, frozen once, 1915 train / 240 val / 240 test.

Invocation used to produce these results. `python run_pilot_v2.py --num-nodes 800 --num-edges 2400 --seeds 10 --epochs 20 --lr 1e-2`.

| Arm      | Mean test MRR | 1.96 x SE | Parameters |
|----------|---------------|-----------|------------|
| T        | 0.0149        | 0.0021    | 34,020     |
| F        | 0.0235        | 0.0030    | 29,730     |
| F-large  | 0.0226        | 0.0039    | 34,022     |
| T+topo   | 0.0166        | 0.0022    | 34,092     |

**Paired Wilcoxon signed-rank, one-sided, exact null distribution, alternative that A ranks above B.**

| Comparison             | Mean diff | p       | Cohen's d |
|------------------------|-----------|---------|-----------|
| T > F                  | -0.0086   | 1.000   | -2.06     |
| T > F-large            | -0.0077   | 0.990   | -1.53     |
| **T+topo > T**         | **+0.0017** | **0.116** | **+0.49** |
| T+topo > F-large       | -0.0060   | 0.986   | -1.18     |

Reading with the preregistered decision rule from `PREREGISTRATION.md` Section 10.

- **T vs F under Condition A.** Primary Wilcoxon (H1: T greater than F) fails to reject at p equal to 1.000. Diagnostic Wilcoxon in the complementary direction (H1: F greater than T) is significant at p close to zero. Under revised `PREREGISTRATION.md` Section 10, this is Case IVa (F significantly greater than T under Condition A). The preregistered follow-up is a `num_bases` ablation, deferred to the Reactome-scale run where a Case IVa outcome would carry more weight than at 800 synthetic nodes.
- **T vs F-large under Condition B.** Primary Wilcoxon fails to reject at p equal to 0.990. Diagnostic Wilcoxon is significant at p approximately 0.010. Because Condition A also reached Case IV, v2's outcome is Case IVa (typing actively harms on this substrate), whereas v1's outcome was Case IVb (capacity beats typing while F did not significantly beat T at Condition A).
- **T+topo vs T (the topological intervention on the typed encoder).** Mean improvement of +0.0017 in MRR, medium effect size (d=+0.49), but not statistically significant at alpha 0.05 (Wilcoxon one-sided p=0.116). This is a **suggestive positive direction that does not reach significance** at pilot scale with pilot data.
- **T+topo vs F-large.** F-large still wins significantly. The topological signature does not lift T+topo above the flat-with-more-capacity baseline on this substrate.

---

## What this tells us and what it does not

**What it tells us.**

1. On synthetic preferential-attachment data, relation typing does not by itself add link-prediction signal. Two pilots at different scales (400 nodes and 800 nodes) both reproduce F > T. This is exactly the outcome the preregistration was built to handle. The likely mechanism is that the synthetic generator has no reason to make relation labels carry information about node pairs, so a model that pretends they do gets worse. The v1 alternatives (undertraining, num_bases instability at the synthetic scale) remain live, though the reproducibility of F > T across 400 and 800 nodes shifts weight toward the substrate explanation.
2. Adding a per-node per-relation degree signature to the typed model produces a **small positive direction** (+0.49 SD, +0.0017 MRR) that fails to reach significance at n=10. This is the sort of trend that gets bigger with more seeds, more epochs, or more capacity. It is also the sort of trend that could vanish under a different seed selection. Do not read this as a real T+topo advantage yet.
3. Neither typing nor typing + topological augmentation reaches the capacity-matched flat baseline on this substrate. That is the substrate telling us it has no relational signal to exploit. Real curated pathway data is where the substrate's answer changes.

**What it does not tell us.**

1. It does not tell us anything about typed R-GCN performance on Reactome. Two orders of magnitude of scale sit between this pilot and the preregistered target.
2. It does not tell us whether persistent path homology or magnitude homology (`DESIGN_NOTE.md`'s intended topological arm) would behave the way the degree signature does. The degree signature is a cheap surrogate. The real topological methods might produce a much larger effect, or a much smaller one, or the opposite direction.
3. It does not tell us whether a soft-manifold embedding (Marinoni et al. 2026 TPAMI) on top of the reified graph would recover the signal. That is the four-arm study preregistration proposed in `../reification_vs_missingness/NOTE.md`.

---

## Deviations and scope extensions from `PREREGISTRATION.md`

Section 12 of the preregistration says every deviation is documented in the same branch, before results are read. v2's deviations, listed below, were fixed before the runs started.

1. **Scale.** Preregistration targets Reactome-scale (~10^4-10^5 nodes). v2 ran at 800 nodes. Two orders of magnitude below the target. Preregistered ceiling of Section 12 (graph must exceed 10^3 arcs) is met at 2400 arcs, but the pilot is still not the preregistered finding.
2. **Data source.** Preregistration requires Reactome SBGN-ML. v2 uses synthetic preferential attachment. Same rationale as v1. The Reactome ContentService SBGN exporter was Cloudflare-gated at run time.
3. **Training length.** Preregistration specifies 500 epochs. v2 uses 20. Compute budget.
4. **Hidden dimension.** Preregistration specifies 128. v2 uses 32. Compute budget.
5. **Compute environment.** Preregistration budgets 20 GPU-hours. v2 ran on CPU.
6. **Filtered ranking.** Preregistration requires filtered ranking as a sensitivity check. v2 reports raw ranking only.
7. **Learning rate.** Preregistration Section 5 specifies `lr = 1e-3`. v2 was invoked with `lr = 1e-2`. Ten-fold larger. Compute budget shortening called for faster convergence at 20 epochs. Logged as a deviation here rather than left implicit in the CLI invocation.
8. **Early stopping.** Preregistration specifies patience 30 with minimum improvement 0.005 on validation MRR. v2 uses a 20-epoch hard cap with periodic checkpointing of the best-val model, no patience criterion.
9. **Sensitivity paired t-test.** Preregistration Section 8 specifies a Shapiro-Wilk normality precheck at alpha 0.10 followed by a paired t-test as sensitivity companion to Wilcoxon. Not reported. Wilcoxon and Cohen's d_z are reported.
10. **Secondary metrics.** Preregistration Section 7 requires Hits@1, Hits@3, and Hits@10 alongside MRR. v2 reports MRR only.
11. **Case IV `num_bases` ablation deferred.** Preregistration Section 10 Case IV asks for a `num_bases` ablation when F significantly beats T. Deferred to the Reactome-scale run, where the outcome would carry more weight than at 800 synthetic nodes.

Scope extension (not a deviation).

12. **New arm not in the preregistration.** T+topo is an addition to the preregistered T / F / F-large three-arm structure. Under Section 12, this is a scope extension. It is documented here and not treated as part of the preregistered finding. If a proper Reactome-scale run were to include the topo arm, that inclusion would go into a separate `PREREGISTRATION-topo.md` first.

All eleven deviations and the one scope extension were known before running.

---

## What comes next under the preregistration

1. Full Reactome-scale run once ContentService SBGN exporter access is unblocked (Cloudflare workaround, or bulk download route). Runs T / F / F-large under the preregistered protocol at ~10^4 nodes, 500 epochs, hidden_dim 128, 20 GPU-hours.
2. Separate `PREREGISTRATION-topo.md` filing before running a topological arm at Reactome scale. The current T+topo pilot is a substrate test, not a preregistered study.
3. Separate `PREREGISTRATION-4arm.md` for the reification-vs-embedding-space four-arm study from `../reification_vs_missingness/NOTE.md`. This is a larger design and belongs in its own preregistration.

---

## Files

- [`run_pilot_v2.py`](run_pilot_v2.py) is the script that produced these results
- [`outputs/pilot_v2_results.json`](outputs/pilot_v2_results.json) has the full per-seed and per-arm output
- [`outputs/pilot_v2_curves.png`](outputs/pilot_v2_curves.png) shows training curves and the bar chart

## References

- Schlichtkrull, Kipf, Bloem, van den Berg, Titov, Welling. *Modeling Relational Data with Graph Convolutional Networks*. ESWC 2018.
- Yang, Yih, He, Gao, Deng. *Embedding Entities and Relations for Learning and Inference in Knowledge Bases*. ICLR 2015.
- Marinoni, Liò, Barp, Girolami. *Improving Embedding of Graphs With Missing Data by Soft Manifolds*. IEEE TPAMI 48(3), 2026.
- Leinster. *The magnitude of metric spaces*. Documenta Mathematica, 2013.
