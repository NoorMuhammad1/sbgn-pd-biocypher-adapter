# Design memo: Directed topological descriptors on reified SBGN-PD graphs

**Author.** Noor Muhammad ([github.com/NoorMuhammad1](https://github.com/NoorMuhammad1))
**Date.** 2026-07-05
**Companion artifacts.** [`sketches/rgcn_link_prediction/`](../rgcn_link_prediction/) (R-GCN sketch) and [`sketches/rgcn_link_prediction/PREREGISTRATION.md`](../rgcn_link_prediction/PREREGISTRATION.md) (typed-vs-flat protocol)
**Status.** Working design memo. Not a preregistration. Committed publicly before any of the experiments discussed have been run.

---

## Why this memo exists

The R-GCN sketch in `sketches/rgcn_link_prediction/` demonstrates that the SBGN-PD BioCypher adapter's output feeds a relation-aware graph neural network end-to-end. The accompanying preregistration compares a typed encoder against a flattened one on link prediction, under a Wilcoxon signed-rank protocol with a parameter-matched control. Alternative encoders (block-diagonal R-GCN, HetGNN, HGT, GAT) are held out of scope by design.

That scope decision leaves a natural next question. If R-GCN with basis decomposition is bounded by the multi-relational 1-WL test (Barceló et al. 2022), and the Vietoris-Rips PH vs 1-WL result (Ballester and Rieck 2023) shows the two are incomparable rather than one dominating, then adding a topological descriptor as a second view should recover signal that neither view has alone. This memo works through which topological descriptor is the right one for SBGN-PD, and what the smallest useful preregistered follow-up would look like.

---

## 1. Which flavor of topology is right for SBGN-PD?

SBGN-PD graphs are **directed and typed**. Each process has ordered participants (inputs, outputs, catalysts, modulators), and the direction encodes reaction semantics that a downstream ML pipeline should not discard. The choice of descriptor family has to respect this.

Three options, arranged from least to most direction-aware.

**Vietoris-Rips persistent homology (VR-PH)** on the underlying graph with the shortest-path metric. Detects H_0 (connected components) and H_1 (undirected 1-cycles). VR-PH takes a symmetric dissimilarity, so applied to a directed graph one first symmetrises the shortest-path quasi-metric. After that step orientation information is gone. A feedback loop A→B→C→A and an undirected triangle then produce identical VR filtrations and identical H_1 barcodes. This is the wrong tool for SBGN-PD if we care about direction.

**Persistent path homology (Chowdhury and Mémoli 2018, extending Grigor'yan, Lin, Muranov, and Yau 2012)** on the directed graph. Uses regular directed paths as chains and takes a boundary that respects orientation. Distinguishes directed cycles from undirected ones. On a directed 3-cycle A→B→C→A the H_1 generator survives because the triangle cannot be filled by a 2-path in the digraph. A linear waterfall A→B→C→D produces no H_1 class at all. This is the natural descriptor family for feedback-loop detection. The Chowdhury-Mémoli contribution is the persistent version on filtered directed networks, with a stability theorem under network-distance perturbations that matters for curated-versus-noisy pathway data.

**Persistent homology of directed flag complexes (Lütgehetmann, Govc, Levi, and Hess 2020).** Build a simplicial complex where a p-simplex is a totally ordered clique with all edges pointing in the total-order direction. On the participant-projected SBGN-PD graph (obtained by contracting each process node to induce directed input-to-output edges), a coherent feedforward loop X→Y, X→Z, Y→Z is a directed 2-simplex. A feedback loop X→Y, Y→Z, Z→X admits no consistent total order and so is not a simplex at all. This separates coherent-feedforward from feedback. It does not, however, separate coherent from incoherent feedforward, since those differ by the sign of the Y→Z modulation rather than by orientation. Signed persistent homology and edge-type-filtered directed flag complexes are the two candidates for bridging that gap.

**Concrete claim.** For SBGN-PD, the empirical study I want to run is not "PH vs 1-WL" but "path homology plus directed flag complex vs R-GCN." Undirected VR-PH is a Ballester-Rieck 2023 baseline for graph learning in general. On this substrate it is the wrong starting point.

**Substrate ambiguity worth flagging.** SBGN-PD reactions are bipartite (process nodes connect participant nodes with no direct participant-to-participant edges). The directed flag complex in the paragraph above is being described on the participant-projected graph, obtained by contracting each process node to induce directed input-to-output edges. Building the same complex directly on the bipartite reified graph gives a different construction whose simplices span process and participant nodes, and whose expressivity relative to R-GCN would need to be studied separately. Neither construction has been implemented here.

---

## 2. What object should we compute the homology on?

SBGN-PD process glyphs are effectively n-ary in the sense that matters here. The generic `process`, `association`, and `dissociation` glyphs bundle several participants under a single reaction event, whether or not the arity is fixed by the glyph type. See the sibling note `../reification_vs_missingness/NOTE.md` for a fuller carve-out of which glyphs are genuinely n-ary and which are fixed-arity or unary.

Given an SBGN process P with participants {A, B, C, D} carrying roles {input, input, output, catalyst}, three representations are on the table.

**Reified graph.** P becomes a `BiologicalProcess` node. A, B, C, D become participant nodes. Four typed edges (`has_input`, `has_input`, `has_output`, `catalyzes`) connect P to the participants. This is what the adapter emits today. Every homology construction above operates natively on it. The cost is that the n-ary structure of P has been unfolded into 1-cells, so the "these participants co-occur in one reaction" invariant is not encoded at the level of the underlying complex.

**Hypergraph.** P becomes a hyperedge over {A, B, C, D}, labeled with the participant roles. The n-arity is preserved. Persistent homology on hypergraphs is well-defined via the associated hypergraph nerve or by lifting to a bipartite representation. Direction and role labeling are harder to encode natively.

**Simplicial complex.** P becomes a 3-simplex spanning {A, B, C, D} (with lower-dimensional faces filled automatically). Higher-order TDA methods operate natively on simplicial complexes and would see P as a filled tetrahedron. Directed variants (directed flag complexes from Section 1) can carry the participant ordering.

**Worked example.** Consider the hexokinase reaction from glycolysis. Glucose plus ATP produces glucose-6-phosphate plus ADP, catalyzed by hexokinase. Participants:

- Substrates: glucose, ATP
- Products: glucose-6-phosphate, ADP
- Catalyst: hexokinase

In the reified representation the adapter emits, this is a 6-node subgraph (5 participants plus 1 process node P) with 5 typed edges. Explicitly, `has_input(P, glucose)`, `has_input(P, ATP)`, `has_output(P, G6P)`, `has_output(P, ADP)`, and `catalyzes(hexokinase, P)` (catalyst-to-process orientation). In the hypergraph representation, it is a 5-node hyperedge with per-participant role labels {substrate, substrate, product, product, catalyst}. In the simplicial representation, it is a 4-simplex on the 5 participants. Carrying role information at that level requires an oriented and role-colored simplex, a construction that goes beyond the plain simplicial complex, because catalyst is topologically distinct from substrate and product but a bare simplex is symmetric on its vertices.

**Concrete claim.** Which representation is best is an empirical question that has to be answered per downstream task. For **link prediction on relations that already exist in Biolink** (`has_input`, `has_output`, `catalyzes`, `regulates`), the reified graph is the natural substrate because the prediction target is literally the reified edges. For **motif-level or process-level classification** (does this reaction subgraph correspond to a known pathway motif?), the simplicial representation is likely a better substrate because the invariant of interest lives at the level of the higher-order simplex, not at the level of individual edges.

---

## 3. Magnitude and magnitude homology

**Graph magnitude (Leinster 2013).** For a finite metric space (V, d), define a matrix Z with Z_ij = exp(-d(i, j)). The magnitude is the sum of all entries of Z^{-1}. For graphs, d is usually the shortest-path metric. Magnitude gives a real-valued invariant that behaves like an "effective size" of the space.

**Magnitude homology (Hepworth and Willerton 2017).** A bigraded homology theory whose graded Euler characteristic recovers the magnitude power series. Adds structure that magnitude alone does not see. On graphs, magnitude homology encodes information about paths and their length distributions.

**What is known and what is open on typed and directed cases.**

- Magnitude of directed graphs uses the asymmetric shortest-path quasi-metric and gives a well-defined invariant.
- Magnitude of typed graphs (where each edge carries a categorical label from a small alphabet) does not have a clean off-the-shelf definition. Naive approaches assign a single edge weight per typed edge and lose the categorical information. Per-relation magnitude, meaning compute magnitude separately on each relation subgraph and concatenate the resulting scalars, is a defensible engineering choice that has to my knowledge not been benchmarked on curated pathway data.

**Concrete claim.** For SBGN-PD, a per-relation magnitude fingerprint concatenated across the seven Level 1 arc classes (consumption, production, catalysis, modulation, stimulation, inhibition, necessary stimulation) gives a 7-vector descriptor per graph. That descriptor is a candidate topological feature to concatenate with an R-GCN embedding. It is the smallest meaningful topological arm I can add to the preregistered T-vs-F comparison without changing the encoder. The seven-versus-six question depends on whether necessary stimulation is collapsed into stimulation for the adapter's emission, which is an adapter-configuration choice rather than a homology question.

---

## 4. Proposed follow-up to the existing preregistration

The current preregistration at `sketches/rgcn_link_prediction/PREREGISTRATION.md` compares a typed R-GCN (T) against a flattened R-GCN (F) on link prediction over the adapter output. Alternative encoders are out of scope by design.

**Extension proposal.** Add one topological-descriptor arm under the same Wilcoxon protocol.

- **T-topo.** Each node's T-encoder output is concatenated with the graph-level per-relation magnitude fingerprint (a 7-vector broadcast to every node) before the DistMult decoder scores candidate triples. Broadcasting is what makes a graph-level readout usable in a per-node link-prediction pipeline. The alternative, a per-node per-relation magnitude localised via a nested-neighbourhood filtration, is left as a separate arm.
- **Same protocol.** Ten seeds, stratified 80/10/10 arc split (frozen once), raw MRR primary metric, exact one-sided Wilcoxon signed-rank paired test at alpha 0.05, Cohen's d_z effect size on paired differences.
- **Matched controls, not optional.** F-topo (flat R-GCN with the same 7-scalar broadcast concatenation) is the matched control and is part of the primary comparison. F-large-topo (F-topo with hidden dim grown to match T-topo parameter count within 5 percent) is the capacity-matched control, mirroring the R-GCN preregistration's Condition B.
- **Pre-committed effect-size prior.** Cohen's d_z in the 0.2 to 0.5 range on paired MRR difference against T is the range where the fingerprint carries measurable signal without being dominant. Anything below 0.15 is not evidence for the fingerprint's utility. Anything above 0.5 is a surprising strong effect that would earn a follow-up ablation on the fingerprint's dimensionality (7 scalars is a narrow channel and a large effect from such a channel would need explaining).
- **Pre-committed outcomes.** T-topo significantly greater than T means the magnitude fingerprint carries link-prediction signal beyond what typed message passing captures on its own. T-topo not significantly different from T means the magnitude fingerprint is redundant with what R-GCN's basis-decomposed encoder already sees at this graph scale. T-topo significantly less than T means the concatenation is destabilising the readout, which is unlikely at this scale but reportable.
- **Compute budget.** Per-relation graph magnitude on a Reactome-scale SBGN-PD graph is O(V^3) via a direct linear solve of `Z x = 1` with `Z_ij = exp(-d(i, j))`. Iterative and low-rank approximation methods exist (see for instance Gimperlein and Goffeng on magnitude of compact spaces) and are out of scope for this note. If the compute cost is prohibitive on the full Reactome graph, the fingerprint is computed on connected components of curated pathway modules and the module-level fingerprints are aggregated. That aggregation is a real design choice and would need to be pinned in the separate preregistration.
- **Fingerprint dimensionality caveat.** A 7-scalar per-relation channel is deliberately minimal. It may be too narrow to move a 128-dim R-GCN embedding meaningfully. A magnitude function evaluated at multiple scales `t_1, ..., t_k` per relation gives `7k` scalars and would be the natural upgrade if the minimal channel underperforms.
- **Preregistration status.** This memo is a design note, not a preregistration. If the T-topo comparison were run, a separate `PREREGISTRATION-topo.md` following the same template as the existing preregistration would be filed before any topo run.

**Out of scope for this memo.**

- Persistent path homology on the SBGN-PD graph would be a separate, larger design study, because the choice of filtration (which is not trivial for directed graphs) has to be justified independently.
- Directed flag complexes at the process-node level require the simplicial representation from Section 2, which changes the substrate the R-GCN sees and is therefore not comparable under a fixed-encoder protocol.
- Signed persistent homology for coherent-vs-incoherent feedforward-loop distinction is a research problem in its own right, not an engineering extension.

---

## 5. What this memo is not

- It is not a preregistration. Any experimental commitment lives in a separate PREREGISTRATION file.
- It is not a claim that I have implemented any of these methods. I have not. Every construction cited above is from a paper I have read.
- It is not a critique of the Ballester and Rieck 2023 paper. That paper's undirected VR-PH vs 1-WL result is the theoretical starting point I care about. This memo is about what happens when the same question is taken onto directed typed data.

---

## References

- Ballester, Rieck. *On the Expressivity of Persistent Homology in Graph Learning*. Third Learning on Graphs Conference (LoG), 2023. [arXiv:2302.09826](https://arxiv.org/abs/2302.09826).
- Barceló, Galkin, Morris, Orth. *Weisfeiler and Leman Go Relational*. Learning on Graphs Conference (LoG), 2022. [arXiv:2211.17113](https://arxiv.org/abs/2211.17113).
- Chowdhury, Mémoli. *Persistent Path Homology of Directed Networks*. Symposium on Discrete Algorithms (SODA), 2018.
- Grigor'yan, Lin, Muranov, Yau. *Homologies of path complexes and digraphs*. arXiv:1207.2834, 2012.
- Lütgehetmann, Govc, Levi, Hess. *Computing persistent homology of directed flag complexes*. Algorithms, 2020.
- Hepworth, Willerton. *Categorifying the magnitude of a graph*. Homology, Homotopy and Applications, 2017.
- Leinster. *The magnitude of metric spaces*. Documenta Mathematica, 2013.
- Schlichtkrull, Kipf, Bloem, van den Berg, Titov, Welling. *Modeling Relational Data with Graph Convolutional Networks*. Extended Semantic Web Conference (ESWC), 2018.
- Yang, Yih, He, Gao, Deng. *Embedding Entities and Relations for Learning and Inference in Knowledge Bases*. International Conference on Learning Representations (ICLR), 2015.
