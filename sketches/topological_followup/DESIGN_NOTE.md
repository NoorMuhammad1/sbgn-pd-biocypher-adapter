# Design memo: Directed topological descriptors on reified SBGN-PD graphs

**Author.** Noor Muhammad ([github.com/NoorMuhammad1](https://github.com/NoorMuhammad1))
**Date.** 2026-07-05
**Companion artifacts.** [`sketches/rgcn_link_prediction/`](../rgcn_link_prediction/) (R-GCN sketch) and [`sketches/rgcn_link_prediction/PREREGISTRATION.md`](../rgcn_link_prediction/PREREGISTRATION.md) (typed-vs-flat protocol)
**Status.** Working design memo. Not a preregistration. Committed publicly before any of the experiments discussed have been run.

---

## Why this memo exists

The R-GCN sketch in `sketches/rgcn_link_prediction/` demonstrates that the SBGN-PD BioCypher adapter's output feeds a relation-aware graph neural network end-to-end. The accompanying preregistration compares a typed encoder against a flattened one on link prediction, under a Wilcoxon signed-rank protocol with a parameter-matched control. Alternative encoders (block-diagonal R-GCN, HetGNN, HGT, GAT) are held out of scope by design.

That scope decision leaves a natural next question. If R-GCN with basis decomposition is bounded by the multi-relational 1-WL test (Barcelo et al. 2022), and the Vietoris-Rips PH vs 1-WL result (Ballester and Rieck 2023) shows the two are incomparable rather than one dominating, then adding a topological descriptor as a second view should recover signal that neither view has alone. This memo works through which topological descriptor is the right one for SBGN-PD, and what the smallest useful preregistered follow-up would look like.

---

## 1. Which flavor of topology is right for SBGN-PD?

SBGN-PD graphs are **directed and typed**. Each process has ordered participants (inputs, outputs, catalysts, modulators), and the direction encodes reaction semantics that a downstream ML pipeline should not discard. The choice of descriptor family has to respect this.

Three options, arranged from least to most direction-aware.

**Vietoris-Rips persistent homology (VR-PH)** on the underlying graph with the shortest-path metric. Detects H_0 (connected components) and H_1 (undirected 1-cycles). Insensitive to edge orientation by construction. A feedback loop A→B→C→A and an undirected triangle produce the same VR-PH signature. This is the wrong tool for SBGN-PD if we care about direction.

**Persistent path homology (Chowdhury and Memoli 2018)** on the directed graph. Uses regular paths as chains and takes a boundary that respects orientation. Distinguishes directed cycles from undirected ones. A feedback loop appears as an H_1 class with a specific orientation signature. A linear waterfall A→B→C→D produces no H_1 class at all. This is the natural descriptor for feedback-loop detection.

**Persistent homology of directed flag complexes.** Build a simplicial complex where a p-simplex is a totally ordered clique with all edges pointing in the total-order direction. On SBGN-PD, a coherent feedforward loop (X→Y, X→Z, Y→Z) is a directed 2-simplex. A feedback loop is not a simplex at all. PH on the directed flag complex then separates these two motif classes cleanly.

**Concrete claim.** For SBGN-PD, the empirical study I want to run is not "PH vs 1-WL" but "path homology plus directed flag complex vs R-GCN." Undirected VR-PH is a Ballester-Rieck 2023 baseline for graph learning in general. On this substrate it is the wrong starting point.

**Open question I would want to open with.** For distinguishing coherent from incoherent feedforward loops on SBGN-PD (both are directed triangles, differing only in whether the modulation on Y→Z is activation or inhibition), the directed flag complex does not natively encode edge sign. Signed persistent homology or an edge-type-filtered directed flag complex are the two candidates I have read about. I have not implemented either.

---

## 2. What object should we compute the homology on?

Given an SBGN process P with participants {A, B, C, D} carrying roles {input, input, output, catalyst}, three representations are on the table.

**Reified graph.** P becomes a `BiologicalProcess` node. A, B, C, D become participant nodes. Four typed edges (`has_input`, `has_input`, `has_output`, `catalyzes`) connect P to the participants. This is what the adapter emits today. Every homology construction above operates natively on it. The cost is that the n-ary structure of P has been unfolded into 1-cells, so the "these participants co-occur in one reaction" invariant is not encoded at the level of the underlying complex.

**Hypergraph.** P becomes a hyperedge over {A, B, C, D}, labeled with the participant roles. The n-arity is preserved. Persistent homology on hypergraphs is well-defined via the associated hypergraph nerve or by lifting to a bipartite representation. Direction and role labeling are harder to encode natively.

**Simplicial complex.** P becomes a 3-simplex spanning {A, B, C, D} (with lower-dimensional faces filled automatically). Higher-order TDA methods operate natively on simplicial complexes and would see P as a filled tetrahedron. Directed variants (directed flag complexes from Section 1) can carry the participant ordering.

**Worked example.** Consider the hexokinase reaction from glycolysis. Glucose plus ATP produces glucose-6-phosphate plus ADP, catalyzed by hexokinase. Participants:

- Substrates: glucose, ATP
- Products: glucose-6-phosphate, ADP
- Catalyst: hexokinase

In the reified representation the adapter emits, this is a 6-node subgraph (5 participants plus 1 process node) with 5 typed edges. In the hypergraph representation, it is a 5-node hyperedge with role labels. In the simplicial representation, it is a 4-simplex over the 5 participants, with internal orientation carried by the process's role assignments.

**Concrete claim.** Which representation is best is an empirical question that has to be answered per downstream task. For **link prediction on relations that already exist in Biolink** (`has_input`, `has_output`, `catalyzes`, `regulates`), the reified graph is the natural substrate because the prediction target is literally the reified edges. For **motif-level or process-level classification** (does this reaction subgraph correspond to a known pathway motif?), the simplicial representation is likely a better substrate because the invariant of interest lives at the level of the higher-order simplex, not at the level of individual edges.

---

## 3. Magnitude and magnitude homology

**Graph magnitude (Leinster 2013).** For a finite metric space (V, d), define a matrix Z with Z_ij = exp(-d(i, j)). The magnitude is the sum of all entries of Z^{-1}. For graphs, d is usually the shortest-path metric. Magnitude gives a real-valued invariant that behaves like an "effective size" of the space.

**Magnitude homology (Hepworth and Willerton 2017).** A bigraded homology theory whose graded Euler characteristic recovers the magnitude power series. Adds structure that magnitude alone does not see. On graphs, magnitude homology encodes information about paths and their length distributions.

**What is known and what is open on typed and directed cases.**

- Magnitude of directed graphs uses the asymmetric shortest-path quasi-metric and gives a well-defined invariant.
- Magnitude of typed graphs (where each edge carries a categorical label from a small alphabet) does not have a clean off-the-shelf definition. Naive approaches assign a single edge weight per typed edge and lose the categorical information. Per-relation magnitude, meaning compute magnitude separately on each relation subgraph and concatenate the resulting scalars, is a defensible engineering choice that has to my knowledge not been benchmarked on curated pathway data.

**Concrete claim.** For SBGN-PD, a per-relation magnitude fingerprint concatenated across the six arc classes gives a 6-vector descriptor per graph. That descriptor is a candidate topological feature to concatenate with an R-GCN embedding at the graph-level readout. It is the smallest meaningful topological arm I can add to the preregistered T-vs-F comparison without changing the encoder.

---

## 4. Proposed follow-up to the existing preregistration

The current preregistration at `sketches/rgcn_link_prediction/PREREGISTRATION.md` compares a typed R-GCN (T) against a flattened R-GCN (F) on link prediction over the adapter output. Alternative encoders are out of scope by design.

**Extension proposal.** Add one topological-descriptor arm under the same Wilcoxon protocol.

- **T-topo.** T's R-GCN embedding concatenated with a per-relation graph-magnitude fingerprint (6 scalars, one per SBGN-PD arc class) at the graph-level readout, followed by the same DistMult decoder for link scoring.
- **Same protocol.** Ten seeds, stratified 80/10/10 arc split (frozen once), raw MRR primary metric, Wilcoxon signed-rank paired test at alpha 0.05, Cohen's d effect size.
- **Parameter matching.** T-topo has 6 additional scalars at the readout compared to T. F-large-topo would be the natural fourth arm to match capacity.
- **Pre-committed interpretation.** T-topo > T significantly → the magnitude fingerprint carries link-prediction signal beyond what typed message passing captures on its own. T-topo not significantly different from T → the magnitude fingerprint is redundant with what R-GCN's basis-decomposed encoder already sees at this graph scale. T-topo < T significantly → concatenation is destabilising the readout (unlikely at this scale but reportable).
- **Compute budget.** Per-relation magnitude on a Reactome-scale SBGN-PD graph is O(V^3) if computed naively via matrix inversion. Sparse approximations exist in the literature and are out of scope for this note.
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
- Barcelo, Galkin, Morris, Orth. *Weisfeiler and Leman Go Relational*. Learning on Graphs Conference (LoG), 2022. [arXiv:2211.17113](https://arxiv.org/abs/2211.17113).
- Chowdhury, Memoli. *Persistent Path Homology of Directed Networks*. Symposium on Discrete Algorithms (SODA), 2018.
- Hepworth, Willerton. *Categorifying the magnitude of a graph*. Homology, Homotopy and Applications, 2017.
- Leinster. *The magnitude of metric spaces*. Documenta Mathematica, 2013.
- Schlichtkrull, Kipf, Bloem, van den Berg, Titov, Welling. *Modeling Relational Data with Graph Convolutional Networks*. Extended Semantic Web Conference (ESWC), 2018.
- Yang, Yih, He, Gao, Deng. *Embedding Entities and Relations for Learning and Inference in Knowledge Bases*. International Conference on Learning Representations (ICLR), 2015.
