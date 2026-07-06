# Coherent versus incoherent feedforward loops under three operator regimes

**Author.** Noor Muhammad ([github.com/NoorMuhammad1](https://github.com/NoorMuhammad1))
**Date.** 2026-07-06
**Companion artifacts.** [`sketches/topological_followup/DESIGN_NOTE.md`](../topological_followup/DESIGN_NOTE.md) (design memo on directed topological descriptors), [`sketches/rgcn_link_prediction/PREREGISTRATION.md`](../rgcn_link_prediction/PREREGISTRATION.md) (typed-versus-flat protocol).
**Status.** Working technical note. Not a preregistration. Committed publicly before any of the operators discussed have been implemented by me. Every operator claim below rests on a paper reading, not on code I have run.

---

## Why this note exists

The design memo on directed topological descriptors argues that the coherent-versus-incoherent feedforward loop is the natural test case for a directed higher-order descriptor on SBGN-PD graphs. This note takes one step down. It writes the operator equations for three regimes side by side (relation-typed R-GCN, sign-aware simplicial neural network, directed semi-simplicial upper Laplacian) and asks, at each level, whether the two configurations remain distinguishable at the graph-signal readout. The point is to isolate where the coherent-versus-incoherent distinction lives, and where it is averaged away.

I have not implemented any of these operators. The note reads as a specification of the empirical study I want to run, with the invariance argument for each regime written down carefully enough that a supervisor can push back on the claims. Where the claim rests on a paper's own proof (Ebli et al. 2020 on sign-aware SNNs, Lecha et al. 2025 on directed semi-simplicial learning), I cite and quote rather than reproduce.

---

## 1. Setup. The two feedforward loops

The feedforward loop is a three-node directed motif from Alon's transcription-network taxonomy (Mangan and Alon 2003 PNAS, and Alon 2007 *Nature Reviews Genetics*). Nodes are `s` (source), `i` (intermediate), and `t` (target). Edges are `s → i`, `s → t`, `i → t`. Each edge carries a sign in {+1, -1} indicating activation or inhibition.

**Coherent type-1 FFL (C1-FFL).** All three edges positive. The two paths from `s` to `t` (direct `s → t` and indirect `s → i → t`) both activate `t`.

**Incoherent type-1 FFL (I1-FFL).** `s → i` positive, `s → t` positive, `i → t` negative. The direct path activates `t`. The indirect path inhibits `t`. The two paths disagree at the sink.

The two configurations share the same underlying directed graph (same node set, same three edges with the same orientation). They differ on a single edge sign, namely the sign of `i → t`. Every operator regime below sees an identical unlabelled directed graph and receives the edge signs as either a per-edge scalar attribute in {+1, -1} or a per-edge relation-type label in {`positively_regulates`, `negatively_regulates`}. The two encoding choices matter and will be flagged explicitly at each step.

**Why this pair matters.** Under a downstream biological task like edge-sign prediction or motif-role classification, C1 and I1 have different phenotypes. C1 acts as a sign-sensitive delay element on the ON step (Mangan and Alon 2003, fig. 2). I1 acts as a pulse-generator when the direct arm activates before the indirect arm represses (Mangan and Alon 2003 shows the mechanism separately from fig. 2). A representation that maps C1 and I1 to the same embedding under a permutation-equivariant readout has thrown away the biological signal at the substrate level.

## 2. Reified representation

Both configurations, as emitted by the SBGN-PD BioCypher adapter, become bipartite subgraphs with a process node reifying the edge and typed participant edges attaching to it. The C1 configuration.

- `has_input(P_si, s)`, `has_output(P_si, i)`, both process edges with a positive sign
- `has_input(P_st, s)`, `has_output(P_st, t)`, positive
- `has_input(P_it, i)`, `has_output(P_it, t)`, positive

The I1 configuration is identical except the last process (`P_it`) carries a negative sign on its `has_output` edge, or equivalently is emitted with the relation label `negatively_regulates(i, t)` where C1 has `positively_regulates(i, t)`.

For the operator arguments below I write the sign as a per-edge scalar attribute `σ_e ∈ {+1, -1}` on the participant-projected graph (each process contracted to induce a directed participant-to-participant edge). The two-encoding-choices question (scalar attribute versus separate relation type) is exactly the necessary-stimulation collapse question flagged in the topological design memo, where SBGN-PD's seven Level 1 arc classes may collapse onto six Biolink predicates depending on whether typed relations subsume signed scalars, and is left open here.

## 3. Regime A. Relation-typed R-GCN with mean readout

**Operator.** A two-layer R-GCN with basis decomposition (Schlichtkrull et al. 2018) computes node embeddings as

    h_v^(l+1) = phi( sum_{r in R} sum_{u in N_r(v)} (1 / c_{v,r}) W_r^(l) h_u^(l) + W_self^(l) h_v^(l) )

where `N_r(v)` is the neighbourhood of `v` under relation `r`, `c_{v,r}` is a normaliser, `phi` is the elementwise nonlinearity, and `W_r^(l)` are per-relation weight matrices factorised through basis decomposition (Schlichtkrull et al. 2018 eq. 3).

Under the signed-scalar encoding, we treat the edge sign as a scalar multiplier on the message. A signed R-GCN then computes

    h_v^(l+1) = phi( sum_{r in R} sum_{u in N_r(v)} (1 / c_{v,r}) sigma_{uv,r} W_r^(l) h_u^(l) + W_self^(l) h_v^(l) )

with `sigma_{uv,r} ∈ {+1, -1}` the sign of the edge from `u` to `v` under `r` (a scalar coefficient distinct from the elementwise nonlinearity `phi`).

**Claim.** On the C1 versus I1 pair, all messages are identical except the message on `i → t`, which enters with `sigma_{it} = +1` under C1 and `sigma_{it} = -1` under I1. At the first layer, before the elementwise nonlinearity `phi`, the target-node pre-activation differs by

    Delta h_t^(1)_pre = 2 (1 / c_{t,r}) W_r^(0) h_i^(0)

a factor-of-2 sign-flipped vector. Under the participant-projected graph, `i` is upstream of `t` and receives no message from `t`, so the other nodes' first-layer embeddings are identical between C1 and I1. At the second layer, that first-layer perturbation on `h_t` propagates outward through `W_self` and any outgoing edges from `t`. In the three-node FFL, `t` has no outgoing edges, so the second-layer perturbation is confined to `h_t^(2)` and to whatever aggregates it.

**Under a graph-level mean readout** across `{s, i, t}` at layer 1, the C1-minus-I1 embedding difference is `(0 + 0 + 2 (1 / c_{t,r}) W_r^(0) h_i^(0)) / 3 = (2 / (3 c_{t,r})) W_r^(0) h_i^(0)`. This is a non-zero sign-flipped vector. It does not average to zero. **Under this construction, the two configurations appear distinguishable at the graph-level readout under a signed R-GCN with mean readout, subject to the caveats below.**

**Where the argument narrows.** If the readout is an equivariant nonlinearity followed by mean pooling, and if `h_i^(0)` and `W_r^(0)` are such that the sign-flipped contribution lives in a symmetry direction of the nonlinearity, the two graph embeddings could align. That happens for specific parameter settings, not generically. On a downstream task like edge-sign prediction, the representation gap should be preserved.

**Where the argument breaks.** Under the separate-relation-type encoding (`positively_regulates` vs `negatively_regulates` as distinct relations, distinct weight matrices), the two configurations receive completely different message-passing computations from the first layer. They are distinguishable, but the distinction is entirely absorbed by the two weight matrices' independence. There is no clean sign-flip statement.

**Reading.** For a standard R-GCN over the reified SBGN-PD graph, coherent and incoherent feedforward loops appear distinguishable at the graph-level readout under the argument above. The failure mode this note originally worried about (embeddings identical up to a sign flip absorbed by a symmetric readout) requires a specific untyped-flat encoder configuration that a preregistered R-GCN would not use. This is a concession relative to the design memo's earlier framing and is corrected here.

## 4. Regime B. Sign-aware simplicial neural network

**Operator.** Ebli, Defferrard, Spreemann 2020 (*Simplicial Neural Networks*, NeurIPS 2020 Workshop on Topological Data Analysis and Beyond) introduced the SNN as Hodge-Laplacian polynomial operators on oriented simplicial complexes. Bodnar, Frasca, Wang, Otter, Montúfar, Liò, Bronstein 2021 (*Weisfeiler and Lehman Go Topological. Message Passing Simplicial Networks*, ICML 2021) generalised the construction to message passing on faces and cofaces (MPSN) and formalised the orientation-equivariance property that the earlier Ebli operators possess implicitly. Bodnar's Section 4 gives the informal Theorem 14 (an MPSN layer with orientation-multiplied messages and odd nonlinearities is orientation equivariant), with the formal treatment in Appendix D (Def. 36 for the orientation transformation `T H = (T_0 H_0, T_1 H_1, ...)`, Prop. 39 for the equivariance of message-passing under `T`, and Lemma 40 for the invariance of orientation-invariant readouts). Bodnar's own Appendix C notes that "the orientation equivariance properties of the convolutional operators from Ebli et al. (2020) and Bunch et al. (2020) were not considered" when those papers were introduced, so the equivariance framing here belongs to Bodnar rather than to Ebli.

**FFL as a 2-simplex.** The three-node FFL sits as a 2-simplex `{s, i, t}` (the flag complex fills every clique in the underlying undirected graph as a simplex) with 1-faces `{s, i}`, `{s, t}`, `{i, t}`. The 2-simplex has two possible orientations, `(s, i, t)` and `(s, t, i)`, differing by a sign in the chain complex. The 1-face `{i, t}` has two possible orientations, `i → t` and `t → i`, again differing by a sign.

**Claim (hedged).** A signed simplicial neural network in the MPSN formulation appears unable to separate C1 from I1, provided the regulatory sign is encoded in the same channel as the chain-complex orientation. The argument runs as follows. Bodnar's Def. 36 defines the orientation transformation `T` acting on both the boundary matrices and the feature vectors, so flipping the orientation of edge `{i, t}` flips the sign of the entire feature vector on that edge. If the regulatory sign in {+1, -1} lives in that feature channel (a single-channel encoding where orientation and regulatory sign share the same scalar), then C1 and I1 fall into the same orbit of `T`, and Bodnar's Lemma 40 tells us an orientation-invariant readout absorbs the difference.

**Where the claim narrows.** Under a two-channel encoding (orientation in channel 1, regulatory sign in channel 2, treated as an independent per-simplex attribute), Bodnar's `T` acts only on channel 1. Channel 2 carries `+1` under C1 and `-1` under I1 unchanged, and a permutation-invariant readout does distinguish them. The claim above therefore holds under the specific encoding where regulatory sign is not carried as an independent feature channel. Under the more natural two-channel encoding, sign-aware SNNs distinguish C1 from I1 too, and Regime B collapses into Regime A.

**Reading.** Under a sign-aware SNN with single-channel sign-in-orientation encoding, C1 and I1 appear to sit in one equivalence class of the simplex-orientation sign group and are not distinguishable at the graph-level readout under this reading. Under two-channel encoding they are distinguishable. Which encoding is "natural" is itself a design decision the preregistered study would need to fix. This is the failure mode the design memo argues against, and it is the one the letter to Isufi highlights.

**Caveat I want to state honestly.** I have not derived this claim from first principles for the specific FFL case. I have read Bodnar's Def. 36, Prop. 39, Lemma 40, and Thm. 14, and inferred that they apply to the FFL orientation sign under the single-channel encoding. The two-channel-encoding sidestep above was pointed out to me in review and is the right hedge to state up front. A supervisor's push-back on this claim is entirely fair. The right way to nail it down is to write the C1 and I1 chain-complex generators explicitly under both encodings and check that they lie in the same or different orbits of the sign group acting on `C_1(K)` (1-chains on the flag complex).

## 5. Regime C. Directed semi-simplicial learning

**Operator.** Lecha, Cavallo, Dominici, Levi, Del Bue, Isufi, Morerio, and Battiloro 2025 (*Directed Semi-Simplicial Learning with Applications to Brain Activity Decoding*, arXiv:2505.17939) define learning on a semi-simplicial set `(S_n)` with face maps `d_i` satisfying the semi-simplicial identity `d_i d_j = d_{j-1} d_i` for `i < j`. Their directed simplicial complexes have `n`-simplices as ordered tuples `(v_0, ..., v_n)` with a directed edge for every `(v_i, v_j)`, `i < j`. Crucially, their Section 2 (last paragraph, page 4) introduces **attributed semi-simplicial sets**, `S_F = (S, F)`, where `F` assigns each simplex `sigma` a feature vector `x_sigma`. Per-simplex feature attributes are first-class citizens of the construction. SSN layers (their eq. 1) update those features through operators `omega_R` that respect the ordered face maps but do not act on the attribute channel through any sign symmetry.

**FFL as a directed 2-semi-simplex.** The directed 2-semi-simplex `(s, i, t)` requires all three ordered edges (`s → i`, `s → t`, `i → t`) to exist. That is the FFL. The two configurations C1 and I1 both realise the same directed 2-semi-simplex. They differ on the *sign attribute* of the 1-face `(i, t)`. That attribute is a separate per-face scalar in the edge-feature channel `x_{(i,t)}`, not the face-map orientation. The semi-simplicial machinery does not identify the flipped face `(t, i)` with `(i, t)` in the first place (ordered face maps are asymmetric), and it does not act on the feature channel through any sign group.

**Claim (hedged).** Under a directed semi-simplicial upper Laplacian in the Lecha et al. 2025 formulation, C1 and I1 look likely to produce distinguishable simplex features because the sign attribute on `(i, t)` enters the operator's edge-feature channel unchanged. Whether the resulting graph-level readout distinguishes them depends on whether the readout is permutation-equivariant across simplices (in which case the distinction is preserved by the different sign attributes on the same face) or averages across a chain-complex sign group (which does not apply here because the face maps are ordered, not signed).

Lecha et al. 2025 Theorem 2 states that SSNs are strictly more expressive than MPSN at distinguishing non-isomorphic directed simplicial complexes. C1 and I1 with a per-face sign attribute in {+1, -1} are the smallest interesting instance of that expressivity gap on a biological substrate.

**Reading.** The directed semi-simplicial regime is the natural learning-side complement to the schema-side reification the adapter already commits to. Semi-simplicial sets replace the signed orientation with ordered face maps and remove the equivalence class that a signed SNN identifies. Whether it also preserves the sign attribute on `(i, t)` in a way that a downstream classifier can pick up is the open empirical question this note flags for a PhD.

**Caveat.** The argument for regime C rests on the Lecha 2025 paper's construction as I have read it, especially the attributed-semi-simplicial-set definition in their Section 2 and the expressivity theorem (Theorem 2). I have not verified the claim by writing the upper Laplacian on the FFL directly. Doing so is the natural first exercise inside a PhD in Isufi's group.

## 6. Where the pilots do and do not test this

Pilot v2 in [`sketches/rgcn_link_prediction/PILOT_V2_RESULTS.md`](../rgcn_link_prediction/PILOT_V2_RESULTS.md) added a T+topo arm (which I would now relabel T+structural in retrospect because the log-degree signature is not a topological invariant in the Leinster sense, and the letter to Isufi uses the relabelled form) carrying a per-node per-relation log-degree signature. That arm is not a test of the regime A vs B vs C question above. The log-degree signature is a first-order structural statistic and lives at the R-GCN readout level, not at the simplicial or semi-simplicial level. Node degrees in `{s, i, t}` are identical between C1 and I1 under one-relation-typed encoding, so log-degree signatures cannot separate the pair at the motif level. Pilot v2's finding (F-large wins on the synthetic substrate) is informative about the synthetic-substrate signal, not about coherent-versus-incoherent separability. Per the pilot's own hedge, "the degree signature is a cheap surrogate" and does not test whether persistent path homology or magnitude homology would behave the way the degree signature does.

The honest test would run three encoders side by side on a curated pathway corpus with FFL motifs identified (Reactome or a subset), with an edge-sign or motif-role prediction task, and pre-committed effect-size priors for each pairwise comparison. The three encoders would be

- **Encoder A.** Signed R-GCN (Schlichtkrull et al. 2018 with basis decomposition), signed-scalar edge attribute encoding. The scalar-versus-two-relation-type encoding choice from Section 3 should be preregistered as an ablation, not a footnote.
- **Encoder B.** Sign-aware simplicial neural network (Bodnar et al. 2021 MPSN, formalisation) built on the flag complex over the participant-projected graph. Directed edges must be mapped to chain-complex orientations, which is non-trivial and is exactly the ambiguity Section 4 wobbles on. Both the single-channel and two-channel encodings from Section 4 should be preregistered arms.
- **Encoder C.** Directed semi-simplicial learning (Lecha et al. 2025) built on the semi-simplicial set over the same graph. Their released code at the arXiv link can be adapted.

At Reactome scale, coherent and incoherent FFLs are extractable from the curated pathway data. FFL enumeration is a subgraph-isomorphism problem against the three-node FFL template and runs in a few hours of code (NetworkX `subgraph_monomorphisms_iter` or a dedicated motif enumeration routine). Alon's TRANSFAC-derived motif catalogue intersects with Reactome's regulatory arcs on gene symbols. Absolute FFL counts on 10^4-node regulatory graphs are typically hundreds to low thousands, and the C1 to I1 ratio in *E. coli* is roughly 3 to 1 (Mangan and Alon 2003 Table 1). At Reactome scale for human, expect the same order. Class imbalance handling (stratified sampling or class weights) must be preregistered.

At synthetic scale the FFLs have to be planted in the graph generator, which changes the substrate question in ways the current pilot pipeline does not address.

The preregistered comparison would need its own protocol file, `PREREGISTRATION-ffl.md`, filed before any run.

## 7. What this memo is not

- It is not a preregistration. Any experimental commitment lives in a separate `PREREGISTRATION-ffl.md`.
- It is not a claim that I have implemented any of the three operators. I have not. Every operator in Sections 3, 4, and 5 rests on paper reading.
- It is not a critique of Ebli 2020, Bodnar 2021, or Lecha 2025. Those papers are the theoretical starting points. This memo is about what happens when the same question is taken onto directed typed pathway data.
- It is not a proof that C1 and I1 are distinguishable in regime C. It is a specification of the empirical test that would settle the question.

## References

- Alon. *Network Motifs. Theory and Experimental Approaches*. Nature Reviews Genetics 8, 450-461, 2007.
- Bodnar, Frasca, Wang, Otter, Montúfar, Liò, Bronstein. *Weisfeiler and Lehman Go Topological. Message Passing Simplicial Networks*. ICML 2021.
- Ebli, Defferrard, Spreemann. *Simplicial Neural Networks*. NeurIPS 2020 TDA Workshop.
- Lecha, Cavallo, Dominici, Levi, Del Bue, Isufi, Morerio, Battiloro. *Directed Semi-Simplicial Learning with Applications to Brain Activity Decoding*. arXiv:2505.17939, 2025.
- Mangan, Alon. *Structure and Function of the Feed-Forward Loop Network Motif*. Proceedings of the National Academy of Sciences 100(21), 11980-11985, 2003.
- Schlichtkrull, Kipf, Bloem, van den Berg, Titov, Welling. *Modeling Relational Data with Graph Convolutional Networks*. Extended Semantic Web Conference (ESWC), 2018.
- Yang, Yih, He, Gao, Deng. *Embedding Entities and Relations for Learning and Inference in Knowledge Bases*. International Conference on Learning Representations (ICLR), 2015.
