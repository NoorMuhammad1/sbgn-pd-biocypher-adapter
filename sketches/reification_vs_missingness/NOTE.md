# Two axes, not one. Reification and missing-data handling on typed biological graphs

**Author.** Noor Muhammad
**Date.** 2026-07-05 (v2 after 4-agent review 2026-07-06)
**Scope.** Companion note to the [SBGN-PD BioCypher adapter](../..) and to the [R-GCN preregistration](../rgcn_link_prediction/PREREGISTRATION.md). Written to disentangle two design axes that an earlier framing conflated, and that a geometric-DL reader would spot immediately.

---

## The confusion I want to clear

An earlier framing claimed that my adapter and the soft-manifolds work of Marinoni, Liò, Barp, Jutten, and Girolami (*IEEE TPAMI* 2026) both refuse to let structural information disappear into a training loss. That is right in spirit and wrong in detail. The two moves are not on the same axis. This note names what each axis actually is, where they meet, where they do not, and where the message-passing encoder can undo what the graph carefully preserved.

## Axis 1. Reification. A schema-level lossy-compression fix.

**Problem.** SBGN Process Description has several process-glyph subtypes (`process`, `omitted process`, `uncertain process`, `association`, `dissociation`, `phenotype`). The generic `process` glyph is **n-ary** in the sense that matters here. A single `process` glyph can carry an arbitrary number of substrates and products through consumption and production arcs, plus modulators (with catalysis, stimulation, inhibition, and necessary stimulation as modulator subtypes). Stoichiometry coefficients attach to consumption and production arcs, not to participant nodes. `association` and `dissociation` are fixed structural arities and `phenotype` is effectively unary. Those cases are out of scope of the argument below. Biolink's surface predicate layer is binary. Predicates like `biolink:has_input` and `biolink:catalyzes` take one subject and one object. Biolink does support reified associations through its `Association` class hierarchy. The naive lowering of an SBGN process into a set of binary predicates (rather than an `Association` instance) drops one or both of two things.

1. The **co-occurrence** invariant. Three participants in one process are not the same as three pairwise-linked participants in three separate processes with a shared enzyme.
2. The **role-and-coefficient** payload. Stoichiometry lives on the process-participant arc, not on the participant.

**My fix.** Reify each process as a `BiologicalProcess` node with typed `has_input`, `has_output`, `catalyzes`, and `regulates` edges. The stoichiometry coefficient goes on the edge as a property.

**What this fix does at the graph level.** It preserves the n-ary structure of the generic `process` glyph as a two-hop pattern through a reification node. Co-occurrence becomes reachability through that node. Roles become edge types. Coefficients become edge properties.

**What this fix does not do.**

- **Encoder recovery is not automatic.** Preserving the invariant in the graph is a necessary but not sufficient condition. A per-relation mean aggregation (default R-GCN) still collapses the participant set into a bag of pairwise contributions at each layer. Barcelo et al. 2022 shows R-GCN and CompGCN are bounded by the multi-relational 1-WL test. Reification widens the class of graphs the encoder sees but does not raise its expressive ceiling. Set-symmetric readouts at the reification node (a set transformer over the neighbourhood, or a hyperedge-aware encoder) are what actually close the loop.
- **Absence is untouched.** Reification is silent about annotation absence. If a substrate is uncurated in the source, the reified graph has one fewer edge, and the training loop cannot tell whether the edge is absent for biological reasons or for curatorial ones. SBGN-PD does mark process-level uncertainty at the glyph type (`uncertain process`, `omitted process`), and the adapter can carry those through as node-level type distinctions, but that is a partial structural handling of uncertainty, not a treatment of edge-level missingness.

## Axis 2. Missing-data handling. An embedding-space geometry fix.

**Problem.** Curated biomedical graphs have **uneven coverage**. Some pathways are heavily annotated. Others are stubs. Some entities appear frequently, others once. A standard Euclidean graph embedding assumes uniform information geometry across the space. Under uneven coverage that assumption fails, and the embedding tends to pull sparsely-annotated regions toward the mean.

**Marinoni et al. 2026 fix.** Represent the embedding space as a **soft manifold** whose tangent spaces at each point are **hypocycloids** rather than flat planes. The hypocycloid shape is set by the local **velocity of information propagation** across the data points (the paper's headline framing in the abstract). Where information propagates slowly (sparse annotation), the local tangent geometry deforms so that distances stretch and the sparsely-annotated node stops being dragged into the mean. The paper's technical treatment also describes the resulting per-node quantities as *conductivity*, *diffusivity*, and *diffusion rate* of an induced material, which change with feature availability. The two framings (information-propagation and material-diffusion) are the paper's own, describing the same construction from two angles.

**What this fix does at the embedding level.** It absorbs uneven annotation into the geometry itself. Missingness becomes a first-class property of the space rather than a hole to be filled by imputation or an assumption to be papered over by regularisation.

**What this fix does not do.** It does not handle schema-level lossy compression of the input graph. If the graph fed to the soft manifold was produced by naive n-ary-to-binary flattening, the manifold has no way to recover the co-occurrence and role information that flattening threw away.

## Where the two axes meet, and where they do not

The two moves happen at **different pipeline stages**.

```
   raw source   ->  schema translation  ->  graph  ->  encoder  ->  embedding  ->  downstream
   (SBGN, BioPAX,      (Axis 1:            (typed)   (R-GCN,    (Axis 2:        (link pred,
    curated dumps)      reification)                  set-      soft manifold)   drug response)
                                                     transformer,
                                                     etc.)
```

Axis 1 (reification) operates at schema translation. Axis 2 (soft manifolds) operates at embedding. Fixing one does not fix the other. Skipping one leaks structural information downstream of it, whether or not the other is applied.

They **do** share an instinct. Both refuse to let the training objective absorb information that should have been kept as structure.

They **do not** interact by any information-theoretic guarantee I can prove. The claim I want to test is weaker: the two losses are **non-overlapping in pipeline stage**, so composing the two fixes should recover more signal than either alone. Nothing in this note proves that. Prediction 3 below tests it directly.

## Arm definitions and mapping to the existing preregistration

The four-arm study below uses the following definitions. **Flat** here means the generic n-ary `process` glyph is lowered to the set of pairwise `participant--participant` edges implied by the reaction without a reification node. This is stricter than the *F* encoder in the [R-GCN preregistration](../rgcn_link_prediction/PREREGISTRATION.md), which reifies but collapses all arc types to one relation. **Reified** matches the adapter's current output. **Euclidean** means a standard vector-space embedding as in vanilla R-GCN. **Soft-manifold** means the embedding space is the Marinoni et al. 2026 construction, with the per-node hypocycloid tangent geometry applied to the final node embedding before scoring.

- **Arm A. Flat + Euclidean.** Pairwise n-ary lowering, no reification node, standard Euclidean R-GCN embedding.
- **Arm B. Reified + Euclidean.** Reification (the adapter's current move), standard Euclidean R-GCN embedding. This is the closest cousin of *T* in the existing preregistration. *T* keeps typed arcs while B could in principle collapse them, so in the four-arm study B is defined to keep typed arcs.
- **Arm C. Flat + soft-manifold.** Pairwise n-ary lowering, soft-manifold embedding at readout.
- **Arm D. Reified + soft-manifold.** Reification and soft-manifold embedding.

Under this mapping the existing preregistration's F sits between A and B (F reifies but flattens typing), and the existing preregistration's T sits at B with typing preserved. Arms C and D require the soft-manifold construction, which no existing preregistration in this repo covers.

The four arms form a 2 (reified vs flat) x 2 (soft-manifold vs Euclidean) factorial. Parameter matching follows the preregistration's Condition B protocol. Each arm's hidden dimension is grown until parameter counts match within 5%. The four-arm test statistic is a 2x2 factorial analysis on per-seed MRR, with main effects for reification and for soft-manifold, and an interaction term.

## Concrete testable predictions

The composition claim commits to four predictions.

1. **B > A**, controlling for parameter count. Reification alone recovers signal that pairwise flattening loses at the schema stage.
2. **C > A**, controlling for parameter count. Soft-manifold embedding alone recovers signal that Euclidean smoothing loses at the embedding stage.
3. **D > B and D > C**, controlling for parameter count. The two fixes address non-overlapping information losses and compose. Expected effect-size prior is Cohen's d in the 0.3 to 0.8 range on log-MRR difference against the stronger of B or C. Anything below 0.2 counts as "does not compose."
4. Under the 2x2 factorial, both main effects (reification and soft-manifold) are positive and the interaction is not significantly negative. A significantly negative interaction would refute the pipeline-locality framing even if D is nominally the best arm.

**Sub-case falsifier.** If **D is not significantly better than max(B, C)**, the composition claim is refuted. Three sub-cases carry different biological readings.

- **(a) D approximately equal to C, and C > B.** The soft-manifold subsumes reification on this substrate. The information-geometry-aware embedding recovers the schema-level information that reification also recovers, without needing the reification step.
- **(b) D approximately equal to B, and B > C.** Reification subsumes the soft manifold. The schema-level fix already recovers what the geometric fix would have added.
- **(c) D approximately equal to B approximately equal to C, all greater than A.** Both fixes address the same information loss on this substrate, and the "two axes" framing is wrong. This is the strongest refutation.

## Nearby data point already on record

The [PILOT_V2 report](../rgcn_link_prediction/PILOT_V2_RESULTS.md) already ran a cheap proxy for composition. Specifically it added a T+topo arm where the typed R-GCN encoder was augmented with a per-node per-relation degree signature at readout. That is not the same as Arm D above (a degree signature is not the soft-manifold construction), but it is a nearby composition test on the same synthetic substrate that Arms B and D would run on if the study were executed at pilot scale.

The T+topo arm produced a suggestive positive direction over the plain T arm (Cohen's d equal to +0.49, one-sided Wilcoxon p equal to 0.116, n equal to 10), which failed to reach significance. On the same substrate, F-large stayed significantly ahead of T and of T+topo. This tells us two things.

1. On synthetic preferential-attachment data with no genuine relation-typed signal, adding a small topological augmentation on top of typing does not produce a large easy effect. The composition signal exists but is small at this scale.
2. The result is not evidence for or against Prediction 3 on the actual soft-manifold axis. The proxy is a cousin, not the same construction.

The Reactome-scale run is what will actually test the four predictions above. That run is scheduled pending the Cloudflare-gated Reactome ContentService SBGN exporter access noted in the [pilot deviations](../rgcn_link_prediction/PILOT_V2_RESULTS.md).

## Missing considerations worth flagging

Three things this note does not attempt but a topological-ML reader would raise.

- **Hypergraph representation as an alternative to reification.** A hypergraph lift preserves n-ary structure natively without introducing an intermediate reification node. Persistent homology and other topological methods handle hypergraphs directly. Whether reification-plus-R-GCN or hypergraph-plus-hypergraph-GNN is the better substrate is an empirical question this note does not address. The [topological design memo](../topological_followup/DESIGN_NOTE.md) discusses the trade-off from the topology side.
- **Where in the R-GCN pipeline the soft-manifold enters.** Marinoni et al. define the geometry on the embedding space. In an R-GCN that produces per-node embeddings across layers, the soft-manifold construction could be applied at every layer, only at the final layer, or only at the readout that feeds the DistMult decoder. Different placements are different arms and would need to be specified in a preregistration before Arms C and D are runnable.
- **Interaction with persistent path homology.** [DESIGN_NOTE.md](../topological_followup/DESIGN_NOTE.md) proposes persistent path homology and directed flag complexes as a third axis of intervention. A full study would extend the 2x2 factorial to include topological descriptors, but a 2x2x2 factorial with three axes needs at least three times as many seeds for the interaction terms to be estimable and is out of scope of this note.

## Scope of what this note is not

- It is not a preregistration. The four-arm study above would need its own `PREREGISTRATION-4arm.md` filing before it is run.
- It is not a claim that I have implemented the soft-manifold arm. I have not. The soft-manifold construction is one of the geometric methods listed in the [topological design memo](../topological_followup/DESIGN_NOTE.md) as "read but not shipped in."
- It is not a critique of Marinoni et al. 2026. It is an attempt to be precise about what problem that paper solves and what problem the adapter's reification solves, so a reader who cares about both does not have to correct my framing.

## References

- Marinoni, Liò, Barp, Jutten, Girolami. *Improving Embedding of Graphs With Missing Data by Soft Manifolds*. IEEE Transactions on Pattern Analysis and Machine Intelligence 48(3), 2221-2235, 2026.
- Schlichtkrull, Kipf, Bloem, van den Berg, Titov, Welling. *Modeling Relational Data with Graph Convolutional Networks*. Extended Semantic Web Conference, 2018.
- Yang, Yih, He, Gao, Deng. *Embedding Entities and Relations for Learning and Inference in Knowledge Bases*. International Conference on Learning Representations, 2015.
- Barcelo, Galkin, Morris, Orth. *Weisfeiler and Leman Go Relational*. Learning on Graphs Conference, 2022.
- Chowdhury, Memoli. *Persistent Path Homology of Directed Networks*. Symposium on Discrete Algorithms, 2018.
- Lobentanzer et al. *Democratizing knowledge representation with BioCypher*. Nature Biotechnology 41, 1056-1059, 2023.
