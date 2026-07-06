# Two axes, not one: reification and missing-data handling on typed biological graphs

**Author.** Noor Muhammad
**Date.** 2026-07-05
**Scope.** Companion note to the [SBGN-PD BioCypher adapter](../..) and to the [R-GCN preregistration](../rgcn_link_prediction/PREREGISTRATION.md). Written to disentangle two design axes that the cover letter I sent to Prof. Pietro Lio conflated in a first pass, and that a Cambridge geometric-DL reader would spot immediately.

---

## The confusion I want to clear

In the letter, I claimed that my adapter and the soft-manifolds work of Marinoni, Liò, Barp, Jutten, and Girolami (*IEEE TPAMI* 2026) both refuse to let structural information disappear into a training loss. The framing is right in spirit and wrong in detail. The two moves are not on the same axis. This note names what each axis actually is, where they meet, and where they do not.

## Axis 1. Reification. A schema-level lossy-compression fix.

**Problem.** SBGN Process Description encodes reactions as **n-ary process glyphs**. A single glyph can carry an arbitrary number of substrates, products, catalysts, and modulators, each with a role and a stoichiometry coefficient. Biolink's predicate layer is **mostly binary**. Naive lowering of an n-ary process into a set of binary Biolink predicates loses one or both of two things.

1. The **co-occurrence** invariant. Two substrates that appear together in the same process are not the same as two substrates that appear separately in two processes involving the same enzyme.
2. The **role-and-coefficient** payload. Stoichiometry lives on the process-participant relationship, not on the participant node.

**My fix.** Reify each process as a `BiologicalProcess` node with typed `has_input`, `has_output`, `catalyzes`, and `regulates` edges. The stoichiometry coefficient and any modifier semantics go on the edge as a property. Every participant that shared the process is now connected to the same `BiologicalProcess` node.

**What this fix does.** It preserves the n-ary structure as a two-hop pattern through a reification node. Co-occurrence survives as reachability through that node. Roles survive as edge types. Coefficients survive as edge properties.

**What this fix does not do.** It does not handle **absence**. If a substrate is missing from a reaction because the source database has not curated it, the reification is silent about that. The graph has one fewer edge, and the training loop cannot tell whether the edge is absent because the reaction genuinely does not involve that substrate or because the curator has not annotated it yet.

## Axis 2. Missing-data handling. An embedding-space geometry fix.

**Problem.** Curated biomedical graphs have **uneven coverage**. Some pathways are heavily annotated. Others are stubs. Some entities appear frequently. Others appear once. A standard Euclidean graph embedding assumes uniform information geometry across the space. Under uneven coverage that assumption fails, and the embedding pulls sparsely-annotated regions toward the mean.

**Marinoni et al. 2026 fix.** Represent the embedding space as a **soft manifold** whose tangent geometry varies with feature availability. Tangent spaces are hypocycloids rather than flat planes. Local material properties (conductivity, diffusivity, diffusion rate) are per-node and are a function of how much information is present at that node. Where information is scarce, the manifold locally deforms so that distances stretch and the sparsely-annotated node stops being dragged into the mean.

**What this fix does.** It absorbs uneven annotation into the embedding geometry itself. Missingness becomes a first-class property of the space rather than a hole to be filled by imputation or an assumption to be papered over by regularisation.

**What this fix does not do.** It does not handle **schema-level lossy compression** of the input graph. If the graph fed to the soft manifold was produced by naive n-ary-to-binary flattening, the manifold cannot recover the co-occurrence and role information that flattening threw away.

## Where the two axes meet, and where they do not

The two moves happen at **different pipeline stages**.

```
   raw source  ->  schema translation  ->  graph  ->  embedding  ->  downstream
   (SBGN, BioPAX,      (Axis 1)            (typed)    (Axis 2)      (link pred,
    curated dumps)                                                   drug response)
```

Axis 1 (reification) operates at schema translation. Axis 2 (soft manifolds) operates at embedding. Fixing one does not fix the other. Skipping one leaks structural information downstream of it, whether or not the other is applied.

They **do** share an instinct. Both refuse to let the training objective absorb information that should have been kept as structure. That is the family resemblance I meant to signal in the cover letter. It is a real resemblance but a shallow one.

They **do not** interact commutatively. Reifying and then embedding on a soft manifold is not the same as reifying only, and it is also not the same as embedding on a soft manifold over a flattened graph. The two fixes address orthogonal information losses. Composing them should recover more signal than either alone. Nothing in this note proves that. It is a testable prediction.

## Concrete testable predictions

Under the [preregistered protocol](../rgcn_link_prediction/PREREGISTRATION.md), four experimental arms would separate the two axes cleanly.

- **Arm A. Flat + Euclidean.** Naive n-ary flattening, standard Euclidean R-GCN embedding.
- **Arm B. Reified + Euclidean.** Reification (my adapter's move), standard Euclidean R-GCN embedding.
- **Arm C. Flat + soft-manifold.** Naive flattening, soft-manifold embedding.
- **Arm D. Reified + soft-manifold.** Reification, soft-manifold embedding.

The predictions the composition claim commits to.

1. **B > A**, controlling for parameter count. Reification alone recovers signal that flattening loses at the schema stage.
2. **C > A**, controlling for parameter count. Soft-manifold embedding alone recovers signal that Euclidean smoothing loses at the embedding stage.
3. **D > B and D > C**, controlling for parameter count. The two fixes address non-overlapping information losses and compose.
4. If **D is not significantly better than max(B, C)**, the two axes are not orthogonal in the way I have claimed, and one of the fixes is subsuming the other on this substrate. That is a real finding and would tell me my mental model needs revision.

Prediction 4 is what makes this note falsifiable. If it fails, the composition claim in this note is wrong.

## Scope of what this note is not

- It is not a preregistration. The four-arm study above would need its own PREREGISTRATION file before it is run, filed under the same discipline as the existing typed-vs-flat preregistration.
- It is not a claim that I have implemented the soft-manifold arm. I have not. The soft-manifold construction is one of the topological / geometric methods listed in the topological design memo as "read but not shipped in."
- It is not a critique of the Marinoni et al. 2026 paper. It is an attempt to be precise about what problem that paper solves and what problem my adapter solves, so that a Cambridge reader who cares about both does not have to correct my framing at interview.

## References

- Marinoni, Liò, Barp, Jutten, Girolami. *Improving Embedding of Graphs With Missing Data by Soft Manifolds*. IEEE Transactions on Pattern Analysis and Machine Intelligence 48(3), 2221-2235, 2026.
- Schlichtkrull, Kipf, Bloem, van den Berg, Titov, Welling. *Modeling Relational Data with Graph Convolutional Networks*. Extended Semantic Web Conference, 2018.
- Yang, Yih, He, Gao, Deng. *Embedding Entities and Relations for Learning and Inference in Knowledge Bases*. International Conference on Learning Representations, 2015.
- BioCypher project team (Lobentanzer, Aloy, Saez-Rodriguez, et al.). *Democratizing knowledge representation with BioCypher*. Nature Biotechnology 41, 1056-1059, 2023.
- Chowdhury and Memoli. *Persistent Path Homology of Directed Networks*. Symposium on Discrete Algorithms, 2018.
