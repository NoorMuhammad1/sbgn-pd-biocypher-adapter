# R-GCN link-prediction sketch on the SBGN-PD adapter output

A small end-to-end learning sketch on the multi-relational graph the parent
SBGN-PD BioCypher adapter emits. Intended for readers of the parent adapter
README who want to see that its output feeds a standard graph-representation-learning
pipeline. Run `python sketches/rgcn_link_prediction/rgcn_sketch.py --epochs 300`
from the repository root, read `outputs/results.json`, done in under a minute.

The sketch trains a Relational Graph Convolutional Network (Schlichtkrull et
al., ESWC 2018) with basis decomposition and scores candidate edges with a
DistMult decoder (Yang et al., ICLR 2015). Everything is implemented by hand
in `rgcn_sketch.py` in around 370 lines. No PyTorch Geometric or DGL dependency.

## What this is

A wiring demonstration. It reads an SBGN-ML file, runs the adapter to
produce typed multi-relational triples, feeds those into the R-GCN encoder
and DistMult decoder, and reports MRR and Hits@K. The point is to show that
the adapter's output is directly consumable by a standard link-prediction
pipeline without special preprocessing. The reification decision (SBGN
process glyphs mapped to Biolink `BiologicalProcess` nodes with typed
`has_input`, `has_output`, `catalyzes`, and `regulates` edges) leaves the
graph in exactly the shape an R-GCN expects. The parent adapter also emits
`positively_regulates`, `negatively_regulates`, and `same_as` predicates on
modulation and equivalence arcs. Six distinct Biolink predicates in total
for the flow-and-modulation subset, mapped from the seven SBGN-PD Level 1
arc classes. See `PREREGISTRATION.md` Section 4 for the exact enumeration.

## Companion documents

- [`PREREGISTRATION.md`](PREREGISTRATION.md) states the frozen typed-versus-flat
  R-GCN comparison protocol this sketch is a wiring demo for. Committed before
  any Reactome-scale results are read.
- [`PILOT_RESULTS.md`](PILOT_RESULTS.md) reports v1 pilot outcomes at 400
  nodes and 1197 arcs on synthetic preferential-attachment data. F-large won
  under the diagnostic Wilcoxon direction at p approximately 0.010, Case IVb.
- [`PILOT_V2_RESULTS.md`](PILOT_V2_RESULTS.md) reports v2 pilot outcomes at
  800 nodes and 2400 arcs. Adds a T+topo arm. Same substrate direction, Case IVa.
- [`../reification_vs_missingness/NOTE.md`](../reification_vs_missingness/NOTE.md)
  disentangles the schema-level lossy-compression fix (reification) from the
  embedding-space geometry fix (soft manifolds after Marinoni et al. 2026 TPAMI).
- [`../topological_followup/DESIGN_NOTE.md`](../topological_followup/DESIGN_NOTE.md)
  works through the case for directed topological descriptors as a
  complementary view to R-GCN.

The T, F, F-large, and T+topo arms referenced in the pilot documents are
defined in `PREREGISTRATION.md` Section 4 and, for T+topo, in
`PILOT_V2_RESULTS.md`.

## What this is not

A benchmark. The default fixtures produce a graph with 14 nodes and 12
arcs across 3 relations. Splitting that gives 10 training triples and 1
validation and 1 test triple, which is far too little to draw any
conclusion about R-GCN's performance on real pathway data. Any numbers you
see on the default run are illustrative of the pipeline, not the model.

## Running it

Run this from the repository root.

```bash
python sketches/rgcn_link_prediction/rgcn_sketch.py --epochs 300
```

Outputs land in `sketches/rgcn_link_prediction/outputs/`.

- `results.json` with graph statistics, split sizes, test metrics, and config.
- `training_curves.png` with training loss and validation MRR over epochs.

To scale to real data, point `--data-dir` at a directory of Reactome or
PANTHER SBGN-ML files.

```bash
python sketches/rgcn_link_prediction/rgcn_sketch.py \
    --data-dir /path/to/reactome/sbgn --epochs 500 --hidden-dim 64
```

## Reproducing the pilots

The exact invocations behind `PILOT_RESULTS.md` and `PILOT_V2_RESULTS.md`.

- v1. `python sketches/rgcn_link_prediction/run_pilot.py --lr 1e-2` (script
  defaults 400 nodes, 1197 arcs after self-loop skips, 25 epochs, 10 seeds,
  `hidden_dim = 32`, `num_bases = 3`).
- v2. `python sketches/rgcn_link_prediction/run_pilot_v2.py --num-nodes 800 --num-edges 2400 --seeds 10 --epochs 20 --lr 1e-2`
  (overrides the argparse defaults of 1500 nodes, 4500 arcs, 30 epochs).

## Implementation notes

- **R-GCN layer.** Follows Schlichtkrull et al. 2018 equation 2. Per-relation
  weight matrices `W_r` are factorised through basis decomposition (equation
  3) with `num_bases` shared bases and per-relation coefficient vectors. In
  the default 3-relation setting `num_bases=2` is a mild compression. On the
  six SBGN-PD relations preserved at Reactome scale (see `PREREGISTRATION.md`
  Section 4) the compression starts to bite.
- **Normalisation.** Messages are averaged by per-target in-degree within
  each relation, matching the eq. 2 `1/c_{i,r}` term.
- **Decoder.** DistMult scoring `f(s, r, o) = sum_k h_s[k] * r_r[k] * h_o[k]`,
  the Hadamard-product trilinear form of Yang et al. 2015. Chosen because it
  is the scorer the original R-GCN link-prediction paper reports on. It
  cannot model asymmetric relations well, which is one of the honest
  limitations of this sketch (see below).
- **Negative sampling.** Standard uniform corruption of either the head or
  the tail, 5 negatives per positive.
- **Loss.** Binary cross-entropy with logits on the concatenated positive
  and negative scores.
- **Evaluation.** For each held-out triple `(s, r, o)`, the true object `o`
  is ranked against every node in the graph. MRR, Hits@1, Hits@3, Hits@10
  are computed on the resulting ranks. This is the raw (non-filtered)
  ranking. On real data a filtered ranking (excluding other known-true
  triples with the same `(s, r, ?)` from the ranking) is standard, and
  adding that is a natural extension.

## Honest limitations

1. **Scale.** The default graph is 14 nodes. Any metric on this is a coin
   flip. The pipeline works, and that is the deliverable.
2. **Symmetric decoder.** DistMult treats `(s, r, o)` and `(o, r, s)` with
   the same score. This is wrong for asymmetric relations like `has_input`
   and `has_output`. ComplEx, TransE, or a small MLP decoder would fix
   this. The default here is DistMult because it is what the R-GCN paper
   reports on.
3. **Random init on node features.** The encoder starts from a learned
   `nn.Embedding` per node. There is no biological feature vector
   (sequence, expression profile, structural embedding) attached to each
   node. Wiring in real features would be the second natural extension.
4. **No filtered evaluation.** See the evaluation note above. The preregistered
   follow-up (`PREREGISTRATION.md` Section 7) treats filtered ranking as a
   required sensitivity check, and both pilots log the omission as a deviation.
5. **CPU only.** The sketch runs in a few seconds on the default data. On
   Reactome scale a GPU is advisable.

## Why this exists

Written as a companion to the SBGN-PD BioCypher adapter in the parent
directory, to demonstrate that the reified property-graph output the
adapter emits sits at the right shape for standard graph-representation
learning. The parent adapter README explains the reification design call.
This sketch shows what you can do with the resulting graph.

## References

- Schlichtkrull, Kipf, Bloem, van den Berg, Titov, Welling. *Modeling
  Relational Data with Graph Convolutional Networks*. ESWC 2018.
- Yang, Yih, He, Gao, Deng. *Embedding Entities and Relations for
  Learning and Inference in Knowledge Bases*. ICLR 2015.
- Bordes, Usunier, Garcia-Duran, Weston, Yakhnenko. *Translating
  Embeddings for Modeling Multi-relational Data*. NeurIPS 2013.
- Trouillon, Welbl, Riedel, Gaussier, Bouchard. *Complex Embeddings for
  Simple Link Prediction*. ICML 2016.
