"""R-GCN link-prediction sketch on the SBGN-PD adapter output.

Trains a Relational Graph Convolutional Network (Schlichtkrull et al., ESWC
2018) with basis decomposition on the multi-relational graph emitted by the
SBGN-PD BioCypher adapter. Scoring is DistMult (Yang et al., ICLR 2015).

Scope. The default fixtures ship 14 nodes and 12 arcs (two glycolysis
sub-pathways). Numbers on that size are illustrative, not benchmark quality.
Point ``--data-dir`` at a directory of real Reactome SBGN-ML files to train
on a production-scale corpus.

Usage
    python rgcn_sketch.py --epochs 300
    python rgcn_sketch.py --data-dir path/to/reactome/sbgn --epochs 500
"""

from __future__ import annotations

import argparse
import json
import logging
import random
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from matplotlib import pyplot as plt

# Make the adapter importable without installing the package.
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
from sbgn_pd_adapter import SBGNPDAdapter  # noqa: E402


# ---------- data preparation ----------

def build_graph(sbgn_files, matcher_threshold=0.7):
    """Run the adapter and reduce its output to (head, relation, tail) triples."""
    adapter = SBGNPDAdapter(sbgn_files, matcher_threshold=matcher_threshold)
    adapter.load()

    nodes = list(adapter.get_nodes())
    edges = list(adapter.get_edges())

    node_to_idx = {n[0]: i for i, n in enumerate(nodes)}
    node_labels = [n[1] for n in nodes]

    relations = sorted({e[3] for e in edges})
    relation_to_idx = {r: i for i, r in enumerate(relations)}

    triples = []
    for _, src, tgt, label, _ in edges:
        if src in node_to_idx and tgt in node_to_idx:
            triples.append(
                (node_to_idx[src], relation_to_idx[label], node_to_idx[tgt])
            )

    per_rel_counts = defaultdict(int)
    for _, r, _ in triples:
        per_rel_counts[relations[r]] += 1

    return {
        "num_nodes": len(nodes),
        "num_relations": len(relations),
        "triples": triples,
        "node_labels": node_labels,
        "relation_names": relations,
        "per_relation_counts": dict(per_rel_counts),
    }


def split_triples(triples, val_frac=0.1, test_frac=0.1, seed=42):
    rng = random.Random(seed)
    triples = list(triples)
    rng.shuffle(triples)
    n = len(triples)
    n_test = max(1, int(round(n * test_frac)))
    n_val = max(1, int(round(n * val_frac)))
    test = triples[:n_test]
    val = triples[n_test:n_test + n_val]
    train = triples[n_test + n_val:]
    return train, val, test


def edge_index_per_relation(triples, num_relations, device):
    per_rel = {}
    by_rel = defaultdict(list)
    for s, r, o in triples:
        by_rel[r].append((s, o))
    for r in range(num_relations):
        pairs = by_rel[r]
        if not pairs:
            continue
        src = torch.tensor([p[0] for p in pairs], device=device, dtype=torch.long)
        tgt = torch.tensor([p[1] for p in pairs], device=device, dtype=torch.long)
        per_rel[r] = (src, tgt)
    return per_rel


# ---------- model ----------

class RGCNLayer(nn.Module):
    """R-GCN layer with basis decomposition, Schlichtkrull et al. 2018 eq. 2-3.

    h_i^{l+1} = sigma( sum_r sum_{j in N_r(i)} 1/c_{i,r} W_r h_j + W_0 h_i )
    W_r = sum_b a_{r,b} V_b (shared bases V_b, per-relation coefficients a_{r,b}).
    """

    def __init__(self, in_dim, out_dim, num_relations, num_bases):
        super().__init__()
        self.in_dim = in_dim
        self.out_dim = out_dim
        self.num_relations = num_relations
        self.num_bases = num_bases
        self.bases = nn.Parameter(torch.empty(num_bases, in_dim, out_dim))
        self.coefs = nn.Parameter(torch.empty(num_relations, num_bases))
        self.self_loop = nn.Parameter(torch.empty(in_dim, out_dim))
        nn.init.xavier_uniform_(self.bases)
        nn.init.xavier_uniform_(self.coefs)
        nn.init.xavier_uniform_(self.self_loop)

    def forward(self, h, edge_index_per_rel):
        num_nodes = h.size(0)
        out = h @ self.self_loop
        for r, (src, tgt) in edge_index_per_rel.items():
            w_r = torch.einsum("b,bij->ij", self.coefs[r], self.bases)
            msg = h[src] @ w_r
            degrees = torch.zeros(num_nodes, device=h.device)
            degrees.scatter_add_(0, tgt, torch.ones_like(tgt, dtype=torch.float))
            norm = 1.0 / degrees.clamp(min=1.0)
            agg = torch.zeros(num_nodes, self.out_dim, device=h.device)
            agg.index_add_(0, tgt, msg * norm[tgt].unsqueeze(-1))
            out = out + agg
        return F.relu(out)


class RGCNLinkPredictor(nn.Module):
    """Two-layer R-GCN encoder with a DistMult decoder."""

    def __init__(self, num_nodes, num_relations, hidden_dim=32, num_bases=2):
        super().__init__()
        self.node_embedding = nn.Embedding(num_nodes, hidden_dim)
        self.rgcn1 = RGCNLayer(hidden_dim, hidden_dim, num_relations, num_bases)
        self.rgcn2 = RGCNLayer(hidden_dim, hidden_dim, num_relations, num_bases)
        self.rel_embedding = nn.Embedding(num_relations, hidden_dim)
        nn.init.xavier_uniform_(self.node_embedding.weight)
        nn.init.xavier_uniform_(self.rel_embedding.weight)

    def encode(self, edge_index_per_rel):
        h = self.node_embedding.weight
        h = self.rgcn1(h, edge_index_per_rel)
        h = self.rgcn2(h, edge_index_per_rel)
        return h

    def score(self, h_all, s, r, o):
        return (h_all[s] * self.rel_embedding.weight[r] * h_all[o]).sum(-1)


# ---------- training and evaluation ----------

def sample_negatives(pos_triples, num_nodes, num_neg_per_pos, rng):
    negs = []
    for s, r, o in pos_triples:
        for _ in range(num_neg_per_pos):
            if rng.random() < 0.5:
                negs.append((rng.randrange(num_nodes), r, o))
            else:
                negs.append((s, r, rng.randrange(num_nodes)))
    return negs


def train_step(model, train_triples, edge_idx_train, num_nodes, optimizer,
               device, rng, num_neg=5):
    model.train()
    optimizer.zero_grad()
    h = model.encode(edge_idx_train)
    neg = sample_negatives(train_triples, num_nodes, num_neg, rng)

    def to_tensors(triples):
        s = torch.tensor([t[0] for t in triples], device=device, dtype=torch.long)
        r = torch.tensor([t[1] for t in triples], device=device, dtype=torch.long)
        o = torch.tensor([t[2] for t in triples], device=device, dtype=torch.long)
        return s, r, o

    s_p, r_p, o_p = to_tensors(train_triples)
    s_n, r_n, o_n = to_tensors(neg)
    pos_scores = model.score(h, s_p, r_p, o_p)
    neg_scores = model.score(h, s_n, r_n, o_n)
    labels = torch.cat([torch.ones_like(pos_scores), torch.zeros_like(neg_scores)])
    scores = torch.cat([pos_scores, neg_scores])
    loss = F.binary_cross_entropy_with_logits(scores, labels)
    loss.backward()
    optimizer.step()
    return loss.item()


def evaluate(model, eval_triples, edge_idx_train, num_nodes, device):
    """Ranking eval. For each true (s, r, o) rank the true o against all nodes."""
    model.eval()
    with torch.no_grad():
        h = model.encode(edge_idx_train)
        ranks = []
        for s, r, o in eval_triples:
            s_t = torch.full((num_nodes,), s, device=device, dtype=torch.long)
            r_t = torch.full((num_nodes,), r, device=device, dtype=torch.long)
            o_t = torch.arange(num_nodes, device=device, dtype=torch.long)
            scores = model.score(h, s_t, r_t, o_t)
            true_score = scores[o]
            rank = int((scores >= true_score).sum().item())
            ranks.append(rank)
    ranks_np = np.array(ranks) if ranks else np.array([np.nan])
    return {
        "mrr": float(np.nanmean(1.0 / ranks_np)),
        "hits@1": float(np.nanmean(ranks_np <= 1)),
        "hits@3": float(np.nanmean(ranks_np <= 3)),
        "hits@10": float(np.nanmean(ranks_np <= 10)),
    }


# ---------- entry point ----------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=_REPO_ROOT / "data",
        help="Directory of SBGN-ML files (recursed). Defaults to repo /data.",
    )
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--hidden-dim", type=int, default=32)
    parser.add_argument("--num-bases", type=int, default=2)
    parser.add_argument("--lr", type=float, default=1e-2)
    parser.add_argument("--num-neg", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "outputs",
    )
    args = parser.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    rng = random.Random(args.seed)
    device = "cpu"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s"
    )

    sbgn_files = sorted(
        list(args.data_dir.rglob("*.sbgn"))
        + list(args.data_dir.rglob("*.sbgnml"))
        + list(args.data_dir.rglob("*.xml"))
    )
    if not sbgn_files:
        print(f"ERROR: no SBGN-ML files found under {args.data_dir}")
        return 2
    print(f"Adapter input: {len(sbgn_files)} SBGN-ML files")

    graph = build_graph(sbgn_files)
    print(
        f"Graph: {graph['num_nodes']} nodes, {len(graph['triples'])} edges, "
        f"{graph['num_relations']} relations"
    )
    for rel, count in graph["per_relation_counts"].items():
        print(f"  {rel}: {count}")

    train, val, test = split_triples(graph["triples"], seed=args.seed)
    print(f"Split: {len(train)} train / {len(val)} val / {len(test)} test")
    if len(train) < 8:
        print(
            "NOTE. Training split is small. Metrics are illustrative, not "
            "benchmark. Use --data-dir on a Reactome corpus for scale."
        )

    edge_idx_train = edge_index_per_relation(train, graph["num_relations"], device)

    model = RGCNLinkPredictor(
        graph["num_nodes"],
        graph["num_relations"],
        hidden_dim=args.hidden_dim,
        num_bases=min(args.num_bases, graph["num_relations"]),
    ).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-5)

    losses = []
    val_history = []
    best_val_mrr = -1.0
    best_state = None
    for epoch in range(1, args.epochs + 1):
        loss = train_step(
            model, train, edge_idx_train, graph["num_nodes"],
            optimizer, device, rng, args.num_neg,
        )
        losses.append(loss)
        if epoch % 10 == 0 or epoch == args.epochs:
            v = evaluate(model, val, edge_idx_train, graph["num_nodes"], device)
            val_history.append((epoch, v["mrr"]))
            print(
                f"epoch {epoch:4d}  loss {loss:.4f}  val_mrr {v['mrr']:.3f}  "
                f"val_h@1 {v['hits@1']:.2f}  val_h@3 {v['hits@3']:.2f}"
            )
            if v["mrr"] > best_val_mrr:
                best_val_mrr = v["mrr"]
                best_state = {
                    k: v_.detach().clone() for k, v_ in model.state_dict().items()
                }

    if best_state is not None:
        model.load_state_dict(best_state)
    test_metrics = evaluate(
        model, test, edge_idx_train, graph["num_nodes"], device
    )
    print(
        f"\nTest  MRR {test_metrics['mrr']:.3f}  H@1 {test_metrics['hits@1']:.2f}  "
        f"H@3 {test_metrics['hits@3']:.2f}  H@10 {test_metrics['hits@10']:.2f}"
    )

    results = {
        "graph": {
            "num_nodes": graph["num_nodes"],
            "num_relations": graph["num_relations"],
            "relation_names": graph["relation_names"],
            "per_relation_counts": graph["per_relation_counts"],
            "num_edges": len(graph["triples"]),
            "splits": {"train": len(train), "val": len(val), "test": len(test)},
        },
        "test_metrics": test_metrics,
        "best_val_mrr": best_val_mrr,
        "config": {
            "epochs": args.epochs,
            "hidden_dim": args.hidden_dim,
            "num_bases": args.num_bases,
            "lr": args.lr,
            "num_neg": args.num_neg,
            "seed": args.seed,
        },
    }
    with open(args.output_dir / "results.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.plot(losses)
    ax1.set_xlabel("epoch")
    ax1.set_ylabel("training BCE loss")
    ax1.set_title("R-GCN training loss")
    if val_history:
        eps, mrrs = zip(*val_history)
        ax2.plot(eps, mrrs, marker="o")
        ax2.set_xlabel("epoch")
        ax2.set_ylabel("validation MRR")
        ax2.set_title(f"validation MRR (best {best_val_mrr:.3f})")
    plt.tight_layout()
    plt.savefig(args.output_dir / "training_curves.png", dpi=120)
    print(f"Saved {args.output_dir / 'results.json'}")
    print(f"Saved {args.output_dir / 'training_curves.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
