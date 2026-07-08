"""Extended pilot with a T+topo arm added to the preregistered comparison.

Runs T, F, F-large, and T+topo at higher synthetic scale (default 1500 nodes,
4500 arcs) as an extension of PILOT_RESULTS.md. Not the preregistered
Reactome-scale finding. Reactome-scale run remains scheduled follow-up
(the Reactome ContentService SBGN exporter was Cloudflare-gated when this
pilot was executed).

Arms
    T        : typed R-GCN with basis decomposition, num_relations=6
    F        : flat R-GCN (relations collapsed), num_relations=1
    F-large  : F parameter-matched to T (Condition B in the preregistration)
    T+topo   : T augmented with a per-node per-relation degree signature
               concatenated to the R-GCN output embedding. The signature is
               a 2 x num_relations = 12-dim vector per node (in-degree and
               out-degree per relation), computed from the training graph.
               Not learned. Feeds a widened relation embedding in DistMult.

The T+topo arm is a defensible cousin of the per-relation magnitude
fingerprint proposed in DESIGN_NOTE.md. True graph magnitude is O(V^3) per
relation and is out of scope for this pilot's compute budget. A per-node
degree signature captures related local structural information cheaply.
DESIGN_NOTE.md remains the reference for the intended full-magnitude arm at
Reactome scale.

Usage
    python run_pilot_v2.py --seeds 10 --epochs 30 --num-nodes 1500 --num-edges 4500
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
from scipy import stats

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
from rgcn_sketch import (
    RGCNLayer,
    RGCNLinkPredictor,
    edge_index_per_relation,
    sample_negatives,
    split_triples,
)
from run_pilot import (
    count_parameters,
    find_matched_hidden_dim,
    make_synthetic_sbgn_graph,
)

# ---------- topological descriptor: per-node per-relation degree signature ----------


def per_node_relation_degrees(triples, num_nodes, num_relations, device):
    """Return a (num_nodes, 2*num_relations) tensor of in/out degrees per relation.

    Computed from training triples only. Not learned. Broadcast to every
    forward pass through the model.
    """
    sig = torch.zeros(num_nodes, 2 * num_relations, device=device)
    for s, r, o in triples:
        # out-degree of s under relation r
        sig[s, r] += 1.0
        # in-degree of o under relation r
        sig[o, num_relations + r] += 1.0
    # Log-scale to keep gradients healthy when we concat with learned embedding.
    sig = torch.log1p(sig)
    return sig


# ---------- topo-augmented model ----------


class RGCNTopoLinkPredictor(nn.Module):
    """R-GCN + DistMult with a topological signature concatenated to the encoder output.

    Relation embedding is widened to match hidden_dim + topo_dim.
    """

    def __init__(self, num_nodes, num_relations, hidden_dim, num_bases, topo_features):
        super().__init__()
        self.hidden_dim = hidden_dim
        self.topo_dim = topo_features.size(1)
        self.aug_dim = hidden_dim + self.topo_dim
        self.node_embedding = nn.Embedding(num_nodes, hidden_dim)
        self.rgcn1 = RGCNLayer(hidden_dim, hidden_dim, num_relations, num_bases)
        self.rgcn2 = RGCNLayer(hidden_dim, hidden_dim, num_relations, num_bases)
        self.rel_embedding = nn.Embedding(num_relations, self.aug_dim)
        nn.init.xavier_uniform_(self.node_embedding.weight)
        nn.init.xavier_uniform_(self.rel_embedding.weight)
        self.register_buffer("topo", topo_features)

    def encode(self, edge_index_per_rel):
        h = self.node_embedding.weight
        h = self.rgcn1(h, edge_index_per_rel)
        h = self.rgcn2(h, edge_index_per_rel)
        return torch.cat([h, self.topo], dim=1)

    def score(self, h_all, s, r, o):
        return (h_all[s] * self.rel_embedding.weight[r] * h_all[o]).sum(-1)


# ---------- training and evaluation ----------


def train_step_generic(model, train_triples, edge_idx_train, num_nodes,
                       optimizer, device, rng, num_neg=5):
    model.train()
    optimizer.zero_grad()
    h = model.encode(edge_idx_train)
    neg = sample_negatives(train_triples, num_nodes, num_neg, rng)

    def tt(triples):
        s = torch.tensor([t[0] for t in triples], device=device, dtype=torch.long)
        r = torch.tensor([t[1] for t in triples], device=device, dtype=torch.long)
        o = torch.tensor([t[2] for t in triples], device=device, dtype=torch.long)
        return s, r, o

    sp, rp, op_ = tt(train_triples)
    sn, rn, on_ = tt(neg)
    pos_scores = model.score(h, sp, rp, op_)
    neg_scores = model.score(h, sn, rn, on_)
    labels = torch.cat([torch.ones_like(pos_scores), torch.zeros_like(neg_scores)])
    scores = torch.cat([pos_scores, neg_scores])
    loss = F.binary_cross_entropy_with_logits(scores, labels)
    loss.backward()
    optimizer.step()
    return loss.item()


def evaluate_generic(model, eval_triples, edge_idx_train, num_nodes, device):
    model.eval()
    with torch.no_grad():
        h = model.encode(edge_idx_train)
        ranks = []
        for s, r, o in eval_triples:
            s_t = torch.full((num_nodes,), s, device=device, dtype=torch.long)
            r_t = torch.full((num_nodes,), r, device=device, dtype=torch.long)
            o_t = torch.arange(num_nodes, device=device, dtype=torch.long)
            scores = model.score(h, s_t, r_t, o_t)
            rank = int((scores >= scores[o]).sum().item())
            ranks.append(rank)
    ranks_np = np.array(ranks) if ranks else np.array([np.nan])
    return {
        "mrr": float(np.nanmean(1.0 / ranks_np)),
        "hits@1": float(np.nanmean(ranks_np <= 1)),
        "hits@3": float(np.nanmean(ranks_np <= 3)),
        "hits@10": float(np.nanmean(ranks_np <= 10)),
    }


def run_one_arm(
    train, val, test, num_nodes, num_relations, num_bases, hidden_dim,
    epochs, lr, num_neg, seed, device, arm_name, topo_features=None,
):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    rng = random.Random(seed)

    if topo_features is not None:
        model = RGCNTopoLinkPredictor(
            num_nodes, num_relations,
            hidden_dim=hidden_dim,
            num_bases=min(num_bases, max(1, num_relations)),
            topo_features=topo_features,
        ).to(device)
    else:
        model = RGCNLinkPredictor(
            num_nodes, num_relations,
            hidden_dim=hidden_dim,
            num_bases=min(num_bases, max(1, num_relations)),
        ).to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)
    edge_idx = edge_index_per_relation(train, num_relations, device)

    losses = []
    best_val, best_state = -1.0, None
    for epoch in range(1, epochs + 1):
        loss = train_step_generic(
            model, train, edge_idx, num_nodes, optimizer, device, rng, num_neg
        )
        losses.append(loss)
        if epoch % 10 == 0 or epoch == epochs:
            m = evaluate_generic(model, val, edge_idx, num_nodes, device)
            if m["mrr"] > best_val:
                best_val = m["mrr"]
                best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}

    if best_state is not None:
        model.load_state_dict(best_state)
    test_metrics = evaluate_generic(model, test, edge_idx, num_nodes, device)

    return {
        "arm": arm_name,
        "seed": seed,
        "hidden_dim": hidden_dim,
        "params": count_parameters(model),
        "final_loss": float(losses[-1]),
        "test_mrr": test_metrics["mrr"],
        "test_hits1": test_metrics["hits@1"],
        "test_hits3": test_metrics["hits@3"],
        "test_hits10": test_metrics["hits@10"],
        "best_val_mrr": best_val,
        "loss_curve": losses,
    }


# ---------- entry point ----------


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--seeds", type=int, default=10)
    p.add_argument("--epochs", type=int, default=30)
    p.add_argument("--num-nodes", type=int, default=1500)
    p.add_argument("--num-edges", type=int, default=4500)
    p.add_argument("--num-relations", type=int, default=6)
    p.add_argument("--hidden-dim", type=int, default=32)
    p.add_argument("--num-bases", type=int, default=3)
    p.add_argument("--lr", type=float, default=1e-2)
    p.add_argument("--num-neg", type=int, default=5)
    p.add_argument("--graph-seed", type=int, default=42)
    p.add_argument("--output-dir", type=Path, default=_HERE / "outputs")
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")

    device = "cpu"
    print(
        f"Pilot v2 (extended). nodes={args.num_nodes} edges={args.num_edges} "
        f"relations={args.num_relations} seeds={args.seeds} epochs={args.epochs}"
    )
    print(
        "SCOPE. Extended synthetic pilot with T+topo arm added. Not the "
        "preregistered Reactome-scale run. See PILOT_V2_RESULTS.md for the "
        "deviations from PREREGISTRATION.md and from PILOT_RESULTS.md."
    )

    graph = make_synthetic_sbgn_graph(
        args.num_nodes, args.num_edges, args.num_relations, args.graph_seed
    )
    print(
        f"Graph. {graph['num_nodes']} nodes, {len(graph['triples'])} edges, "
        f"per-relation counts {graph['per_relation_counts']}"
    )

    typed = graph["triples"]
    train_t, val_t, test_t = split_triples(typed, seed=args.graph_seed)
    train_f = [(s, 0, o) for s, _, o in train_t]
    val_f = [(s, 0, o) for s, _, o in val_t]
    test_f = [(s, 0, o) for s, _, o in test_t]
    print(f"Split. train={len(train_t)} val={len(val_t)} test={len(test_t)}")

    # Topo features are computed ONCE from the training triples only, matched
    # to typed relation labels. F/F-large arms do not see them (they are the
    # baseline comparison arm; T+topo is the intervention).
    topo = per_node_relation_degrees(
        train_t, args.num_nodes, args.num_relations, device
    )
    print(f"Topo signature. {topo.shape} (per-node in/out degree per relation)")

    ref_t = RGCNLinkPredictor(
        args.num_nodes, args.num_relations,
        hidden_dim=args.hidden_dim, num_bases=args.num_bases,
    )
    t_params = count_parameters(ref_t)
    f_large_dim, _ = find_matched_hidden_dim(
        t_params, args.num_nodes, num_relations=1, num_bases=1
    )
    print(f"T params {t_params:,}. F-large hidden_dim={f_large_dim}.")

    all_runs = []
    for seed in range(args.seeds):
        for arm, num_rel, num_bases, hd, train, val, test, topo_feat in [
            ("T",       args.num_relations, args.num_bases, args.hidden_dim, train_t, val_t, test_t, None),
            ("F",       1, 1, args.hidden_dim, train_f, val_f, test_f, None),
            ("F_large", 1, 1, f_large_dim, train_f, val_f, test_f, None),
            ("T_topo",  args.num_relations, args.num_bases, args.hidden_dim, train_t, val_t, test_t, topo),
        ]:
            r = run_one_arm(
                train, val, test, args.num_nodes, num_rel, num_bases, hd,
                args.epochs, args.lr, args.num_neg, seed, device, arm, topo_feat,
            )
            all_runs.append(r)
            print(
                f"  seed={seed:2d} arm={arm:<8s} params={r['params']:>7d} "
                f"mrr={r['test_mrr']:.3f} h@10={r['test_hits10']:.2f}"
            )

    per_arm = defaultdict(list)
    for r in all_runs:
        per_arm[r["arm"]].append(r["test_mrr"])

    def summary(arm):
        v = np.array(per_arm[arm])
        return {
            "mean": float(v.mean()),
            "std": float(v.std(ddof=1)) if len(v) > 1 else 0.0,
            "se": float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else 0.0,
            "values": v.tolist(),
        }

    stats_sum = {a: summary(a) for a in ("T", "F", "F_large", "T_topo")}

    def paired(a, b, alt="greater"):
        va = np.array(per_arm[a])
        vb = np.array(per_arm[b])
        d = va - vb
        try:
            w = stats.wilcoxon(va, vb, alternative=alt)
            wp = float(w.pvalue)
        except Exception:
            wp = None
        pooled = np.sqrt((va.var(ddof=1) + vb.var(ddof=1)) / 2.0) if len(va) > 1 else 1.0
        cd = float(d.mean() / pooled) if pooled > 0 else 0.0
        return {"mean_diff": float(d.mean()), "wilcoxon_p": wp, "cohen_d": cd, "n": len(d)}

    tests = {
        "T_gt_F":         paired("T", "F", "greater"),
        "T_gt_F_large":   paired("T", "F_large", "greater"),
        "T_topo_gt_T":    paired("T_topo", "T", "greater"),
        "T_topo_gt_F_large": paired("T_topo", "F_large", "greater"),
    }

    print("\n=== Aggregate ===")
    for arm, s in stats_sum.items():
        print(
            f"  {arm:<8s} mean_MRR={s['mean']:.4f} +/- {1.96 * s['se']:.4f} "
            f"(n={len(s['values'])})"
        )
    print("\n=== Paired Wilcoxon (one-sided) ===")
    for k, v in tests.items():
        print(
            f"  {k:<20s} mean_diff={v['mean_diff']:+.4f} "
            f"p={v['wilcoxon_p']} d={v['cohen_d']:+.3f}"
        )

    with open(args.output_dir / "pilot_v2_results.json", "w") as f:
        json.dump({
            "config": vars(args),
            "graph": graph,
            "split_sizes": {"train": len(train_t), "val": len(val_t), "test": len(test_t)},
            "topo_dim": int(topo.size(1)),
            "runs": [{k: v for k, v in r.items() if k != "loss_curve"} for r in all_runs],
            "per_arm_summary": stats_sum,
            "paired_tests": tests,
        }, f, indent=2, default=str)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    for arm in ("T", "F", "F_large", "T_topo"):
        curves = [r["loss_curve"] for r in all_runs if r["arm"] == arm]
        ax1.plot(np.mean(curves, axis=0), label=arm)
    ax1.set_xlabel("epoch")
    ax1.set_ylabel("training BCE loss (mean over seeds)")
    ax1.set_title("training loss (pilot v2)")
    ax1.legend()

    arms = ["T", "F", "F_large", "T_topo"]
    means = [stats_sum[a]["mean"] for a in arms]
    errs = [1.96 * stats_sum[a]["se"] for a in arms]
    ax2.bar(arms, means, yerr=errs, capsize=6,
            color=["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd"])
    ax2.set_ylabel("test MRR (mean +/- 1.96 SE)")
    p_ttopo_t = tests["T_topo_gt_T"]["wilcoxon_p"]
    p_ttopo_fl = tests["T_topo_gt_F_large"]["wilcoxon_p"]
    ax2.set_title(
        f"test MRR by arm  |  T_topo>T p={p_ttopo_t}  T_topo>F-large p={p_ttopo_fl}"
    )
    plt.tight_layout()
    plt.savefig(args.output_dir / "pilot_v2_curves.png", dpi=120)
    print(f"\nSaved {args.output_dir / 'pilot_v2_results.json'}")
    print(f"Saved {args.output_dir / 'pilot_v2_curves.png'}")


if __name__ == "__main__":
    raise SystemExit(main())
