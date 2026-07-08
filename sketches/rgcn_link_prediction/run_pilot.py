"""Pilot execution of the preregistered typed-vs-flat R-GCN comparison.

Runs the three-arm ablation from ``PREREGISTRATION.md`` at reduced scale
on a synthetic SBGN-PD-like graph. The result is illustrative, not the
preregistered Reactome-scale finding. Every deviation from the
preregistration is listed in ``PILOT_RESULTS.md`` alongside these numbers.

Model arms
    T        : typed R-GCN with basis decomposition, num_relations=6
    F        : flat R-GCN (all relations collapsed to one), num_relations=1
    F-large  : F with hidden dim increased until parameter count is within
               5% of T (Condition B in the preregistration)

Decoder is DistMult in all three arms. Metric is raw MRR on held-out arcs.
Statistical test is Wilcoxon signed-rank on paired per-seed MRR differences.

Usage
    python run_pilot.py --seeds 10 --epochs 50 --num-nodes 500 --num-edges 1500

Outputs
    outputs/pilot_results.json       per-seed metrics for T, F, F-large
    outputs/pilot_curves.png         training loss and MRR bar chart
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
from matplotlib import pyplot as plt
from scipy import stats

# Reuse the hand-written R-GCN and DistMult classes from the sketch.
_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
from rgcn_sketch import (
    RGCNLinkPredictor,
    edge_index_per_relation,
    evaluate,
    split_triples,
    train_step,
)

# ---------- synthetic SBGN-PD-like graph ----------


def make_synthetic_sbgn_graph(
    num_nodes: int, num_edges: int, num_relations: int, seed: int
):
    """Preferential-attachment-style multi-relational graph.

    Relation distribution follows a rough SBGN-PD prior (has_input and
    has_output dominant, catalyzes/regulates smaller, containment rare).
    """
    rng = np.random.default_rng(seed)
    rel_weights = np.array([0.35, 0.30, 0.15, 0.10, 0.06, 0.04])[:num_relations]
    rel_weights = rel_weights / rel_weights.sum()

    edges = []
    degree = np.zeros(num_nodes)
    for _ in range(num_edges):
        p = (degree + 1.0) / (degree + 1.0).sum()
        src = int(rng.choice(num_nodes, p=p))
        tgt = int(rng.choice(num_nodes, p=p))
        if src == tgt:
            continue
        r = int(rng.choice(num_relations, p=rel_weights))
        edges.append((src, r, tgt))
        degree[src] += 1
        degree[tgt] += 1

    return {
        "num_nodes": num_nodes,
        "num_relations": num_relations,
        "triples": edges,
        "per_relation_counts": {
            f"rel_{i}": sum(1 for _, r, _ in edges if r == i)
            for i in range(num_relations)
        },
    }


# ---------- flat encoder view of the same graph ----------


def flatten_triples(triples, num_relations):
    """Collapse all relations to a single one (relation index 0)."""
    return [(s, 0, o) for s, _, o in triples]


# ---------- parameter-count utilities ----------


def count_parameters(model):
    return sum(p.numel() for p in model.parameters())


def find_matched_hidden_dim(target_params, num_nodes, num_relations, num_bases):
    """Binary search for a hidden dim whose flat R-GCN parameter count
    matches target_params to within 5%.
    """
    lo, hi = 32, 1024
    best_dim, best_diff = None, None
    while lo <= hi:
        mid = (lo + hi) // 2
        try:
            m = RGCNLinkPredictor(
                num_nodes, num_relations, hidden_dim=mid, num_bases=num_bases
            )
            pc = count_parameters(m)
        except Exception:
            hi = mid - 1
            continue
        diff = abs(pc - target_params)
        if best_diff is None or diff < best_diff:
            best_diff = diff
            best_dim = mid
        if pc < target_params:
            lo = mid + 1
        else:
            hi = mid - 1
    return best_dim, best_diff


# ---------- one full T-vs-F-vs-F-large run at a single seed ----------


def run_one_arm(
    train, val, test, num_nodes, num_relations, num_bases,
    hidden_dim, epochs, lr, num_neg, seed, device, arm_name,
):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    rng = random.Random(seed)

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
        loss = train_step(
            model, train, edge_idx, num_nodes, optimizer, device, rng, num_neg
        )
        losses.append(loss)
        if epoch % 10 == 0 or epoch == epochs:
            m = evaluate(model, val, edge_idx, num_nodes, device)
            if m["mrr"] > best_val:
                best_val = m["mrr"]
                best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}

    if best_state is not None:
        model.load_state_dict(best_state)
    test_metrics = evaluate(model, test, edge_idx, num_nodes, device)

    return {
        "arm": arm_name,
        "seed": seed,
        "hidden_dim": hidden_dim,
        "num_relations": num_relations,
        "num_bases": min(num_bases, max(1, num_relations)),
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
    p.add_argument("--epochs", type=int, default=50)
    p.add_argument("--num-nodes", type=int, default=500)
    p.add_argument("--num-edges", type=int, default=1500)
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
        f"Pilot config. nodes={args.num_nodes} edges={args.num_edges} "
        f"relations={args.num_relations} seeds={args.seeds} epochs={args.epochs}"
    )
    print(
        "SCOPE. This is a pilot at synthetic scale. See PILOT_RESULTS.md for "
        "the list of deviations from PREREGISTRATION.md."
    )

    graph = make_synthetic_sbgn_graph(
        args.num_nodes, args.num_edges, args.num_relations, args.graph_seed
    )
    print(
        f"Graph. {graph['num_nodes']} nodes, {len(graph['triples'])} edges, "
        f"per-relation counts {graph['per_relation_counts']}"
    )

    typed_triples = graph["triples"]
    train_t, val_t, test_t = split_triples(typed_triples, seed=args.graph_seed)
    train_f = [(s, 0, o) for s, _, o in train_t]
    val_f = [(s, 0, o) for s, _, o in val_t]
    test_f = [(s, 0, o) for s, _, o in test_t]
    print(f"Split. train={len(train_t)} val={len(val_t)} test={len(test_t)}")

    ref_model = RGCNLinkPredictor(
        args.num_nodes, args.num_relations,
        hidden_dim=args.hidden_dim, num_bases=args.num_bases,
    )
    t_params = count_parameters(ref_model)
    print(f"T reference param count. {t_params:,}")

    f_large_dim, f_large_diff = find_matched_hidden_dim(
        t_params, args.num_nodes, num_relations=1, num_bases=1
    )
    f_large_model = RGCNLinkPredictor(
        args.num_nodes, 1, hidden_dim=f_large_dim, num_bases=1
    )
    print(
        f"F-large hidden dim {f_large_dim} (param count {count_parameters(f_large_model):,}, "
        f"diff {f_large_diff:,} vs T)"
    )

    all_results = []
    for seed in range(args.seeds):
        for arm, num_rel, num_bases, hd, train, val, test in [
            ("T", args.num_relations, args.num_bases, args.hidden_dim,
             train_t, val_t, test_t),
            ("F", 1, 1, args.hidden_dim, train_f, val_f, test_f),
            ("F_large", 1, 1, f_large_dim, train_f, val_f, test_f),
        ]:
            r = run_one_arm(
                train, val, test,
                num_nodes=args.num_nodes, num_relations=num_rel,
                num_bases=num_bases, hidden_dim=hd, epochs=args.epochs,
                lr=args.lr, num_neg=args.num_neg, seed=seed, device=device,
                arm_name=arm,
            )
            all_results.append(r)
            print(
                f"  seed={seed:2d} arm={arm:<8s} params={r['params']:>7d} "
                f"mrr={r['test_mrr']:.3f} h@10={r['test_hits10']:.2f}"
            )

    per_arm = defaultdict(list)
    for r in all_results:
        per_arm[r["arm"]].append(r["test_mrr"])

    def summary(arm):
        vals = np.array(per_arm[arm])
        return {
            "mean": float(vals.mean()),
            "std": float(vals.std(ddof=1)) if len(vals) > 1 else 0.0,
            "se": float(vals.std(ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else 0.0,
            "values": vals.tolist(),
        }

    stats_summary = {arm: summary(arm) for arm in ("T", "F", "F_large")}

    def paired_test(a, b):
        va = np.array(per_arm[a])
        vb = np.array(per_arm[b])
        d = va - vb
        try:
            w = stats.wilcoxon(va, vb, alternative="greater")
            wilc_stat, wilc_p = float(w.statistic), float(w.pvalue)
        except Exception:
            wilc_stat, wilc_p = None, None
        try:
            t = stats.ttest_rel(va, vb, alternative="greater")
            t_stat, t_p = float(t.statistic), float(t.pvalue)
        except Exception:
            t_stat, t_p = None, None
        pooled = np.sqrt((va.var(ddof=1) + vb.var(ddof=1)) / 2.0) if len(va) > 1 else 1.0
        cohen_d = float(d.mean() / pooled) if pooled > 0 else 0.0
        return {
            "mean_diff": float(d.mean()),
            "n_pairs": len(d),
            "wilcoxon_stat": wilc_stat,
            "wilcoxon_p_one_sided": wilc_p,
            "t_stat": t_stat,
            "t_p_one_sided": t_p,
            "cohen_d": cohen_d,
        }

    tests = {
        "T_vs_F": paired_test("T", "F"),
        "T_vs_F_large": paired_test("T", "F_large"),
    }

    print("\n=== Aggregate ===")
    for arm, s in stats_summary.items():
        print(
            f"  {arm:<8s} mean_MRR={s['mean']:.4f} +/- {1.96 * s['se']:.4f} "
            f"(n={len(s['values'])})"
        )
    print("\n=== Paired tests (H1: T > baseline) ===")
    for k, v in tests.items():
        print(
            f"  {k:<15s} mean_diff={v['mean_diff']:+.4f} "
            f"Wilcoxon p={v['wilcoxon_p_one_sided']} "
            f"t p={v['t_p_one_sided']} d={v['cohen_d']:+.3f}"
        )

    with open(args.output_dir / "pilot_results.json", "w") as f:
        json.dump({
            "config": vars(args),
            "graph": graph,
            "split_sizes": {
                "train": len(train_t), "val": len(val_t), "test": len(test_t)
            },
            "t_params": t_params,
            "f_large_hidden_dim": f_large_dim,
            "f_large_params": count_parameters(f_large_model),
            "runs": [{k: v for k, v in r.items() if k != "loss_curve"}
                     for r in all_results],
            "per_arm_summary": stats_summary,
            "paired_tests": tests,
        }, f, indent=2, default=str)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    for arm in ("T", "F", "F_large"):
        loss_curves = [r["loss_curve"] for r in all_results if r["arm"] == arm]
        mean_loss = np.mean(loss_curves, axis=0)
        ax1.plot(mean_loss, label=arm)
    ax1.set_xlabel("epoch")
    ax1.set_ylabel("training BCE loss (mean over seeds)")
    ax1.set_title("training loss")
    ax1.legend()

    arms = ["T", "F", "F_large"]
    means = [stats_summary[a]["mean"] for a in arms]
    errs = [1.96 * stats_summary[a]["se"] for a in arms]
    ax2.bar(arms, means, yerr=errs, capsize=6,
            color=["#1f77b4", "#ff7f0e", "#2ca02c"])
    ax2.set_ylabel("test MRR (mean +/- 1.96 SE)")
    tvf_p = tests["T_vs_F"]["wilcoxon_p_one_sided"]
    tvfl_p = tests["T_vs_F_large"]["wilcoxon_p_one_sided"]
    p1 = f"{tvf_p:.3f}" if tvf_p is not None else "N/A"
    p2 = f"{tvfl_p:.3f}" if tvfl_p is not None else "N/A"
    ax2.set_title(
        f"test MRR by arm\nWilcoxon T>F p={p1}, T>F-large p={p2}"
    )
    plt.tight_layout()
    plt.savefig(args.output_dir / "pilot_curves.png", dpi=120)
    print(f"\nSaved {args.output_dir / 'pilot_results.json'}")
    print(f"Saved {args.output_dir / 'pilot_curves.png'}")


if __name__ == "__main__":
    raise SystemExit(main())
