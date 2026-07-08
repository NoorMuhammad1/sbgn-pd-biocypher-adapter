"""FFL operator-regime demo. Regime A (signed R-GCN with mean readout).

Verifies empirically the Regime A argument in
[`NOTE.md`](NOTE.md) Section 3. Coherent-1 (C1) and Incoherent-1 (I1)
feedforward loops on the same three nodes with the same three directed
edges differ only on the sign of the (i -> t) edge. Under a signed R-GCN
with mean readout, the pre-nonlinearity difference at the target node
`t` should be

    Delta h_t^(1) = 2 * (1 / c_{t,r}) * W_r^(0) * h_i^(0)

a factor-of-2 sign-flipped vector, and the graph-level mean-readout
difference should be `Delta h_t^(1) / 3`. That is what this script
checks.

Regimes B (sign-aware MPSN, Bodnar et al. 2021) and C (directed
semi-simplicial, Lecha et al. 2025) remain as reading in `NOTE.md`
rather than as code. This demo covers only Regime A.

Run: `python rgcn_ffl_demo.py`. No external dependencies beyond PyTorch.
"""

import torch


def build_ffl_signs(kind):
    """Return the signed adjacency for a type-1 FFL.

    Nodes are labelled `s` (source), `i` (intermediate), `t` (target).
    Edges are always (s->i), (s->t), (i->t). The regulatory sign on
    (i->t) is the only difference between C1 and I1.
    """
    if kind == "C1":
        return {("s", "i"): +1.0, ("s", "t"): +1.0, ("i", "t"): +1.0}
    if kind == "I1":
        return {("s", "i"): +1.0, ("s", "t"): +1.0, ("i", "t"): -1.0}
    raise ValueError(f"unknown FFL kind {kind}")


def signed_rgcn_layer(H, W_r, signs, node_index):
    """One signed R-GCN layer with a single relation type.

    H:          (N, D) node features
    W_r:        (D, D) relation weight matrix
    signs:      dict[(str, str) -> float] with per-edge signs in {+1, -1}
    node_index: dict[str -> int] mapping node label to row in H

    Returns pre-nonlinearity aggregate at each target node. Skips the
    self-loop term because it is identical between C1 and I1 and
    therefore cancels in the C1-minus-I1 difference we care about.
    """
    edges = list(signs.keys())
    N, D = H.shape
    out = torch.zeros_like(H)

    in_degree = {}
    for (_u, v) in edges:
        in_degree[v] = in_degree.get(v, 0) + 1

    for (u, v) in edges:
        c_vr = in_degree[v]
        sigma_uv = signs[(u, v)]
        u_idx = node_index[u]
        v_idx = node_index[v]
        message = (sigma_uv / c_vr) * (H[u_idx] @ W_r.T)
        out[v_idx] = out[v_idx] + message

    return out


def run_config(kind, H0, W_r, node_index):
    signs = build_ffl_signs(kind)
    H1_pre = signed_rgcn_layer(H0, W_r, signs, node_index)
    graph_readout = H1_pre.mean(dim=0)
    return H1_pre, graph_readout


def main(seed=0, hidden=8):
    torch.manual_seed(seed)
    node_index = {"s": 0, "i": 1, "t": 2}

    H0 = torch.randn(3, hidden)
    W_r = torch.randn(hidden, hidden) * 0.5

    H_c1, graph_c1 = run_config("C1", H0, W_r, node_index)
    H_i1, graph_i1 = run_config("I1", H0, W_r, node_index)

    delta = H_c1 - H_i1

    print("Per-node pre-nonlinearity difference (C1 minus I1):")
    for label, idx in node_index.items():
        norm = delta[idx].norm().item()
        print(f"  node {label}: L2 norm = {norm:.6f}")

    c_tr = 2  # in-degree of t under the single relation (edges s->t and i->t)
    predicted_delta_t = 2.0 * (1.0 / c_tr) * (H0[node_index["i"]] @ W_r.T)
    observed_delta_t = delta[node_index["t"]]
    delta_t_err = (predicted_delta_t - observed_delta_t).norm().item()

    print()
    print("Delta at target node t:")
    print(f"  predicted 2/(c_t,r) * W_r * h_i : L2 norm = {predicted_delta_t.norm().item():.6f}")
    print(f"  observed                         : L2 norm = {observed_delta_t.norm().item():.6f}")
    print(f"  |predicted minus observed|       : L2 norm = {delta_t_err:.3e}")

    graph_delta = graph_c1 - graph_i1
    predicted_graph_delta = predicted_delta_t / 3.0
    graph_err = (predicted_graph_delta - graph_delta).norm().item()

    print()
    print("Graph-level mean-readout difference:")
    print(f"  predicted delta_t / 3     : L2 norm = {predicted_graph_delta.norm().item():.6f}")
    print(f"  observed                  : L2 norm = {graph_delta.norm().item():.6f}")
    print(f"  |predicted minus observed|: L2 norm = {graph_err:.3e}")

    print()
    print("=" * 60)
    print("Regime A (signed R-GCN with mean readout).")
    print(f"  Graph embedding L2 distance C1 vs I1 = {graph_delta.norm().item():.6f}")
    if graph_delta.norm().item() > 1e-4:
        print("  DISTINGUISHABLE. The mean readout does not average the")
        print("  sign flip away. Regime A does separate C1 from I1 at the")
        print("  graph-level readout under a signed-scalar encoding.")
    else:
        print("  INDISTINGUISHABLE. The mean readout absorbed the sign flip.")
        print("  This would contradict the argument in NOTE.md Section 3.")
    print("=" * 60)

    return {
        "delta_t_err": delta_t_err,
        "graph_err": graph_err,
        "graph_l2": graph_delta.norm().item(),
    }


if __name__ == "__main__":
    result = main()
    assert result["delta_t_err"] < 1e-5, (
        f"delta_t prediction failed by {result['delta_t_err']:.3e}; "
        "check the memo Section 3 arithmetic against the code above"
    )
    assert result["graph_err"] < 1e-5, (
        f"graph-readout prediction failed by {result['graph_err']:.3e}"
    )
    assert result["graph_l2"] > 1e-4, (
        "C1 and I1 collapsed to the same graph embedding, which would "
        "contradict the Regime A argument in NOTE.md"
    )
    print()
    print("All three checks passed. Regime A argument holds empirically.")
