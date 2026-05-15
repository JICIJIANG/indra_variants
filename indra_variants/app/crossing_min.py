"""LNS crossing minimisation (Wilson et al., IEEE TVCG 2025).

Pure layout algorithm — no Dash / pandas dependencies.  Used by the network
builders in :mod:`indra_variants.app.graph_builders` to compute layer
orderings that minimise edge crossings.
"""
import itertools as _it
import logging
import random as _rand
import time as _time
from collections import defaultdict

from scipy.optimize import linprog
from scipy.sparse import lil_matrix


_log = logging.getLogger(__name__)


def _n_comb3(n):
    if n < 3:
        return 0
    return n * (n - 1) * (n - 2) // 6


def _add_anchors(layer_nodes, cross_layer_edges, node_layer):
    """Insert dummy anchor nodes for edges spanning >1 layer.

    Returns (aug_layers, adj_edges, node_layer) – the augmented graph
    where every edge connects adjacent layers.
    """
    aug = {li: list(ns) for li, ns in layer_nodes.items()}
    adj_edges: list[tuple] = []
    _aid = 0
    for src, tgt, w in cross_layer_edges:
        ls, lt = node_layer[src], node_layer[tgt]
        if ls > lt:
            src, tgt = tgt, src
            ls, lt = lt, ls
        if lt - ls == 1:
            adj_edges.append((src, tgt, w))
        else:
            prev = src
            for mid in range(ls + 1, lt):
                aname = f"__anch_{_aid}"
                _aid += 1
                aug.setdefault(mid, []).append(aname)
                node_layer[aname] = mid
                adj_edges.append((prev, aname, w))
                prev = aname
            adj_edges.append((prev, tgt, w))
    return aug, adj_edges


def _count_adj_crossings(order, adj_edges, node_layer, pos):
    """Count crossings between adjacent-layer edges using position dict."""
    ebl: dict[int, list] = defaultdict(list)
    for s, t, _w in adj_edges:
        ls = node_layer[s]
        ebl[ls].append((pos.get(s, 0.0), pos.get(t, 0.0)))
    total = 0
    for _li, edges in ebl.items():
        for a in range(len(edges)):
            for b in range(a + 1, len(edges)):
                if (edges[a][0] - edges[b][0]) * \
                   (edges[a][1] - edges[b][1]) < 0:
                    total += 1
    return total


def _optimize_layer_ordering(layer_nodes, cross_layer_edges,
                             time_budget=3.0, sub_time=0.3,
                             neighbourhood_k=12, fixed_layers=None):
    """Minimise edge crossings with Large Neighbourhood Search (LNS).

    1. Add anchor (dummy) nodes so every edge spans exactly one layer.
    2. Build a barycentric initial solution  (fast, O(n·iter)).
    3. Repeatedly pick a random candidate node, collect a small
       neighbourhood (≤ *neighbourhood_k* nodes per layer), solve
       that sub-problem optimally via ILP, and splice the improved
       ordering back into the global solution.
    4. Stop when *time_budget* seconds have elapsed.
    5. Strip anchor nodes from the final output.
    """
    fixed_layers = set(fixed_layers or [])
    node_layer: dict = {}
    for li, nodes in layer_nodes.items():
        for n in nodes:
            node_layer[n] = li

    aug, adj_edges = _add_anchors(layer_nodes, cross_layer_edges, node_layer)

    # Adjacency (including anchors) for barycentric
    adj: dict = defaultdict(set)
    for s, t, _w in adj_edges:
        adj[s].add(t)
        adj[t].add(s)

    # Edge weight lookup keyed on (lower-layer-node, upper-layer-node)
    edge_w: dict[tuple, float] = {}
    for s, t, w in adj_edges:
        edge_w[(s, t)] = edge_w.get((s, t), 0) + w

    sorted_layers = sorted(aug.keys())
    movable_nodes = [
        n for li in sorted_layers if li not in fixed_layers for n in aug[li]
    ]

    # ---------- Phase 1: barycentric initial solution ---------------------
    order: dict[int, list] = {li: list(ns) for li, ns in aug.items()}
    pos: dict = {}

    def _assign_pos(li):
        for i, n in enumerate(order[li]):
            pos[n] = float(i)

    for li in sorted_layers:
        _assign_pos(li)

    def _bary(n):
        vals = [pos[nb] for nb in adj[n] if nb in pos]
        return sum(vals) / len(vals) if vals else pos.get(n, 0.0)

    for _ in range(20):
        for li in sorted_layers:
            if li in fixed_layers:
                continue
            order[li].sort(key=lambda n: (_bary(n), n))
            _assign_pos(li)
        for li in reversed(sorted_layers):
            if li in fixed_layers:
                continue
            order[li].sort(key=lambda n: (_bary(n), n))
            _assign_pos(li)

    # ---------- Phase 2: LNS – local ILP improvements --------------------
    best_crossings = _count_adj_crossings(order, adj_edges, node_layer, pos)
    t0 = _time.time()
    n_improvements = 0
    n_iters = 0

    while movable_nodes and _time.time() - t0 < time_budget and best_crossings > 0:
        n_iters += 1
        candidate = _rand.choice(movable_nodes)
        c_layer = node_layer[candidate]

        # Neighbourhood: candidate + 1-hop + 2-hop, capped per layer
        sub_set: dict[int, set] = defaultdict(set)
        sub_set[c_layer].add(candidate)
        for nb in adj[candidate]:
            sub_set[node_layer[nb]].add(nb)
        for nb in list(adj[candidate]):
            for nb2 in adj[nb]:
                sub_set[node_layer[nb2]].add(nb2)

        sub_nodes: dict[int, list] = {}
        for li, ns in sub_set.items():
            if li in fixed_layers:
                continue
            lst = [n for n in order[li] if n in ns][:neighbourhood_k]
            if lst:
                sub_nodes[li] = lst

        if all(len(ns) < 2 for ns in sub_nodes.values()):
            continue

        # Expand each touched layer to include ALL its nodes (up to the cap).
        # Without this, two source nodes in the same layer (e.g. PTGS2 and
        # glucocorticoid) whose anchor chains only share a fixed ancestor
        # (layer 0 variants) are never simultaneously visible in a single
        # sub-problem, so the ILP cannot detect or fix their crossing edges.
        for li in list(sub_nodes.keys()):
            sub_nodes[li] = order[li][:neighbourhood_k]

        # Collect sub-edges (only between adjacent-layer sub-node pairs)
        sub_node_all = {n for ns in sub_nodes.values() for n in ns}
        sub_edges = [(s, t, w) for (s, t), w in edge_w.items()
                     if s in sub_node_all and t in sub_node_all]
        if not sub_edges:
            continue

        sub_result = _solve_sub_ilp(sub_nodes, sub_edges, node_layer,
                                    sub_time)
        if sub_result is None:
            continue

        # Splice improved sub-ordering back into global order
        old_order = {li: list(order[li]) for li in sub_result}
        for li, sub_ordered in sub_result.items():
            sset = set(sub_ordered)
            idx_map = sorted(i for i, n in enumerate(order[li]) if n in sset)
            base = [n for n in order[li] if n not in sset]
            for slot, n in zip(idx_map, sub_ordered):
                base.insert(slot, n)
            order[li] = base
            _assign_pos(li)

        new_crossings = _count_adj_crossings(order, adj_edges, node_layer, pos)
        if new_crossings < best_crossings:
            best_crossings = new_crossings
            n_improvements += 1
        else:
            for li, old in old_order.items():
                order[li] = old
                _assign_pos(li)

    elapsed = _time.time() - t0
    _log.info("LNS crossing minimisation: %d crossings, %d improvements, "
              "%d iters in %.2fs", best_crossings, n_improvements,
              n_iters, elapsed)

    # Strip anchors, return only original-node orderings
    result: dict[int, list] = {}
    for li in sorted(layer_nodes.keys()):
        result[li] = [n for n in order[li]
                      if not str(n).startswith("__anch_")]
    return result


def _solve_sub_ilp(layer_nodes, edges, node_layer, cutoff):
    """Solve a small sub-problem exactly via ILP.  Returns optimised
    orderings or *None* on failure."""
    all_n = []
    for li in sorted(layer_nodes):
        all_n.extend(layer_nodes[li])
    if len(all_n) < 2:
        return None
    iid = {n: i for i, n in enumerate(all_n)}

    adj_e: list[tuple] = []
    for s, t, w in edges:
        ls, lt = node_layer[s], node_layer[t]
        if ls == lt:
            continue
        if ls > lt:
            s, t = t, s
        adj_e.append((iid[s], iid[t], w))

    x_vars: dict[tuple, int] = {}
    nv = 0
    for li in sorted(layer_nodes):
        ids = sorted(iid[n] for n in layer_nodes[li])
        for a, b in _it.combinations(ids, 2):
            x_vars[(a, b)] = nv
            nv += 1

    # Group edges by (src_layer, tgt_layer) so only same-span pairs interact
    ebl: dict[tuple, list] = defaultdict(list)
    for s, t, w in adj_e:
        ls, lt = node_layer[all_n[s]], node_layer[all_n[t]]
        ebl[(min(ls, lt), max(ls, lt))].append((s, t, w))

    c_vars: dict = {}
    c_wt: dict = {}
    for _lpair, elist in ebl.items():
        for i in range(len(elist)):
            u1, v1, w1 = elist[i]
            for j in range(i + 1, len(elist)):
                u2, v2, w2 = elist[j]
                if u1 != u2 and v1 != v2:
                    c_vars[((u1, v1), (u2, v2))] = nv
                    c_wt[nv] = w1 * w2
                    nv += 1

    if not c_vars:
        return None

    n_trans = sum(
        2 * _n_comb3(len(layer_nodes[li]))
        for li in sorted(layer_nodes)
    )
    n_cons = n_trans + 2 * len(c_vars)
    A = lil_matrix((n_cons, nv))
    b_ub = [0.0] * n_cons
    r = 0

    for li in sorted(layer_nodes):
        ids = sorted(iid[n] for n in layer_nodes[li])
        for i, j, k in _it.combinations(ids, 3):
            A[r, x_vars[(i, j)]] = -1
            A[r, x_vars[(j, k)]] = -1
            A[r, x_vars[(i, k)]] = 1
            r += 1
            A[r, x_vars[(i, j)]] = 1
            A[r, x_vars[(j, k)]] = 1
            A[r, x_vars[(i, k)]] = -1
            b_ub[r] = 1
            r += 1

    def _xdir(a, b_id):
        if (a, b_id) in x_vars:
            return x_vars[(a, b_id)], 1, 0
        return x_vars[(b_id, a)], -1, 1

    for ((u1, v1), (u2, v2)), ci in c_vars.items():
        xi_u, du, fu = _xdir(u1, u2)
        xi_v, dv, fv = _xdir(v1, v2)
        A[r, xi_u] = du;  A[r, xi_v] = -dv;  A[r, ci] = -1
        b_ub[r] = fv - fu
        r += 1
        A[r, xi_u] = -du; A[r, xi_v] = dv;  A[r, ci] = -1
        b_ub[r] = fu - fv
        r += 1

    obj = [0.0] * nv
    for vi, w in c_wt.items():
        obj[vi] = w

    try:
        res = linprog(
            obj, method="highs",
            A_ub=A.tocsc(), b_ub=b_ub,
            bounds=[(0, 1)] * nv,
            integrality=[1] * nv,
            options={"time_limit": cutoff, "disp": False},
        )
    except Exception:
        return None
    if res.x is None:
        return None

    result: dict[int, list] = {}
    for li in sorted(layer_nodes):
        ids = sorted(iid[n] for n in layer_nodes[li])
        if len(ids) < 2:
            result[li] = [all_n[ids[0]]] if ids else []
            continue
        rank = {n: 0 for n in ids}
        for a, b_id in _it.combinations(ids, 2):
            if round(res.x[x_vars[(a, b_id)]]) == 1:
                rank[b_id] += 1
            else:
                rank[a] += 1
        ordered = sorted(ids, key=lambda n: rank[n])
        result[li] = [all_n[n] for n in ordered]
    return result
