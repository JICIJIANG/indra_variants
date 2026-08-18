"""Network graph builders."""

import copy
import math
import re
from collections import defaultdict, deque
from functools import lru_cache
from typing import Optional

import networkx as nx
import pandas as pd

from indra_variants.app.crossing_min import _optimize_layer_ordering
from indra_variants.app.data_index import ENDPOINT_INDEX, TSV_FILES
from indra_variants.app.utils import (
    _protein_href,
    _sort_text,
    _variant_aa_position,
    format_star_rating,
)


# ----------------------Build Graph--------------------------–
def build_elements(prot: str, variant_aa_range: Optional[tuple[float, float]] = None):
    df_path = TSV_FILES[prot]
    df = pd.read_csv(df_path, sep="\t").fillna('')

    def choose_best_clinvar(existing: Optional[dict], candidate: Optional[dict]):
        if not candidate:
            return existing
        if not existing:
            return candidate
        # Prefer richer condition lists, then keep existing.
        prev_len = len(existing.get("conditions", ""))
        cand_len = len(candidate.get("conditions", ""))
        if cand_len > prev_len:
            return candidate
        return existing

    G = nx.MultiDiGraph()
    G.add_node(prot)
    variants = set(df["variant_info"])
    endpoints = set(df["biological_process/disease"])
    endpoint_freq = df["biological_process/disease"].value_counts().to_dict()
    protein_uniprot_id = ""

    variant_meta = defaultdict(lambda: {
        "domains": set(),
        "domain_notes": set(),
        "clinvar_data": None,
        "protein_pos": None,
        "allele_id": None,
        "dbsnp_rs": None,
    })
    edge_bucket = defaultdict(lambda: {
        "pmids": set(),
        "notes": set(),
        "clinvar_data": None,
        "count": 0
    })
    chain_pos: dict[str, int] = {}

    for _, row in df.iterrows():
        var = row["variant_info"]
        name_label = row.get("Name", "")
        if not protein_uniprot_id:
            domain_protein_id = str(row.get("DomainProteinID", "")).strip()
            if domain_protein_id and domain_protein_id.lower() != "nan":
                protein_uniprot_id = domain_protein_id.split(";")[0].strip()
        protein_pos = _variant_aa_position(var, name_label)
        if protein_pos is not None:
            prev = variant_meta[var]["protein_pos"]
            variant_meta[var]["protein_pos"] = protein_pos if prev is None else min(prev, protein_pos)

        all_conditions = []
        for i in range(1, 11):
            disease = row.get(f"disease_{i}", "")
            if disease and 'not provided' not in disease.lower():
                all_conditions.append(disease)
        clinvar_data = None
        if all_conditions:
            clinvar_data = {
                "pathogenicity": row.get("significance_1", "N/A"),
                "review": format_star_rating(row.get("star_1", 0.0)),
                "conditions": "; ".join(all_conditions)
            }
        variant_meta[var]["clinvar_data"] = choose_best_clinvar(variant_meta[var]["clinvar_data"], clinvar_data)
        allele_id = str(row.get("#AlleleID", "")).strip()
        if allele_id and allele_id.lower() != "nan":
            try:
                variant_meta[var]["allele_id"] = str(int(float(allele_id)))
            except (TypeError, ValueError):
                variant_meta[var]["allele_id"] = allele_id
        dbsnp_rs = str(row.get("RS# (dbSNP)", "")).strip()
        if dbsnp_rs and dbsnp_rs.lower() != "nan":
            try:
                variant_meta[var]["dbsnp_rs"] = str(int(float(dbsnp_rs)))
            except (TypeError, ValueError):
                variant_meta[var]["dbsnp_rs"] = dbsnp_rs

        # Keep domain information on variants, but hide domain nodes in layout.
        features = [f.strip() for f in (row.get("DomainFeature", "") or "").split(';') if f.strip()]
        notes = [n.strip() for n in (row.get("DomainNote", "") or "").split(';') if n.strip()]
        for d in features:
            if d != "CHAIN":
                variant_meta[var]["domains"].add(d)
        for note in notes:
            if note and note != "1":
                variant_meta[var]["domain_notes"].add(note)

        # Baseline edge: protein -> variant
        pv_key = (prot, var, "PV")
        edge_bucket[pv_key]["count"] += 1
        pmid_val = str(row.get("pmid", "")).strip()
        if pmid_val:
            edge_bucket[pv_key]["pmids"].add(pmid_val)
        edge_bucket[pv_key]["clinvar_data"] = choose_best_clinvar(
            edge_bucket[pv_key]["clinvar_data"], variant_meta[var]["clinvar_data"]
        )

        # Causal chain from variant onwards.
        src = var
        _hop = 0
        for seg in str(row.get("chain", "")).split(" -[")[1:]:
            if "]->" not in seg:
                continue
            rel, tgt = seg.split("]->", 1)
            rel = rel.strip()
            tgt = tgt.strip()
            if not tgt:
                continue
            _hop += 1
            chain_pos[tgt] = max(chain_pos.get(tgt, 0), _hop)
            key = (src, tgt, rel)
            edge_bucket[key]["count"] += 1
            pmid = str(row.get("pmid", "")).strip()
            if pmid:
                edge_bucket[key]["pmids"].add(pmid)
            edge_bucket[key]["clinvar_data"] = choose_best_clinvar(edge_bucket[key]["clinvar_data"], clinvar_data)
            src = tgt

    for (u, v, rel), payload in edge_bucket.items():
        edge_attrs = {"relation": rel, "weight": payload["count"]}
        if payload["pmids"]:
            edge_attrs["pmid"] = "; ".join(sorted(payload["pmids"]))
        if payload["notes"]:
            edge_attrs["note"] = "; ".join(sorted(payload["notes"]))
        if payload["clinvar_data"]:
            edge_attrs["clinvar_data"] = payload["clinvar_data"]
        G.add_edge(u, v, **edge_attrs)

    subgraph_data: dict = {}


    if variant_aa_range is not None:
        lo, hi = float(variant_aa_range[0]), float(variant_aa_range[1])
        exclude = set()
        for v in variants:
            if not G.has_node(v):
                continue
            vm = variant_meta.get(v, {})
            p = vm.get("protein_pos")
            if p is not None:
                try:
                    pi = int(p)
                except (TypeError, ValueError):
                    pi = None
            else:
                pi = None
            if pi is None:
                po = _variant_aa_position(v, "")
                pi = int(po) if po is not None else None
            if pi is None or pi < lo or pi > hi:
                exclude.add(v)
        for v in exclude:
            if G.has_node(v):
                G.remove_node(v)
        if G.has_node(prot):
            seen = {prot}
            dq = deque([prot])
            while dq:
                u = dq.popleft()
                for _, v, _ in G.out_edges(u, data=True):
                    if v not in seen:
                        seen.add(v)
                        dq.append(v)
            for n in list(G.nodes()):
                if n not in seen:
                    G.remove_node(n)
        variants = {v for v in variants if G.has_node(v)}
        endpoints = {e for e in endpoints if G.has_node(e)}
        for k in list(chain_pos.keys()):
            if not G.has_node(k):
                del chain_pos[k]

    # ---------- Kind-aware layered layout with pseudo nodes ----------
    # Layer 0: protein + variants (fixed)
    # Layer -1: endpoints reachable ONLY directly from variants
    # Layers 1..N: longest-path depth plus extra spacing for repeated kinds
    # (e.g. protein -> protein or endpoint -> endpoint).

    def _node_kind(n: str) -> str:
        if n in variants:
            return "variant"
        if n in endpoints:
            return "endpoint"
        return "intermediate"

    node_kind = {n: _node_kind(n) for n in G.nodes()}

    # -- Step 1: identify "direct-only" endpoints → layer -1 ------------
    # An endpoint qualifies as "direct-only" when ALL its graph neighbours
    # (both predecessors and successors, excluding the root protein) are
    # variants.  If it also feeds into intermediates or other endpoints it
    # participates in longer chains and must go through normal layering.
    _direct_only_nodes: dict[str, int] = {}

    def _is_direct_only(node: str) -> bool:
        preds = {u for u, _, _ in G.in_edges(node, data=True) if u != prot}
        succs = {v for _, v, _ in G.out_edges(node, data=True)}
        if not preds:
            return False
        all_variant_preds = all(p in variants for p in preds)
        has_non_variant_succs = any(
            s not in variants and s != prot for s in succs
        )
        return all_variant_preds and not has_non_variant_succs

    for ep in endpoints:
        if _is_direct_only(ep):
            _direct_only_nodes[ep] = -1
    _direct_only_eps = {n for n in _direct_only_nodes if n in endpoints}

    # -- Step 2: layer seeds from chain positions -------------------------
    node_depth: dict[str, int] = {}
    for v in variants:
        node_depth[v] = 0
    node_depth[prot] = 0

    _init_depth: dict[str, int] = {}
    for n in G.nodes():
        if n == prot or n in variants or n in _direct_only_nodes:
            continue
        if n in chain_pos:
            _init_depth[n] = chain_pos[n]
        else:
            _init_depth[n] = 1

    # -- Step 3: condensation DAG propagation (cycle-safe) ---------------
    #    Collapse SCCs so cycles don't cascade layer depths, then align
    #    later kinds after the deepest preceding layer of the source kind.
    _non_vp = [
        n for n in G.nodes()
        if n != prot and n not in variants and n not in _direct_only_nodes
    ]
    if _non_vp:
        _H = nx.DiGraph()
        for n in _non_vp:
            _H.add_node(n)
        for u, v, _ in G.edges(data=True):
            if u in _H and v in _H:
                _H.add_edge(u, v)
        _C = nx.condensation(_H)
        _scc_map = _C.graph["mapping"]
        _scc_depth: dict[int, int] = {}
        for n in _non_vp:
            sid = _scc_map[n]
            _scc_depth[sid] = max(_scc_depth.get(sid, 0),
                                  _init_depth.get(n, 1))
        _topo_sids = list(nx.topological_sort(_C))
        for sid in _topo_sids:
            for succ_sid in _C.successors(sid):
                if _scc_depth.get(succ_sid, 0) <= _scc_depth.get(sid, 0):
                    _scc_depth[succ_sid] = _scc_depth[sid] + 1

        for n in _non_vp:
            node_depth[n] = _scc_depth[_scc_map[n]]

        # Rewrite endpoint depths: use only within-endpoint distance
        # so that all endpoints directly reachable from intermediates
        # land on the same (first) endpoint layer, and only
        # endpoint-to-endpoint chains create additional endpoint layers.
        _ep_nodes = {n for n in _non_vp if node_kind[n] == "endpoint"}
        if _ep_nodes:
            _ep_depth: dict[str, int] = {}
            for n in _ep_nodes:
                has_non_ep_pred = any(
                    node_kind.get(u, "variant") != "endpoint"
                    for u, _, _ in G.in_edges(n, data=True)
                )
                _ep_depth[n] = 0 if has_non_ep_pred else 1

            _ep_sub = nx.DiGraph()
            for n in _ep_nodes:
                _ep_sub.add_node(n)
            for u, v, _ in G.edges(data=True):
                if u in _ep_nodes and v in _ep_nodes:
                    _ep_sub.add_edge(u, v)

            if _ep_sub.edges:
                _ep_C = nx.condensation(_ep_sub)
                _ep_scc = _ep_C.graph["mapping"]
                _ep_scc_d: dict[int, int] = {}
                for n in _ep_nodes:
                    sid = _ep_scc[n]
                    _ep_scc_d[sid] = min(
                        _ep_scc_d.get(sid, 999), _ep_depth[n])
                for sid in nx.topological_sort(_ep_C):
                    for succ in _ep_C.successors(sid):
                        if _ep_scc_d.get(succ, 0) <= _ep_scc_d[sid]:
                            _ep_scc_d[succ] = _ep_scc_d[sid] + 1
                for n in _ep_nodes:
                    _ep_depth[n] = _ep_scc_d[_ep_scc[n]]

            for n in _ep_nodes:
                node_depth[n] = _ep_depth[n]

        # Kind-zone separation: shift each kind so zones don't overlap.
        # Canonical order: intermediate first, then endpoint.
        # When the kind-feed graph has cycles (e.g. endpoint <-> intermediate),
        # fall back to the canonical order instead of skipping the shift.
        _kind_layers: dict[str, set] = defaultdict(set)
        for n in _non_vp:
            _kind_layers[node_kind[n]].add(node_depth[n])

        _kind_feeds: set = set()
        for u, v, _ in G.edges(data=True):
            ku, kv = node_kind.get(u, ""), node_kind.get(v, "")
            if ku and kv and ku != kv and ku != "variant" and kv != "variant":
                _kind_feeds.add((ku, kv))

        _CANONICAL_KIND_ORDER = ["intermediate", "endpoint"]

        if _kind_feeds:
            _kind_dag = nx.DiGraph(list(_kind_feeds))
            if _kind_dag.nodes and nx.is_directed_acyclic_graph(_kind_dag):
                _ordered_kinds = list(nx.topological_sort(_kind_dag))
            else:
                _ordered_kinds = [k for k in _CANONICAL_KIND_ORDER
                                  if k in _kind_layers]

            _kind_shift: dict[str, int] = defaultdict(int)
            for i in range(len(_ordered_kinds) - 1):
                src_k = _ordered_kinds[i]
                dst_k = _ordered_kinds[i + 1]
                if not _kind_layers.get(src_k) or not _kind_layers.get(dst_k):
                    continue
                src_max = max(_kind_layers[src_k]) + _kind_shift[src_k]
                dst_min = min(_kind_layers[dst_k]) + _kind_shift[dst_k]
                if dst_min <= src_max:
                    _kind_shift[dst_k] += src_max + 1 - dst_min

            for n in _non_vp:
                node_depth[n] += _kind_shift.get(node_kind[n], 0)
    else:
        for n, d in _init_depth.items():
            node_depth[n] = d

    node_depth.update(_direct_only_nodes)

    # Safety net: no endpoint (except direct-only at -1) may share layer 0
    # with variants, and no endpoint may share a layer with intermediates.
    _max_int_depth = max(
        (node_depth[n] for n in node_depth
         if node_kind.get(n) == "intermediate"),
        default=0,
    )
    _min_ep_floor = max(_max_int_depth + 1, 1)
    for n in list(node_depth):
        if node_kind.get(n) == "endpoint" and n not in _direct_only_nodes:
            if node_depth[n] < _min_ep_floor:
                node_depth[n] = _min_ep_floor

    def get_layer(n):
        if n == prot or n in variants:
            return 0
        return node_depth.get(n, 1)

    # -- Step 4: collect all layers -------------------------------------
    layers: dict[int, list] = defaultdict(list)
    for n in G.nodes():
        layers[get_layer(n)].append(n)

    def _variant_sort_key(variant_name: str):
        protein_pos = variant_meta[variant_name]["protein_pos"]
        return (
            protein_pos is None,
            protein_pos if protein_pos is not None else math.inf,
            variant_name,
        )

    if 0 in layers:
        ordered_variants = sorted(
            [n for n in layers[0] if n in variants],
            key=_variant_sort_key,
        )
        layer_zero_other = [n for n in layers[0] if n not in variants and n != prot]
        layers[0] = ordered_variants + layer_zero_other + ([prot] if prot in layers[0] else [])

    # -- Step 5: crossing minimisation (LNS) ---------------------------
    _cross = {}
    for u, v, _d in G.edges(data=True):
        lu, lv = get_layer(u), get_layer(v)
        if lu == lv:
            continue
        pair = (u, v) if lu < lv else (v, u)
        _cross[pair] = _cross.get(pair, 0) + 1
    weighted_cross = [(u, v, w) for (u, v), w in _cross.items()]

    optimized = _optimize_layer_ordering(
        dict(layers),
        weighted_cross,
        fixed_layers={0},
    )

    # -- Step 6: assign positions from optimised orderings --------------
    x_pos: dict[str, float] = {}
    sorted_layer_keys = sorted(optimized.keys())
    max_layer = max(sorted_layer_keys) if sorted_layer_keys else 0
    min_layer = min(sorted_layer_keys) if sorted_layer_keys else 0
    layer_gap = 190.0
    y_pos = {
        li: li * layer_gap
        for li in range(min_layer, max_layer + 1)
    }

    for li, ordered in optimized.items():
        if li == 0 and prot in ordered:
            ordered = [n for n in ordered if n != prot]
            ordered.append(prot)
        n_nodes = len(ordered)
        sp = 160.0 if li == 0 else 210.0
        start = -((n_nodes - 1) * sp) / 2.0
        for i, n in enumerate(ordered):
            x_pos[n] = start + i * sp

    for n in G.nodes():
        if n not in x_pos:
            x_pos[n] = 0.0

    # -- Step 6.5: Barycentric repositioning for non-variant layers ----------
    # The uniform even-spacing above looks artificially symmetric.  Reposition
    # each node to the centroid of its "toward-layer-0" parents so the layout
    # reflects actual connection structure, while preserving the crossing-
    # minimised ordering the LNS found.  Process layers in order of increasing
    # distance from layer 0 so parent positions are finalised first.
    _max_half_px = max(
        (len(nodes) - 1) * 170.0 / 2.0
        for nodes in optimized.values() if len(nodes) > 1
    ) if any(len(v) > 1 for v in optimized.values()) else 170.0

    # For each node, collect the adjacent-layer neighbors that are closer to
    # layer 0 (the variant layer).  Those are the "parents" that anchor x.
    _parent_adj: dict[str, list] = defaultdict(list)
    for _u, _v, _d in G.edges(data=True):
        _lu, _lv = get_layer(_u), get_layer(_v)
        if _lu == _lv:
            continue
        if _lu == 0:
            _parent_adj[_v].append(_u)
        elif _lv == 0:
            _parent_adj[_u].append(_v)
        elif abs(_lu) < abs(_lv):
            _parent_adj[_v].append(_u)
        else:
            _parent_adj[_u].append(_v)

    _bary_min_sep = 80.0
    for _li in sorted((li for li in optimized if li != 0), key=abs):
        _ordered = list(optimized[_li])
        if not _ordered:
            continue
        # Target x = centroid of parents already placed in x_pos
        _bx: dict[str, float] = {}
        for _nd in _ordered:
            _nbrs = [_nb for _nb in _parent_adj.get(_nd, []) if _nb in x_pos]
            _bx[_nd] = (
                sum(x_pos[_nb] for _nb in _nbrs) / len(_nbrs)
                if _nbrs else x_pos.get(_nd, 0.0)
            )
        # Keep LNS order; just slide each node toward its centroid x
        _xs = [max(-_max_half_px, min(_max_half_px, _bx[_nd])) for _nd in _ordered]
        # Forward pass: enforce minimum separation
        for _idx in range(1, len(_xs)):
            if _xs[_idx] < _xs[_idx - 1] + _bary_min_sep:
                _xs[_idx] = _xs[_idx - 1] + _bary_min_sep
        # Shift left if right boundary exceeded
        if _xs and _xs[-1] > _max_half_px:
            _shift = _xs[-1] - _max_half_px
            _xs = [_x - _shift for _x in _xs]
        # Backward pass: fix left-side violations after the shift
        for _idx in range(len(_xs) - 2, -1, -1):
            if _xs[_idx] > _xs[_idx + 1] - _bary_min_sep:
                _xs[_idx] = _xs[_idx + 1] - _bary_min_sep
        for _nd, _x in zip(_ordered, _xs):
            x_pos[_nd] = _x

    pos = {n: (x_pos[n], y_pos.get(get_layer(n), 0.0)) for n in G.nodes()}

    def _css_safe(name: str) -> str:
        return re.sub(r'[^A-Za-z0-9_-]', '-', name)

    def _short_label(text, max_len=28):
        if len(text) <= max_len:
            return text
        cut = text[:max_len].rfind(' ')
        if cut < max_len // 2:
            cut = max_len
        return text[:cut].rstrip() + "..."

    raw_rel_types = sorted({
        d['relation'] for _, _, d in G.edges(data=True)
        if d['relation'] not in {'PV', 'DV', 'has_domain'}
    })
    palette = ["#e74c3c", "#2ecc71", "#3498db", "#f39c12", "#9b59b6"]
    rel_color_safe = {_css_safe(r): palette[i % len(palette)] for i, r in enumerate(raw_rel_types)}
    rel_display = {_css_safe(r): r for r in raw_rel_types}

    els = []
    for n, (x, y) in pos.items():
        layer = get_layer(n)

        if n == prot:
            size = 90
        elif n in endpoints or n in _direct_only_eps:
            size = 52 + min(24, 5 * endpoint_freq.get(n, 1))
        else:
            size = 54

        if n == prot:
            label = n
            role_class = "role-protein"
        elif n in variants:
            label = n
            role_class = "role-variant"
        elif n in endpoints or n in _direct_only_eps:
            label = _short_label(n)
            role_class = "role-endpoint"
        else:
            label = _short_label(n)
            role_class = "role-intermediate"

        node_el = {
            "data": {
                "id": n,
                "label": label,
                "real": n,
                "role": role_class.replace("role-", ""),
                "motif_type": "",
                "motif_members": "",
                "domain_notes": "; ".join(sorted(variant_meta[n]["domain_notes"])) if n in variant_meta else ""
            },
            "classes": f"L{layer} {role_class}",
            "style": {"width": size, "height": size}
        }
        if n == prot:
            node_el["data"]["uniprot_id"] = protein_uniprot_id or prot
        if n in variants and n in variant_meta:
            vm = variant_meta[n]
            node_el["data"]["gene_symbol"] = prot
            if vm.get("protein_pos") is not None:
                node_el["data"]["protein_pos"] = vm["protein_pos"]
            if vm["allele_id"]:
                node_el["data"]["clinvar_allele"] = vm["allele_id"]
            if vm["dbsnp_rs"]:
                node_el["data"]["dbsnp_rs"] = vm["dbsnp_rs"]
            if vm["clinvar_data"]:
                node_el["data"]["clinvar_data"] = vm["clinvar_data"]

        if n in pos:
            node_el["position"] = {"x": pos[n][0], "y": pos[n][1]}
        els.append(node_el)

    for u, v, d in G.edges(data=True):
        relation = d.get('relation', '')
        if relation == 'PV': cls = 'edge-PV'
        elif relation == 'DV': cls = 'edge-DV'
        elif relation == 'has_domain': cls = 'edge-has_domain'
        else: cls = f"edge-{_css_safe(relation)}"

        src4indra = prot if u in variants else u
        edge_data = {
            "id": f"{u}->{v}_{d.get('pmid', '')}_{d.get('note', '')}",
            "source": u,
            "target": v,
            "rel": relation,
            "src4indra": src4indra,
        }
        if 'pmid' in d and d['pmid']:
            edge_data['pmid'] = d['pmid']
        if 'note' in d and d['note']:
            edge_data['note'] = d['note']
        if 'clinvar_data' in d and d['clinvar_data']:
            edge_data['clinvar_data'] = d['clinvar_data']
        if 'weight' in d:
            edge_data['evidence_count'] = d['weight']

        els.append({"data": edge_data, "classes": cls})

    edge_set = {(u, v) for u, v, _ in G.edges(data=True)}

    legend_rels = ['Gene–variant'] + raw_rel_types
    legend_colors = {
        'Gene–variant': '#c9c4bf',
        **{rel_display.get(k, k): v for k, v in rel_color_safe.items()}
    }

    return els, legend_rels, legend_colors, rel_color_safe, list(edge_set), subgraph_data


def build_endpoint_elements(endpoint: str):
    return copy.deepcopy(_build_endpoint_elements_cached(endpoint))


@lru_cache(maxsize=32)
def _build_endpoint_elements_cached(endpoint: str):
    prot_stats = ENDPOINT_INDEX.get(endpoint, {})
    if not prot_stats:
        return [], [], {}, {}, [], {}

    def _css_safe(name: str) -> str:
        return re.sub(r'[^A-Za-z0-9_-]', '-', name)

    def _short_label(text, max_len=28):
        if len(text) <= max_len:
            return text
        cut = text[:max_len].rfind(' ')
        if cut < max_len // 2:
            cut = max_len
        return text[:cut].rstrip() + "..."

    def choose_best_clinvar(existing: Optional[dict], candidate: Optional[dict]):
        if not candidate:
            return existing
        if not existing:
            return candidate
        prev_len = len(existing.get("conditions", ""))
        cand_len = len(candidate.get("conditions", ""))
        if cand_len > prev_len:
            return candidate
        return existing

    protein_nodes: dict[str, dict] = {}
    variant_nodes: dict[str, dict] = {}
    intermediate_nodes: dict[str, dict] = {}
    edge_bucket: dict[tuple, dict] = defaultdict(lambda: {
        "count": 0, "pmids": set(), "notes": set(),
    })
    edge_set: set[tuple[str, str]] = set()

    read_columns = [
        "biological_process/disease",
        "variant_info",
        "chain",
        "pmid",
        "#AlleleID",
        "RS# (dbSNP)",
        "significance_1",
        "star_1",
        *[f"disease_{i}" for i in range(1, 11)],
    ]

    for prot in sorted(prot_stats, key=_sort_text):
        tsv_path = TSV_FILES[prot]
        try:
            df = pd.read_csv(tsv_path, sep="\t", usecols=read_columns).fillna('')
        except ValueError:
            df = pd.read_csv(tsv_path, sep="\t").fillna('')

        if "biological_process/disease" not in df.columns or "variant_info" not in df.columns:
            continue

        endpoint_series = df["biological_process/disease"].astype(str).str.strip()
        sub = df[endpoint_series.eq(endpoint)]
        if sub.empty:
            continue

        protein_entry = protein_nodes.setdefault(prot, {
            "row_count": 0,
            "variant_ids": set(),
        })

        for _, row in sub.iterrows():
            variant_label = str(row.get("variant_info", "")).strip()
            if not variant_label:
                continue

            variant_id = f"{prot}::{variant_label}"
            protein_entry["row_count"] += 1
            protein_entry["variant_ids"].add(variant_id)

            variant_entry = variant_nodes.setdefault(variant_id, {
                "label": variant_label,
                "real": f"{prot} {variant_label}",
                "protein": prot,
                "row_count": 0,
                "pmids": set(),
                "clinvar_data": None,
                "allele_id": None,
                "dbsnp_rs": None,
            })
            variant_entry["row_count"] += 1

            pmid_val = str(row.get("pmid", "")).strip()
            if pmid_val:
                variant_entry["pmids"].add(pmid_val)

            all_conditions = []
            for i in range(1, 11):
                disease = str(row.get(f"disease_{i}", "")).strip()
                if disease and 'not provided' not in disease.lower():
                    all_conditions.append(disease)
            clinvar_data = None
            if all_conditions:
                clinvar_data = {
                    "pathogenicity": row.get("significance_1", "N/A"),
                    "review": format_star_rating(row.get("star_1", 0.0)),
                    "conditions": "; ".join(all_conditions),
                }
            variant_entry["clinvar_data"] = choose_best_clinvar(
                variant_entry["clinvar_data"], clinvar_data
            )

            allele_id = str(row.get("#AlleleID", "")).strip()
            if allele_id and allele_id.lower() != "nan":
                try:
                    variant_entry["allele_id"] = str(int(float(allele_id)))
                except (TypeError, ValueError):
                    variant_entry["allele_id"] = allele_id

            dbsnp_rs = str(row.get("RS# (dbSNP)", "")).strip()
            if dbsnp_rs and dbsnp_rs.lower() != "nan":
                try:
                    variant_entry["dbsnp_rs"] = str(int(float(dbsnp_rs)))
                except (TypeError, ValueError):
                    variant_entry["dbsnp_rs"] = dbsnp_rs

            pv_key = (prot, variant_id, "PV")
            edge_bucket[pv_key]["count"] += 1
            if pmid_val:
                edge_bucket[pv_key]["pmids"].add(pmid_val)
            edge_set.add((prot, variant_id))

            chain_str = str(row.get("chain", ""))
            segs = chain_str.split(" -[")[1:]
            src = variant_id
            for seg in segs:
                if "]->" not in seg:
                    continue
                rel, tgt = seg.split("]->", 1)
                rel = rel.strip()
                tgt = tgt.strip()
                if not tgt:
                    continue

                if tgt == endpoint:
                    tgt_id = endpoint
                elif tgt == prot:
                    tgt_id = prot
                elif tgt in prot_stats:
                    # This target also has its own variants for this disease,
                    # so it already exists as a protein node.  Give it a
                    # distinct intermediate-role ID so both can coexist.
                    tgt_id = f"{tgt}::inode"
                    if tgt_id not in intermediate_nodes:
                        intermediate_nodes[tgt_id] = {
                            "label": _short_label(tgt),
                            "real": tgt,
                            "row_count": 0,
                            "is_known_protein": True,
                        }
                    intermediate_nodes[tgt_id]["row_count"] += 1
                else:
                    tgt_id = tgt
                    if tgt_id not in intermediate_nodes:
                        is_known_prot = tgt_id in prot_stats or tgt_id in TSV_FILES
                        intermediate_nodes[tgt_id] = {
                            "label": _short_label(tgt_id),
                            "real": tgt_id,
                            "row_count": 0,
                            "is_known_protein": is_known_prot,
                        }
                    intermediate_nodes[tgt_id]["row_count"] += 1

                ekey = (src, tgt_id, rel)
                edge_bucket[ekey]["count"] += 1
                if pmid_val:
                    edge_bucket[ekey]["pmids"].add(pmid_val)
                edge_set.add((src, tgt_id))
                src = tgt_id

    ordered_proteins = sorted(
        protein_nodes.items(),
        key=lambda kv: (-len(kv[1]["variant_ids"]), -kv[1]["row_count"], kv[0].casefold())
    )
    ordered_variants = sorted(
        variant_nodes.items(),
        key=lambda kv: (-kv[1]["row_count"], kv[1]["protein"].casefold(), kv[1]["label"].casefold())
    )
    # Map variant_id → gene symbol so INDRA links use the gene name, not the variant label.
    _variant_gene: dict[str, str] = {vid: info["protein"] for vid, info in ordered_variants}

    total_proteins = len(ordered_proteins)
    total_variants = len(ordered_variants)
    total_records = sum(info["row_count"] for info in protein_nodes.values())

    raw_rel_types = sorted({
        rel for (_s, _t, rel) in edge_bucket if rel != "PV"
    })
    palette = ["#e74c3c", "#2ecc71", "#3498db", "#f39c12", "#9b59b6"]
    rel_color_safe = {_css_safe(r): palette[i % len(palette)] for i, r in enumerate(raw_rel_types)}
    rel_display = {_css_safe(r): r for r in raw_rel_types}

    # --- Layered layout --------------------------------------------------
    # Indirect paths (gene→variant→intermediates→disease) go BELOW disease.
    # Direct paths (gene→variant→disease, no intermediates) go ABOVE disease.
    #
    # Bottom section:  layer 0  = indirect genes
    #                  layer 1  = indirect variants
    #                  layer 2+ = intermediate nodes (longest-path depth)
    # Middle:          layer N  = disease endpoint
    # Top section:     layer N+1 = direct genes
    #                  layer N+2 = direct variants

    _G_ep = nx.DiGraph()
    for _s, _t in edge_set:
        _G_ep.add_edge(_s, _t)

    # Classify genes: direct-only if every variant edge leads straight to endpoint
    _direct_genes: set[str] = set()
    for _prot in protein_nodes:
        _my_vars = {vid for vid in variant_nodes
                    if variant_nodes[vid]["protein"] == _prot}
        _out = {tgt for (src, tgt) in edge_set if src in _my_vars}
        if _out and _out <= {endpoint}:
            _direct_genes.add(_prot)
    _indirect_genes = set(protein_nodes.keys()) - _direct_genes

    # Nodes whose layer is set from the start (will be excluded from
    # the condensation-DAG propagation for intermediates)
    _direct_nodes: set[str] = _direct_genes | {
        vid for vid in variant_nodes
        if variant_nodes[vid]["protein"] in _direct_genes
    }

    node_layer: dict[str, int] = {}
    for _prot in _indirect_genes:
        node_layer[_prot] = 0
    for _vid in variant_nodes:
        if variant_nodes[_vid]["protein"] not in _direct_genes:
            node_layer[_vid] = 1

    # Longest-path depth for intermediate nodes (cycle-safe via condensation)
    _non_fixed = [n for n in _G_ep.nodes()
                  if n not in node_layer and n not in _direct_nodes and n != endpoint]
    if _non_fixed:
        _H = nx.DiGraph()
        _H.add_nodes_from(_non_fixed)
        for _u, _v in _G_ep.edges():
            if _u in _H and _v in _H:
                _H.add_edge(_u, _v)
        _init: dict[str, int] = {}
        for n in _non_fixed:
            _kp = [p for p in _G_ep.predecessors(n) if p in node_layer]
            _init[n] = (max(node_layer[p] for p in _kp) + 1) if _kp else 2
        _C = nx.condensation(_H)
        _scc_map = _C.graph["mapping"]
        _scc_d: dict[int, int] = {}
        for n in _non_fixed:
            sid = _scc_map[n]
            _scc_d[sid] = max(_scc_d.get(sid, 0), _init.get(n, 2))
        for sid in nx.topological_sort(_C):
            for succ_sid in _C.successors(sid):
                if _scc_d.get(succ_sid, 0) <= _scc_d.get(sid, 0):
                    _scc_d[succ_sid] = _scc_d[sid] + 1
        for n in _non_fixed:
            node_layer[n] = _scc_d[_scc_map[n]]

    # Disease endpoint sits above all indirect/intermediate layers
    _ep_layer = max(node_layer.values(), default=1) + 1
    node_layer[endpoint] = _ep_layer

    # Direct paths placed above disease; variant layer is closer to disease,
    # gene layer is one step further – mirroring the indirect section structure.
    for _prot in _direct_genes:
        node_layer[_prot] = _ep_layer + 2
    for _vid in variant_nodes:
        if variant_nodes[_vid]["protein"] in _direct_genes:
            node_layer[_vid] = _ep_layer + 1

    # --- Build layer groups with crossing-minimising orderings -----------
    _layers: dict[int, list] = defaultdict(list)
    for n, li in node_layer.items():
        _layers[li].append(n)

    def _order_gene_variant_layers(gene_layer: int, var_layer: int,
                                   genes: set, var_rank_source: list):
        """Sort variants grouped by gene; sort genes by variant centroid."""
        vl = sorted(
            _layers.get(var_layer, []),
            key=lambda vid: (
                variant_nodes[vid]["protein"].casefold() if vid in variant_nodes else "",
                vid.casefold(),
            ),
        )
        _layers[var_layer] = vl
        vrank = {vid: i for i, vid in enumerate(vl)}
        def _centroid(g):
            mv = [v for v in vl if v in variant_nodes and variant_nodes[v]["protein"] == g]
            return (sum(vrank[v] for v in mv) / len(mv)) if mv else 0.0
        _layers[gene_layer] = sorted(_layers.get(gene_layer, []), key=_centroid)
        return vl

    _indirect_var_list = _order_gene_variant_layers(0, 1, _indirect_genes, [])
    _direct_var_list   = _order_gene_variant_layers(
        _ep_layer + 2, _ep_layer + 1, _direct_genes, []
    )

    # Intermediate and endpoint layers: alphabetical baseline
    for li in _layers:
        if li not in {0, 1, _ep_layer + 1, _ep_layer + 2}:
            _layers[li].sort(key=lambda n: n.casefold())

    # LNS crossing optimisation – fix all gene/variant layers, optimise the rest
    _fixed_layers = {0, 1}
    if _direct_genes:
        _fixed_layers |= {_ep_layer + 1, _ep_layer + 2}

    _cross_w: dict[tuple, int] = {}
    for _s, _t in edge_set:
        ls, lt = node_layer.get(_s, 0), node_layer.get(_t, 0)
        if ls == lt:
            continue
        pair = (_s, _t) if ls < lt else (_t, _s)
        _cross_w[pair] = _cross_w.get(pair, 0) + 1
    _optimized = _optimize_layer_ordering(
        dict(_layers),
        [(u, v, w) for (u, v), w in _cross_w.items()],
        fixed_layers=_fixed_layers,
    )

    _layer_gap = 190.0
    _sp = 210.0
    pos: dict[str, tuple[float, float]] = {}
    for li, ordered in _optimized.items():
        n_nodes = len(ordered)
        start = -((n_nodes - 1) * _sp) / 2.0
        y = li * _layer_gap
        for i, n in enumerate(ordered):
            pos[n] = (start + i * _sp, y)

    # Align gene x to exact centroid of its variants' x-positions (both sections)
    for _var_list_section in (_indirect_var_list, _direct_var_list):
        for prot in protein_nodes:
            mv = [v for v in _var_list_section
                  if v in variant_nodes and variant_nodes[v]["protein"] == prot]
            xs = [pos[v][0] for v in mv if v in pos]
            if xs:
                pos[prot] = (sum(xs) / len(xs), pos.get(prot, (0, 0))[1])

    # --- Barycentric x for intermediate layers ----------------------------
    # Sparse intermediate layers tend to cluster in the centre; instead,
    # position each node at the average x of its lower-layer neighbours,
    # bounded by the span of the widest layer.
    _max_half = max(
        (len(nodes) - 1) * _sp / 2.0
        for nodes in _optimized.values()
        if len(nodes) > 1
    ) if any(len(v) > 1 for v in _optimized.values()) else _sp

    _adj_lo: dict[str, list] = defaultdict(list)
    for _s, _t in edge_set:
        if node_layer.get(_s, 0) < node_layer.get(_t, 0):
            _adj_lo[_t].append(_s)

    _fixed_li = {0, 1, _ep_layer, _ep_layer + 1, _ep_layer + 2}
    _inter_layers = sorted(li for li in _optimized if li not in _fixed_li)
    _min_sep = 80.0

    for li in _inter_layers:
        ordered = list(_optimized[li])
        if not ordered:
            continue
        y = li * _layer_gap

        # Barycentric x from lower-layer neighbours; fall back to current x
        bary: dict[str, float] = {}
        for nd in ordered:
            nbrs = [nb for nb in _adj_lo.get(nd, []) if nb in pos]
            bary[nd] = (
                sum(pos[nb][0] for nb in nbrs) / len(nbrs)
                if nbrs else pos.get(nd, (0.0, 0.0))[0]
            )

        # Sort by barycentric x, then enforce minimum spacing within bounds
        sorted_nodes = sorted(ordered, key=lambda nd: bary[nd])
        xs = [max(-_max_half, min(_max_half, bary[nd])) for nd in sorted_nodes]

        # Forward pass: minimum separation
        for idx in range(1, len(xs)):
            if xs[idx] < xs[idx - 1] + _min_sep:
                xs[idx] = xs[idx - 1] + _min_sep
        # Shift left if right boundary exceeded
        if xs and xs[-1] > _max_half:
            shift = xs[-1] - _max_half
            xs = [x - shift for x in xs]
        # Backward pass: fix any left-side violations after shift
        for idx in range(len(xs) - 2, -1, -1):
            if xs[idx] > xs[idx + 1] - _min_sep:
                xs[idx] = xs[idx + 1] - _min_sep

        for nd, x in zip(sorted_nodes, xs):
            pos[nd] = (x, y)

    def _ep_el(nid, data, role_cls, size):
        li = node_layer.get(nid, 0)
        el = {
            "data": data,
            "classes": f"L{li} {role_cls}",
            "style": {"width": size, "height": size},
        }
        if nid in pos:
            el["position"] = {"x": pos[nid][0], "y": pos[nid][1]}
        return el

    els = [_ep_el(endpoint, {
        "id": endpoint,
        "label": _short_label(endpoint, max_len=36),
        "real": endpoint,
        "role": "endpoint",
        "n_proteins": total_proteins,
        "n_variants": total_variants,
        "n_records": total_records,
    }, "role-endpoint", 110)]

    for prot, info in ordered_proteins:
        size = 48 + min(26, 5 * math.sqrt(max(len(info["variant_ids"]), 1)))
        els.append(_ep_el(prot, {
            "id": prot,
            "label": prot,
            "real": prot,
            "role": "protein",
            "n_variants": len(info["variant_ids"]),
            "n_records": info["row_count"],
            "protein_page": _protein_href(prot),
        }, "role-protein", size))

    for variant_id, info in ordered_variants:
        size = 42 + min(18, 4 * math.sqrt(max(info["row_count"], 1)))
        node_data = {
            "id": variant_id,
            "label": info["label"],
            "real": info["real"],
            "role": "variant",
            "gene_symbol": info["protein"],
            "n_records": info["row_count"],
            "protein_page": _protein_href(info["protein"]),
        }
        if info["allele_id"]:
            node_data["clinvar_allele"] = info["allele_id"]
        if info["dbsnp_rs"]:
            node_data["dbsnp_rs"] = info["dbsnp_rs"]
        if info["clinvar_data"]:
            node_data["clinvar_data"] = info["clinvar_data"]
        els.append(_ep_el(variant_id, node_data, "role-variant", size))

    for nid, info in sorted(intermediate_nodes.items(), key=lambda kv: kv[0].casefold()):
        size = 44 + min(20, 3 * math.sqrt(max(info["row_count"], 1)))
        node_data = {
            "id": nid,
            "label": info["label"],
            "real": info["real"],
            "role": "intermediate",
            "n_records": info["row_count"],
        }
        if info["is_known_protein"]:
            node_data["protein_page"] = _protein_href(info["real"])
        els.append(_ep_el(nid, node_data, "role-intermediate", size))

    for (src, tgt, rel), payload in sorted(edge_bucket.items()):
        if rel == "PV":
            cls = "edge-PV"
        else:
            cls = f"edge-{_css_safe(rel)}"
        edge_data = {
            "id": f"{src}->{tgt}::{rel}",
            "source": src,
            "target": tgt,
            "rel": rel,
            "src4indra": _variant_gene.get(src, src.split("::")[0]),
            "evidence_count": payload["count"],
        }
        if payload["pmids"]:
            pmids = sorted(payload["pmids"])
            edge_data["pmid"] = "; ".join(pmids[:20])
        els.append({
            "data": edge_data,
            "classes": cls,
        })

    legend_rels = ["Gene–variant"] + raw_rel_types
    legend_colors = {
        "Gene–variant": "#c9c4bf",
        **{rel_display.get(k, k): v for k, v in rel_color_safe.items()},
    }
    return els, legend_rels, legend_colors, rel_color_safe, list(edge_set), {}


def _cy_net_layout_preset() -> dict:
    return {"name": "preset"}


def _adjacency_from_elements(els: list) -> tuple[dict, dict]:
    """Directed forward and reverse adjacency from edge elements only."""
    fwd: dict = defaultdict(set)
    rev: dict = defaultdict(set)
    for el in els:
        d = el.get("data") or {}
        if "source" not in d:
            continue
        s, t = d["source"], d["target"]
        fwd[s].add(t)
        rev[t].add(s)
    return fwd, rev
