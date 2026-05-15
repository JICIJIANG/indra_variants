import copy
import re

import dash
from dash import Input, MATCH, Output, State, ctx, dcc, html

from indra_variants.app.data_index import (
    ENDPOINT_INDEX,
    ENDPOINTS,
    PROTS,
    TSV_FILES,
)
from indra_variants.app.figures import (
    _norm_map_range,
    _protein_lollipop_figure,
    _relayout_xaxis_range,
)
from indra_variants.app.graph_builders import (
    _adjacency_from_elements,
    _cy_net_layout_preset,
    build_elements,
)
from indra_variants.app.pages import (
    endpoint_network_page,
    network_page,
    search_page,
    statistics_page,
)
from indra_variants.app.sidebar import _build_edge_info, _build_node_info
from indra_variants.app.stats import (
    STATS_FIG_HEIGHT_9,
    STATS_FIG_HEIGHT_10,
    STATS_FIG_HEIGHT_20,
    _STATS_BP_ROWS,
    _STATS_DISEASE_ROWS,
    _STATS_GENE_ROWS,
    _build_alpha_directory,
    _stats_bar_figure,
    _stats_point_href,
    _stats_value_bp_disease,
    _stats_value_gene,
)
from indra_variants.app.utils import (
    GRAPH_ENDPOINT_BG,
    GRAPH_ENDPOINT_FG,
    GRAPH_INTERMEDIATE_BG,
    GRAPH_INTERMEDIATE_FG,
    GRAPH_NODE_BG_OPACITY,
    GRAPH_PROTEIN_BG,
    GRAPH_PROTEIN_FG,
    GRAPH_VARIANT_BG,
    GRAPH_VARIANT_FG,
    U,
    _css_safe_global,
    _decode_route_value,
    _endpoint_href,
    _protein_href,
)


def register(app):
    """Register every @app.callback on the provided Dash ``app`` instance."""

    @app.callback(Output("page", "children"), Input("url", "pathname"))
    def router(path):
        if path in (None, "/", "/search"):
            return search_page()
        if path == "/statistics":
            return statistics_page()
        if path.startswith("/protein/"):
            prot = _decode_route_value(path.split("/protein/", 1)[1])
            if prot in PROTS:
                return network_page(prot)
        if path.startswith("/endpoint/"):
            endpoint = _decode_route_value(path.split("/endpoint/", 1)[1])
            if endpoint in ENDPOINT_INDEX:
                return endpoint_network_page(endpoint)
        return html.Div(
            [
                html.H3("404 Not found",
                        style={'color': U['ink'], 'fontFamily': U['font_display'],
                               'fontWeight': 600}),
                html.P(
                    [
                        dcc.Link("Browser", href="/",
                                 style={'marginRight': 16, 'color': U['link'],
                                        'fontWeight': 600}),
                        dcc.Link("Statistics", href="/statistics",
                                 style={'color': U['link'], 'fontWeight': 600}),
                    ],
                    style={'fontSize': 15, 'fontFamily': U['font_ui']},
                ),
            ],
            style={'maxWidth': 560, 'margin': '80px auto',
                   'fontFamily': U['font_ui'], 'background': U['paper'],
                   'minHeight': '100vh', 'padding': '0 20px'},
        )

    @app.callback(
        Output({'type': 'edge-info', 'prot': MATCH}, 'children'),
        Input({'type': 'cy-net', 'prot': MATCH}, 'tapNodeData'),
        Input({'type': 'cy-net', 'prot': MATCH}, 'tapEdgeData'),
        prevent_initial_call=True)
    def show_sidebar_info(node, edge):
        if not dash.ctx.triggered:
            return dash.no_update

        prop = dash.ctx.triggered[0]['prop_id'].split('.')[-1]
        if prop == 'tapNodeData' and node:
            return _build_node_info(node)
        if prop == 'tapEdgeData' and edge:
            return _build_edge_info(edge)
        return dash.no_update

    @app.callback(
        Output({'type': 'gene-map', 'prot': MATCH}, 'figure'),
        Input({'type': 'gene-map-reset', 'prot': MATCH}, 'n_clicks'),
        Input({'type': 'gene-map-clin-filter', 'prot': MATCH}, 'value'),
        State({'type': 'gene-map', 'prot': MATCH}, 'id'),
        prevent_initial_call=True)
    def update_gene_variant_map(_n_clicks, clin_filter, gid):
        """Regenerate lollipop on reset or ClinVar filter change."""
        rid = (gid or {}).get('prot', '')
        if not rid.startswith('protein::'):
            return dash.no_update
        gene = rid.split('::', 1)[1]
        fig = _protein_lollipop_figure(gene, clin_filter=clin_filter)
        return fig if fig is not None else dash.no_update

    @app.callback(
        Output({'type': 'store-map-range', 'prot': MATCH}, 'data'),
        Input({'type': 'gene-map', 'prot': MATCH}, 'relayoutData'),
        Input({'type': 'gene-map-reset', 'prot': MATCH}, 'n_clicks'),
        prevent_initial_call=True)
    def sync_gene_map_x_range(relayout, _n_reset):
        trig = ctx.triggered_id
        if isinstance(trig, dict) and trig.get('type') == 'gene-map-reset':
            return None
        xr = _relayout_xaxis_range(relayout if isinstance(relayout, dict) else None)
        if xr is None:
            return dash.no_update
        return [xr[0], xr[1]]

    # ---------------------- subgraph modal callback ----------------------------
    @app.callback(
        [Output({'type': 'subgraph-modal', 'prot': MATCH}, 'is_open'),
         Output({'type': 'subgraph-title', 'prot': MATCH}, 'children'),
         Output({'type': 'cy-subgraph',    'prot': MATCH}, 'elements'),
         Output({'type': 'cy-subgraph',    'prot': MATCH}, 'stylesheet')],
        Input({'type': 'cy-net', 'prot': MATCH}, 'tapNodeData'),
        [State({'type': 'store-subgraphs',  'prot': MATCH}, 'data'),
         State({'type': 'store-relcolors',  'prot': MATCH}, 'data')],
        prevent_initial_call=True)
    def open_subgraph_modal(node, subgraphs, rel_colors):
        empty = (False, "", [], [])
        if not node or not subgraphs:
            return empty
        nid = node.get("id", "")
        if nid not in subgraphs:
            return empty

        sg = subgraphs[nid]
        parent_name = sg["parent_name"]
        nodes_map = sg["nodes"]
        edges_list = sg["edges"]
        rel_colors = rel_colors or {}

        role_to_layer = {"protein": 0, "variant": 0, "intermediate": 1, "endpoint": 2}
        layer_y = {0: 400, 1: 220, 2: 40}
        role_cls = {
            "protein": "role-protein", "variant": "role-variant",
            "intermediate": "role-intermediate", "endpoint": "role-endpoint",
        }
        role_size = {"protein": 70, "variant": 50, "intermediate": 46, "endpoint": 50}

        buckets = {0: [], 1: [], 2: []}
        for nid_key, info in nodes_map.items():
            ly = role_to_layer.get(info["role"], 2)
            buckets[ly].append(nid_key)

        x_pos = {}
        for ly, ids in buckets.items():
            spacing = 160
            start = -((len(ids) - 1) * spacing) / 2.0
            for i, n in enumerate(sorted(ids)):
                x_pos[n] = start + i * spacing

        sub_els = []
        n_endpoints = len(buckets[2])
        for nid_key, info in nodes_map.items():
            ly = role_to_layer.get(info["role"], 2)
            cls = role_cls.get(info["role"], "role-endpoint")
            sz = role_size.get(info["role"], 46)
            sub_els.append({
                "data": {"id": nid_key, "label": info["label"]},
                "position": {"x": x_pos.get(nid_key, 0), "y": layer_y[ly]},
                "classes": cls,
                "style": {"width": sz, "height": sz}
            })

        for e in edges_list:
            rel = e.get("rel", "")
            if rel == "PV":
                cls = "edge-PV"
            else:
                cls = f"edge-{_css_safe_global(rel)}"
            edge_data = {
                "source": e["source"], "target": e["target"],
                "rel": rel, "label": rel,
            }
            for k in ("pmid", "note", "clinvar_data", "evidence_count"):
                if k in e and e[k]:
                    edge_data[k] = e[k]
            sub_els.append({"data": edge_data, "classes": cls})

        stylesheet = [
            {"selector": "node", "style": {
                "label": "data(label)",
                "text-valign": "center", "text-halign": "center",
                "font-size": 11, "font-weight": "bold",
                "background-opacity": GRAPH_NODE_BG_OPACITY,
                "text-wrap": "wrap", "text-max-width": 110}},
            {"selector": ".role-protein", "style": {
                "background-color": GRAPH_PROTEIN_BG, "color": GRAPH_PROTEIN_FG}},
            {"selector": ".role-variant", "style": {
                "background-color": GRAPH_VARIANT_BG, "color": GRAPH_VARIANT_FG}},
            {"selector": ".role-intermediate", "style": {
                "background-color": GRAPH_INTERMEDIATE_BG,
                "color": GRAPH_INTERMEDIATE_FG}},
            {"selector": ".role-endpoint", "style": {
                "background-color": GRAPH_ENDPOINT_BG,
                "color": GRAPH_ENDPOINT_FG}},
            {"selector": ".edge-PV", "style": {
                "line-color": "#c9c4bf", "target-arrow-color": "#c9c4bf",
                "target-arrow-shape": "triangle", "curve-style": "bezier",
                "width": 1.5, "label": "data(label)",
                "font-size": 9, "color": "#a8a29e", "text-rotation": "autorotate"}},
        ]
        for css_cls, color in rel_colors.items():
            stylesheet.append({
                "selector": f".edge-{css_cls}", "style": {
                    "line-color": color, "target-arrow-color": color,
                    "target-arrow-shape": "triangle", "curve-style": "bezier",
                    "width": 2, "label": "data(label)",
                    "font-size": 9, "color": color, "text-rotation": "autorotate"}
            })

        title = f"{parent_name} — {n_endpoints} endpoints"
        return True, title, sub_els, stylesheet

    # ---------------------- subgraph edge-info callback ------------------------
    @app.callback(
        Output({'type': 'subgraph-edge-info', 'prot': MATCH}, 'children'),
        Input({'type': 'cy-subgraph', 'prot': MATCH}, 'tapEdgeData'),
        prevent_initial_call=True)
    def show_subgraph_edge_info(edge):
        if not edge:
            return [html.Div("Select an edge for details.",
                             style={'color': U['muted'], 'fontSize': 14, 'padding': 16})]

        rel = edge.get('rel', '')
        source = edge.get('source', '')
        target = edge.get('target', '')

        card_style = {'background': U['card'], 'padding': 12, 'borderRadius': 4,
                      'boxShadow': U['shadow'],
                      'marginBottom': 12, 'border': f'1px solid {U["rule"]}',
                      'fontFamily': U['font_ui']}

        content = [
            html.Div("Edge",
                     style={'fontSize': 16, 'fontWeight': 'bold', 'color': U['ink'],
                            'borderBottom': f'1px solid {U["rule"]}', 'paddingBottom': 8,
                            'marginBottom': 12, 'padding': '12px 16px 8px'}),
            html.Div([
                html.Div(f"{source} → {target}",
                         style={'fontSize': 13, 'fontWeight': 'bold',
                                'color': U['ink'], 'marginBottom': 6}),
                *([] if rel in ('PV',) else [
                    html.Div(f"Type: {rel}",
                             style={'fontSize': 13, 'color': U['muted']})]),
                *([] if not edge.get('evidence_count') else [
                    html.Div(f"Statements: {edge['evidence_count']}",
                             style={'fontSize': 12, 'color': U['muted'], 'marginTop': 4})]),
            ], style={**card_style, 'margin': '0 12px 12px'}),
        ]

        if edge.get('note'):
            content.append(html.Div([
                html.Div("Description",
                         style={'fontSize': 13, 'fontWeight': 'bold',
                                'color': U['ink_soft'], 'marginBottom': 6}),
                html.Div(edge['note'],
                         style={'fontSize': 13, 'color': U['ink'], 'lineHeight': '1.4'})
            ], style={**card_style, 'margin': '0 12px 12px'}))

        if edge.get('clinvar_data'):
            data = edge['clinvar_data']
            content.append(html.Div([
                html.Div("ClinVar",
                         style={'fontSize': 13, 'fontWeight': 'bold',
                                'color': U['ink_soft'], 'marginBottom': 6}),
                html.Div([
                    html.Span("Pathogenicity: ", style={'fontWeight': 'bold'}),
                    html.Span(data.get('pathogenicity', 'N/A'))
                ], style={'fontSize': 12, 'marginBottom': 4}),
                html.Div([
                    html.Span("Review: ", style={'fontWeight': 'bold'}),
                    html.Span(data.get('review', ''))
                ], style={'fontSize': 12, 'marginBottom': 4}),
                html.Div([
                    html.Span("Condition: ", style={'fontWeight': 'bold'}),
                    html.Span(data.get('conditions', ''))
                ], style={'fontSize': 12, 'lineHeight': '1.4'}),
            ], style={**card_style, 'margin': '0 12px 12px'}))

        if edge.get('pmid'):
            raw = str(edge['pmid'])
            pmids = [p.strip() for p in re.split(r'[;,]', raw) if p.strip()]
            links = []
            for p in pmids:
                links.append(html.A(
                    f"PMID {p}",
                    href=f"https://pubmed.ncbi.nlm.nih.gov/{p}/",
                    target="_blank",
                    style={'display': 'block', 'color': U['link'],
                           'textDecoration': 'none', 'fontSize': 12,
                           'marginBottom': 3}))
            content.append(html.Div([
                html.Div("PubMed",
                         style={'fontSize': 13, 'fontWeight': 'bold',
                                'color': U['ink_soft'], 'marginBottom': 6}),
                *links
            ], style={**card_style, 'margin': '0 12px 12px'}))

        return content

    # ---------------------- highlight callback ---------------------------------
    @app.callback(
        [Output({'type': 'cy-net', 'prot': MATCH}, 'elements'),
         Output({'type': 'cy-net', 'prot': MATCH}, 'layout')],
        Input({'type': 'cy-net', 'prot': MATCH}, 'tapNodeData'),
        Input({'type': 'store-map-range', 'prot': MATCH}, 'data'),
        [State({'type': 'store-els', 'prot': MATCH}, 'data'),
         State({'type': 'store-root', 'prot': MATCH}, 'data')],
        prevent_initial_call=True)
    def highlight(node, map_range, elements, root_prot):
        base = copy.deepcopy(elements)
        xr = _norm_map_range(map_range)
        trig = ctx.triggered_id
        range_changed = (
            isinstance(trig, dict) and trig.get("type") == "store-map-range"
        )
        # When the variant map is zoomed, clear any node selection so the graph
        # rebuilds in the default (un-highlighted) state.
        if range_changed:
            node = None
        if xr and root_prot in TSV_FILES:
            els = build_elements(root_prot, variant_aa_range=xr)[0]
            # Always preset so node `position` from the same layered algorithm as the full graph is applied.
            layout_out = _cy_net_layout_preset()
        else:
            els = base
            layout_out = (
                _cy_net_layout_preset() if range_changed else dash.no_update
            )

        fwd, rev = _adjacency_from_elements(els)

        def _strip_faded():
            for el in els:
                c = el.get('classes') or ''
                el['classes'] = c.replace(' faded', '')

        if not node:
            _strip_faded()
            return els, layout_out

        if node['id'] == root_prot:
            _strip_faded()
            return els, layout_out

        sel = node['id']
        keep_nodes = {sel}
        keep_edges = set()

        stack = [sel]
        while stack:
            cur = stack.pop()
            for t in fwd.get(cur, ()):
                if (cur, t) not in keep_edges:
                    keep_edges.add((cur, t))
                    keep_nodes.add(t)
                    stack.append(t)

        stack = [sel]
        while stack:
            cur = stack.pop()
            for s in rev.get(cur, ()):
                if (s, cur) not in keep_edges:
                    keep_edges.add((s, cur))
                    keep_nodes.add(s)
                    stack.append(s)

        for el in els:
            d = el.get('data') or {}
            if 'source' in d:
                keep = ((d['source'], d['target']) in keep_edges
                        or d.get('rel') == 'PV')
            else:
                keep = d.get('id') in keep_nodes

            c = el.get('classes') or ''
            if keep:
                el['classes'] = c.replace(' faded', '')
            else:
                if 'faded' not in c:
                    el['classes'] = c + ' faded'

        return els, layout_out

    # ---------------------- Search -------------------------------
    @app.callback(Output('prot-directory', 'children'),
                  Input('prot-search', 'value'))
    def filter_directory(query):
        return _build_alpha_directory(PROTS, query, _protein_href, columns=3)

    @app.callback(Output('endpoint-directory', 'children'),
                  Input('endpoint-search', 'value'))
    def filter_endpoint_directory(query):
        return _build_alpha_directory(ENDPOINTS, query, _endpoint_href, columns=2)

    @app.callback(Output('url', 'href', allow_duplicate=True),
                  Input('submit-prot', 'n_clicks'),
                  State('prot-search', 'value'),
                  prevent_initial_call=True)
    def jump_to_protein(_, value):
        return _protein_href(value) if value else dash.no_update

    @app.callback(Output('url', 'href', allow_duplicate=True),
                  Input('submit-endpoint', 'n_clicks'),
                  State('endpoint-search', 'value'),
                  prevent_initial_call=True)
    def jump_to_endpoint(_, value):
        return _endpoint_href(value) if value else dash.no_update

    @app.callback(
        Output("url", "href", allow_duplicate=True),
        Input("stats-fig-bp", "clickData"),
        Input("stats-fig-disease", "clickData"),
        Input("stats-fig-genes", "clickData"),
        prevent_initial_call=True,
    )
    def stats_bar_open_network(bp_click, dis_click, genes_click):
        trig = ctx.triggered_id
        if trig == "stats-fig-bp":
            payload = bp_click
        elif trig == "stats-fig-disease":
            payload = dis_click
        elif trig == "stats-fig-genes":
            payload = genes_click
        else:
            return dash.no_update
        if not payload or not payload.get("points"):
            return dash.no_update
        href = _stats_point_href(payload["points"][0])
        return href if href else dash.no_update

    @app.callback(
        Output("stats-metric", "data"),
        Input("stats-btn-path", "n_clicks"),
        Input("stats-btn-gene", "n_clicks"),
        Input("stats-btn-variant", "n_clicks"),
        Input("stats-btn-pmid", "n_clicks"),
        prevent_initial_call=True,
    )
    def stats_set_metric(_p, _g, _v, _m):
        key = ctx.triggered_id
        return {
            "stats-btn-path": "path",
            "stats-btn-gene": "gene",
            "stats-btn-variant": "variant",
            "stats-btn-pmid": "pmid",
        }.get(key, "path")

    @app.callback(
        Output("stats-fig-bp", "figure"),
        Output("stats-fig-disease", "figure"),
        Output("stats-fig-genes", "figure"),
        Output("stats-btn-path", "color"),
        Output("stats-btn-path", "outline"),
        Output("stats-btn-gene", "color"),
        Output("stats-btn-gene", "outline"),
        Output("stats-btn-variant", "color"),
        Output("stats-btn-variant", "outline"),
        Output("stats-btn-pmid", "color"),
        Output("stats-btn-pmid", "outline"),
        Input("stats-metric", "data"),
    )
    def stats_render(metric):
        metric = metric or "path"
        fig_bp = _stats_bar_figure(
            _STATS_BP_ROWS, metric,
            "Biological processes",
            for_gene_chart=False,
            value_fn=_stats_value_bp_disease,
            bar_color=U["chart_bp"],
            plot_bg=U["plot_bp"],
            paper_bg="rgba(0,0,0,0)",
            height=STATS_FIG_HEIGHT_9,
        )
        fig_dis = _stats_bar_figure(
            _STATS_DISEASE_ROWS, metric,
            "Diseases",
            for_gene_chart=False,
            value_fn=_stats_value_bp_disease,
            bar_color=U["chart_dis"],
            plot_bg=U["plot_dis"],
            paper_bg="rgba(0,0,0,0)",
            height=STATS_FIG_HEIGHT_10,
        )
        fig_genes = _stats_bar_figure(
            _STATS_GENE_ROWS, metric,
            "Genes",
            for_gene_chart=True,
            value_fn=_stats_value_gene,
            bar_color=U["chart_gene"],
            plot_bg=U["plot_gene"],
            paper_bg="rgba(0,0,0,0)",
            height=STATS_FIG_HEIGHT_20,
        )
        btn = []
        for m in ("path", "gene", "variant", "pmid"):
            active = metric == m
            btn.append("primary" if active else "secondary")
            btn.append(not active)
        return (fig_bp, fig_dis, fig_genes, *btn)
