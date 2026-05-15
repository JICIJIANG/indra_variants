"""Page renderers: search/browse, statistics landing, and the network view
shared by protein-centric and phenotype-centric routes."""
from typing import Optional

import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import plotly.graph_objects as go
from dash import dcc, html

from indra_variants.app.data_index import (
    ENDPOINTS,
    ENDPOINT_OPTIONS,
    PROTS,
    PROT_OPTIONS,
)
from indra_variants.app.figures import _protein_lollipop_figure
from indra_variants.app.graph_builders import build_elements, build_endpoint_elements
from indra_variants.app.stats import (
    STATS_FIG_HEIGHT_9,
    STATS_FIG_HEIGHT_10,
    STATS_FIG_HEIGHT_20,
    _STATS_BP_ROWS,
    _STATS_DISEASE_ROWS,
    _STATS_GENE_ROWS,
    _stats_network_href,
)
from indra_variants.app.utils import (
    GRAPH_ENDPOINT_BG,
    GRAPH_ENDPOINT_FG,
    GRAPH_INTERMEDIATE_BG,
    GRAPH_INTERMEDIATE_FG,
    GRAPH_PROTEIN_BG,
    GRAPH_PROTEIN_FG,
    GRAPH_VARIANT_BG,
    GRAPH_VARIANT_FG,
    U,
    _CLINVAR_DOT_COLORS,
    _GENE_STRIP_GRAPH_INNER_H,
    _GENE_STRIP_VH,
    _NETWORK_HEADER_PX,
    _network_frame_style,
    _variant_strip_wrap_style,
)


def _browse_panel(title: str, helper_text: str, dropdown_id: str, options: list,
                  placeholder: str, button_id: str, directory_id: str,
                  summary_text: str, note_text: Optional[str] = None):
    return html.Div([
        html.H4(title, style={'marginTop': 0, 'marginBottom': 10,
                              'color': U['ink'],
                              'fontFamily': U['font_display'],
                              'fontWeight': 600}),
        html.P(helper_text, style={'fontSize': 16, 'margin': '0 0 16px 0',
                                   'color': U['ink_soft'],
                                   'fontFamily': U['font_ui'],
                                   'lineHeight': 1.45}),
        dcc.Dropdown(id=dropdown_id, options=options,
                     placeholder=placeholder,
                     style={'fontSize': 16, 'fontFamily': U['font_ui']},
                     clearable=True,
                     searchable=True),
        dbc.Button("Search", id=button_id, n_clicks=0,
                   color="primary", style={'marginTop': 18,
                                           'fontSize': 16,
                                           'fontFamily': U['font_ui'],
                                           'fontWeight': 600,
                                           'padding': '8px 22px'}),
        html.Div(summary_text,
                 style={'marginTop': 16, 'marginBottom': 10,
                        'fontSize': 14, 'color': U['muted'],
                        'fontFamily': U['font_ui']}),
        *([] if not note_text else [
            html.Div(note_text,
                     style={'marginBottom': 14, 'fontSize': 14,
                            'color': U['muted'], 'fontStyle': 'italic',
                            'fontFamily': U['font_ui'], 'lineHeight': 1.45})
        ]),
        html.Div(id=directory_id,
                 style={'fontFamily': U['font_ui']})
    ], style={'padding': '12px 6px 6px'})


def _render_network_page(view_key: str, root_node_id: str, title: str,
                         graph_tuple: tuple, layout: Optional[dict] = None,
                         lollipop_figure: Optional[go.Figure] = None):
    els, legend_rels, legend_colors, rel_color_safe, edge_set, subgraph_data = graph_tuple
    layout = layout or {'name': 'preset'}

    def rel_style(css_cls, c):
        return {'selector': f'.edge-{css_cls}',
                'style': {
                    'line-color': c,
                    'target-arrow-color': c,
                    'target-arrow-shape': 'triangle',
                    'curve-style': 'bezier',
                    # Thicker edges for higher evidence counts.
                    'width': 'mapData(evidence_count, 1, 30, 1, 6.5)',
                }}

    sidebar = html.Div(
        id={'type': 'edge-info', 'prot': view_key},
        children=[
            html.Div("Details",
                    style={'fontSize': 17, 'fontWeight': 600,
                           'fontFamily': U['font_display'],
                           'marginBottom': 15, 'color': U['ink'],
                           'borderBottom': f'1px solid {U["rule"]}',
                           'paddingBottom': 10}),
            html.Div("Select a node or edge for attributes, evidence, and external links.",
                    style={'color': U['muted'], 'fontSize': 14,
                           'lineHeight': '1.45', 'fontFamily': U['font_ui']})
        ],
        style={
            'position': 'fixed',
            'left': 0,
            'top': 0,
            'width': 350,
            'height': '100vh',
            'background': U['graph_bg'],
            'padding': '22px 20px',
            'boxShadow': U['shadow_strong'],
            'borderRight': f'1px solid {U["rule"]}',
            'fontSize': 15,
            'fontFamily': U['font_ui'],
            'zIndex': 1000,
            'overflowY': 'auto'
        }
    )

    main_content = html.Div([
        html.Div([
            dcc.Link("← Browser", href="/",
                    style={'color': U['link'], 'textDecoration': 'none',
                           'fontSize': 14, 'fontWeight': 600,
                           'fontFamily': U['font_ui']}),
            html.H4(title,
                   style={'textAlign': 'center', 'margin': '6px 0 2px',
                          'color': U['ink'],
                          'fontFamily': U['font_display'],
                          'fontWeight': 600}),
            html.P("Select the root node to clear highlighting.",
                   style={'textAlign': 'center', 'marginTop': 0,
                          'marginBottom': 15, 'color': U['muted'],
                          'fontFamily': U['font_ui'], 'fontSize': 13})
        ], style={'padding': '14px 16px', 'background': U['panel'],
                  'borderBottom': f'1px solid {U["rule"]}'}),

        dcc.Store(id={'type': 'store-els',  'prot': view_key},  data=els),
        dcc.Store(id={'type': 'store-edges', 'prot': view_key},  data=edge_set),
        dcc.Store(id={'type': 'store-root', 'prot': view_key},  data=root_node_id),
        dcc.Store(id={'type': 'store-subgraphs', 'prot': view_key}, data=subgraph_data),
        dcc.Store(id={'type': 'store-relcolors', 'prot': view_key}, data=rel_color_safe),
        dcc.Store(id={'type': 'store-map-range', 'prot': view_key}, data=None),
        *([] if lollipop_figure is not None else [
            html.Div([
                dcc.Graph(
                    id={'type': 'gene-map', 'prot': view_key},
                    figure={'data': [], 'layout': {'margin': {'t': 0, 'b': 0, 'l': 0, 'r': 0}}},
                    config={'displayModeBar': False, 'staticPlot': True},
                    style={'display': 'none', 'height': 0, 'width': 0},
                ),
                html.Button(
                    id={'type': 'gene-map-reset', 'prot': view_key},
                    n_clicks=0,
                    style={'display': 'none'},
                ),
            ], style={'display': 'none'}),
        ]),

        dbc.Modal([
            dbc.ModalHeader(dbc.ModalTitle(
                id={'type': 'subgraph-title', 'prot': view_key})),
            dbc.ModalBody(
                html.Div([
                    cyto.Cytoscape(
                        id={'type': 'cy-subgraph', 'prot': view_key},
                        elements=[], layout={'name': 'preset'},
                        style={'width': '70%', 'height': '100%',
                               'backgroundColor': U['paper']},
                        stylesheet=[]),
                    html.Div(
                        id={'type': 'subgraph-edge-info', 'prot': view_key},
                        children=[
                            html.Div("Select an edge for details.",
                                     style={'color': U['muted'], 'fontSize': 14,
                                            'padding': 16,
                                            'fontFamily': U['font_ui']})
                        ],
                        style={'width': '30%', 'height': '100%',
                               'overflowY': 'auto',
                               'borderLeft': f'1px solid {U["rule"]}',
                               'background': U['paper']})
                ], style={'display': 'flex', 'height': '70vh'}),
                style={'padding': 0}),
        ], id={'type': 'subgraph-modal', 'prot': view_key},
           size="xl", is_open=False),

        html.Div(
            [
                html.Div(
                    [
                        cyto.Cytoscape(
                            id={'type': 'cy-net', 'prot': view_key},
                            elements=els,
                            layout=layout,
                            style={
                                'width': '100%',
                                'height': '100%',
                                'flex': '1 1 auto',
                                'minHeight': 0,
                                'backgroundColor': U['card'],
                            },
                            stylesheet=[
                {'selector': 'node', 'style': {
                    'shape': 'ellipse', 'background-opacity': 0.92,
                    'font-size': 15, 'font-weight': '600',
                    'label': 'data(label)',
                    'text-wrap': 'wrap',
                    'text-max-width': 100,
                    'text-valign': 'center',
                    'text-halign': 'center'}},
                {'selector': '.role-protein',
                 'style': {'background-color': GRAPH_PROTEIN_BG,
                           'color': GRAPH_PROTEIN_FG}},
                {'selector': '.role-variant',
                 'style': {'background-color': GRAPH_VARIANT_BG,
                           'color': GRAPH_VARIANT_FG}},
                {'selector': '.role-intermediate',
                 'style': {'background-color': GRAPH_INTERMEDIATE_BG,
                           'color': GRAPH_INTERMEDIATE_FG}},
                {'selector': '.role-endpoint',
                 'style': {'background-color': GRAPH_ENDPOINT_BG,
                           'color': GRAPH_ENDPOINT_FG}},
                {'selector': '.edge-PV',
                 'style': {
                     'line-color': '#c9c4bf',
                     'target-arrow-shape': 'none',
                     'width': 1.5,
                     'opacity': 0.6,
                 }},
                *[rel_style(css_cls, c) for css_cls, c in rel_color_safe.items()],
                {'selector': '.faded', 'style': {'opacity': 0.15}}
            ],
                        ),
                    ],
                    style={
                        **_network_frame_style(),
                        'display': 'flex',
                        'flexDirection': 'column',
                        'flex': '1 1 auto',
                        'minHeight': 0,
                    },
                ),
                *([] if lollipop_figure is None else [
                    html.Div(
                        [
                            html.Div([
                                # ── Title row ───────────────────────────────
                                html.Div([
                                    html.Span(
                                        "Variant map",
                                        style={
                                            'fontSize': 12, 'color': U['ink_soft'],
                                            'fontFamily': U['font_ui'],
                                            'fontWeight': 600,
                                            'letterSpacing': '0.06em',
                                            'textTransform': 'uppercase',
                                            'flex': '1',
                                        }),
                                    dbc.Button(
                                        "Reset view",
                                        id={'type': 'gene-map-reset', 'prot': view_key},
                                        n_clicks=0, size="sm", outline=True,
                                        style={
                                            'fontSize': 11, 'padding': '2px 10px',
                                            'fontFamily': U['font_ui'],
                                        }),
                                ], style={
                                    'display': 'flex', 'alignItems': 'center',
                                    'justifyContent': 'space-between',
                                    'padding': '6px 10px 2px',
                                }),
                                # ── ClinVar pathogenicity filter row ─────────
                                html.Div([
                                    html.Span(
                                        "Filter",
                                        style={
                                            'fontSize': 10, 'color': U['muted'],
                                            'fontFamily': U['font_ui'],
                                            'fontWeight': 600,
                                            'letterSpacing': '0.05em',
                                            'textTransform': 'uppercase',
                                            'marginRight': 8,
                                        }),
                                    dcc.Checklist(
                                        id={'type': 'gene-map-clin-filter', 'prot': view_key},
                                        options=[
                                            {
                                                'label': html.Span([
                                                    html.Span(
                                                        '●',
                                                        style={'color': _CLINVAR_DOT_COLORS['pathogenic'],
                                                               'marginRight': 3, 'fontSize': 11}),
                                                    'Pathogenic / likely pathogenic',
                                                ], style={'fontSize': 11, 'fontFamily': U['font_ui']}),
                                                'value': 'pathogenic',
                                            },
                                            {
                                                'label': html.Span([
                                                    html.Span(
                                                        '●',
                                                        style={'color': _CLINVAR_DOT_COLORS['vus'],
                                                               'marginRight': 3, 'fontSize': 11}),
                                                    'VUS / conflicting',
                                                ], style={'fontSize': 11, 'fontFamily': U['font_ui']}),
                                                'value': 'vus',
                                            },
                                            {
                                                'label': html.Span([
                                                    html.Span(
                                                        '●',
                                                        style={'color': _CLINVAR_DOT_COLORS['benign'],
                                                               'marginRight': 3, 'fontSize': 11}),
                                                    'Benign / likely benign',
                                                ], style={'fontSize': 11, 'fontFamily': U['font_ui']}),
                                                'value': 'benign',
                                            },
                                        ],
                                        value=['pathogenic', 'vus', 'benign'],
                                        inline=True,
                                        labelStyle={
                                            'display': 'inline-flex',
                                            'alignItems': 'center',
                                            'cursor': 'pointer',
                                            'marginRight': 10,
                                        },
                                        inputStyle={
                                            'marginRight': 4,
                                            'cursor': 'pointer',
                                            # White box + app-background khaki when checked
                                            'accentColor': '#e6e1d8',
                                            'backgroundColor': '#ffffff',
                                            'border': f"1px solid {U['wash']}",
                                            'borderRadius': 3,
                                            'width': 15,
                                            'height': 15,
                                            'verticalAlign': 'middle',
                                            'boxSizing': 'border-box',
                                        },
                                        style={'display': 'inline-flex', 'alignItems': 'center',
                                               'gap': 0},
                                    ),
                                ], style={
                                    'display': 'flex', 'alignItems': 'center',
                                    'padding': '3px 10px 5px',
                                    'borderTop': f'1px solid {U["rule"]}',
                                }),
                            ], style={
                                'background': U['card'],
                            }),
                            dcc.Graph(
                                id={'type': 'gene-map', 'prot': view_key},
                                figure=lollipop_figure,
                                config={
                                    'displayModeBar': False,
                                    'scrollZoom': True,
                                    'doubleClick': False,
                                },
                                style={
                                    'height': _GENE_STRIP_GRAPH_INNER_H,
                                    'width': '100%', 'margin': 0,
                                    'paddingLeft': 6, 'paddingRight': 6,
                                    'boxSizing': 'border-box',
                                },
                            ),
                        ],
                        style={
                            **_variant_strip_wrap_style(),
                            'flex': f'0 0 {_GENE_STRIP_VH}',
                            'minHeight': 0,
                            'display': 'flex',
                            'flexDirection': 'column',
                        },
                    ),
                ]),
            ],
            style={
                'display': 'flex', 'flexDirection': 'column',
                'flex': '1 1 auto', 'minHeight': 0,
                'height': f'calc(100vh - {_NETWORK_HEADER_PX}px)',
                'gap': 10,
                'padding': '12px 14px 14px',
                'boxSizing': 'border-box',
            },
        ),

        html.Div([
            html.H4("Edge types",
                    style={'margin': 0, 'fontSize': 14,
                           'fontWeight': 600,
                           'fontFamily': U['font_display'],
                           'color': U['ink'],
                           'letterSpacing': '0.02em'}),
            html.Ul([
                html.Li([html.Span('→',
                                  style={'color': legend_colors.get(r, '#c9c4bf'),
                                         'marginRight': 8,
                                         'fontSize': 14}), r],
                        style={'fontSize': 13, 'listStyle': 'none',
                               'margin': '6px 0',
                               'color': U['ink_soft'],
                               'fontFamily': U['font_ui']})
                for r in legend_rels
            ], style={'paddingLeft': 0, 'margin': '10px 0 0 0'})
        ], style={'position': 'absolute',
                  'top': _NETWORK_HEADER_PX + 47, 'right': 22,
                  'background': U['card'],
                  'padding': '14px 18px',
                  'borderRadius': 4,
                  'border': f'1px solid {U["rule"]}',
                  'boxShadow': U['shadow'],
                  'fontFamily': U['font_ui'],
                  'maxHeight': '60vh',
                  'overflowY': 'auto',
                  'zIndex': 10})

    ], style={
        'marginLeft': 350,
        'position': 'relative',
        'height': '100vh',
        'backgroundColor': U['wash'],
        'display': 'flex',
        'flexDirection': 'column',
    })

    return html.Div([sidebar, main_content])


def _app_footer():
    return html.Div(
        [
            html.Span("Maintained by the "),
            html.A("Gyori Lab", href="https://gyorilab.github.io",
                   target="_blank",
                   style={'color': U['link'], 'fontWeight': 500}),
            html.Span(", Northeastern University."),
            html.Br(),
            html.Span("Supported by DARPA ASKEM / ARPA-H BDF (HR00112220036).")
        ],
        style={'background': U['wash'],
               'padding': '16px 24px',
               'textAlign': 'center', 'fontSize': 13,
               'fontFamily': U['font_ui'],
               'color': U['muted'],
               'borderTop': f'1px solid {U["rule"]}'},
    )


# ---Search / browse (gene & endpoint pickers)---
def _quick_stats_panel():
    """Compact right-sidebar panel: top-5 genes, diseases, and pathways as links."""
    # Reuse pre-computed TOP_BP_ROWS / TOP_DISEASE_ROWS / TOP_GENE_ROWS
    def _top5_links(rows, for_gene):
        items = []
        for row in rows[:5]:
            label = row[0]
            href = _stats_network_href(label, for_gene)
            if href:
                items.append(html.Li(
                    dcc.Link(label, href=href,
                             style={'color': U['link'], 'textDecoration': 'none',
                                    'fontSize': 13, 'fontFamily': U['font_ui']}),
                    style={'marginBottom': 5, 'lineHeight': 1.4}))
            else:
                items.append(html.Li(label,
                    style={'color': U['ink_soft'], 'fontSize': 13,
                           'marginBottom': 5, 'fontFamily': U['font_ui']}))
        return html.Ul(items, style={'paddingLeft': 16, 'margin': '6px 0 0'})

    def _section(title, rows, for_gene, accent):
        return html.Div([
            html.Div(title, style={
                'fontSize': 11, 'fontWeight': 700, 'letterSpacing': '0.10em',
                'textTransform': 'uppercase', 'color': U['ink_soft'],
                'fontFamily': U['font_ui'], 'marginBottom': 2,
                'borderLeft': f'3px solid {accent}', 'paddingLeft': 8,
            }),
            _top5_links(rows, for_gene),
        ], style={'marginBottom': 20})

    return html.Div([
        html.Div("Quick Stats", style={
            'fontSize': 13, 'fontWeight': 700, 'color': U['ink'],
            'fontFamily': U['font_display'], 'marginBottom': 16,
            'paddingBottom': 8, 'borderBottom': f'1px solid {U["rule"]}',
        }),
        _section("Genes",    _STATS_GENE_ROWS,    True,  U['accent_card_gene']),
        _section("Diseases", _STATS_DISEASE_ROWS, False, U['accent_card_dis']),
        _section("Biological Processes", _STATS_BP_ROWS,      False, U['accent_card_bp']),
        dcc.Link("Full statistics →", href="/statistics",
                 style={'fontSize': 12, 'color': U['link'],
                        'fontFamily': U['font_ui'], 'fontWeight': 600}),
    ], style={
        'background': U['panel'], 'padding': '22px 20px',
        'borderRadius': 4, 'border': f'1px solid {U["rule"]}',
        'boxShadow': U['shadow'], 'position': 'sticky', 'top': 24,
        'minWidth': 200, 'maxWidth': 240,
    })


def search_page():
    search_card = html.Div(
        [
            html.Div([
                html.Span("INDRA Variant", style={
                    'fontSize': 11, 'letterSpacing': '0.14em',
                    'textTransform': 'uppercase', 'color': U['hero_muted'],
                    'fontWeight': 700, 'fontFamily': U['font_ui']}),
            ], style={'marginBottom': 6}),
            html.H1("Variant Network",
                    style={'marginTop': 0, 'marginBottom': 10,
                           'color': U['ink'],
                           'fontFamily': U['font_display'],
                           'fontWeight': 600,
                           'fontSize': '1.9rem',
                           'letterSpacing': '-0.01em'}),
            html.P(
                "Protein-centric graphs are keyed by gene symbol; "
                "phenotype-centric graphs by disease or biological process term.",
                style={'fontSize': 14, 'margin': '0 0 22px 0',
                       'color': U['ink_soft'], 'fontFamily': U['font_ui'],
                       'lineHeight': 1.5}),
            dcc.Tabs([
                dcc.Tab(
                    label="Protein-centric",
                    children=[
                        _browse_panel(
                            title="Protein",
                            helper_text="Search or browse A–Z",
                            dropdown_id='prot-search',
                            options=PROT_OPTIONS,
                            placeholder="Protein symbol…",
                            button_id='submit-prot',
                            directory_id='prot-directory',
                            summary_text=f"{len(PROTS)} graphs available.",
                        )
                    ]
                ),
                dcc.Tab(
                    label="Phenotype-centric",
                    children=[
                        _browse_panel(
                            title="Phenotype",
                            helper_text="Search or browse A–Z",
                            dropdown_id='endpoint-search',
                            options=ENDPOINT_OPTIONS,
                            placeholder="Biological process or disease term…",
                            button_id='submit-endpoint',
                            directory_id='endpoint-directory',
                            summary_text=f"{len(ENDPOINTS)} indexed phenotype nodes.",
                        )
                    ]
                )
            ])
        ],
        style={'flex': '1 1 0', 'minWidth': 0,
               'background': U['panel'], 'padding': '36px 40px',
               'borderRadius': 4, 'boxShadow': U['shadow'],
               'fontFamily': U['font_ui'],
               'border': f'1px solid {U["rule"]}'}
    )

    body = html.Div([
        html.Div([search_card, _quick_stats_panel()],
                 style={'display': 'flex', 'gap': 24, 'alignItems': 'flex-start',
                        'maxWidth': 1100, 'margin': '48px auto', 'padding': '0 22px'}),
        _app_footer(),
    ], style={'minHeight': '100vh', 'background': U['paper']})

    return body


def statistics_page():
    """Landing page: static reference tables for top BPs / diseases / genes."""
    _cta_wrap = {
        'textDecoration': 'none',
        'display': 'block',
        'width': '100%',
    }
    hero = html.Div(
        [
            html.Div(
                style={'height': 3, 'background': '#b8a06e', 'opacity': 0.85}),
            html.Div(
                [
                    html.Div(
                        [
                            html.Span(
                                "INDRA variant networks",
                                style={'fontSize': 11, 'letterSpacing': '0.14em',
                                       'textTransform': 'uppercase',
                                       'color': U['hero_muted'],
                                       'fontWeight': 600,
                                       'fontFamily': U['font_ui']}),
                            html.H1(
                                "Variant networks overview",
                                style={'fontSize': '2.05rem', 'fontWeight': 600,
                                       'fontFamily': U['font_display'],
                                       'color': U['hero_text'], 'margin': '12px 0 10px',
                                       'lineHeight': 1.2,
                                       'letterSpacing': '-0.02em'}),
                            html.P(
                                "Ranked biological processes, diseases, and genes "
                                "from the reference path statistics (static summary).",
                                style={'fontSize': 15, 'color': U['hero_muted'],
                                       'maxWidth': 620, 'margin': 0, 'lineHeight': 1.55,
                                       'fontFamily': U['font_ui']}),
                        ],
                        style={'flex': '1 1 320px', 'minWidth': 0}),
                    html.Div(
                        [
                            html.Div("Graphs", style={
                                'fontSize': 11, 'textTransform': 'uppercase',
                                'letterSpacing': '0.12em', 'color': U['hero_muted'],
                                'marginBottom': 10, 'fontWeight': 600,
                                'fontFamily': U['font_ui']}),
                            dcc.Link(
                                dbc.Button(
                                    "Open browser",
                                    className="w-100 mb-2",
                                    style={
                                        'fontWeight': 600,
                                        'fontFamily': U['font_ui'],
                                        'background': '#b8a06e',
                                        'color': '#1c1b18',
                                        'border': 'none',
                                        'letterSpacing': '0.02em',
                                        'boxShadow': '0 2px 8px rgba(184,160,110,0.35)',
                                    }),
                                href="/",
                                style=_cta_wrap),
                            dcc.Link(
                                dbc.Button(
                                    "Protein-centric",
                                    className="w-100 mb-2",
                                    style={
                                        'fontWeight': 500,
                                        'fontFamily': U['font_ui'],
                                        'background': '#4d7c8a',
                                        'color': '#f7f4ec',
                                        'border': 'none',
                                        'letterSpacing': '0.02em',
                                    }),
                                href="/",
                                style=_cta_wrap),
                            dcc.Link(
                                dbc.Button(
                                    "Phenotype-centric",
                                    className="w-100 mb-2",
                                    style={
                                        'fontWeight': 500,
                                        'fontFamily': U['font_ui'],
                                        'background': '#a1665f',
                                        'color': '#f7f4ec',
                                        'border': 'none',
                                        'letterSpacing': '0.02em',
                                    }),
                                href="/",
                                style=_cta_wrap),
                            html.P(
                                "Opens /search — use tabs to switch view.",
                                style={'fontSize': 12, 'color': U['hero_muted'],
                                       'marginTop': 12, 'marginBottom': 0,
                                       'lineHeight': 1.45,
                                       'fontFamily': U['font_ui']}),
                        ],
                        style={'flex': '0 0 260px', 'maxWidth': '100%'}),
                ],
                style={'display': 'flex', 'flexWrap': 'wrap',
                       'gap': '32px 40px', 'alignItems': 'flex-start',
                       'justifyContent': 'space-between',
                       'padding': '32px 40px 36px'},
            ),
        ],
        style={
            'background': f'linear-gradient(165deg, {U["hero_hi"]} 0%, {U["hero"]} 55%, {U["hero_deep"]} 100%)',
            'borderBottom': f'1px solid {U["rule"]}',
            'boxShadow': U['shadow_strong'],
        },
    )

    metric_row = html.Div(
        [
            html.Span(
                "Metric",
                style={'fontWeight': 600, 'marginRight': 14, 'alignSelf': 'center',
                       'color': U['ink_soft'], 'fontSize': 13,
                       'fontFamily': U['font_ui'],
                       'letterSpacing': '0.04em', 'textTransform': 'uppercase'}),
            dbc.ButtonGroup(
                [
                    dbc.Button("Paths", id="stats-btn-path", n_clicks=0,
                               color="primary", outline=False,
                               style={'fontFamily': U['font_ui']}),
                    dbc.Button("Genes", id="stats-btn-gene", n_clicks=0,
                               color="secondary", outline=True,
                               style={'fontFamily': U['font_ui']}),
                    dbc.Button("Variants", id="stats-btn-variant", n_clicks=0,
                               color="secondary", outline=True,
                               style={'fontFamily': U['font_ui']}),
                    dbc.Button("PMIDs", id="stats-btn-pmid", n_clicks=0,
                               color="secondary", outline=True,
                               style={'fontFamily': U['font_ui']}),
                ],
                size="md",
            ),
            html.Span(
                "Same selection for all charts.",
                style={'marginLeft': 14, 'color': U['muted'], 'fontSize': 13,
                       'alignSelf': 'center', 'fontFamily': U['font_ui']}),
        ],
        style={'display': 'flex', 'flexWrap': 'wrap', 'alignItems': 'center',
               'marginBottom': 8, 'padding': '14px 18px',
               'background': U['panel'],
               'borderRadius': 4,
               'border': f'1px solid {U["rule"]}',
               'boxShadow': U['shadow'],
               'width': '100%', 'boxSizing': 'border-box'},
    )

    explain = html.P(
        "Bars sort by the active metric (descending). "
        "Genes: gene count on process/disease charts; BP/disease term count on the gene chart. "
        "Click a label name on the y-axis to open its network.",
        style={'fontSize': 14, 'color': U['muted'], 'maxWidth': 920,
               'margin': '12px 0 0', 'lineHeight': 1.55,
               'fontFamily': U['font_ui']},
    )

    _card = lambda gid, accent, body_pad, fig_h: dbc.Card(
        dbc.CardBody(
            dcc.Graph(
                id=gid,
                config={'displayModeBar': False},
                style={'height': fig_h, 'width': '100%'},
            ),
            style={'paddingTop': body_pad, 'paddingBottom': 12},
        ),
        className="mb-4",
        style={
            'width': '100%',
            'boxSizing': 'border-box',
            'border': f'1px solid {U["rule"]}',
            'borderRadius': 4,
            'overflow': 'hidden',
            'boxShadow': U['shadow'],
            'borderTop': f'3px solid {accent}',
            'background': U['card'],
        },
    )

    charts = html.Div(
        [
            _card("stats-fig-bp", U['accent_card_bp'], 4, STATS_FIG_HEIGHT_9),
            _card("stats-fig-disease", U['accent_card_dis'], 4, STATS_FIG_HEIGHT_10),
            _card("stats-fig-genes", U['accent_card_gene'], 4, STATS_FIG_HEIGHT_20),
        ],
        style={
            'marginTop': 22,
            'display': 'flex',
            'flexDirection': 'column',
            'gap': 16,
            'alignItems': 'stretch',
            'width': '100%',
            'maxWidth': 920,
        },
    )

    bottom_cta = html.Div(
        [
            html.Hr(style={'border': 'none', 'borderTop': f'1px solid {U["rule"]}',
                           'margin': '32px 0 22px'}),
            html.Div(
                [
                    html.Span(
                        "Open a gene or endpoint graph",
                        style={'fontSize': 15, 'color': U['ink_soft'],
                               'fontWeight': 600, 'marginRight': 16,
                               'fontFamily': U['font_display']}),
                    dcc.Link(
                        dbc.Button(
                            "Browser →",
                            color="primary",
                            style={'fontWeight': 600,
                                   'fontFamily': U['font_ui'],
                                   'padding': '8px 20px'}),
                        href="/",
                        style={'textDecoration': 'none'}),
                ],
                style={'display': 'flex', 'flexWrap': 'wrap', 'alignItems': 'center',
                       'gap': 12, 'justifyContent': 'center'},
            ),
        ],
    )

    body = html.Div(
        [
            hero,
            html.Div(
                [metric_row, explain, charts, bottom_cta, dcc.Store(id="stats-metric", data="path")],
                style={'maxWidth': 920, 'margin': '0 auto', 'padding': '28px 22px 48px',
                       'fontFamily': U['font_ui']},
            ),
            _app_footer(),
        ],
        style={'minHeight': '100vh', 'background': U['paper']},
    )
    return body


# ---Network Page---
def network_page(prot: str):
    return _render_network_page(
        view_key=f"protein::{prot}",
        root_node_id=prot,
        title=f"{prot} — variant network",
        graph_tuple=build_elements(prot),
        lollipop_figure=_protein_lollipop_figure(prot),
    )


def endpoint_network_page(endpoint: str):
    return _render_network_page(
        view_key=f"endpoint::{endpoint}",
        root_node_id=endpoint,
        title=f"{endpoint} — phenotypecentric network",
        graph_tuple=build_endpoint_elements(endpoint),
        layout={'name': 'preset'},
    )
