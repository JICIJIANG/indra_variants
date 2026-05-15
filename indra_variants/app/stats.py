"""Reference statistics: hard-coded top-rows tables, helper accessors, the
horizontal bar-chart figure builder, and the A–Z alpha directory used by the
search page."""
from typing import Optional

import plotly.graph_objects as go
from dash import html

from indra_variants.app.utils import (
    U,
    _alpha_bucket,
    _endpoint_href,
    _protein_href,
    _sort_text,
)
from indra_variants.app.data_index import ENDPOINT_INDEX, PROTS


# Reference coverage tables (publication / supplementary-style summary).
_STATS_BP_ROWS = [
    ("Cell population proliferation", 2206, 749, 1370, 1742),
    ("Apoptotic process", 1878, 667, 1252, 1495),
    ("Cell death", 1105, 413, 763, 849),
    ("Neoplasm invasiveness", 812, 339, 567, 623),
    ("Localization", 776, 418, 679, 572),
    ("Cell differentiation", 749, 350, 538, 630),
    ("Cell growth", 699, 332, 526, 532),
    ("Cell survival", 695, 309, 528, 529),
    ("Cell migration", 570, 277, 472, 408),
]
_STATS_DISEASE_ROWS = [
    ("Neoplasms", 503, 231, 359, 382),
    ("Parkinson disease", 210, 37, 79, 156),
    ("Alzheimer disease", 143, 55, 100, 116),
    ("Melanoma", 75, 24, 37, 61),
    ("Breast neoplasms", 55, 33, 48, 42),
    ("Syndrome", 46, 29, 46, 28),
    ("Infections", 38, 31, 38, 35),
    ("Amyotrophic lateral sclerosis", 34, 4, 32, 6),
    ("Disease", 33, 28, 33, 26),
    ("Frontotemporal lobar degeneration", 32, 4, 32, 5),
]
# (gene, paths, variants, bps_or_diseases, pmids)
_STATS_GENE_ROWS = [
    ("TP53", 1298, 87, 110, 261),
    ("BRAF", 1218, 30, 130, 562),
    ("KRAS", 1152, 27, 162, 365),
    ("TARDBP", 909, 35, 54, 30),
    ("LRRK2", 398, 17, 77, 166),
    ("SOD1", 305, 16, 53, 141),
    ("JAK2", 231, 16, 50, 125),
    ("HRAS", 225, 14, 56, 89),
    ("MAPT", 187, 19, 77, 55),
    ("RPS6KB1", 183, 4, 58, 100),
    ("PIK3CA", 180, 9, 47, 55),
    ("EGFR", 179, 29, 38, 78),
    ("SNCA", 178, 10, 45, 52),
    ("NRAS", 173, 14, 47, 63),
    ("CTNNB1", 156, 20, 36, 59),
    ("DNM1L", 156, 16, 41, 48),
    ("RAC1", 151, 16, 51, 54),
    ("MDM2", 136, 22, 41, 37),
    ("GSK3B", 128, 15, 39, 40),
    ("PTEN", 127, 23, 39, 40),
]


def _stats_value_bp_disease(row: tuple, metric: str) -> int:
    _name, paths, genes, variants, pmids = row
    return {"path": paths, "gene": genes, "variant": variants,
            "pmid": pmids}[metric]


def _stats_value_gene(row: tuple, metric: str) -> int:
    _gene, paths, variants, bp_dis, pmids = row
    return {"path": paths, "gene": bp_dis, "variant": variants,
            "pmid": pmids}[metric]


def _stats_axis_label(metric: str, for_gene_chart: bool) -> str:
    if metric == "path":
        return "Paths"
    if metric == "variant":
        return "Variants"
    if metric == "pmid":
        return "PMIDs"
    return "BPs / diseases" if for_gene_chart else "Genes"


def _stats_network_href(label: str, for_gene_chart: bool) -> str:
    """Resolve bar label to an in-app network URL, or '' if not indexed."""
    key = (label or "").strip()
    if not key:
        return ""
    if for_gene_chart:
        for p in PROTS:
            if p.casefold() == key.casefold():
                return _protein_href(p)
        return ""
    for ep in ENDPOINT_INDEX:
        if ep.casefold() == key.casefold():
            return _endpoint_href(ep)
    return ""


def _stats_point_href(point: dict) -> Optional[str]:
    cd = point.get("customdata")
    if isinstance(cd, (list, tuple)) and len(cd) >= 1:
        href = cd[0]
        if isinstance(href, str) and href.startswith("/"):
            return href
    if isinstance(cd, str) and cd.startswith("/"):
        return cd
    return None


# Per-bar pixel height – kept constant so all three charts have the same bar size.
STATS_BAR_PX = 34
STATS_FIG_MARGIN_PX = 118   # title + x-axis + padding
# Convenience heights (used by both the figure and the card container)
STATS_FIG_HEIGHT_9 = STATS_BAR_PX * 9 + STATS_FIG_MARGIN_PX     # BP (9 rows)
STATS_FIG_HEIGHT_10 = STATS_BAR_PX * 10 + STATS_FIG_MARGIN_PX   # disease
STATS_FIG_HEIGHT_20 = STATS_BAR_PX * 20 + STATS_FIG_MARGIN_PX   # gene
STATS_FIG_HEIGHT = STATS_FIG_HEIGHT_10   # backward-compat alias


def _stats_bar_figure(
    rows: list[tuple],
    metric: str,
    title: str,
    for_gene_chart: bool,
    value_fn,
    *,
    bar_color: str,
    plot_bg: str = "#f8fafc",
    paper_bg: str = "rgba(0,0,0,0)",
    height: int = STATS_FIG_HEIGHT_10,
) -> go.Figure:
    scored = [(value_fn(r, metric), r) for r in rows]
    scored.sort(key=lambda t: t[0], reverse=True)
    values = [v for v, _ in scored]
    labels = [r[0] for _, r in scored]
    axis_title = _stats_axis_label(metric, for_gene_chart)
    hrefs = [_stats_network_href(lab, for_gene_chart) for lab in labels]
    hover_hint = [
        "Click to open graph" if h else "No matching graph in this build"
        for h in hrefs
    ]
    bar_height = 22
    fig = go.Figure(
        data=[go.Bar(
            x=values,
            y=labels,
            orientation="h",
            marker=dict(
                color=bar_color,
                line=dict(color="rgba(255,255,255,0.45)", width=1),
            ),
            text=values,
            textposition="outside",
            textfont=dict(color=U["ink_soft"], size=11),
            cliponaxis=False,
            # No customdata on bars — navigation is via y-axis label annotations.
            hovertemplate="<b>%{y}</b><br>%{x} " + axis_title + "<extra></extra>",
        )]
    )
    fig.update_layout(
        title=dict(text=title, font=dict(size=15, color=U["ink"],
                                          family=U["font_display"])),
        xaxis_title=axis_title,
        yaxis=dict(autorange="reversed", title="", showticklabels=False),
        margin=dict(l=10, r=88, t=52, b=44),
        height=height,
        font=dict(family=U["font_ui"], size=12, color=U["ink_soft"]),
        paper_bgcolor=paper_bg,
        plot_bgcolor=plot_bg,
        hoverlabel=dict(bgcolor=U["panel"], font_size=12,
                        font_family=U["font_ui"]),
    )
    fig.update_xaxes(
        gridcolor="rgba(45, 42, 36, 0.08)",
        zeroline=False,
        title_font=dict(color=U["muted"], size=12, family=U["font_ui"]),
    )
    # Replace y-axis tick labels with annotations that carry <a href> links so
    # clicking the label name navigates to the network, not clicking the bar.
    # Estimate left margin: ~6.5px per character at size-11 font.
    max_label_len = max((len(lb) for lb in labels), default=1)
    left_margin = max(100, min(int(max_label_len * 6.8) + 12, 280))
    fig.update_layout(margin=dict(l=left_margin, r=88, t=52, b=44))
    for lb, href in zip(labels, hrefs):
        link_text = (
            f'<a href="{href}">{lb}</a>'
            if href else lb
        )
        fig.add_annotation(
            x=0, y=lb,
            xref="paper", yref="y",
            text=link_text,
            showarrow=False,
            xanchor="right",
            font=dict(size=11, color=U["link"] if href else U["ink_soft"],
                      family=U["font_ui"]),
            align="right",
        )
    return fig


def _build_alpha_directory(items, query: str, href_builder, columns: int = 3):
    query_norm = (query or "").strip().casefold()
    blocks = []
    bucket_order = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ") + ["#"]

    for letter in bucket_order:
        group = [item for item in items if _alpha_bucket(item) == letter]
        if query_norm:
            group = [item for item in group if query_norm in item.casefold()]
        if not group:
            continue

        summary_label = "Other" if letter == "#" else letter
        blocks.append(
            html.Details([
                html.Summary(f"{summary_label} ({len(group)})",
                             style={'cursor': 'pointer', 'fontSize': 20}),
                html.Ul([
                    html.Li(html.A(
                        item,
                        href=href_builder(item),
                        style={'textDecoration': 'none',
                               'color': U['link'],
                               'fontWeight': 'bold'
                               if query_norm and query_norm in item.casefold()
                               else 'normal'}))
                    for item in sorted(group, key=_sort_text)
                ], style={'columnCount': columns, 'listStyle': 'none',
                          'padding': 0, 'margin': '6px 0'})
            ], open=bool(query_norm))
        )

    if blocks:
        return blocks
    return html.Div("No results.",
                    style={'color': U['muted'], 'fontSize': 15,
                           'padding': '12px 0', 'fontFamily': U['font_ui']})
