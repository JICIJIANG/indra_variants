"""Protein lollipop figure and small relayout helpers."""
from collections import defaultdict
from typing import Optional

import pandas as pd
import plotly.graph_objects as go

from indra_variants.app.data_index import TSV_FILES
from indra_variants.app.utils import (
    GRAPH_PROTEIN_BG,
    GRAPH_PROTEIN_FG,
    GRAPH_VARIANT_BG,
    GRAPH_VARIANT_FG,
    U,
    _CLINVAR_DOT_COLORS,
    _clinsig_category,
    _graph_node_fill_hex,
    _sort_text,
    _variant_aa_position,
)


def _protein_lollipop_figure(
    prot: str,
    clin_filter: Optional[list] = None,
) -> Optional[go.Figure]:
    """Lollipop plot: wide protein bar as backbone, variant discs stacked vertically
    above with stems.  Variants at the same amino-acid position are stacked along
    the y-axis.  A small ClinVar pathogenicity dot is drawn below each disc.

    Args:
        prot:        Gene/protein symbol.
        clin_filter: List of categories to show ('pathogenic', 'vus', 'benign').
                     None = show all.  Unknown/unannotated variants always shown.
    """
    if prot not in TSV_FILES:
        return None
    df = pd.read_csv(TSV_FILES[prot], sep="\t").fillna("")

    # ── Collect variant positions and ClinVar significance ───────────────────
    pos_to_vars: dict[int, set[str]] = defaultdict(set)
    var_to_sig:  dict[str, str]      = {}
    for _, row in df.iterrows():
        var = str(row.get("variant_info", "")).strip()
        if not var:
            continue
        pos = _variant_aa_position(var, str(row.get("Name", "")).strip())
        if pos is None:
            continue
        pos_to_vars[pos].add(var)
        sig_raw = str(row.get("significance_1", "")).strip()
        cat = _clinsig_category(sig_raw)
        # Keep the most informative (non-unknown) annotation per variant.
        if var not in var_to_sig or cat != "unknown":
            var_to_sig[var] = cat

    if not pos_to_vars:
        return None

    ordered_all = {p: sorted(pos_to_vars[p], key=_sort_text) for p in sorted(pos_to_vars)}
    max_pos = max(ordered_all)

    # ── Apply ClinVar filter (unknown/unannotated variants are always shown) ─
    if clin_filter is None:
        clin_filter = list(_CLINVAR_DOT_COLORS.keys())
    active_cats = set(clin_filter) | {"unknown"}
    ordered: dict[int, list[str]] = {}
    for pos, labs in ordered_all.items():
        kept = [lb for lb in labs if var_to_sig.get(lb, "unknown") in active_cats]
        if kept:
            ordered[pos] = kept
    if not ordered:
        ordered = {}  # empty – figure will show protein bar only

    # ── Layout constants ─────────────────────────────────────────────────────
    y_track  = 0.0     # vertical centre of the protein bar
    bar_half = 0.18    # half-height of the backbone bar
    y_var    = 0.74    # vertical centre of the lowest (first) disc per position
    y_step   = 0.72    # spacing between stacked discs at the same position
    dot_off  = -0.32   # dot centre offset below disc centre (data units)

    max_stack = max((len(v) for v in ordered.values()), default=1)

    # ── Collect per-disc data (x, y, label, sig) ─────────────────────────────
    vx:   list = []
    vy:   list = []
    vtxt: list = []
    vhov: list = []
    vsig: list = []
    for pos, labs in ordered.items():
        for i, lab in enumerate(labs):
            vx.append(float(pos))
            vy.append(y_var + i * y_step)
            vtxt.append(lab)
            sig = var_to_sig.get(lab, "unknown")
            cat_label = {
                "pathogenic": "Pathogenic / likely pathogenic",
                "benign": "Benign / likely benign",
                "vus": "VUS / conflicting",
            }.get(sig, "")
            suffix = f" · {cat_label}" if cat_label else ""
            vhov.append(f"{lab} · aa {pos}{suffix}<extra></extra>")
            vsig.append(sig)

    fig = go.Figure()
    fill_protein = _graph_node_fill_hex(GRAPH_PROTEIN_BG, U["card"])
    fill_variant = _graph_node_fill_hex(GRAPH_VARIANT_BG, U["card"])

    # ── Protein backbone bar ──────────────────────────────────────────────────
    fig.add_shape(
        type="rect",
        xref="paper", yref="y",
        x0=0.0, x1=0.93,
        y0=y_track - bar_half, y1=y_track + bar_half,
        fillcolor=fill_protein, line=dict(width=0),
        layer="below",
    )

    # ── Stems: vertical line from disc bottom to protein bar top ─────────────
    stem_x: list = []
    stem_y: list = []
    # Group stems by position: the lowest disc drives the stem bottom
    pos_to_lowest_y: dict[float, float] = {}
    for px, py in zip(vx, vy):
        if px not in pos_to_lowest_y or py < pos_to_lowest_y[px]:
            pos_to_lowest_y[px] = py
    for px, lowest_y in pos_to_lowest_y.items():
        stem_x.extend([px, px, None])
        stem_y.extend([lowest_y, y_track + bar_half, None])
    fig.add_trace(go.Scatter(
        x=stem_x, y=stem_y, mode="lines",
        line=dict(color=fill_variant, width=1.4),
        showlegend=False, hoverinfo="skip",
    ))

    # ── Variant discs ─────────────────────────────────────────────────────────
    fig.add_trace(go.Scatter(
        x=vx, y=vy, mode="markers+text",
        marker=dict(size=28, color=fill_variant, line=dict(width=0)),
        text=vtxt, textposition="middle center",
        textfont=dict(size=10, color=GRAPH_VARIANT_FG, family=U["font_ui"], weight=600),
        hovertemplate=vhov,
        showlegend=False,
    ))

    # ── ClinVar pathogenicity dots (below each disc) ──────────────────────────
    dot_x:     list = []
    dot_y:     list = []
    dot_color: list = []
    dot_hover: list = []
    for px, py, sig in zip(vx, vy, vsig):
        col = _CLINVAR_DOT_COLORS.get(sig)
        if col:
            dot_x.append(px)
            dot_y.append(py + dot_off)
            dot_color.append(col)
            cat_lbl = {
                "pathogenic": "Pathogenic / likely pathogenic",
                "benign": "Benign / likely benign",
                "vus": "VUS / conflicting",
            }[sig]
            dot_hover.append(f"ClinVar: {cat_lbl}<extra></extra>")
    if dot_x:
        fig.add_trace(go.Scatter(
            x=dot_x, y=dot_y, mode="markers",
            marker=dict(
                size=8,
                color=dot_color,
                line=dict(width=1.2, color=U["card"]),
                symbol="circle",
            ),
            hovertemplate=dot_hover,
            showlegend=False,
        ))

    # ── Protein cap node (secondary x-axis, always fixed at right) ────────────
    fig.add_trace(go.Scatter(
        x=[0.5], y=[y_track],
        mode="markers+text",
        marker=dict(size=32, color=fill_protein, line=dict(width=0)),
        text=[prot], textposition="middle center",
        textfont=dict(size=12, color=GRAPH_PROTEIN_FG, family=U["font_ui"], weight=600),
        hovertemplate=f"<b>{prot}</b><extra></extra>",
        showlegend=False,
        xaxis="x2",
    ))

    # ── Axis ranges & figure sizing ───────────────────────────────────────────
    _pad   = float(max_pos) * 0.03 + 15
    x_lo   = -_pad
    x_hi   = float(max_pos) * 1.06 + 30
    y_top  = y_var + (max_stack - 1) * y_step + 0.55
    fig_h  = max(120, int(60 + (y_top + 0.52) * 52))

    fig.update_layout(
        xaxis2=dict(
            range=[0, 1], domain=[0.87, 1.0],
            visible=False, fixedrange=True,
        ),
        margin=dict(l=54, r=14, t=14, b=28),
        plot_bgcolor=U["card"],
        paper_bgcolor=U["card"],
        font=dict(family=U["font_ui"], size=10, color=U["ink_soft"]),
        xaxis=dict(
            title="Amino acid position",
            range=[x_lo, x_hi],
            domain=[0, 0.87],
            showgrid=True, gridcolor="rgba(45,42,36,0.06)",
            zeroline=False,
            tickfont=dict(size=9, color=U["ink_soft"], family=U["font_ui"]),
            fixedrange=False,
        ),
        yaxis=dict(visible=False, range=[-0.52, y_top], fixedrange=True),
        height=fig_h,
        hovermode="closest",
        dragmode="zoom",
    )
    return fig


def _relayout_xaxis_range(relayout: Optional[dict]) -> Optional[tuple[float, float]]:
    """Parse Plotly relayout dict for x-axis range; None means no zoom box / full span."""
    if not relayout:
        return None
    if relayout.get("xaxis.autorange") is True:
        return None
    k0, k1 = "xaxis.range[0]", "xaxis.range[1]"
    if k0 not in relayout or k1 not in relayout:
        return None
    try:
        a = float(relayout[k0])
        b = float(relayout[k1])
    except (TypeError, ValueError):
        return None
    return (a, b) if a <= b else (b, a)


def _norm_map_range(mr) -> Optional[tuple[float, float]]:
    if not mr or not isinstance(mr, (list, tuple)) or len(mr) != 2:
        return None
    try:
        a, b = float(mr[0]), float(mr[1])
    except (TypeError, ValueError):
        return None
    return (a, b) if a <= b else (b, a)
