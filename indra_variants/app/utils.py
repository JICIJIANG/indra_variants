"""Small foundational utilities: routing helpers, tiny pure functions, and
shared UI theme constants/styles.  Imported by almost every other module
in ``indra_variants.app``."""
import re
import urllib.parse as _url
from typing import Optional


# ---------- Variant-name regular expressions ---------------------------------
AA_SUB_RE = re.compile(r"^[A-Za-z\*]+(?P<pos>\d+)[A-Za-z\*=]+$")
P_DOT_RE = re.compile(r"p\.[A-Za-z\*]+(?P<pos>\d+)[A-Za-z\*=]+", re.I)


# ---------- Sorting / bucketing ---------------------------------------------
def _sort_text(value: str) -> str:
    return value.casefold()


def _alpha_bucket(value: str) -> str:
    if not value:
        return "#"
    first = value[0].upper()
    return first if "A" <= first <= "Z" else "#"


# ---------- Route helpers ----------------------------------------------------
def _encode_route_value(value: str) -> str:
    return _url.quote(value, safe="")


def _decode_route_value(value: str) -> str:
    return _url.unquote(value)


def _protein_href(prot: str) -> str:
    return f"/protein/{_encode_route_value(prot)}"


def _endpoint_href(endpoint: str) -> str:
    return f"/endpoint/{_encode_route_value(endpoint)}"


def _css_safe_global(name: str) -> str:
    return re.sub(r'[^A-Za-z0-9_-]', '-', name)


# ---------- ClinVar / variant helpers ----------------------------------------
def format_star_rating(star_val):
    review_map = {
        4.0: "practice guideline",
        3.0: "expert panel",
        2.0: "multiple submitters, no conflicts",
        1.0: "single submitter",
        0.0: "no assertion criteria provided" # Updated for clarity
    }
    try:
        s_val = float(star_val)
        num_stars = int(s_val)
        review_text = review_map.get(s_val, "review status not specified")

        if num_stars > 0:
            return f"{'★ ' * num_stars} ({review_text})"
        else:
            return f"({review_text})"

    except (ValueError, TypeError):
        return "(no review info)"


def _variant_aa_position(variant_label: str, name_label: str = "") -> Optional[int]:
    m = AA_SUB_RE.match((variant_label or "").strip())
    if m:
        return int(m.group("pos"))
    m = P_DOT_RE.search((name_label or "").strip())
    if m:
        return int(m.group("pos"))
    return None


def _clinsig_category(sig: str) -> str:
    """Map a raw ClinVar significance string to one of: pathogenic, benign, vus, unknown."""
    s = (sig or "").lower().strip()
    if not s or s in ("n/a", "nan", "not provided"):
        return "unknown"
    if "conflicting" in s:
        return "vus"
    if "pathogenic" in s:
        return "pathogenic"
    if "benign" in s:
        return "benign"
    if "uncertain" in s or "vus" in s:
        return "vus"
    return "unknown"


# --- Shared UI theme (warm paper / ink; avoids generic “AI gradient” look) ---
U = {
    "font_ui": (
        "system-ui, -apple-system, 'Segoe UI', Roboto, 'Helvetica Neue', "
        "Arial, sans-serif"
    ),
    "font_display": "Georgia, Cambria, 'Times New Roman', Times, serif",
    "ink": "#1c1b18",
    "ink_soft": "#454036",
    "muted": "#6f6b63",
    "paper": "#f2efe6",
    "panel": "#fdfbf7",
    "card": "#fffcf7",
    "wash": "#ebe6dc",
    "rule": "#e0dbd2",
    "link": "#2a4a66",
    "hero": "#262422",
    "hero_hi": "#32302c",
    "hero_deep": "#1a1917",
    "hero_text": "#f7f4ec",
    "hero_muted": "#c4beb4",
    "shadow": "0 1px 2px rgba(28, 27, 24, 0.05), 0 6px 20px rgba(28, 27, 24, 0.06)",
    "shadow_strong": "0 2px 8px rgba(28, 27, 24, 0.1)",
    "chart_bp": "#4d5f52",
    "chart_dis": "#8b534c",
    "chart_gene": "#4d5866",
    "plot_bp": "#eef1ec",
    "plot_dis": "#f3eceb",
    "plot_gene": "#eceff2",
    "accent_card_bp": "#6d7f72",
    "accent_card_dis": "#a1665f",
    "accent_card_gene": "#6a7484",
    "graph_bg": "#f7f5f0",
    "legend_bg": "rgba(253, 251, 247, 0.96)",
}

GRAPH_PROTEIN_BG = "#c5d2ce"
GRAPH_PROTEIN_FG = "#2a3d38"
GRAPH_VARIANT_BG = "#c9c0d4"
GRAPH_VARIANT_FG = "#3d324d"
GRAPH_INTERMEDIATE_BG = "#cfd9c3"
GRAPH_INTERMEDIATE_FG = "#35422e"
GRAPH_ENDPOINT_BG = "#e8d4bc"
GRAPH_ENDPOINT_FG = "#5c3f24"

# Cytoscape main graph `node` uses this opacity over `backgroundColor` (U["paper"]).
GRAPH_NODE_BG_OPACITY = 0.92


def _hex_to_rgb(h: str) -> tuple[int, int, int]:
    h = (h or "").strip().lstrip("#")
    if len(h) != 6:
        return 0, 0, 0
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


def _rgb_to_hex(r: int, g: int, b: int) -> str:
    return f"#{r:02x}{g:02x}{b:02x}"


def _graph_node_fill_hex(fg_hex: str, blend_on: Optional[str] = None) -> str:
    """Blend fg onto canvas (default paper, e.g. card for the variant strip)."""
    fg = _hex_to_rgb(fg_hex)
    bg = _hex_to_rgb(blend_on or U["paper"])
    a = GRAPH_NODE_BG_OPACITY
    return _rgb_to_hex(
        int(fg[0] * a + bg[0] * (1 - a)),
        int(fg[1] * a + bg[1] * (1 - a)),
        int(fg[2] * a + bg[2] * (1 - a)),
    )


def _network_frame_style() -> dict:
    """White functional panel like stats chart cards."""
    return {
        "background": U["card"],
        "border": f"1px solid {U['rule']}",
        "borderRadius": 4,
        "boxShadow": U["shadow"],
        "overflow": "hidden",
    }


def _variant_strip_wrap_style() -> dict:
    """Variant strip: white surface, no frame (per design), allow markers to paint fully."""
    return {
        "background": U["card"],
        "border": "none",
        "boxShadow": "none",
        "borderRadius": 4,
        "overflow": "visible",
    }


# Header row under URL bar (approx.) — remainder is graph + bottom strip
_NETWORK_HEADER_PX = 96
_GENE_STRIP_VH_PCT = 38
_GENE_STRIP_VH = f"{_GENE_STRIP_VH_PCT}vh"
_GENE_STRIP_GRAPH_INNER_H = f"calc({_GENE_STRIP_VH_PCT}vh - 56px)"

# ── ClinVar pathogenicity colour palette ─────────────────────────────────────
_CLINVAR_DOT_COLORS: dict[str, str] = {
    "pathogenic": "#e03131",   # red
    "benign":     "#74c69d",   # light green
    "vus":        "#2f9e44",   # dark green
}
