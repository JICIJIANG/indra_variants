"""Entry point for the INDRA variant-network Dash app.

 Run with
``python variant_network.py`` to start a local development server.

Module breakdown:

* :mod:`indra_variants.app.config`         — environment-driven settings
* :mod:`indra_variants.app.utils`          — routing helpers, theme, tiny utils
* :mod:`indra_variants.app.data_index`     — TSV discovery, PROTS/ENDPOINTS index
* :mod:`indra_variants.app.stats`          — reference tables + bar charts
* :mod:`indra_variants.app.crossing_min`   — LNS layer-ordering algorithm
* :mod:`indra_variants.app.graph_builders` — protein/endpoint network builders
* :mod:`indra_variants.app.figures`        — protein lollipop + relayout helpers
* :mod:`indra_variants.app.sidebar`        — node/edge info panels
* :mod:`indra_variants.app.pages`          — page renderers (search, stats, network)
* :mod:`indra_variants.app.callbacks`      — all @app.callback definitions
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import dash
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
from dash import dcc, html

from indra_variants.app import callbacks
from indra_variants.app.config import DEBUG, PORT


cyto.load_extra_layouts()


# ------------------------Dash App------------------------–
app = dash.Dash(__name__,
                title="INDRA Variant",
                suppress_callback_exceptions=True,
                external_stylesheets=[dbc.themes.SANDSTONE])

server = app.server
app.layout = html.Div([dcc.Location(id="url"), html.Div(id="page")])


callbacks.register(app)


# --------------------Run App----------------------------–
if __name__ == "__main__":
    app.run(debug=DEBUG, port=PORT, host="0.0.0.0")
