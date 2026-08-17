"""Sidebar info builders shown on the network page when a node or edge is
selected (or by default, before any selection)."""
import urllib.parse as _url

from dash import dcc, html

from indra_variants.app.utils import U


def _sidebar_default():
    return [
        html.Div("Details",
                 style={'fontSize': 17, 'fontWeight': 600,
                        'fontFamily': U['font_display'],
                        'marginBottom': 15, 'color': U['ink'],
                        'borderBottom': f'1px solid {U["rule"]}',
                        'paddingBottom': 10}),
        html.Div("Select a node or edge for attributes, evidence, and external links.",
                 style={'color': U['muted'], 'fontSize': 14, 'lineHeight': '1.45',
                        'fontFamily': U['font_ui']})
    ]


def _sidebar_card(children):
    return html.Div(children, style={
        'background': U['card'],
        'padding': 12,
        'borderRadius': 4,
        'boxShadow': U['shadow'],
        'marginBottom': 15,
        'border': f'1px solid {U["rule"]}',
        'fontFamily': U['font_ui'],
    })


def _build_edge_info(edge):
    if not edge:
        return _sidebar_default()

    content = [
        html.Div("Edge",
                 style={'fontSize': 18, 'fontWeight': 'bold',
                        'marginBottom': 15, 'color': U['ink'],
                        'borderBottom': f'1px solid {U["rule"]}',
                        'paddingBottom': 10})
    ]

    rel = edge.get('rel', 'N/A')
    source = edge.get('source', 'N/A')
    target = edge.get('target', 'N/A')

    content.append(_sidebar_card([
        html.Div("Relation",
                 style={'fontSize': 14, 'fontWeight': 'bold',
                        'color': U['ink_soft'], 'marginBottom': 5}),
        html.Div(f"{source} → {target}",
                 style={'fontSize': 14, 'color': U['ink'],
                        'marginBottom': 8, 'fontWeight': 'bold'}),
        *([] if rel in ['DV', 'PV', 'has_domain'] else [
            html.Div(f"Type: {rel}",
                     style={'fontSize': 14, 'color': U['muted']})
        ]),
        *([] if not edge.get('evidence_count') else [
            html.Div(f"Statements: {edge['evidence_count']}",
                     style={'fontSize': 13, 'color': U['muted'], 'marginTop': 4})
        ])
    ]))

    if rel == 'DV':
        if edge.get('note'):
            content.append(_sidebar_card([
                html.Div("Domain",
                         style={'fontSize': 14, 'fontWeight': 'bold',
                                'color': U['ink_soft'], 'marginBottom': 8}),
                html.Div(edge['note'],
                         style={'fontSize': 14, 'color': U['ink'],
                                'lineHeight': '1.4'})
            ]))
    else:
        if edge.get('note'):
            content.append(_sidebar_card([
                html.Div("Description",
                         style={'fontSize': 14, 'fontWeight': 'bold',
                                'color': U['ink_soft'], 'marginBottom': 8}),
                html.Div(edge['note'],
                         style={'fontSize': 14, 'color': U['ink'],
                                'lineHeight': '1.4'})
            ]))

    if edge.get('clinvar_data'):
        data = edge['clinvar_data']
        content.append(_sidebar_card([
            html.Div("ClinVar",
                     style={'fontSize': 14, 'fontWeight': 'bold',
                            'color': U['ink_soft'], 'marginBottom': 10}),
            html.Div([
                html.Div([
                    html.Span("Pathogenicity: ", style={'fontWeight': 'bold'}),
                    html.Span(data.get('pathogenicity', 'N/A'))
                ], style={'marginBottom': 6}),
                html.Div([
                    html.Span("Review: ", style={'fontWeight': 'bold'}),
                    html.Span(data.get('review', 'N/A'))
                ], style={'marginBottom': 6}),
                html.Div([
                    html.Span("Condition: ", style={'fontWeight': 'bold'}),
                    html.Div(data.get('conditions', 'N/A'),
                             style={'marginTop': 4, 'fontStyle': 'italic'})
                ])
            ], style={'fontSize': 13, 'color': U['ink'], 'lineHeight': '1.4'})
        ]))

    if edge.get('pmid'):
        pmid_raw = edge['pmid']
        pmid_list = [p.strip() for p in pmid_raw.replace(";", ",").split(",") if p.strip()]
        pubmed_links = []
        for pid in pmid_list:
            pubmed_links.append(
                html.A(f"PMID {pid}",
                       href=f"https://pubmed.ncbi.nlm.nih.gov/{pid}/",
                       target="_blank",
                       style={'display': 'inline-block', 'marginRight': 10,
                              'marginBottom': 4, 'color': U['link'],
                              'textDecoration': 'none', 'fontSize': 13})
            )
        content.append(_sidebar_card([
            html.Div("Links",
                     style={'fontSize': 14, 'fontWeight': 'bold',
                            'color': U['ink_soft'], 'marginBottom': 10}),
            html.Div(pubmed_links,
                     style={'marginBottom': 8, 'lineHeight': '1.8'}),
            html.A("INDRA Discovery",
                   href=(f"https://discovery.indra.bio/search/"
                         f"?agent={_url.quote_plus(edge['src4indra'])}"
                         f"&other_agent={_url.quote_plus(edge['target'].split('::')[0])}"
                         "&agent_role=subject&other_role=object"),
                   target="_blank",
                   style={'display': 'block', 'color': U['link'],
                          'textDecoration': 'none', 'fontSize': 13,
                          'fontWeight': 'bold'})
        ]))

    return content


def _build_node_info(node):
    if not node:
        return _sidebar_default()

    real_name = node.get('real', node.get('label', 'N/A'))
    role = node.get('role', 'unknown').replace('-', ' ').replace('_', ' ').title()

    content = [
        html.Div("Node",
                 style={'fontSize': 18, 'fontWeight': 'bold',
                        'marginBottom': 15, 'color': U['ink'],
                        'borderBottom': f'1px solid {U["rule"]}',
                        'paddingBottom': 10})
    ]

    content.append(_sidebar_card([
        html.Div("Name",
                 style={'fontSize': 14, 'fontWeight': 'bold',
                        'color': U['ink_soft'], 'marginBottom': 5}),
        html.Div(real_name,
                 style={'fontSize': 14, 'color': U['ink'],
                        'marginBottom': 8, 'fontWeight': 'bold',
                        'lineHeight': '1.4'}),
        html.Div(f"Role: {role}",
                 style={'fontSize': 13, 'color': U['muted']})
    ]))

    stat_lines = []
    if node.get('n_proteins') is not None:
        stat_lines.append(html.Div([
            html.Span("Proteins: ", style={'fontWeight': 'bold'}),
            html.Span(str(node['n_proteins']))
        ], style={'marginBottom': 6}))
    if node.get('n_variants') is not None:
        stat_lines.append(html.Div([
            html.Span("Variants: ", style={'fontWeight': 'bold'}),
            html.Span(str(node['n_variants']))
        ], style={'marginBottom': 6}))
    if node.get('n_records') is not None:
        stat_lines.append(html.Div([
            html.Span("Source rows: ", style={'fontWeight': 'bold'}),
            html.Span(str(node['n_records']))
        ]))
    if stat_lines:
        content.append(_sidebar_card([
            html.Div("Counts",
                     style={'fontSize': 14, 'fontWeight': 'bold',
                            'color': U['ink_soft'], 'marginBottom': 10}),
            *stat_lines
        ]))

    if node.get('domain_notes'):
        content.append(_sidebar_card([
            html.Div("Domain notes",
                     style={'fontSize': 14, 'fontWeight': 'bold',
                            'color': U['ink_soft'], 'marginBottom': 8}),
            html.Div(node['domain_notes'],
                     style={'fontSize': 14, 'color': U['ink'],
                            'lineHeight': '1.4'})
        ]))

    if node.get('clinvar_data'):
        data = node['clinvar_data']
        content.append(_sidebar_card([
            html.Div("ClinVar",
                     style={'fontSize': 14, 'fontWeight': 'bold',
                            'color': U['ink_soft'], 'marginBottom': 10}),
            html.Div([
                html.Div([
                    html.Span("Pathogenicity: ", style={'fontWeight': 'bold'}),
                    html.Span(data.get('pathogenicity', 'N/A'))
                ], style={'marginBottom': 6}),
                html.Div([
                    html.Span("Review: ", style={'fontWeight': 'bold'}),
                    html.Span(data.get('review', 'N/A'))
                ], style={'marginBottom': 6}),
                html.Div([
                    html.Span("Condition: ", style={'fontWeight': 'bold'}),
                    html.Div(data.get('conditions', 'N/A'),
                             style={'marginTop': 4, 'fontStyle': 'italic'})
                ])
            ], style={'fontSize': 13, 'color': U['ink'], 'lineHeight': '1.4'})
        ]))

    # For variant nodes use the gene symbol so INDRA resolves the query correctly.
    _indra_agent = node.get('gene_symbol') or real_name
    external_links = [
        html.A("INDRA Discovery",
               href=f"https://discovery.indra.bio/search/?agent={_url.quote_plus(_indra_agent)}",
               target="_blank",
               style={'display': 'block', 'color': U['link'],
                      'textDecoration': 'none', 'fontSize': 13,
                      'fontWeight': 'bold', 'marginBottom': 6})
    ]
    if node.get('protein_page'):
        external_links.insert(
            0,
            dcc.Link("Protein-centric graph",
                     href=node['protein_page'],
                     style={'display': 'block', 'color': U['link'],
                            'textDecoration': 'none', 'fontSize': 13,
                            'fontWeight': 'bold', 'marginBottom': 6})
        )
    if node.get('uniprot_id'):
        uid = node['uniprot_id']
        if '_' in uid:
            uniprot_href = (
                f"https://www.uniprot.org/uniprotkb/"
                f"{_url.quote_plus(uid)}/entry"
            )
        else:
            uniprot_href = (
                "https://www.uniprot.org/uniprotkb?query="
                f"{_url.quote_plus(f'gene_exact:{uid} AND organism_id:9606')}"
            )
        external_links.append(
            html.A("UniProt",
                   href=uniprot_href,
                   target="_blank",
                   style={'display': 'block', 'color': U['link'],
                          'textDecoration': 'none', 'fontSize': 13,
                          'marginBottom': 6})
        )
    if node.get('clinvar_allele'):
        external_links.append(
            html.A("ClinVar",
                   href=("https://www.ncbi.nlm.nih.gov/clinvar/?term="
                         f"{_url.quote_plus(str(node['clinvar_allele']) + '[alleleid]')}"),
                   target="_blank",
                   style={'display': 'block', 'color': U['link'],
                          'textDecoration': 'none', 'fontSize': 13,
                          'marginBottom': 6})
        )
    elif node.get('gene_symbol'):
        clinvar_query = f"{node['gene_symbol']}[gene] AND {real_name}"
        external_links.append(
            html.A("ClinVar (search)",
                   href=("https://www.ncbi.nlm.nih.gov/clinvar/?term="
                         f"{_url.quote_plus(clinvar_query)}"),
                   target="_blank",
                   style={'display': 'block', 'color': U['link'],
                          'textDecoration': 'none', 'fontSize': 13,
                          'marginBottom': 6})
        )
    if node.get('dbsnp_rs'):
        external_links.append(
            html.A("dbSNP",
                   href=f"https://www.ncbi.nlm.nih.gov/snp/rs{node['dbsnp_rs']}",
                   target="_blank",
                   style={'display': 'block', 'color': U['link'],
                          'textDecoration': 'none', 'fontSize': 13})
        )

    content.append(_sidebar_card([
        html.Div("Links",
                 style={'fontSize': 14, 'fontWeight': 'bold',
                        'color': U['ink_soft'], 'marginBottom': 10}),
        *external_links
    ]))

    return content
