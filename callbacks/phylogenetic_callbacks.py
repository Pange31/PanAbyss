#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 22:08:19 2025

@author: fgraziani
"""

import base64
import dash
from dash import html, Input, Output, callback, State, callback_context, dcc, ctx, exceptions, no_update
from dash.exceptions import PreventUpdate



import os
import sys
import io
import math

from Bio import Phylo
root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if root_path not in sys.path:
    sys.path.append(root_path)
from app import *
from neo4j_requests import *
import logging

last_tree = "./export/phylo/last_tree.nwk"
logger = logging.getLogger("panabyss_logger")

EXPORT_DIR = "./export/phylo/"

def generate_elements(newick_str, xlen=30, ylen=30, grabbable=False):
    tree = Phylo.read(io.StringIO(newick_str), "newick")
    def get_col_positions(tree, column_width=80):
        taxa = tree.get_terminals()

        # Some constants for the drawing calculations
        max_label_width = max(len(str(taxon)) for taxon in taxa)
        drawing_width = column_width - max_label_width - 1
    
        """Create a mapping of each clade to its column position."""
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
        fudge_margin = int(math.ceil(math.log(len(taxa), 2)))
        cols_per_branch_unit = ((drawing_width - fudge_margin) /
                                float(max(depths.values())))
        return dict((clade, int(blen * cols_per_branch_unit + 1.0))
                    for clade, blen in depths.items())

    def get_row_positions(tree):
        taxa = tree.get_terminals()
        positions = dict((taxon, 2 * idx) for idx, taxon in enumerate(taxa))
    
        def calc_row(clade):
            for subclade in clade:
                if subclade not in positions:
                    calc_row(subclade)
            positions[clade] = ((positions[clade.clades[0]] +
                                 positions[clade.clades[-1]]) // 2)
    
        calc_row(tree.root)
        return positions
    
    def add_to_elements(clade, clade_id):
        children = clade.clades
    
        pos_x = col_positions[clade] * xlen
        pos_y = row_positions[clade] * ylen
    
        cy_source = {
            "data": {"id": clade_id},
            'position': {'x': pos_x, 'y': pos_y},
            'classes': 'nonterminal',
            'grabbable': grabbable
        }
        nodes.append(cy_source)
    
        if clade.is_terminal():
            cy_source['data']['name'] = clade.name
            cy_source['classes'] = 'terminal'
    
        for n, child in enumerate(children):
            support_id = clade_id + 's' + str(n)
            child_id = clade_id + 'c' + str(n)
            pos_y_child = row_positions[child] * ylen
    
            cy_support_node = {
                'data': {'id': support_id},
                'position': {'x': pos_x, 'y': pos_y_child},
                'grabbable': grabbable,
                'classes': 'support'
            }
    
            cy_support_edge = {
                'data': {
                    'source': clade_id,
                    'target': support_id,
                    'sourceCladeId': clade_id
                },
            }
    
            cy_edge = {
                'data': {
                    'source': support_id,
                    'target': child_id,
                    'length': clade.branch_length,
                    'sourceCladeId': clade_id
                },
            }
    
            if clade.confidence and clade.confidence.value:
                cy_source['data']['confidence'] = clade.confidence.value
    
            nodes.append(cy_support_node)
            edges.extend([cy_support_edge, cy_edge])
    
            add_to_elements(child, child_id)

    col_positions = get_col_positions(tree)
    row_positions = get_row_positions(tree)
    
    nodes = []
    edges = []
    
    add_to_elements(tree.clade, 'r')
    
    return nodes+edges




stylesheet = [
    {
        'selector': '.nonterminal',
        'style': {
            'label': 'data(confidence)',
            'background-opacity': 0,
            "text-halign": "left",
            "text-valign": "top",
        }
    },
    {
        'selector': 'node',
        'style': {
            'label': 'data(name)',
            'font-size':'50px',
            "text-halign": "left",
            "text-valign": "center"
        }
    },
    {
        'selector': '.support',
        'style': {'background-opacity': 0}
    },
    {
        'selector': 'edge',
        'style': {
            "source-endpoint": "inside-to-node",
            "target-endpoint": "inside-to-node",
        }
    },
    {
        'selector': '.terminal',
        'style': {
            'label': 'data(name)',
            'width': 10,
            'height': 10,
            "text-valign": "center",
            "text-halign": "right",
            'background-color': '#222222'
        }
    }
]


#Populate chromosome droplist
@app.callback(
    Output('phylogenetic_chromosomes_dropdown', 'options'),
    Input('shared_storage', 'data')
)
def update_dropdown(data):
    if not data or "chromosomes" not in data:
        return []

    chromosomes = data["chromosomes"]
    options = [{"label": str(chrom), "value": str(chrom)} for chrom in chromosomes]
    return options



@app.callback(
    Output("phylo-job-trigger", "data", allow_duplicate=True),
    Output("phylo-spinner-container", "style", allow_duplicate=True),
    Output("phylo-job-status", "data", allow_duplicate=True),
    Output("phylogenetic-page-store", "data", allow_duplicate=True),
    Input('upload-newick', 'contents'),
    Input('btn-plot-global-tree', 'n_clicks'),
    Input("btn-load-last-tree", "n_clicks"),
    State("phylogenetic-page-store", "data"),
    prevent_initial_call=True
)
def trigger_job(contents, n_clicks_global_tree, n_clicks_last_tree, phylo_data):

    triggered_id = ctx.triggered_id
    if triggered_id not in ["upload-newick","btn-plot-global-tree","btn-load-last-tree"]:
        raise exceptions.PreventUpdate
    if triggered_id=="upload-newick" and contents is None :
        raise exceptions.PreventUpdate
    if phylo_data is not None :
        phylo_data.pop("newick_global", None)
    return {"run": True, "contents":contents, "triggered_id":triggered_id}, {"display": "block", "marginTop": "20px"}, {"status":"running"}, phylo_data


#Callback to plot the global tree
@app.callback(
    Output("phylogenetic-page-store", "data", allow_duplicate=True),
    Output("global-notification", "data", allow_duplicate=True),
    Output("phylo-job-status", "data", allow_duplicate=True),
    Input("phylo-job-trigger", "data"),
    State('method-dropdown', 'value'),
    State('upload-newick', 'filename'),
    State("phylogenetic_chromosomes_dropdown", 'value'),
    State("phylogenetic-page-store", "data"),
    prevent_initial_call=True,
    background=True,
    # running=[
    #     (Output("btn-plot-global-tree", "disabled"), True, False),
    #     #(Output("phylo-spinner-container", "style"), {"display": "block", "marginTop": "20px"}, {"display": "none", "marginTop": "20px"}),
    #     (Output("btn-load-last-tree", "disabled"), True, False),
    #     (Output("btn-cancel-plot-global-tree", "disabled"), False, True)
    # ],
    cancel=[Input("btn-cancel-plot-global-tree", "n_clicks")],
    suppress_callback_exceptions=True
)
def update_phylo_graph(trigger_data, method, filename, chromosome, phylo_data):
    triggered_id = trigger_data.get("triggered_id", None)
    if not triggered_id or not trigger_data:
        raise exceptions.PreventUpdate
    else:
        contents = trigger_data.get("contents", None)

    if phylo_data is None:
        phylo_data = {}

    if triggered_id == "btn-plot-global-tree":
        if chromosome == None or chromosome == "":
            c = None
        else:
            c = chromosome
        logger.debug(f"Compute global phylo tree using {method} method.")
        newick_str = compute_global_phylo_tree_from_nodes(method=method, chromosome=c)
        phylo_data["message"] = f"Tree successfully computed."
    elif triggered_id == "btn-load-last-tree" :
        if os.path.exists(last_tree):
            with open(last_tree, "r") as f:
                newick_str = f.read()
                phylo_data["message"] = f"File '{filename}' successfully load."

        else:
            toast_message = {
                "title": "Global tree construction error.",
                "message": f"No existing last tree file : {last_tree}",
                "type": "danger",
            }
            phylo_data["message"] = f"No existing last tree file : {last_tree}"
            status_data = {"status":"done"}
            return  phylo_data, toast_message, status_data
    else:
        if contents is None:
            toast_message = {
                "title": "Global tree construction error.",
                "message": "No file loaded.",
                "type": "danger",
            }
            phylo_data["message"] = "No file loaded."
            status_data = {"status":"done"}
            return  phylo_data, toast_message, status_data
    
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        
        try:
            newick_str = decoded.decode('utf-8')
        except Exception as e:
            toast_message = {
                "title": "Global tree construction error.",
                "message": f"Parsing error : {str(e)}",
                "type": "danger",
            }
            phylo_data["message"] = f"Parsing error : {str(e)}"
            status_data = {"status":"done"}
            return  phylo_data, toast_message, status_data
    if newick_str is not None and newick_str != "":
        phylo_data["newick_global"] = newick_str
        toast_message = {
            "title": "Global tree construction.",
            "message": "Process terminated.",
            "type": "success"
        }
        status_data = {"status":"done"}
        return  phylo_data, toast_message, status_data
    else:
        toast_message = {
            "title": "Global tree construction error.",
            "message": f"Unable to compute global tree.",
            "type": "danger",
        }
        phylo_data["message"] = f"Unable to compute global tree."
        status_data = {"status":"done"}
        return phylo_data, toast_message, status_data


@app.callback(
    Output('phylo-job-status', 'data', allow_duplicate=True),
    Input('btn-cancel-plot-global-tree', 'n_clicks'),
    prevent_initial_call=True,
)
def handle_cancel_click(n_clicks):
    if not n_clicks:
        return no_update
    status_data = {"status": "done"}
    return status_data


@app.callback(
    Output('cytoscape-phylo', 'stylesheet', allow_duplicate=True),
    Input('cytoscape-phylo', 'mouseoverEdgeData'),
    prevent_initial_call=True
    )
def color_children(edgeData):
    if edgeData is None:
        return stylesheet

    if 's' in edgeData['source']:
        val = edgeData['source'].split('s')[0]
    else:
        val = edgeData['source']

    children_style = [{
        'selector': 'edge[source *= "{}"]'.format(val),
        'style': {
            'line-color': 'blue'
        }
    }]

    return stylesheet + children_style

#Callback to plot tree of the current displayed region
@app.callback(
    #Output('cytoscape-phylo-region', 'elements', allow_duplicate=True),
    Output("phylogenetic-message", "children", allow_duplicate=True),
    Output("phylogenetic-page-store", "data", allow_duplicate=True),
    Output("phylo-local-tree-job-status", "data"),
    Input('btn-plot-region', 'n_clicks'),
    State('shared_storage_nodes', 'data'),
    State("phylogenetic-page-store", "data"),
    prevent_initial_call=True
)
def plot_region(n_clicks, stored_data, phylo_data):
    if not n_clicks:
        raise exceptions.PreventUpdate

    ctx = dash.callback_context
    if ctx.triggered_id == "btn-plot-region" and n_clicks == 0:
        raise PreventUpdate
    phylo_local_data = {"status": "done"}
    if not stored_data:
        phylo_local_data = {"status": "done"}
        return html.Div(html.P([
        "❌ No data to compute tree. Select a region to visualise on the ",
        dcc.Link("home page", href="/", style={'color': 'blue', 'textDecoration': 'underline'}),
        " or on the ",
        dcc.Link("gwas page", href="/gwas", style={'color': 'blue', 'textDecoration': 'underline'})
        ], style=error_style)), phylo_data, phylo_local_data

    try:
        # Step 1 : compute tree of the region
        newick_str = compute_phylo_tree_from_nodes(stored_data)
        if phylo_data is None:
            phylo_data = {"newick_region":newick_str}
        else:
            phylo_data["newick_region"] = newick_str

        # Step2 : draw tree
        phylo_local_data = {"status": "done"}
        return "", phylo_data, phylo_local_data

    except Exception as e:
        logger.error(f"Error while computing tree : {e}")
        phylo_local_data = {"status": "done"}
        return html.Div(f"❌ Error while computing tree : {e}", style=error_style), phylo_data, phylo_local_data



@app.callback(
    Output("phylogenetic-message", "children", allow_duplicate=True),
    Output("download-tree", "data", allow_duplicate=True),
    Input('btn-save-tree', 'n_clicks'),
    State("phylogenetic-page-store", "data"),
    prevent_initial_call=True
)
def save_tree(n_clicks, phylo_data):
    if not phylo_data or "newick_region" not in phylo_data or not phylo_data["newick_region"]:
        return "No data to save", None

    newick_content = phylo_data["newick_region"]
    if not SERVER_MODE:
        save_path = os.path.join(os.getcwd(), EXPORT_DIR, "tree.nwk")
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        with open(save_path, "w", encoding="utf-8") as f:
            f.write(newick_content)
        return f"File saved : {save_path}", None
    else:
        return "File downloaded.", dcc.send_string(newick_content, "tree.nwk")

    
@app.callback(
    Output('upload-status', 'children'),
    Output('cytoscape-phylo', 'elements'),
    Output('cytoscape-phylo-region', 'elements'),
    Output("phylo-spinner-container", "style"),
    Output("btn-plot-global-tree", "disabled"),
    Output("btn-load-last-tree", "disabled"),
    Output("btn-cancel-plot-global-tree", "disabled"),
    Input('url', 'pathname'),
    Input("phylo-job-status", "data"),
    Input("phylo-local-tree-job-status", "data"),
    State('phylogenetic-page-store', 'data'),
    #prevent_initial_call='initial_duplicate'
    prevent_initial_call=False

)
def update_graph_on_page_load(pathname, status_data, local_tree_status, phylo_data ):
    if pathname != "/phylogenetic":
        raise PreventUpdate
    button_search = False
    button_load = False
    button_cancel = True
    spinner_container = {"display": "none", "marginTop": "20px"}
    if status_data and "status" in status_data and status_data["status"] == "running":
        button_search = True
        button_load = True
        button_cancel = False
        spinner_container = {"display": "block", "marginTop": "20px"}
    elements_region = []
    elements_global = []
    if phylo_data is None:
        phylo_data = {}
    message = phylo_data.get("message", "")
    if "newick_region" in phylo_data and phylo_data["newick_region"] is not None:
        elements_region = generate_elements(phylo_data["newick_region"])
    if "newick_global" in phylo_data:
        elements_global = generate_elements(phylo_data["newick_global"])
    #Display last tree button if a tree has already been computed

    return message, elements_global, elements_region, spinner_container, button_search, button_load, button_cancel
