#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 22:08:19 2025

@author: fgraziani
"""

import base64
import dash
from dash import html, Input, Output, callback, State, callback_context, dcc, ctx, exceptions, no_update, ALL
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from cache_manager import *

import os
import sys
import io
import math

from Bio import Phylo
root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if root_path not in sys.path:
    sys.path.append(root_path)
from app import *
from database.services.neo4j_requests import *
import logging

from pages.phylo_styles import stylesheet

logger = logging.getLogger("panabyss_logger")

EXPORT_DIR = "./export/phylo/"
MAX_NODES_FROM_DB = get_max_nodes_from_db()

def generate_elements(
    newick_str,
    xlen=30,
    ylen=30,
    grabbable=False,
    colors=None,
    root=None,
    branch_depth_ratio=1
):
    if colors is None:
        colors = {}

    tree = Phylo.read(io.StringIO(newick_str), "newick")


    if root == "midpoint":
        tree.root_at_midpoint()

    elif root:
        outgroup = next(
            (
                terminal
                for terminal in tree.get_terminals()
                if terminal.name == root
            ),
            None
        )

        if outgroup is not None:
            tree.root_with_outgroup(outgroup)

    def get_col_positions(tree, column_width=80):
        taxa = tree.get_terminals()

        max_label_width = max(len(str(taxon)) for taxon in taxa)
        drawing_width = (column_width - max_label_width - 1)*branch_depth_ratio

        depths = tree.depths()

        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)

        fudge_margin = int(math.ceil(math.log(len(taxa), 2)))

        cols_per_branch_unit = (
            (drawing_width - fudge_margin) /
            float(max(depths.values()))
        )

        return {
            clade: int(blen * cols_per_branch_unit + 1.0)
            for clade, blen in depths.items()
        }

    def get_row_positions(tree):
        taxa = tree.get_terminals()

        positions = {
            taxon: 2 * idx
            for idx, taxon in enumerate(taxa)
        }

        def calc_row(clade):
            for subclade in clade:
                if subclade not in positions:
                    calc_row(subclade)

            positions[clade] = (
                positions[clade.clades[0]] +
                positions[clade.clades[-1]]
            ) // 2

        calc_row(tree.root)

        return positions

    def add_to_elements(clade, clade_id):
        children = clade.clades

        pos_x = col_positions[clade] * xlen
        pos_y = row_positions[clade] * ylen

        cy_source = {
            "data": {
                "id": clade_id
            },
            "position": {
                "x": pos_x,
                "y": pos_y
            },
            "classes": "nonterminal",
            "grabbable": grabbable
        }

        nodes.append(cy_source)

        if clade.is_terminal():
            cy_source['data']['name'] = clade.name
            cy_source['classes'] = 'terminal'

            if clade.name in colors:
                cy_source['data']['individualColor'] = colors[clade.name]
                cy_source['classes'] = 'colored-terminal'
            else:
                cy_source['classes'] = 'terminal'

        for n, child in enumerate(children):

            support_id = clade_id + "s" + str(n)
            child_id = clade_id + "c" + str(n)

            pos_y_child = row_positions[child] * ylen

            cy_support_node = {
                "data": {
                    "id": support_id
                },
                "position": {
                    "x": pos_x,
                    "y": pos_y_child
                },
                "grabbable": grabbable,
                "classes": "support"
            }

            cy_support_edge = {
                'data': {
                    'source': clade_id,
                    'target': support_id,
                    'sourceCladeId': clade_id,
                    'targetCladeId': child_id,
                    'parentCladeId': clade_id,
                },
                'classes': ''
            }

            cy_edge = {
                'data': {
                    'source': support_id,
                    'target': child_id,
                    'length': clade.branch_length,
                    'sourceCladeId': clade_id,
                    'targetCladeId': child_id,
                    'parentCladeId': clade_id,
                },
                'classes': ''
            }

            if clade.confidence and clade.confidence.value:
                cy_source["data"]["confidence"] = clade.confidence.value

            nodes.append(cy_support_node)
            edges.extend([cy_support_edge, cy_edge])

            add_to_elements(child, child_id)

    col_positions = get_col_positions(tree)
    row_positions = get_row_positions(tree)

    nodes = []
    edges = []

    add_to_elements(tree.clade, "r")

    return nodes + edges


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

###################################################################
@app.callback(
    Output("phylo-job-status", "data", allow_duplicate=True),
    Output("phylo-spinner-container", "style", allow_duplicate=True),
    Output("phylo-job-poller", "disabled", allow_duplicate=True),
    Output("phylogenetic-page-store", "data", allow_duplicate=True),
    Input('upload-newick', 'contents'),
    Input('btn-plot-global-tree', 'n_clicks'),
    State("use-cache-checkbox", "value"),
    State('method-dropdown', 'value'),
    State("phylogenetic_chromosomes_dropdown", 'value'),
    State("phylogenetic-page-store", "data"),
    State('upload-newick', 'filename'),
    State('upload-newick', 'last_modified'),

    prevent_initial_call=True
)
def start_phylo_job(
        contents,
        n_clicks,
        use_cache_checkbox_value,
        method,
        chromosome,
        phylo_data,
    filename,
    last_modified
):

    triggered_id = ctx.triggered_id
    if triggered_id not in ["upload-newick", "btn-plot-global-tree"]:
        raise exceptions.PreventUpdate
    if phylo_data is None:
        phylo_data = {}

    use_cache = True
    if "use_cache" not in use_cache_checkbox_value:
        use_cache = False

    phylo_data.pop("newick_global", None)

    c = chromosome if chromosome else None
    # Upload of a newick file
    if triggered_id == "upload-newick":


        if contents is None or filename is None or last_modified is None:
            raise exceptions.PreventUpdate

        if (phylo_data.get("_upload_filename") == filename and phylo_data.get("_upload_last_modified") == last_modified):
            raise exceptions.PreventUpdate

        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        try:
            newick_str = decoded.decode('utf-8')
        except Exception as e:
            return (
                {"status": "done"},
                {"display": "none"},
                True,
                {**phylo_data, "message": str(e)}
            )

        phylo_data["newick_global"] = newick_str
        return (
            {"status": "done"},
            {"display": "none"},
            True,
            phylo_data
        )
    # Launch the tree computation
    results = compute_global_phylo_tree_from_nodes_wrapper(
        method=method,
        chromosome=c,
        force_reload=(not use_cache)
    )

    if results["status"] == "SUCCESS":
        logger.debug(f"Global tree already exists, load it.")

        phylo_data["newick_global"] = results["newick_tree"]
        return (
            {"status": "done"},
            {"display": "none"},
            True,
            phylo_data
        )

    # Status running => polling + spinner
    logger.debug(f"Computation of global tree already in progress.")
    return (
        {"status": "running", "job_id": results["job_id"]},
        {"display": "block", "marginTop": "20px"},
        False,
        phylo_data
    )


@app.callback(
    Output("phylogenetic-page-store", "data", allow_duplicate=True),
    Output("global-notification", "data", allow_duplicate=True),
    Output("phylo-job-status", "data", allow_duplicate=True),
    Output("phylo-spinner-container", "style", allow_duplicate=True),
    Output("phylo-job-poller", "disabled", allow_duplicate=True),
    Input("phylo-job-poller", "n_intervals"),
    State("phylo-job-status", "data"),
    State("phylogenetic-page-store", "data"),
    prevent_initial_call=True
)
def poll_phylo_job(n, status_data, phylo_data):

    if not status_data or status_data.get("status") != "running":
        raise exceptions.PreventUpdate

    job_id = status_data["job_id"]

    status = get_phylo_status(job_id)

    if status == "RUNNING":
        raise exceptions.PreventUpdate

    if status == "SUCCESS":
        job = get_phylo_job(job_id)
        phylo_data["newick_global"] = job["global_tree"]
        phylo_data["message"] = "Tree successfully computed"
        return (
            phylo_data,
            {"title": "Global tree", "message": "Done", "type": "success"},
            {"status": "done"},
            {"display": "none"},
            True  # stop polling
        )

    if status == "ERROR":
        job = get_phylo_job(job_id)
        phylo_data["message"]=job["error_message"]
        return (
            phylo_data,
            {"title": "Error", "message": job["error_message"], "type": "danger"},
            {"status": "done"},
            {"display": "none"},
            True
        )

    if status == "CANCEL":
        job = get_phylo_job(job_id)
        phylo_data["message"]="Job canceled"
        return (
            phylo_data,
            {"title": "Warning", "message": "Job canceled", "type": "warning"},
            {"status": "done"},
            {"display": "none"},
            True
        )


@app.callback(
    Output('phylo-job-status', 'data', allow_duplicate=True),
    Output("phylo-job-poller", "disabled", allow_duplicate=True),
    Output("phylo-spinner-container", "style", allow_duplicate=True),
    Input('btn-cancel-plot-global-tree', 'n_clicks'),
    State("phylo-job-status", "data"),
    prevent_initial_call=True,
)
def handle_cancel_click(n_clicks, status_data):
    if not n_clicks:
        return no_update
    job_id = None
    if status_data:
        job_id = status_data["job_id"]
    if job_id:
        set_phylo_job_cancel(job_id)
    return (
        {"status": "done"},
        True,  # stop polling
        {"display": "none"}
    )



# @app.callback(
#     Output('cytoscape-phylo','stylesheet',allow_duplicate=True),
#     Input('cytoscape-phylo','selectedEdgeData'),
#     prevent_initial_call=True
# )
# def color_children(selected_edges):
#
#     if not selected_edges:
#         return stylesheet
#
#     # On ne veut qu'un seul lien sélectionné
#     edgeData = selected_edges[-1]
#
#     source = edgeData['source']
#
#     if 's' in source:
#         val = source.split('s')[0]
#     else:
#         val = source
#
#     children_style = {
#         'selector': f'edge[source *= "{val}"]',
#         'style': {
#             'line-color': 'blue'
#         }
#     }
#
#     return stylesheet + [children_style]
from dash import ctx

@app.callback(
    Output('cytoscape-phylo','stylesheet',allow_duplicate=True),
    Output('cytoscape-phylo-region','stylesheet',allow_duplicate=True),
    Input('cytoscape-phylo','selectedEdgeData'),
    Input('cytoscape-phylo-region','selectedEdgeData'),
    prevent_initial_call=True
)
def color_children(selected_edges_global, selected_edges_region):

    triggered_id = ctx.triggered_id

    # ------------------------------------------------------------
    # Global tree
    # ------------------------------------------------------------

    if triggered_id == 'cytoscape-phylo':

        if not selected_edges_global:
            return stylesheet, no_update

        edge_data = selected_edges_global[-1]

        source = edge_data['source']

        if 's' in source:
            val = source.split('s')[0]
        else:
            val = source

        children_style = {
            'selector': f'edge[source *= "{val}"]',
            'style': {
                'line-color': 'blue'
            }
        }

        return (
            stylesheet + [children_style],
            no_update
        )

    # ------------------------------------------------------------
    # Region tree
    # ------------------------------------------------------------

    if triggered_id == 'cytoscape-phylo-region':

        if not selected_edges_region:
            return no_update, stylesheet

        edge_data = selected_edges_region[-1]

        source = edge_data['source']

        if 's' in source:
            val = source.split('s')[0]
        else:
            val = source

        children_style = {
            'selector': f'edge[source *= "{val}"]',
            'style': {
                'line-color': 'blue'
            }
        }

        return (
            no_update,
            stylesheet + [children_style]
        )

    return no_update, no_update



#Callback to plot tree of the current displayed region
@app.callback(
    #Output('cytoscape-phylo-region', 'elements', allow_duplicate=True),
    Output("phylogenetic-message", "children", allow_duplicate=True),
    Output("phylogenetic-page-store", "data", allow_duplicate=True),
    Output("phylo-local-tree-job-status", "data"),
    #Output('shared_storage_nodes', 'data',allow_duplicate=True),
    Input('btn-plot-region', 'n_clicks'),
    State('shared_storage_nodes', 'data'),
    State("phylogenetic-page-store", "data"),
    State("checkbox-weight-node-size", "value"),
    State('home-page-store', 'data'),
    State('global_parameters', 'data'),
    prevent_initial_call=True
)
def plot_region(n_clicks, stored_data,
                phylo_data, weighted_checkbox_value,
                home_data_storage, global_parameters):
    if not n_clicks:
        raise exceptions.PreventUpdate



    ctx = dash.callback_context
    if ctx.triggered_id == "btn-plot-region" and n_clicks == 0:
        raise PreventUpdate
    nodes = no_update
    message = ""
    weighted = False
    if "weight_by_size" in weighted_checkbox_value:
        weighted = True
    phylo_local_data = {"status": "done"}

    if "nodes_cache_id" not in stored_data:
        raise PreventUpdate
    nodes_cache_id = stored_data["nodes_cache_id"]
    cached = get_session_cache(nodes_cache_id)
    min_node_size = cached.get("min_node_size", 10)
    nodes = cached.get("nodes", {})

    if not stored_data or len(nodes) == 0:
        phylo_local_data = {"status": "done"}
        return html.Div(html.P([
        "❌ No data to compute tree. Select a region to visualise on the ",
        dcc.Link("home page", href="/", style={'color': 'blue', 'textDecoration': 'underline'}),
        " or on the ",
        dcc.Link("gwas page", href="/gwas", style={'color': 'blue', 'textDecoration': 'underline'})
        ], style=error_style)), phylo_data, phylo_local_data

    try:
        # Step 1: Check if all nodes are in the region (since min node size can be set greater than 1)
        if min_node_size > 1:
            # Get all the nodes from the region
            max_nodes_from_db = MAX_NODES_FROM_DB
            if global_parameters and "max_nodes_from_db" in global_parameters:
                max_nodes_from_db = global_parameters["max_nodes_from_db"]

            genome = home_data_storage.get("selected_genome", None)
            if "genome_zoom" in home_data_storage and home_data_storage["genome_zoom"]:
                genome = home_data_storage["genome_zoom"]

            chromosome = home_data_storage.get("selected_chromosome", None)
            start = home_data_storage.get("start", None)
            end = home_data_storage.get("end", None)
            logger.debug(f"Phylo tree construction: getting all the nodes for the region chr {chromosome} start {start} end {end} on genome {genome}")
            nodes, return_metadata = get_nodes_by_region(
                genome, chromosome=chromosome, start=start, end=end, use_anchor=False, max_nodes_number=max_nodes_from_db)
            if not return_metadata or not "return_code" in return_metadata or return_metadata["return_code"].lower() not in ["ok", "partial"]:
                if not return_metadata or not "return_code" in return_metadata:
                    message = "Unknown error"
                else:
                    match return_metadata["return_code"].lower():
                        case "filter" | "partial":
                            if "removed_genomes" in return_metadata and len(return_metadata["removed_genomes"]) > 0:
                                message = html.Div(f"Warning: incomplete region for these genomes: {return_metadata['removed_genomes']}", style=warning_style)
                            else:
                                message = "Region too wide"
                                return html.Div(f"❌ {message}", style=error_style), phylo_data, phylo_local_data
                        case "wide" | "zoom":
                            message = "region too wide"
                            return html.Div(f"❌ {message}", style=error_style), phylo_data, phylo_local_data
                        case _:
                            message = "Unknown error"
                            return html.Div(f"❌ {message}", style=error_style), phylo_data, phylo_local_data

            cached["min_node_size"] = 1
            cached["nodes"] = nodes
            nodes_cache.set(nodes_cache_id, cached, expire=8 * 3600)
            logger.debug(f"Number of nodes in the region: {len(nodes)}")
            stored_data = nodes
        else:
            logger.debug("Min node size : 1 => constructing local tree")

        if nodes and len(nodes) > 0:
            # Step 2: compute tree of the region
            newick_str = compute_phylo_tree_from_nodes(nodes, weighted=weighted)
            if phylo_data is None:
                phylo_data = {"newick_region":newick_str}
            else:
                phylo_data["newick_region"] = newick_str

            # Step 3: draw tree
            phylo_local_data = {"status": "done"}
            return message, phylo_data, phylo_local_data
        else:
            phylo_local_data = {"status": "done"}
            return html.Div(html.P([
                "❌ Too much nodes to compute tree. Increase Limit nodes number from database in About Page or Select a smaller region to visualise on the ",
                dcc.Link("home page", href="/", style={'color': 'blue', 'textDecoration': 'underline'}),
                " or on the ",
                dcc.Link("gwas page", href="/gwas", style={'color': 'blue', 'textDecoration': 'underline'})
            ], style=error_style)), phylo_data, phylo_local_data

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
    if not n_clicks:
        raise PreventUpdate
    if not phylo_data or "newick_region" not in phylo_data or not phylo_data["newick_region"]:
        return "No data to save", None

    newick_content = phylo_data["newick_region"]
    file_name = "local_tree_"+datetime.now().strftime("%Y_%m_%d_%H_%M_%S")+".nwk"
    if not SERVER_MODE:
        save_path = os.path.join(os.getcwd(), EXPORT_DIR, file_name)
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        with open(save_path, "w", encoding="utf-8") as f:
            f.write(newick_content)
        return f"File saved : {save_path}", None
    else:
        return "File downloaded.", dcc.send_string(newick_content, file_name)


@app.callback(
    Output('upload-status', 'children', allow_duplicate=True),
    Output("download-tree", "data", allow_duplicate=True),
    Input('btn-save-global-tree', 'n_clicks'),
    State("phylogenetic-page-store", "data"),
    State('method-dropdown', 'value'),
    State("phylogenetic_chromosomes_dropdown", 'value'),
    prevent_initial_call=True
)
def save_global_tree(n_clicks, phylo_data, method, chromosome):
    if not n_clicks:
        raise PreventUpdate
    if not phylo_data or "newick_global" not in phylo_data or not phylo_data["newick_global"]:
        global_newick = get_existing_global_tree(method=method,  chromosome=chromosome)
        if global_newick:
            if not phylo_data:
                phylo_data = {}
            phylo_data["newick_global"] = global_newick
        else:
            return "No data to save", None

    newick_content = phylo_data["newick_global"]
    if not SERVER_MODE:
        save_path = os.path.join(os.getcwd(), EXPORT_DIR, "global_tree_"+datetime.now().strftime("%Y_%m_%d_%H_%M_%S")+".nwk")
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        with open(save_path, "w", encoding="utf-8") as f:
            f.write(newick_content)
        return f"File saved : {save_path}", None
    else:
        return "File downloaded.", dcc.send_string(newick_content, "global_tree_"+datetime.now().strftime("%Y_%m_%d_%H_%M_%S")+".nwk")



"""
This function returns the div color picker option
"""
def build_color_picker(shared_store):
    genomes = []
    color_pickers = no_update
    if shared_store:
        genomes = shared_store.get("genomes", []) or []

        color_pickers = [
            html.Div(
                [
                    dbc.Input(
                        id={
                            'type': 'color-picker',
                            'index': genome
                        },
                        type='color',
                        value="#000000",
                        style={
                            'width': '22px',
                            'height': '22px',
                            'minWidth': '22px',
                            'padding': '0',
                            'border': 'none',
                            'borderRadius': '5px',
                            'overflow': 'hidden',
                            'cursor': 'pointer'
                        }
                    ),

                    html.Label(
                        genome,
                        style={
                            "marginLeft": "5px",
                            "fontSize": "12px",
                            "whiteSpace": "nowrap"
                        }
                    )
                ],

                style={
                    'display': 'flex',
                    'alignItems': 'center',
                    'width': '165px',
                    'padding': '1px 2px',
                    "marginBottom": "5px",
                    "overflow": "visible"
                }
            )

            for i, genome in enumerate(genomes)
        ]

    return color_pickers


"""
Callback launched on url change, data change or option selection
"""
@app.callback(
    Output('upload-status', 'children'),
    Output('cytoscape-phylo', 'elements'),
    Output('cytoscape-phylo-region', 'elements'),
    Output("phylo-spinner-container", "style"),
    Output("btn-plot-global-tree", "disabled"),
    Output("btn-cancel-plot-global-tree", "disabled"),
    Output("individual-color-pickers", "children"),
    Output('root-tree-dropdown', 'options'),
    Input('url', 'pathname'),
    Input("phylo-job-status", "data"),
    Input("phylo-local-tree-job-status", "data"),
    Input({'type': 'color-picker', 'index': ALL},'value'),
    Input('root-tree-dropdown', 'value'),
    State({'type': 'color-picker', 'index': ALL},'id'),
    State('phylogenetic-page-store', 'data'),
    State('method-dropdown', 'value'),
    State("phylogenetic_chromosomes_dropdown", 'value'),
    State("shared_storage", "data"),
    prevent_initial_call=False
)
def update_graph(
    pathname,
    status_data,
    local_tree_status,
    selected_colors,
    root,
    color_picker_ids,
    phylo_data,
    method,
    chromosome,
    shared_store
):

    if pathname != "/phylogenetic":
        raise PreventUpdate

    button_search = False
    button_cancel = True

    spinner_container = {
        "display": "none",
        "marginTop": "20px"
    }

    if (
        status_data
        and "status" in status_data
        and status_data["status"] == "running"
    ):
        button_cancel = False
        spinner_container = {
            "display": "block",
            "marginTop": "20px"
        }

    elements_region = []
    elements_global = []

    # ============================================================
    # COLOR PICKERS
    # ============================================================

    genomes = []
    root_options = no_update
    if shared_store:
        genomes = shared_store.get("genomes", []) or []
        root_options = [
            {
                "label": "No rooting",
                "value": None
            },
            {
                "label": "Midpoint",
                "value": "midpoint"
            }
        ]

        root_options += [
            {
                "label": genome,
                "value": genome
            }
            for genome in genomes
        ]

    existing_genomes = [
        picker_id["index"]
        for picker_id in color_picker_ids
        if isinstance(picker_id, dict)
        and picker_id.get("type") == "color-picker"
    ]

    if existing_genomes != genomes:

        color_pickers = build_color_picker(shared_store)

    else:

        color_pickers = no_update
        root_options = no_update

    # ============================================================
    # BUILD COLOR DICTIONARY
    # ============================================================

    colors = {}

    for picker_id, color in zip(color_picker_ids, selected_colors):

        if (isinstance(picker_id, dict) and picker_id.get("type") == "color-picker"):
            genome = picker_id.get("index")

            if (
                    genome is not None
                    and color
                    and color.lower() != "#000000"
            ):
                colors[genome] = color

    # ============================================================
    # TREE DATA
    # ============================================================
    if phylo_data is None:
        phylo_data = {}

    if (
        "newick_global" not in phylo_data
        and (
            status_data is None
            or "status" not in status_data
            or status_data["status"] != "running"
        )
    ):
        global_newick = get_existing_global_tree(method=method,chromosome=chromosome)

        if global_newick:
            phylo_data["newick_global"] = global_newick

    message = phylo_data.get("message", "")
    if ("newick_region" in phylo_data and phylo_data["newick_region"] is not None):
        elements_region = generate_elements(phylo_data["newick_region"],colors=colors,root=root)

    if "newick_global" in phylo_data:
        elements_global = generate_elements(phylo_data["newick_global"],colors=colors,root=root)

    return (
        message,
        elements_global,
        elements_region,
        spinner_container,
        button_search,
        button_cancel,
        color_pickers,
        root_options
    )
