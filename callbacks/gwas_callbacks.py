#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 22:03:10 2025

@author: fgraziani
"""

from dash import html, Output, Input, State, no_update, dcc, ctx, callback_context, exceptions
import dash_bootstrap_components as dbc

import os
import sys
root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if root_path not in sys.path:
    sys.path.append(root_path)

from app import *
from neo4j_requests import *
import base64
import io
import logging
import plotly.graph_objects as go

logger = logging.getLogger("panabyss_logger")

EXPORT_DIR = "./export/gwas/"

NB_POINTS_WEBGL = 1000

MAX_GWAS_STORE, MAX_RUNNING_INACTIVITY_HOURS, MAX_GWAS_REGIONS = get_gwas_conf()

def compute_gwas_file_name(selected_genomes_list, min_node_size, max_node_size, selected_hap_percent, tolerance_percent,
                           max_gap, deletion, unselected_hap_percent):
    header = ""
    file_name = ""
    return header, file_name


#populates genomes checkboxes
@app.callback(
    Output('genome-list', 'options'),
    Input('shared_storage', 'data')
)
def update_genome_checkboxes(data):
    if not data:
        data  = {}
    if "genomes" not in data:
        data['genomes'] = get_genomes()
    #logger.info("update data : " + str(data))
    return [{'label': genome, 'value': genome} for genome in data['genomes']]


#Populate chromosome droplist
@app.callback(
    Output('gwas_chromosomes_dropdown', 'options'),
    Input('shared_storage', 'data')
)
def update_dropdown(data):
    if not data:
        data  = {}
    if "chromosomes" not in data:
        data["chromosomes"] = get_chromosomes()

    chromosomes = data["chromosomes"]
    options = [{"label": "All", "value": "All"}] + [
        {"label": str(chrom), "value": str(chrom)} for chrom in chromosomes
    ]
    return options

#Populate genome droplist
# @app.callback(
#     Output('gwas_ref_genome_dropdown', 'options'),
#     Input('shared_storage', 'data')
# )
# def update_dropdown(data):
#     if not data or "genomes" not in data:
#         return []
#
#     genomes = data["genomes"]
#     options =  [{"label": str(g), "value": str(g)} for g in genomes]
#     return options


#Sort chromosome by name (ex : chr1, chr2, chr3, etc.)
def chromosome_sort_key(chrom):

    chrom_str = str(chrom)

    if chrom_str.isdigit():
        return (0, int(chrom_str))

    match = re.search(r'(\d+)$', chrom_str)
    if match:
        return (1, int(match.group(1)))

    return (2, chrom_str.lower())

#Function used to create the chromosome plot
def build_chromosome_figure(data):
    chromosome_stats = get_chromosomes_stats()
    if not data:
        return go.Figure()

    fig = go.Figure()

    #Check empty results for a chromosome
    chrom_data = {
        chrom: points
        for chrom, points in data.items()
        if points
    }

    if not chrom_data:
        return go.Figure()

    #Sort chromosome
    chromosomes = sorted(chrom_data.keys(), key=chromosome_sort_key)

    chrom_max_lengths = {}
    #Global length for each chromosome
    for chrom in chromosomes:
        if chromosome_stats and chromosome_stats.get(f"{chrom}_max_position_mean") is not None:
            chrom_max_lengths[chrom] = chromosome_stats[f"{chrom}_max_position_mean"]
        else:
            chrom_max_lengths[chrom] = max(p[0] for p in chrom_data[chrom])
    uniform_chrom_length = max(chrom_max_lengths.values())

    all_y = [p[1] for points in chrom_data.values() for p in points]
    y_min, y_max = min(all_y), max(all_y)
    y_pad = 0.1 * max(abs(y_min), abs(y_max), 1)

    colors = ["rgba(200,200,200,0.25)", "rgba(150,150,255,0.25)"]

    x_offset = 0
    chromosome_centers = []
    chromosome_labels = []
    total_points = sum(len(points) for points in chrom_data.values())
    use_webgl = total_points > NB_POINTS_WEBGL
    ScatterClass = go.Scattergl if use_webgl else go.Scatter
    if use_webgl :
        logger.debug("WebGL Scatter class used.")
    else:
        logger.debug("SVG Scatter class used.")
    for i, chrom in enumerate(chromosomes):
        #Sort by genomic coordinates
        points_sorted = sorted(chrom_data[chrom], key=lambda p: p[0])

        x_local = [p[0] for p in points_sorted]
        y = [p[1] for p in points_sorted]

        real_max = chrom_max_lengths[chrom]

        x_start = x_offset
        x_end = x_offset + real_max

        #global coordinates
        x_global = [x + x_offset for x in x_local]

        #chromosome area
        fig.add_shape(
            type="rect",
            x0=x_start,
            x1=x_end,
            y0=y_min - y_pad,
            y1=y_max + y_pad,
            fillcolor=colors[i % len(colors)],
            line_width=0,
            layer="below"
        )

        #Line on beginning ending of chromosome
        for x in (x_start, x_end):
            fig.add_shape(
                type="line",
                x0=x, x1=x,
                y0=y_min - y_pad,
                y1=y_max + y_pad,
                line=dict(color="black", width=1)
            )

        #Points
        fig.add_trace(ScatterClass(
            x=x_global,
            y=y,
            mode="markers",
            name=str(chrom),
            customdata=x_local,
            hovertemplate=(
                f"{chrom}<br>"
                "Position: %{customdata}<br>"
                "Value: %{y}<extra></extra>"
            )
        ))

        chromosome_centers.append(x_offset + real_max / 2)
        chromosome_labels.append(str(chrom))

        x_offset += real_max + 1

    fig.add_hline(
        y=0,
        line_dash="dash",
        line_color="black"
    )

    fig.update_layout(
        title=dict(
            text="Distribution of Shared Regions",
            x=0.5,
            xanchor="center",
            y=0.95,
            yanchor="top"
        ),
        xaxis=dict(
            tickmode="array",
            tickvals=chromosome_centers,
            ticktext=chromosome_labels,
            title="Chromosomes and mean position."
        ),
        yaxis=dict(
            title="Size of the shared region (negative = deletions)",
            range=[y_min - y_pad, y_max + y_pad]
        ),
        showlegend=False,
        plot_bgcolor="white",
        margin=dict(l=40, r=20, t=90, b=40)
    )

    return fig

#Callback to get selected genomes and put them in the genome dropdown
@app.callback(
    Output("gwas_ref_genome_dropdown", "options", allow_duplicate=True),
    Output("gwas_ref_genome_dropdown", "value", allow_duplicate=True),
    Input("genome-list", "value"),
    State("gwas_ref_genome_dropdown", "value"),
    State('parameters-gwas-page-store', 'data'),
    prevent_initial_call=True,
)
def update_ref_genome_dropdown(selected_genomes, current_value, parameters_data):
    if not selected_genomes:
        return [], None

    options = [{"label": g, "value": g} for g in selected_genomes]
    if "ref_genome" in parameters_data:
        value = parameters_data["ref_genome"]
    elif current_value in selected_genomes:
        value = current_value
    else:
        value = None

    return options, value


#Callback triggered by clicking on the launch button
#This callback is required to set the parameters in the store to allow further navigation
@app.callback(
    Output('shared-status', 'children', allow_duplicate=True),
    Output('parameters-gwas-page-store', 'data', allow_duplicate=True),
    Output("gwas-page-store", "data", allow_duplicate=True),
    Output("gwas-poll-interval", "disabled", allow_duplicate=True),
    Output("btn-find-shared", "disabled", allow_duplicate=True),
    Output("btn-cancel-find-shared", "disabled", allow_duplicate=True),
    Input('btn-find-shared', 'n_clicks'),
    State('genome-list', 'value'),
    State("parameters-gwas-page-store", "data"),
    State("gwas-min-node-size-int", 'value'),
    State("gwas-max-node-size-int", 'value'),
    State("gwas-min-percent_selected", 'value'),
    State("tolerance_percentage", 'value'),
    State("gwas-region-gap", 'value'),
    State('gwas-toggle-deletion', 'value'),
    State("gwas_chromosomes_dropdown", 'value'),
    State("gwas_ref_genome_dropdown", 'value'),
    State("deletion-percentage", 'value'),
    State("gwas-page-store", "data"),
    prevent_initial_call=True,

)
def handle_shared_region_search_click(n_clicks, selected_genomes, data, min_node_size, max_node_size,
                                min_percent_selected, tolerance_percentage, region_gap,
                                deletion_checkbox, chromosome, ref_genome, deletion_percentage, gwas_data):
    poll_enabled = True
    find_button = False
    cancel_button = True
    if not n_clicks:
        return no_update, no_update, no_update, no_update, no_update, no_update
    min_size = 10
    if min_node_size is not None and min_node_size != "" and isinstance(min_node_size, int):
        min_size = min_node_size
    if min_node_size is None:
        min_size = 0
    data["min_node_size"] = min_size
    if chromosome == None or chromosome == "All":
        c = None
    else:
        c = [chromosome]
    if data is None:
        data = {}
    if gwas_data is None:
        gwas_data = {}
    data["chromosomes"]=c
    if max_node_size is None or max_node_size == "" or max_node_size == 0:
        max_node_size = 0
    data["max_node_size"] = max_node_size
    data["checkboxes"] = selected_genomes

    if min_percent_selected is not None:
        data["min_percent_selected"] = min_percent_selected
    if tolerance_percentage is not None:
        data["tolerance_percentage"] = tolerance_percentage
    if region_gap is not None:
        data["region_gap"] = region_gap
    if deletion_checkbox is not None:
        data["deletion_checkbox"] = deletion_checkbox
    if deletion_percentage is not None:
        data["deletion_percentage"] = deletion_percentage
    if ref_genome is not None:
        data["ref_genome"] = ref_genome
    if not selected_genomes:
        return "❌ Choose at least one genome.", no_update, no_update, poll_enabled, find_button, cancel_button
    data["launch_ts"] = time.time()
    gwas_data.update({
        "analyse": [],
        "gwas_graph_points": {},
        "message": "",
        "status": "running"
    })
    deletion = 'show' in deletion_checkbox if deletion_checkbox else False

    params = {
        "genomes_list": selected_genomes,
        "genome_ref": ref_genome,
        "chromosomes": c,
        "node_min_size": min_node_size,
        "node_max_size": max_node_size,
        "nodes_max_gap": region_gap,
        "deletion": deletion,
        "min_percent_selected_genomes": min_percent_selected,
        "tolerance_percentage": tolerance_percentage,
        "min_deletion_percentage": deletion_percentage
    }

    job_id = submit_job_gwas(params)
    gwas_data["job_id"] = job_id
    # Enable polling
    poll_enabled = False if job_id else True
    find_button = True
    cancel_button = False
    return "Processing...", data, gwas_data, poll_enabled, find_button, cancel_button


#Callback to poll the current process
@app.callback(
    Output('shared-status', 'children', allow_duplicate=True),
    Output("gwas-page-store", "data"),
    Output("global-notification", "data"),
    Output("gwas-progress-circle", "style", allow_duplicate=True),
    Output("gwas-progress-text", "children", allow_duplicate=True),
    Output("gwas-poll-interval", "disabled", allow_duplicate=True),
    Output("btn-find-shared", "disabled", allow_duplicate=True),
    Output("btn-cancel-find-shared", "disabled", allow_duplicate=True),
    Input("gwas-poll-interval", "n_intervals"),
    State("parameters-gwas-page-store", "data"),
    State("gwas-page-store", "data"),
    Input("gwas-poll-interval", "disabled"),
    prevent_initial_call=True,
    suppress_callback_exceptions=True
)
def poll_gwas_job(n_intervals, parameters_data, gwas_data, poll_disabled):
    if not parameters_data or "job_id" not in gwas_data or poll_disabled:
        raise exceptions.PreventUpdate
    msg = "Processing..."
    job_id = gwas_data["job_id"]
    enable_poll = False
    search_button = True
    cancel_button = False
    progress_style = {"display": "none"}
    progress_value = "0"
    progress_background = {"background": "conic-gradient(#e0e0e0 0deg 360deg)"}
    #job_data = get_job(job_id)  # retrieve job info from SQLite
    status = get_gwas_status(job_id)
    if not status:
        enable_poll = True
        search_button = False
        cancel_button = True
        return "", no_update, no_update, progress_style, progress_value, enable_poll, search_button, cancel_button
        return (no_update,) * 5, True,

    if status == "SUCCESS":
        job_data = get_gwas_job(job_id)
        #Job finished successfully => fetch results
        progress_value = "100"
        msg = ""
        enable_poll = True
        search_button = False
        cancel_button = True
        analyse = job_data.get("result_gwas_regions", "{}")
        gwas_points = job_data.get("result_gwas_points", "{}")
        job_params = job_data.get("params", {})
        selected_genomes = job_params.get("genomes_list", [])
        ref_genome = job_params.get("genome_ref", None)
        if ref_genome is None or ref_genome == "":
            ref_genome = selected_genomes[0]
        analyse_to_plot = analyse[ref_genome]

        # logger.info("analyse to plot : " + str(analyse_to_plot))

        for r in range(len(analyse_to_plot)):
            set_gene_name = set()

            # Get all the genes
            for annot in analyse_to_plot[r].get("annotations", []):
                gene = annot.get("gene_name")
                if gene:
                    set_gene_name.add(gene)

            if set_gene_name:
                genes = sorted(set_gene_name)
                full_text = " ".join(genes)
                annotation = full_text

                # short visualization (60 chars)
                # short_text = full_text[:60]
                # if len(full_text) > 60:
                #     short_text += "..."
                # annotation = f'<details><summary style="cursor:pointer; white-space: pre-wrap;">{short_text}</summary><div style="white-space: pre-wrap;">{full_text}</div></details>'
            else:
                annotation = ""

            analyse_to_plot[r]["annotations"] = annotation

            annot_before = ""
            annot_data = analyse_to_plot[r].get("annotation_before")
            if annot_data:
                gene_name = annot_data.get("gene_name")
                distance = annot_data.get("distance")
                annot_before = f"Gene_name: {gene_name}<br>Distance: {distance}"
            analyse_to_plot[r]["annotation_before"] = annot_before

            annot_after = ""
            annot_data = analyse_to_plot[r].get("annotation_after")
            if annot_data:
                gene_name = annot_data.get("gene_name")
                distance = annot_data.get("distance")

                annot_after = f"Gene_name: {gene_name}<br>Distance: {distance}"
            analyse_to_plot[r]["annotation_after"] = annot_after

        message = f"{len(analyse_to_plot)} shared regions found."
        if MAX_GWAS_REGIONS is not None and MAX_GWAS_REGIONS > 0 and len(analyse_to_plot) >= MAX_GWAS_REGIONS:
            message = f"Too much regions found, limited to the {MAX_GWAS_REGIONS} first regions."

        gwas_data.update({
            "analyse": analyse_to_plot,
            "gwas_graph_points": gwas_points,
            "message": message,
            "status": "done"
        })
        toast = {"title": "Shared regions discovery", "message": "Process terminated.", "type": "success"}
        return msg, gwas_data, toast, progress_style, progress_value, enable_poll, search_button, cancel_button

    elif status == "ERROR":
        job_data = get_gwas_job(job_id)
        enable_poll = True
        msg = f"❌ Error: {job_data.get('error_message', 'Unknown')}"
        #Job finished with error
        gwas_data.update({
            "analyse": [],
            "gwas_graph_points": {},
            "message": f"❌ Error: {job_data.get('error_message', 'Unknown')}",
            "status": "done"
        })
        toast = {"title": "Shared region discovery error", "message": job_data.get("error_message", "Unknown error"), "type": "danger"}
        return msg, gwas_data, toast,  progress_style, progress_value, enable_poll, search_button, cancel_button

    elif status == "RUNNING":
        current_progress = int(get_gwas_progress(job_id))

        deg = current_progress * 3.6
        progress_style = {
            "display": "flex",
            "align-items": "center",
            "justify-content": "center",
            "width": "100px",
            "height": "100px",
            "border-radius": "50%",
            "background": f"conic-gradient(#9d4edd 0deg {deg}deg, #e9ecef {deg}deg 360deg)",
            "margin": "auto",
            "position": "relative"
        }

        progress_value = f"{current_progress}%"

        return msg, no_update, no_update, progress_style, progress_value, enable_poll, search_button, cancel_button
    else:
        enable_poll = True
        msg = f"❌ Error: status {status} unknown"
        # Job finished with error
        gwas_data.update({
            "analyse": [],
            "gwas_graph_points": {},
            "message": f"❌ Error: status {status} unknown",
            "status": "done"
        })
        toast = {"title": "Shared region discovery error", "message": "❌ Error: status "+str(status)+" unknown",
                 "type": "danger"}
        return msg, gwas_data, toast, progress_style, progress_value, enable_poll, search_button, cancel_button


#Cancel clic callback
@app.callback(
    Output('shared-status', 'children', allow_duplicate=True),
    Output("gwas-page-store", "data", allow_duplicate=True),
    Output("global-notification", "data", allow_duplicate=True),
    #Output("spinner-container", "style", allow_duplicate=True),
    Output("gwas-progress-circle", "style", allow_duplicate=True),
    Output("gwas-progress-text", "children", allow_duplicate=True),
    Output("gwas-poll-interval", "disabled", allow_duplicate=True),
    Output("btn-find-shared", "disabled", allow_duplicate=True),
    Output("btn-cancel-find-shared", "disabled", allow_duplicate=True),

    Input('btn-cancel-find-shared', 'n_clicks'),
    State("gwas-page-store", "data"),
    prevent_initial_call=True,

)
def handle_cancel_click(n_clicks, gwas_data):
    if not n_clicks or "job_id" not in gwas_data:
        return (no_update,) * 8

    progress_value = "100"
    progress_style = {"display": "none"}
    logger.debug("Cancel job")
    job_id = gwas_data["job_id"]
    set_gwas_job_cancel(job_id)
    gwas_data["job_id"] = None
    gwas_data.update({
        "analyse": [],
        "gwas_graph_points": {},
        "message": "❌ Job canceled by user.",
        "status": "canceled"
    })
    return ("❌ Job canceled by user.", gwas_data,
            {"title": "Job canceled", "message": "User stopped the job", "type": "warning"},
            progress_style, progress_value, True, False, True)



#Callback to load the graph of the selected region
@app.callback(
    Output('selected-region-output', 'children', allow_duplicate=True),
    Output('shared_storage_nodes', 'data',allow_duplicate=True),
    Output("url", "pathname",allow_duplicate=True),
    Output("home-page-store", "data", allow_duplicate=True),
    #Output("load_spinner_zone", "children", allow_duplicate=True),
    Output('tabs-navigation', 'value'),
    Input('shared-region-table', 'selected_rows'),
    State('shared-region-table', 'data'),
    State('shared_storage_nodes', 'data'),
    State('home-page-store', 'data'),
    running=[
            (Output("spinner-container","style"),{"display": "block", "marginTop": "20px"},{"display": "none", "marginTop": "20px"}),
        ],
    prevent_initial_call=True
)
def handle_row_selection(selected_rows, table_data, data, home_page_data):
    redirect = "/gwas"
    if home_page_data is None:
        home_page_data = {}
    if not selected_rows:
        return no_update, data, redirect, home_page_data, redirect
    logger.debug(table_data[selected_rows[0]])
    row = table_data[selected_rows[0]]
    logger.debug("selected row to plot : " +str(row))
    start = int(row['start'])
    stop = int(row['stop'])
    chromosome = row['chromosome']
    genome = row['genome']
    logger.debug("search region genome " +str(genome) + " chromosome " + str(chromosome) + " start " + str(start) + " stop " + str(stop))
    try:
        nodes, return_metadata = get_nodes_by_region(genome, str(chromosome), start, stop)
        home_page_data["selected_genome"]=genome
        home_page_data["selected_chromosome"]=chromosome
        home_page_data["start"]=start
        home_page_data["end"]=stop
        home_page_data["gene_name"] = ""
        home_page_data["gene_id"] = ""
        home_page_data["search_return_metadata"] = return_metadata
        redirect = "/"
        return html.Div([
            html.P(f"Found nodes into the region : {len(nodes)}")
        ]), nodes,redirect,home_page_data,redirect
    except Exception as e:
        return f"Erreur : {e}", data,redirect,home_page_data,redirect
    
#Update data when navigating or when the process is terminated
@app.callback(
    Output('shared-status', 'children',allow_duplicate=True),
    Output('shared-region-table', 'data',allow_duplicate=True),
    Output("genome-list", "value"),
    Output("gwas-min-node-size-int", 'value'),
    Output("gwas-max-node-size-int", 'value'),
    Output("gwas-min-percent_selected", 'value'),
    Output("tolerance_percentage", 'value'),
    Output("gwas-region-gap", 'value'),
    Output('gwas-toggle-deletion', 'value'),
    Output("deletion-percentage", 'value'),
    Output("chromosome-graph", "figure"),
    Output("chromosome-graph", "style"),
    Output("gwas_chromosomes_dropdown", 'value'),
    Output("gwas_ref_genome_dropdown", 'value'),
    Output("btn-find-shared","disabled", allow_duplicate=True),
    Output("btn-cancel-find-shared","disabled", allow_duplicate=True),
    Output("gwas-poll-interval", "disabled", allow_duplicate=True),
    Input('url', 'pathname'),
    #Input("gwas-page-store", "modified_timestamp"),
    Input("gwas-page-store", "data"),
    State('parameters-gwas-page-store', 'data'),
    prevent_initial_call=True
)
def update_data(path, data, parameters_data):
    #analyse = table_data
    if path != "/gwas":
        raise exceptions.PreventUpdate
    #logger.debug(f"#####################################GWAS Update data")
    analyse = []
    len_analyse = 0
    checkbox = []
    max_node_size = None
    min_node_size = 10
    min_percent_selected = 100
    tolerance_percentage = 0
    region_gap = 10000
    deletion_checkbox = ['show']
    deletion_percentage = 100
    chromosome_figure = {}
    figure_display = {"display": "None"}
    progress_style = {"display": "none"}
    message = ""
    message_analyse = ""
    ref_genome = None
    chromosome = None
    search_button = False
    cancel_button = True
    enable_poll = no_update

    if data is not None:
        if "analyse" in data:
            analyse = data["analyse"]

        job_id = data.get("job_id", None)
        if job_id:
            job_data = get_gwas_job(job_id)  # retrieve job info from SQLite
            if job_data and job_data["status"] == "RUNNING":
                message_analyse = "Processing..."
                search_button = True
                cancel_button = False
                enable_poll = False
            else:
                if "message" in data:
                    message_analyse = data["message"]
                else:
                    if "analyse" in data:
                        len_analyse = len(analyse)
                        if MAX_GWAS_REGIONS is not None and MAX_GWAS_REGIONS > 0 and len_analyse >= MAX_GWAS_REGIONS:
                            message_analyse = f"Too much regions found, limited to the {MAX_GWAS_REGIONS} first regions."
                        else:
                            message_analyse = f"{len_analyse} shared regions found."
        # if data.get("status") == "running":
        #     progress_style = {"display": "block", "marginTop": "20px", "marginBottom": "20px"}
        #     message_analyse = "Progressing..."
        #     search_button = True
        #     cancel_button = False
        else:
            if "analyse" in data:
                if "message" in data :
                    message_analyse = data["message"]
                else:
                    len_analyse = len(analyse)
                    if MAX_GWAS_REGIONS is not None and MAX_GWAS_REGIONS > 0 and len_analyse >= MAX_GWAS_REGIONS:
                        message_analyse = f"Too much regions found, limited to the {MAX_GWAS_REGIONS} first regions."
                    else:
                        message_analyse = f"{len_analyse} shared regions found."
        if "gwas_graph_points" in data and len(data["gwas_graph_points"]) > 0:
            chromosome_figure = build_chromosome_figure(data["gwas_graph_points"])
            figure_display = {"display": "block"}
    if parameters_data is not None:
        if "checkboxes" in parameters_data:
            checkbox = parameters_data["checkboxes"]
        if "min_node_size" in parameters_data:
            min_node_size = parameters_data["min_node_size"]
        if "max_node_size" in parameters_data:
            if parameters_data["max_node_size"] == 0:
                max_node_size = None
            else:
                max_node_size = parameters_data["max_node_size"]
        if "min_percent_selected" in parameters_data:
            min_percent_selected = parameters_data["min_percent_selected"]
        if "tolerance_percentage" in parameters_data:
            tolerance_percentage = parameters_data["tolerance_percentage"]
        if "region_gap" in parameters_data:
            region_gap = parameters_data["region_gap"]
        if "deletion_checkbox" in parameters_data:
            deletion_checkbox = parameters_data["deletion_checkbox"]
        if "deletion_percentage" in parameters_data:
            deletion_percentage = parameters_data["deletion_percentage"]
        if "ref_genome" in parameters_data:
            ref_genome = parameters_data["ref_genome"]
        if "chromosomes" in parameters_data:
            stored = parameters_data["chromosomes"]
            if isinstance(stored, list) and len(stored) > 0:
                chromosome = stored[0]
            else:
                chromosome = stored
        if analyse is not None:
            for i, row in enumerate(analyse):
                row['get_sequence'] = "Get sequence"

    return (message_analyse,analyse, checkbox, min_node_size, max_node_size,
            min_percent_selected,tolerance_percentage,region_gap, deletion_checkbox, deletion_percentage,
            chromosome_figure, figure_display, chromosome, ref_genome,search_button,cancel_button, enable_poll)

#Callback to save the gwas data table into csv file
@app.callback(
    Output('save-feedback', 'children'),
    #Output("load_spinner_zone", "children", allow_duplicate=True),
    Output("download-csv", "data", allow_duplicate=True),
    Input('save-csv-button', 'n_clicks'),
    Input('save-csv-with_seq-button', 'n_clicks'),
    State('shared-region-table', 'data'),
    State('parameters-gwas-page-store', 'data'),
    State("gwas-page-store", "data"),
    running=[
            (Output("spinner-container","style"),{"display": "block", "marginTop": "20px"},{"display": "none", "marginTop": "20px"}),
        ],
    prevent_initial_call=True
)
def save_csv(n_clicks, n_clicks_seq, table_data, parameters_data, gwas_data):
    #logger.info(f"Callback triggered: n_clicks={n_clicks}, table_data={table_data}")
    if not n_clicks and not n_clicks_seq:
        return (no_update,) * 2
    if len(table_data) > 0:
        triggered_id = ctx.triggered_id
        export_sequences = False
        if triggered_id == 'save-csv-with_seq-button':
            export_sequences = True
        if not table_data:
            return "No data.",no_update
        params = parameters_data if parameters_data else {}
        params_str = json.dumps(params)
        df = pd.DataFrame(table_data).drop(columns=['get_sequence'], errors='ignore')
        file_name = "shared_regions_"+datetime.now().strftime("%Y_%m_%d_%H_%M_%S")+".csv"
        if export_sequences :
            sequences = []
            for row in tqdm(table_data):
                sequences.append(get_sequence_from_position(row['genome'], row['chromosome'], row['start'], row['stop']))
            df["sequence"] = sequences
        if not SERVER_MODE:
            save_path = os.path.join(os.getcwd(), EXPORT_DIR, file_name)
            logger.info("save path : " + str(save_path))
            with open(save_path, "w") as f:
                f.write("#PARAMETERS=" + params_str + "\n")
                df.to_csv(f, index=False)
            #df.to_csv(save_path, index=False)
            return f"File saved : {save_path}",no_update
        else:
            logger.info("🌐 Server mode active — file will be downloaded by user.")
            csv_string = "#PARAMETERS=" + params_str + "\n" + df.to_csv(index=False)
            return "File ready for download.", dict(content=csv_string, filename=file_name)
            #return "File ready for download.", dcc.send_data_frame(df.to_csv, "shared_regions.csv", index=False)
    else:
        return "No data to download.", no_update

#Callback to loads csv file into data table
@app.callback(
    Output('shared-status', 'children',allow_duplicate=True),
    Output('shared-region-table', 'data',allow_duplicate=True),
    Output("gwas-page-store", "data", allow_duplicate=True),
    Output("parameters-gwas-page-store", "data", allow_duplicate=True),
    Input('upload-csv', 'contents'),
    State('upload-csv', 'filename'),
    State("gwas-page-store", "data"),
    prevent_initial_call=True
)
def load_csv(contents, filename, gwas_page_store):
    if not contents:
        raise exceptions.PreventUpdate
    if gwas_page_store is None:
        gwas_page_store = {}

    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        content = decoded.decode('utf-8')

        # Extract parameters
        params = None
        lines = content.splitlines()
        if lines and lines[0].startswith("#PARAMETERS="):
            params_line = lines[0]
            try:
                params = json.loads(params_line.replace("#PARAMETERS=", ""))
            except Exception as e:
                logger.info(f"Failed to parse parameters line: {e}")
            content = "\n".join(lines[1:])

        df = pd.read_csv(io.StringIO(content))
        df = df.fillna("")
        analyse = df[[c for c in df.columns if c != "sequence"]].to_dict('records')
        for row in analyse:
            row['get_sequence'] = "Get sequence"

        gwas_page_store["analyse"] = analyse

        if params:
            parameters_store = params
        else:
            parameters_store = no_update
        logger.info("csv file loaded")
        return f"{len(analyse)} shared regions found.", analyse, gwas_page_store, parameters_store

    except Exception as e:
        logger.info(f"Error while loading file: {e}")
        return None, None, gwas_page_store, no_update


    

@app.callback(
    Output('upload-csv', 'style'),
    Input('load-csv-button', 'n_clicks'),
    prevent_initial_call=True
)
def show_upload_area(n_clicks):
    if not n_clicks:
        return no_update
    return {
        'display': 'block',
        'borderWidth': '1px',
        'borderStyle': 'dashed',
        'padding': '10px',
        'marginTop': '10px'
    }


#Callback to display sequence on click on the get sequence column
@app.callback(
    Output('sequence-zone', 'children'),
    Input("shared-region-table", "active_cell"),
    State('shared-region-table', 'data'),
    prevent_initial_call=True,
)
def display_sequence_on_button_click(active_cell, table_data):
    if active_cell and active_cell['column_id'] == 'get_sequence':
        row_index = active_cell["row"]
        row = table_data[row_index]
        sequence = get_sequence_from_position(row['genome'], row['chromosome'], row['start'], row['stop'])
        return html.Div([
            html.P(f"Sequence : {sequence}")
        ])
    return None



