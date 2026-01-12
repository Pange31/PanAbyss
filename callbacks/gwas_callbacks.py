#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 22:03:10 2025

@author: fgraziani
"""

from dash import html, Output, Input, State, no_update, dcc, ctx, callback_context


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


logger = logging.getLogger("panabyss_logger")

EXPORT_DIR = "./export/gwas/"

#populates genomes checkboxes
@app.callback(
    Output('genome-list', 'options'),
    Input('shared_storage', 'data')
)
def update_genome_checkboxes(data):
    
    if not data or "genomes" not in data:
        return []
    #logger.info("update data : " + str(data))
    return [{'label': genome, 'value': genome} for genome in data['genomes']]


#Populate chromosome droplist
@app.callback(
    Output('gwas_chromosomes_dropdown', 'options'),
    Input('shared_storage', 'data')
)
def update_dropdown(data):
    if not data or "chromosomes" not in data:
        return []
    
    chromosomes = data["chromosomes"]
    options = [{"label": "All", "value": "All"}] + [
        {"label": str(chrom), "value": str(chrom)} for chrom in chromosomes
    ]
    return options

#Populate genome droplist
@app.callback(
    Output('gwas_ref_genome_dropdown', 'options'),
    Input('shared_storage', 'data')
)
def update_dropdown(data):
    if not data or "genomes" not in data:
        return []
    
    genomes = data["genomes"]
    options =  [{"label": str(g), "value": str(g)} for g in genomes]
    return options


@app.callback(
    Output('shared-status', 'children'),
    Output("gwas-page-store", "data", allow_duplicate=True),
    Output("load_spinner_zone", "children", allow_duplicate=True),
    Input('btn-find-shared', 'n_clicks'),
    State('genome-list', 'value'),
    State("gwas-page-store", "data"),
    State("gwas-min-node-size-int", 'value'),
    State("gwas-max-node-size-int", 'value'),
    State("gwas-min-percent_selected", 'value'),
    State("tolerance_percentage", 'value'),
    State("gwas-region-gap", 'value'),
    State('gwas-toggle-deletion', 'value'),
    State("gwas_chromosomes_dropdown", 'value'),
    State("gwas_ref_genome_dropdown", 'value'),
    State("deletion-percentage", 'value'),
    prevent_initial_call=True
)
def handle_shared_region_search(n_clicks, selected_genomes, data, min_node_size, max_node_size,
                                min_percent_selected, tolerance_percentage, region_gap,
                                deletion_checkbox, chromosome, ref_genome, deletion_percentage):
    if min_node_size is not None and min_node_size != "" and isinstance(min_node_size, int):
        min_size = min_node_size
    else:
        min_size = 10
    if chromosome == None or chromosome == "All" :
        c = None
    else:
        c = [chromosome]
    if data is None:
        data = {}
    if max_node_size is None or max_node_size == "" or max_node_size == 0:
        max_node_size = 0
        data["max_node_size"] = None
    else:
        data["max_node_size"] = max_node_size

    data["checkboxes"]= selected_genomes
    if min_node_size is not None:
        data["min_node_size"] = min_node_size
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
    if not selected_genomes:
        return "❌ Choose at least one genome.",no_update, ""
    deletion = False
    if 'show' in deletion_checkbox : 
        deletion = True
    try:
        #take an annotated genome if no reference genome selected

        dic_region, analyse = find_shared_regions(selected_genomes, genome_ref = ref_genome, chromosomes = c,
                                                  node_min_size = min_size, node_max_size = max_node_size,
                                                  nodes_max_gap=region_gap, deletion = deletion,
                                                  min_percent_selected_genomes=min_percent_selected,
                                                  tolerance_percentage = tolerance_percentage,
                                                  min_deletion_percentage=deletion_percentage)
        if ref_genome is None or ref_genome == "":
            ref_genome = selected_genomes[0]
        analyse_to_plot = analyse[ref_genome]
                
        #logger.info("analyse to plot : " + str(analyse_to_plot))

        for r in range(len(analyse_to_plot)):
            annotation = ""
            set_gene_name = set()
            if len(analyse_to_plot[r]["annotations"]) > 0:
                for annot in analyse_to_plot[r]["annotations"]:
                    if "gene_name" in annot and annot["gene_name"] is not None:
                        set_gene_name.add(annot["gene_name"])
            if set_gene_name is not None and len(list(set_gene_name)) > 0:
                for gene in list(set_gene_name):
                    annotation += gene + "\n"
            analyse_to_plot[r]["annotations"] = annotation
            annot_before = ""
            if "annotation_before" in analyse_to_plot[r] and "gene_name" in analyse_to_plot[r]["annotation_before"]:
                annot_before = "gene_name : " + analyse_to_plot[r]["annotation_before"]["gene_name"] \
                                +"\nDistance : " + str(analyse_to_plot[r]["annotation_before"]["distance"])

            analyse_to_plot[r]["annotation_before"] = annot_before
            annot_after = ""
            if "annotation_after" in analyse_to_plot[r] and "gene_name" in analyse_to_plot[r]["annotation_after"]:
                annot_after = "gene_name : " + analyse_to_plot[r]["annotation_after"]["gene_name"] \
                                +"\nDistance : " + str(analyse_to_plot[r]["annotation_after"]["distance"])
            
            analyse_to_plot[r]["annotation_after"] = annot_after

        data["analyse"] = analyse_to_plot
        return f"{len(analyse_to_plot)} shared regions found.",data, ""
    
    except Exception as e:
        return f"❌ Error : {e}",no_update, ""
    
@app.callback(
    Output('selected-region-output', 'children', allow_duplicate=True),
    Output('shared_storage_nodes', 'data',allow_duplicate=True),
    Output("url", "pathname",allow_duplicate=True),
    Output("home-page-store", "data", allow_duplicate=True),
    Output("load_spinner_zone", "children", allow_duplicate=True),
    Output('tabs-navigation', 'value'),
    Input('shared-region-table', 'selected_rows'),
    State('shared-region-table', 'data'),
    State('shared_storage_nodes', 'data'),
    State('home-page-store', 'data'),
    prevent_initial_call=True
)
def handle_row_selection(selected_rows, table_data, data, home_page_data):
    redirect = "/gwas"
    if home_page_data is None:
        home_page_data = {}
    if not selected_rows:
        return no_update, data, redirect, home_page_data, "",redirect
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
        ]), nodes,redirect,home_page_data,"",redirect
    except Exception as e:
        return f"Erreur : {e}", data,redirect,home_page_data,"",redirect
    
#Restore checklist
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
    Input('url', 'pathname'),
    Input("gwas-page-store", "modified_timestamp"),
    Input("gwas-page-store", "data"),
    State('shared-region-table', 'data'),
    
    prevent_initial_call=True
)
def restore_checklist_state(path,ts, data, table_data):
    analyse = table_data
    checkbox = []
    max_node_size = None
    min_node_size = 10
    min_percent_selected = 100
    tolerance_percentage = 0
    region_gap = 10000
    deletion_checkbox = ['show']
    deletion_percentage = 100
    if data is not None: 
        if "analyse" in data:
            analyse = data["analyse"]            
        if "checkboxes" in data:
            checkbox = data["checkboxes"]
        if "min_node_size" in data:
            min_node_size = data["min_node_size"]
        if "max_node_size" in data:
            max_node_size = data["max_node_size"]
        if "min_percent_selected" in data:
            min_percent_selected = data["min_percent_selected"]
        if "tolerance_percentage" in data:
            tolerance_percentage = data["tolerance_percentage"]
        if "region_gap" in data:
            region_gap = data["region_gap"]
        if "deletion_checkbox" in data:
            deletion_checkbox = data["deletion_checkbox"]
        if "deletion_percentage" in data:
            deletion_percentage = data["deletion_percentage"]
        if analyse is not None:
            for i, row in enumerate(analyse):
                row['get_sequence'] = "Get it"
    return (f"{len(analyse)} shared regions found.",analyse, checkbox, min_node_size, max_node_size,
            min_percent_selected,tolerance_percentage,region_gap, deletion_checkbox, deletion_percentage)

#Callback to save the gwas data table into csv file
@app.callback(
    Output('save-feedback', 'children'),
    Output("load_spinner_zone", "children", allow_duplicate=True),
    Output("download-csv", "data", allow_duplicate=True),
    Input('save-csv-button', 'n_clicks'),
    Input('save-csv-with_seq-button', 'n_clicks'),
    State('shared-region-table', 'data'),
    prevent_initial_call=True
)
def save_csv(n_clicks, n_clicks_seq, table_data):
    #logger.info(f"Callback triggered: n_clicks={n_clicks}, table_data={table_data}")
    if len(table_data) > 0:
        triggered_id = ctx.triggered_id
        export_sequences = False
        if triggered_id == 'save-csv-with_seq-button':
            export_sequences = True
        if not table_data:
            return "No data.",""
        df = pd.DataFrame(table_data).drop(columns=['get_sequence'], errors='ignore')
        if export_sequences :
            sequences = []
            for row in tqdm(table_data):
                sequences.append(get_sequence_from_position(row['genome'], row['chromosome'], row['start'], row['stop']))
            df["sequence"] = sequences
        if not SERVER_MODE:
            save_path = os.path.join(os.getcwd(), EXPORT_DIR, "shared_regions.csv")
            logger.info("save path : " + str(save_path))
            df.to_csv(save_path, index=False)

            return f"File saved : {save_path}","", no_update
        else:
            logger.info("🌐 Server mode active — file will be downloaded by user.")
            return "File ready for download.", "", dcc.send_data_frame(df.to_csv, "shared_regions.csv", index=False)
    else:
        return "No data to download.", "", no_update

#Callback to loads csv file into data table
@app.callback(
    Output('shared-status', 'children',allow_duplicate=True),
    Output('shared-region-table', 'data',allow_duplicate=True),
    Output("gwas-page-store", "data", allow_duplicate=True),
    Input('upload-csv', 'contents'),
    State('upload-csv', 'filename'),
    State("gwas-page-store", "data"),
    prevent_initial_call=True
)
def load_csv(contents, filename, gwas_page_store):
    logger.info("load csv file")
    if contents is None:
        return None, None, gwas_page_store
    if gwas_page_store is None:
        gwas_page_store = {}
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
        analyse = df[[c for c in df.columns if c!= "sequence"]].to_dict('records')
        if analyse is not None:
            for i, row in enumerate(analyse):
                row['get_sequence'] = "Get it"
        gwas_page_store["analyse"] = analyse
        logger.info("csv file loaded")
        return f"{len(analyse)} shared regions found.", analyse, gwas_page_store
    except Exception as e:
        logger.info(f"Error while loading file : {e}")
        return None, None, gwas_page_store
    

@app.callback(
    Output('upload-csv', 'style'),
    Input('load-csv-button', 'n_clicks'),
    prevent_initial_call=True
)
def show_upload_area(n_clicks):
    return {
        'display': 'block',
        'borderWidth': '1px',
        'borderStyle': 'dashed',
        'padding': '10px',
        'marginTop': '10px'
    }


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



