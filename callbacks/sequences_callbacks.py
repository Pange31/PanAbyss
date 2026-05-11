#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 13 12:36:54 2025

@author: fgraziani
"""
import dash
from dash import html, Input, Output, callback, State, dcc, callback_context, ctx, no_update
from dash.exceptions import PreventUpdate
from Bio.Seq import Seq
from app import *
from neo4j_requests import *
import io


MAX_NODES_FROM_DB = get_max_nodes_from_db()

def generate_sequence_table(sequences_dic):
    return html.Table([
    html.Thead(html.Tr([
        html.Th("Name", style={
            'border': '1px solid black',
            'padding': '6px',
            'width': '20%'
        }),
        html.Th("Sequence", style={
            'border': '1px solid black',
            'padding': '6px',
            'width': '80%',
            'whiteSpace': 'pre-wrap',
            'wordBreak': 'break-word'
        })
    ])),
    html.Tbody([
        html.Tr([
            html.Td(name, style={
                'border': '1px solid black',
                'padding': '6px',
                'width': '20%'
            }),
            html.Td(seq, style={
                'border': '1px solid black',
                'padding': '6px',
                'width': '80%',
                'whiteSpace': 'pre-wrap',
                'wordBreak': 'break-word'
            })
        ]) for name, seq in sequences_dic.items()
        ])
    ], style={
    'border': '1px solid black',
    'borderCollapse': 'collapse',
    'width': '100%',
    'tableLayout': 'fixed'
    })


@app.callback(
    Output('sequences-output', 'children'),
    Input('sequences-page-store', 'data'),
    prevent_initial_call=False
)
def load_sequences_on_page_load(sequences_dic):
    if sequences_dic and isinstance(sequences_dic, dict) and len(sequences_dic) > 0:
        return generate_sequence_table(sequences_dic)
    return html.Div("No sequences to display.")






@app.callback(
    Output('sequences-page-store', 'data', allow_duplicate=True),
    Output("sequences-message", "children", allow_duplicate=True),
    Output('shared_storage_nodes', 'data',allow_duplicate=True),
    Input('get-sequences-btn', 'n_clicks'),
    State('shared_storage_nodes', 'data'),
    State('home-page-store', 'data'),
    State('global_parameters', 'data'),
    prevent_initial_call=True
)
def display_sequences(n_clicks, nodes_data, home_data_storage,global_parameters):
    ctx = dash.callback_context
    if ctx.triggered_id == "get-sequences-btn" and n_clicks == 0:
        raise PreventUpdate
    nodes = no_update
    if not nodes_data:
        return {}, html.Div(html.P([
        "❌ No data to compute sequences. Select a region to visualise on the ",
        dcc.Link("home page", href="/", style={'color': 'blue', 'textDecoration': 'underline'}),
        " or on the ",
        dcc.Link("Shared regions discovery page", href="/gwas", style={'color': 'blue', 'textDecoration': 'underline'})
        ], style=error_style)), nodes
    else:
        # Step 1: Check if all nodes are in the region (since min node size can be set greater than 1)
        if "current_size" not in home_data_storage or home_data_storage["current_size"] > 1:
            # Get all the nodes from the region
            max_nodes_from_db = MAX_NODES_FROM_DB
            if global_parameters and "max_nodes_from_db" in global_parameters:
                max_nodes_from_db = global_parameters["max_nodes_from_db"]
            genome = home_data_storage.get("selected_genome", None)
            if "genome_zoom" in home_data_storage and home_data_storage["genome_zoom"]:
                genome = home_data_storage["genome_zoom"]
            use_anchor = not home_data_storage.get("zoom", False)

            chromosome = home_data_storage.get("selected_chromosome", None)
            start = home_data_storage.get("start", None)
            end = home_data_storage.get("end", None)
            logger.debug(f"Sequences construction: getting all the nodes for the region chr {chromosome} start {start} end {end} on genome {genome}")
            nodes, return_metadata = get_nodes_by_region(
                genome, chromosome=chromosome, start=start, end=end, use_anchor=use_anchor)
            logger.debug(f"Number of nodes in the region: {len(nodes)}")
            nodes_data = nodes

        sequences = []
        genomes_nodes_dic = {}
        set_names = set()
        for n in nodes_data:
            node = nodes_data[n]
            set_names.add(node["ref_node"])
            for g in node["genomes"]:
                if g not in genomes_nodes_dic:
                    genomes_nodes_dic[g] = []
                strand = "P"
                if "strandM" in node and g in node["strandM"]:
                    strand = "M"
                genomes_nodes_dic[g].append({"start":node[g+"_position"], "node_name":node["ref_node"], "strand":strand})

        sorted_names_by_genome = {
            genome: {
                "names": [item["node_name"] for item in sorted(nodes, key=lambda x: x["start"])],
                "strands": [item["strand"] for item in sorted(nodes, key=lambda x: x["start"])]
            }
            for genome, nodes in genomes_nodes_dic.items()
        }
        sequences_list = get_sequence_from_names(list(set_names))
        sequences_dic = {}
        for g in sorted_names_by_genome:
            sequence = ""
            for i in range(len(sorted_names_by_genome[g]["names"])):
                if sorted_names_by_genome[g]["strands"][i] == "M":
                    sequence += Seq(sequences_list[sorted_names_by_genome[g]["names"][i]]).reverse_complement()
                else:
                    sequence += sequences_list[sorted_names_by_genome[g]["names"][i]]
            sequences_dic[g] = str(sequence)
        return sequences_dic, "", nodes


@app.callback(
    Output("download-fasta", "data"),
    Output("sequences-message", "children", allow_duplicate=True),
    Input("get-fasta-btn", "n_clicks"),
    State("sequences-page-store", "data"),
    prevent_initial_call=True
)
def download_fasta(n_clicks, sequences_store):
    if not sequences_store or len(sequences_store) == 0:
        return no_update, html.Div(html.P([
        "❌ No data to compute sequences. Select a region to visualise on the ",
        dcc.Link("home page", href="/", style={'color': 'blue', 'textDecoration': 'underline'}),
        " or on the ",
        dcc.Link("Shared regions discovery page", href="/gwas", style={'color': 'blue', 'textDecoration': 'underline'}), " and compute the sequences."
        ], style=error_style))

    fasta_lines = []

    for name, seq in sequences_store.items():
        if not seq:
            continue

        fasta_lines.append(f">{name}")
        fasta_lines.append(seq)

    fasta_content = "\n".join(fasta_lines)

    return dcc.send_string(
        fasta_content,
        filename="sequences.fasta"
    ), html.Div(html.P(["✅ Fasta successfully downloaded."], style=success_style))