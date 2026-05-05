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
    Input('get-sequences-btn', 'n_clicks'),
    State('shared_storage_nodes', 'data'),
    prevent_initial_call=True
)
def display_sequences(n_clicks, nodes_data):
    ctx = dash.callback_context
    if ctx.triggered_id == "get-sequences-btn" and n_clicks == 0:
        raise PreventUpdate

    if not nodes_data:
        return {}, html.Div(html.P([
            "❌ No data to compute sequences. Select a region to visualise on the ",
            dcc.Link("home page", href="/", style={'color': 'blue', 'textDecoration': 'underline'}),
            " or on the ",
            dcc.Link("Shared regions discovery page", href="/gwas", style={'color': 'blue', 'textDecoration': 'underline'})
        ]))

    genomes_positions = {}
    chromosome = None
    for n in nodes_data:
        node = nodes_data[n]
        if not chromosome and "chromosome" in node:
            chromosome = node["chromosome"]
        for g in node["genomes"]:
            pos = node.get(f"{g}_position")

            if pos is None:
                continue

            if g not in genomes_positions:
                genomes_positions[g] = {"min": pos, "max": pos+node["size"]}
            else:
                genomes_positions[g]["min"] = min(genomes_positions[g]["min"], pos)
                genomes_positions[g]["max"] = max(genomes_positions[g]["max"], pos+node["size"])


    sequences_dic = {}

    for g, bounds in genomes_positions.items():
        start = bounds["min"]
        end = bounds["max"]

        seq = get_sequence_from_position(g, chromosome, start, end)
        sequences_dic[g] = str(seq)

    return sequences_dic, ""




# @app.callback(
#     Output('sequences-page-store', 'data', allow_duplicate=True),
#     Output("sequences-message", "children", allow_duplicate=True),
#     Input('get-sequences-btn', 'n_clicks'),
#     State('shared_storage_nodes', 'data'),
#     prevent_initial_call=True
# )
# def display_sequences(n_clicks, nodes_data):
#     ctx = dash.callback_context
#     if ctx.triggered_id == "get-sequences-btn" and n_clicks == 0:
#         raise PreventUpdate
#
#     if not nodes_data:
#         return {}, html.Div(html.P([
#         "❌ No data to compute sequences. Select a region to visualise on the ",
#         dcc.Link("home page", href="/", style={'color': 'blue', 'textDecoration': 'underline'}),
#         " or on the ",
#         dcc.Link("Shared regions discovery page", href="/gwas", style={'color': 'blue', 'textDecoration': 'underline'})
#         ], style=error_style))
#     else:
#         sequences = []
#         genomes_nodes_dic = {}
#         set_names = set()
#         for n in nodes_data:
#             node = nodes_data[n]
#             set_names.add(node["ref_node"])
#             for g in node["genomes"]:
#                 if g not in genomes_nodes_dic:
#                     genomes_nodes_dic[g] = []
#                 strand = "P"
#                 if "strandM" in node and g in node["strandM"]:
#                     strand = "M"
#                 genomes_nodes_dic[g].append({"start":node[g+"_position"], "node_name":node["ref_node"], "strand":strand})
#
#         sorted_names_by_genome = {
#             genome: {
#                 "names": [item["node_name"] for item in sorted(nodes, key=lambda x: x["start"])],
#                 "strands": [item["strand"] for item in sorted(nodes, key=lambda x: x["start"])]
#             }
#             for genome, nodes in genomes_nodes_dic.items()
#         }
#         sequences_list = get_sequence_from_names(list(set_names))
#         sequences_dic = {}
#         for g in sorted_names_by_genome:
#             sequence = ""
#             for i in range(len(sorted_names_by_genome[g]["names"])):
#                 if sorted_names_by_genome[g]["strands"][i] == "M":
#                     sequence += Seq(sequences_list[sorted_names_by_genome[g]["names"][i]]).reverse_complement()
#                 else:
#                     sequence += sequences_list[sorted_names_by_genome[g]["names"][i]]
#             sequences_dic[g] = str(sequence)
#         return sequences_dic, ""


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