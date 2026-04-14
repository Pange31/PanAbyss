#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 13:52:04 2025

@author: fgraziani
"""

from dash import Dash, html, callback, dcc, dash_table

def layout():
    return html.Div([
        html.H2("Phylogeny Store Management"),

        html.Details([
            html.Summary("ℹ️ Click here to display help"),
            html.Ul([
                html.Li("This page allows to manage database of Phylogeny results."),
            ]),
        ]),

        html.Hr(),

        html.H3("Phylogeny result storage"),
        dash_table.DataTable(
            id="phylo-jobs-table",

            style_table={
                "overflowX": "auto",
                "borderRadius": "10px",
                "border": "1px solid #ddd"
            },

            style_header={
                "backgroundColor": "#f4f4f4",
                "fontWeight": "bold",
                "border": "1px solid #ddd"
            },

            style_cell={
                "textAlign": "left",
                "padding": "8px",
                "border": "1px solid #eee",
                "whiteSpace": "normal",
                "height": "auto"
            },

            style_data_conditional=[
                {
                    "if": {"column_id": "action"},
                    "backgroundColor": "#f8f9fa",
                    "cursor": "pointer",
                    "fontWeight": "bold",
                },
                {
                    "if": {"filter_query": '{status} = "RUNNING"'},
                    "backgroundColor": "#fff8e1",
                },
            ],
        )
    ])

