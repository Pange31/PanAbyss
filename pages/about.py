#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 13:52:04 2025

@author: fgraziani
"""


import dash_cytoscape as cyto
from dash import Dash, html,callback, dcc

from app import DB_VERSION, BLOCK_ADMIN_FUNCTIONNALITIES


import re


def get_current_version(changelog_path="CHANGELOG.md"):
    with open(changelog_path, "r", encoding="utf-8") as f:
        content = f.read()

    match = re.search(r"## \[(\d+\.\d+\.\d+)\]", content)
    if match:
        return match.group(1)

    return "unknown"

panabyss_version = get_current_version()

def layout():
    return html.Div([
        html.H2("About"),
        html.P("PanAbyss is a tool for managing pangenome graph based on a local graph database. It offers a development framework (back office functions) and an IHM based on this data modelisation."),
        #Help section
        html.Div([
            html.Ul([
                html.Li("Information about this tool :"),
                    html.Ul([
                        html.Li(f"PanAbyss version : {panabyss_version}"),
                        html.Li(f"Database version : {DB_VERSION}"),
                        html.Li("PanAbyss is developed and maintained at INRAE in the MIA-T laboratory, Genotoul-Bioinfo team, Toulouse, France." ),
                        html.Li([
                            "Project page: ",
                            html.A("https://github.com/Pange31/PanAbyss",
                                   href="https://github.com/Pange31/PanAbyss",
                                   target="_blank")
                        ]),
                    ]),
                    html.Li("Main functionnalities :"),
                        html.Ul([
                            html.Li("Generate a new database and load GFA / annotations files through the 'DB management' page (not accessible in server mode)."),
                            html.Li("Search and visualize a region of the pangenome graph through the 'Home' page." ), 
                            html.Li("Get the sequences of the displayed region for each haplotype through the 'Sequences' page." ),
                            html.Li("Plot the phylogenetic of the displayed region for each haplotype through the 'Phylogenetic' page." ),
                            html.Li("Find shared regions between selected haplotypes in the pangenome graph through the 'Shared regions finder' page." ),
                            html.Li("The database can be used directly from the 'http://localhost:7474' URL and data can be manipulated with cypher langages." ),
                            html.Li("The database can be used from back office python functions." ),
                        ])
            ]),
            None if BLOCK_ADMIN_FUNCTIONNALITIES else html.Button("Update PanAbyss", id='update-panabyss-btn', n_clicks=0, style={'margin': '15px 0'}),
            dcc.Loading(
                id="update-panabyss",
                #type="default",
                children=html.Div(id='update-panabyss-output'),
            )
            ], style={"marginBottom": "100px"}),
        html.Hr(),
        html.Div([
            html.Img(src='/assets/images/logo_genotoul_bioinfo.jpg',
                     style={'height': '150px', 'marginRight': '100px'}),
            html.Img(src='/assets/images/logo_miat.png',
                     style={'height': '100px', 'marginRight': '100px'}),
            html.Img(src='/assets/images/logo_INRAE.png',
                     style={'height': '75px', 'marginRight': '100px'})
        ],
        style={
            'display': 'flex',
            'alignItems': 'center',
            'justifyContent': 'start',  # ou 'center', selon le besoin
            'marginBottom': '20px'
        })
        
        
        ], style={'width': '48%', 'display': 'inline-block', 'verticalAlign': 'top', 'marginLeft':'4%'}),
        


