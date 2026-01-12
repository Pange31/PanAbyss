#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 12:06:35 2025

@author: fgraziani
"""

from dash import dcc, html, Input, Output, State, ctx
import dash_auth
from app import *
from sidebar import sidebar
import signal
import sys
import argparse
import pages.home as home
import pages.phylogenetic as phylogenetic
import pages.sequences as sequences
import pages.gwas as gwas
import pages.db_management as db_management
import pages.about as about
import callbacks.phylogenetic_callbacks
import callbacks.gwas_callbacks
import callbacks.sequences_callbacks
import callbacks.db_management_callbacks
import callbacks.about_callbacks

from neo4j_requests import *
from neo4j_container_management import *
from config import *
import logging


logger = logging.getLogger("panabyss_logger")

#Limit upload size for gfa / annotations files to 10 Go
MAX_UPLOAD_SIZE = 10 * 1024 * 1024 * 1024



#Limit upload size
app.server.config['MAX_CONTENT_LENGTH'] = MAX_UPLOAD_SIZE
app.server.secret_key = "KEY_PANABYSS_96598421_CDEYUJH"

def clean_exit(signum, frame):
    logger.info("\nStopping docker")
    stop_container()
    time.sleep(10)
    logger.info("\nPanAbyss stopped")
    sys.exit(0)

# Close docker when quitting app
signal.signal(signal.SIGINT, clean_exit)
signal.signal(signal.SIGTERM, clean_exit)

logger.info(f"Server mode : {SERVER_MODE} - Admin mode : {ADMIN_MODE}")

if SERVER_MODE and ADMIN_MODE:
    USERS = get_users()
    auth = dash_auth.BasicAuth(app, USERS)
    logger.info("🔒 Admin functionnalities and server mode = True → credentials needed")
else:
    logger.info("🔓 No Admin functionnalities or no server mode → public access")


tabs = [
    dcc.Tab(label='Home', value='/', className='custom-tab', selected_className='custom-tab--selected'),
    dcc.Tab(label='Shared regions discovery', value='/gwas', className='custom-tab', selected_className='custom-tab--selected'),
    dcc.Tab(label='Phylogenetic', value='/phylogenetic', className='custom-tab', selected_className='custom-tab--selected'),
    dcc.Tab(label='Sequences', value='/sequences', className='custom-tab', selected_className='custom-tab--selected'),
    dcc.Tab(label='About', value='/about', className='custom-tab', selected_className='custom-tab--selected'),
]

if not BLOCK_ADMIN_FUNCTIONNALITIES:
    tabs.insert(-1, dcc.Tab(label='DB management', value='/db_management', className='custom-tab', selected_className='custom-tab--selected'))


app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    dcc.Store(id='shared_storage_nodes', data=[], storage_type='memory'),
    dcc.Store(id='shared_storage', data={'genomes':[], 'chromosomes':[]}, storage_type='session'),
    dcc.Store(id="home-page-store", storage_type='session'),
    dcc.Store(id="phylogenetic-page-store",storage_type="memory"),
    dcc.Store(id="sequences-page-store", data={'sequences':[]},storage_type="memory"),


    html.Div(
        children=[
            dcc.Tabs(
                id="tabs-navigation",
                value=None,
                children=tabs,
            )
        ],
        style={"marginBottom": "20px"}
    ),


# Menu button
# html.Div([
#     dcc.Store(id="sidebar-state", data={"visible":True}),
#     html.Button("☰", id="btn-toggle-sidebar", n_clicks=0),
# ], style={"padding": "10px"}),

# Sidebar
#html.Div(id="sidebar-container", children=sidebar, style={"display": "none", "width": "250px", "backgroundColor": "#f4f4f4", "position": "fixed", "height": "100%", "zIndex": 1}),

# main content
html.Div(id="page-content", style={"marginLeft": "10px", "padding": "20px"})
])

# @app.callback(
#     Output("sidebar-state", "data"),
#     Input("btn-toggle-sidebar", "n_clicks"),
#     Input("url", "pathname"),
#     State("sidebar-state", "data"),
#     prevent_initial_call=True
# )
# def toggle_sidebar(n_clicks, pathname, data):
#     triggered = ctx.triggered_id

#     if triggered == "btn-toggle-sidebar":
#         return {"visible": not data["visible"]}
#     else:
#         # Hide sidebar when changing page
#         return {"visible": False}


# #Toggle menu
# @app.callback(
#     Output("sidebar-container", "style"),
#     Input("sidebar-state", "data")
# )
# def toggle_sidebar(data):
#     if data["visible"]:
#         return {"display": "block"}
#     else:
#         return {"display": "none"}

#Getting chromosomes and genomes
@app.callback(
    Output('shared_storage', 'data'),
    Input('url', 'pathname'),
    prevent_initial_call=False 
)
def init_data(pathname):
    new_data = {}
    all_genomes = get_genomes()
    all_genomes.sort()
    #logger.info("all genomes : " + str(all_genomes))
    new_data["genomes"] = all_genomes
    new_data["chromosomes"]  = get_chromosomes()
    new_data["features"] = get_annotations_features()
    return new_data


@app.callback(
    Output('tabs-navigation', 'value',allow_duplicate=True),
    Input('url', 'pathname'),
    prevent_initial_call=True,
)
def sync_tabs_with_url(pathname):
    """
    Synchronise la tab sélectionnée avec la route actuelle.
    """
    if pathname is None:
        return '/'
    return pathname


@app.callback(
    Output('url', 'pathname'),
    Input('tabs-navigation', 'value'),
    prevent_initial_call=True,
    allow_duplicate=True
)
def update_url_from_tab(tab_value):
    return tab_value

#Routing
@app.callback(
    Output("page-content", "children"),
    Input("url", "pathname")
)
def display_page(pathname):
    #logger.info("callback routing " + str(pathname))
    if pathname in ["/", "", None]:
        return home.layout()
    elif pathname == "/phylogenetic":
        return phylogenetic.layout()
    elif pathname == "/gwas":
        return gwas.layout()
    elif pathname == "/db_management":
        if BLOCK_ADMIN_FUNCTIONNALITIES:
            return html.H3("🚫 Access denied — administration is disabled.", style={"color": "red"})
        return db_management.layout()
    elif pathname == "/sequences":
        return sequences.layout()
    elif pathname == "/about":
        return about.layout()
    else:
        return html.H1("Page non trouvée")




def run():
    parser = argparse.ArgumentParser(description="Launch server.")
    parser.add_argument("--port", type=int, help="HTTP port to use (default : 8050)")
    args = parser.parse_args()
    
    port = args.port or int(8050)
    start_container()

    app.run(debug=True, port = port)
    
    
if __name__ == "__main__":

    run()
