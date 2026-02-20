#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 12:06:35 2025

@author: fgraziani
"""

from dash import Dash, dcc, html, Input, Output, State, ctx,no_update
import dash_bootstrap_components as dbc
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


from neo4j_requests import *
from neo4j_container_management import *
from config import *
import logging

app.config.suppress_callback_exceptions = True
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

start_container()

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
    dcc.Store(id="db-management-page-store", data={}, storage_type="memory"),
    #A bug in dash requires to set memory for gwas-page-store because this storage is used by a background treatment
    #Else it will raise a "Maximum depth" error in the gwas process
    dcc.Store(id="gwas-page-store", storage_type="memory"),
    dcc.Store(id="parameters-gwas-page-store", data={}, storage_type="memory"),
    dcc.Store(id="phylogenetic-page-store",storage_type="memory"),
    dcc.Store(id="phylo-job-trigger",storage_type="memory"),
    dcc.Store(id="phylo-job-status",storage_type="memory"),
    dcc.Store(id="phylo-local-tree-job-status",storage_type="memory"),
    dcc.Store(id="sequences-page-store", data={'sequences':[]},storage_type="memory"),
    dcc.Store(id="db-management-create-db-trigger", storage_type="memory"),
    dcc.Store(id="db-management-load-annotations-trigger", storage_type="memory"),
    dcc.Store(id="db-management-job-trigger",storage_type="memory"),

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

    #Toast and spinner store for asynchron process
    dcc.Store(id="global-notification"),
    dcc.Store(id="global-spinner"),
    #Toast used for asynchron treatments
    html.Div(
        id="global-toast",
        children=[
            html.Div(id="toast-header", style={
                "fontWeight": "bold",
                "marginBottom": "6px",
                "display": "flex",
                "alignItems": "center",
                "gap": "8px"
            }),
            html.Div(id="toast-body", style={
                "fontSize": "14px"
            }),
        ],
        style={
            "position": "fixed",
            "top": "20px",
            "right": "20px",
            "width": "350px",
            "padding": "14px 18px",
            "borderRadius": "8px",
            "boxShadow": "0px 4px 12px rgba(0,0,0,0.2)",
            "zIndex": 9999,
            "display": "none",
            "color": "white",
        },
    ),
    #interval to close the toast
    dcc.Interval(id="toast-auto-close", interval=4000, disabled=True),


# main content
html.Div(id="page-content", style={"marginLeft": "10px", "padding": "20px"}),
])



# app.validation_layout = html.Div([
#     app.layout,
#     home.layout(),
#     phylogenetic.layout(),
#     gwas.layout(),
#     sequences.layout(),
#     db_management.layout(),
#     about.layout()
# ])


app.validation_layout = html.Div([
    dcc.Store(id='shared_storage'),
    dcc.Store(id='shared_storage_nodes'),
    dcc.Store(id='home-page-store'),
    dcc.Store(id='gwas-page-store'),
    dcc.Store(id='parameters-gwas-page-store'),
    dcc.Store(id='phylogenetic-page-store'),
    dcc.Store(id='phylo-job-trigger'),
    dcc.Store(id='phylo-job-status'),
    dcc.Store(id='sequences-page-store'),
    dcc.Store(id='db-management-job-trigger'),
    dcc.Store(id="db-management-page-store", data={}, storage_type="memory"),
])
import callbacks.phylogenetic_callbacks
import callbacks.gwas_callbacks
import callbacks.sequences_callbacks
import callbacks.db_management_callbacks
import callbacks.about_callbacks

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

#callback to display toast
@app.callback(
    Output("global-toast", "style"),
    Output("toast-header", "children"),
    Output("toast-body", "children"),
    Output("toast-auto-close", "disabled"),
    Input("global-notification", "data"),
    prevent_initial_call=True
)
def display_toast(notification):

    if not notification:
        return no_update, no_update, no_update, True

    toast_type = notification.get("type", "success")

    colors = {
        "success": "#28a745",
        "danger": "#dc3545",
        "warning": "#f0ad4e",
        "info": "#17a2b8",
        "primary": "#007bff",
    }

    icons = {
        "success": "✔",
        "danger": "✖",
        "warning": "⚠",
        "info": "ℹ",
        "primary": "🔔",
    }

    background = colors.get(toast_type, "#28a745")
    icon = icons.get(toast_type, "✔")

    style = {
        "position": "fixed",
        "top": "20px",
        "right": "20px",
        "width": "350px",
        "padding": "14px 18px",
        "borderRadius": "8px",
        "boxShadow": "0px 4px 12px rgba(0,0,0,0.2)",
        "zIndex": 9999,
        "color": "white",
        "backgroundColor": background,
        "display": "block",
    }

    header_content = [
        html.Span(icon),
        html.Span(notification.get("title", "Notification"))
    ]

    return (
        style,
        header_content,
        notification.get("message", ""),
        False,  # active l'auto close
    )

#callback to close toast
@app.callback(
    Output("global-toast", "style", allow_duplicate=True),
    Output("toast-auto-close", "disabled", allow_duplicate=True),
    Input("toast-auto-close", "n_intervals"),
    prevent_initial_call=True
)
def close_toast(_):

    style = {
        "display": "none"
    }

    return style, True

@app.callback(
    Output('tabs-navigation', 'value',allow_duplicate=True),
    Input('url', 'pathname'),
    prevent_initial_call=True,
)
def sync_tabs_with_url(pathname):
    """
    Synchronis current selected tab with the current path.
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
    # if pathname in ["/", "", None]:
    #     return home_layout
    # elif pathname == "/phylogenetic":
    #     return phylo_layout
    # elif pathname == "/gwas":
    #     return gwas_layout
    # elif pathname == "/db_management":
    #     if BLOCK_ADMIN_FUNCTIONNALITIES:
    #         return html.H3("🚫 Access denied — administration is disabled.", style={"color": "red"})
    #     return db_mgmt_layout
    # elif pathname == "/sequences":
    #     return sequences_layout
    # elif pathname == "/about":
    #     return about_layout
    # else:
    #     return html.H1("Page non trouvée")

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
    print("SERVER START")
    app.run(debug=True, port = port)
    
    
if __name__ == "__main__":

    run()
