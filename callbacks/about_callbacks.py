#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 11:16:02 2025

@author: fgraziani
"""

import os
import requests
import zipfile
import io
import shutil
from dash import html, Input, Output, callback, State, dcc, no_update
from dash.exceptions import PreventUpdate
from Bio.Seq import Seq
from app import *
from neo4j_requests import *
from auth_utils import require_authorization
import logging

from config import get_max_nodes_from_db, get_max_nodes_to_visualize

MAX_NODES_TO_VISUALIZE = get_max_nodes_to_visualize()
MAX_NODES_FROM_DB = get_max_nodes_from_db()

logger = logging.getLogger("panabyss_logger")


# Dépôt cible (format : "owner/repo")
repo = "Pange31/PanAbyss"
update_dir = "./update"
panabyss_dir = "."

# URL de l'API pour la dernière release
url = f"https://api.github.com/repos/{repo}/releases/latest"

@require_authorization
def copy_update_to_root():
    for root, dirs, files in os.walk(update_dir):
        rel_path = os.path.relpath(root, update_dir)
        target_dir = os.path.join(panabyss_dir, rel_path)

        os.makedirs(target_dir, exist_ok=True)

        for file in files:
            src_file = os.path.join(root, file)
            dest_file = os.path.join(target_dir, file)
            shutil.copy2(src_file, dest_file)

    if os.path.exists(update_dir):
        shutil.rmtree(update_dir)
    logger.info("✅ The PanAbyss update has been successfully completed. It is required to restart the server.")


@app.callback(
    Output('update-panabyss-output', 'children'),
    Input('update-panabyss-btn', 'n_clicks'),
    prevent_initial_call=True
)
@require_authorization
def update_panabyss(n_clicks):

    # Requête GET
    response = requests.get(url)
    
    # Vérification du code de statut
    if response.status_code == 200:
        release = response.json()
        version = release["tag_name"]
        major_version = version.split(".")[0]
        if major_version == DB_VERSION.split(".")[0]:

            logger.info(f"Latest release found : {release['name']}")
            logger.info(f"Tag : {release['tag_name']}")
            logger.info(f"Publication date : {release['published_at']}")
            zip_url = release.get("zipball_url")
            logger.info(f"Download zip file : {zip_url}")
            
            r = requests.get(zip_url)
            if r.status_code != 200:
                logger.error(f"Error {r.status_code} while downloading")
                exit(1)
            
            # Créer dossier update s'il n'existe pas
            os.makedirs(update_dir, exist_ok=True)
            
            # Dézipper le contenu de l'archive dans /update en écrasant
            with zipfile.ZipFile(io.BytesIO(r.content)) as z:
                # Attention : l'archive GitHub zipball a un dossier racine unique
                root_folder = z.namelist()[0].split('/')[0]
            
                for member in z.namelist():
                    # retirer le dossier racine unique
                    filename = member[len(root_folder)+1:]
                    if not filename:
                        continue  # passer le dossier racine
            
                    dest_path = os.path.join(update_dir, filename)
                    if member.endswith('/'):
                        os.makedirs(dest_path, exist_ok=True)
                    else:
                        with open(dest_path, 'wb') as f:
                            f.write(z.read(member))
            copy_update_to_root()
            return html.Div(f"✅ The PanAbyss update has been successfully completed. It is required to restart the server.")
        else:
            logger.info(f"Latest version {version} is not compatible with the current data. To use latest release it is required to regenerate data.")
            return html.Div(f"❌ Latest version {version} is not compatible with the current data. To use latest release it is required to regenerate data.")
            
    else:
        logger.error(f"Error : {response.status_code}")
        return html.Div(f"❌ Error : {response.status_code}")


@app.callback(
    Output('global_parameters', 'data'),
    Input('btn-update-global-parameters', 'n_clicks'),
    State('limit-nodes-from-db-input', 'value'),
    State('limit-nodes-to-visualize', 'value'),
    State('circular-nodes-check', 'value'),
    prevent_initial_call=True
)

def update_parameters(n_clicks, limit_from_db, limit_to_visualize, circle):
    if n_clicks > 0:
        global_parameters = {}
        if limit_from_db and limit_from_db >= 100:
            global_parameters["max_nodes_from_db"] = limit_from_db
            logger.debug(f"Update the nodes limit from db to {limit_from_db}")
        if limit_to_visualize and limit_to_visualize >= 100:
            global_parameters["max_nodes_to_visualize"] = limit_to_visualize
            logger.debug(f"Update the nodes limit to visualize to {limit_to_visualize}")
        if 'circle' in circle:
            global_parameters["circle"] = True
            logger.debug(f"Update the nodes shape to circle")
        return global_parameters
    return no_update

@app.callback(
    Output('global_parameters', 'data', allow_duplicate=True),
    Output('limit-nodes-from-db-input', 'value', allow_duplicate=True),
    Output('limit-nodes-to-visualize', 'value', allow_duplicate=True),
    Output('circular-nodes-check', 'value', allow_duplicate=True),
    Input('btn-reset-global-parameters', 'n_clicks'),
    prevent_initial_call=True
)
def reset_parameters(n_clicks):
    if n_clicks > 0:
        global_parameters = {}

        return global_parameters, MAX_NODES_FROM_DB, MAX_NODES_TO_VISUALIZE, []
    return no_update, no_update, no_update, no_update

@app.callback(
    Output('limit-nodes-from-db-input', 'value'),
    Output('limit-nodes-to-visualize', 'value'),
    Output('circular-nodes-check', 'value'),
    Input('url', 'pathname'),
    State('global_parameters', 'data'),
    prevent_initial_call=False
)
def update_parameters_on_page_load(pathname, global_parameters ):
    if pathname != "/about":
        raise PreventUpdate
    limit_nodes_to_visualize = MAX_NODES_TO_VISUALIZE
    limit_nodes_from_db = MAX_NODES_FROM_DB
    circle_node = []
    if global_parameters:
        print(global_parameters)
        if "max_nodes_from_db" in global_parameters:
            limit_nodes_from_db = global_parameters["max_nodes_from_db"]
        if  "max_nodes_to_visualize" in global_parameters:
            limit_nodes_to_visualize = global_parameters["max_nodes_to_visualize"]
        if 'circle' in global_parameters and global_parameters['circle'] == True:
            circle_node = ['circle']
    return limit_nodes_from_db, limit_nodes_to_visualize, circle_node
