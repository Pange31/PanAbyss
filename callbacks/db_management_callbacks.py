#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 22:03:10 2025

@author: fgraziani
"""

from dash import html, Output, Input, State, no_update, dcc, ctx, no_update, exceptions, background_callback as bg
from dash.dependencies import ALL, MATCH


import os
import sys
root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if root_path not in sys.path:
    sys.path.append(root_path)

from app import *
from neo4j_requests import *
from neo4j_DB_construction import *
from neo4j_container_management import *
from config import *
import base64
import io
import shutil
import logging
import json


logger = logging.getLogger("panabyss_logger")


PREFIX_CONTAINER_NAME = "DB_"+ DB_VERSION + "_"


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
GFA_FOLDER = os.path.join(PROJECT_ROOT, "data", "gfa")
IMPORT_FOLDER = os.path.join(PROJECT_ROOT, "data", "import")
CONF_FILE = os.path.join(PROJECT_ROOT, "conf.json")
DATA_FOLDER = os.path.join(PROJECT_ROOT, "data", "data")
ANNOTATIONS_FOLDER = os.path.join(PROJECT_ROOT, "data", "annotations")
INSTALL_CONF_FILE = os.path.join(PROJECT_ROOT, "install", "data")
DUMP_FILE = os.path.join(PROJECT_ROOT, "data", "import", "neo4j.dump")


def get_container_name_no_prefix(container_name):
    return re.sub(r'^DB_.[^_]+_', '',container_name)


############# Data callbacks#################



@app.callback(
        Output('upload-gfa-output', 'children'),
        Input('upload-gfa-data', 'contents'),
        State('upload-gfa-data', 'filename'),
        State('upload-gfa-data', 'last_modified'),
        prevent_initial_call=True
    )
@require_authorization
def save_uploaded_files(list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        saved_files = []
        for content, filename in zip(list_of_contents, list_of_names):
            if not filename.endswith(".gfa"):
                continue
            try:
                data = content.encode("utf8").split(b";base64,")[1]
                file_path = os.path.join(GFA_FOLDER, filename)
                with open(file_path, "wb") as f:
                    f.write(base64.b64decode(data))
                saved_files.append(html.Li(f"File saved : {filename}"))
            except Exception as e:
                saved_files.append(html.Li(f"Error for saving file {filename} : {str(e)}"))

        return html.Ul(saved_files)
    return html.Div("No file.")



@app.callback(
    Output("add-gfa-message", "children", allow_duplicate=True),
    Output({'type': 'gfa-checkbox', 'index': ALL}, 'value', allow_duplicate=True),
    Input("btn-load-gfa", "n_clicks"),
    State({'type': 'gfa-checkbox', 'index': ALL}, 'value'),
    State({'type': 'gfa-checkbox', 'index': ALL}, 'id'),
    State({'type': 'gfa-input', 'index': ALL}, 'value'),
    State({'type': 'gfa-input', 'index': ALL}, 'id'),
    running=[
        (Output("btn-load-gfa", "disabled"), True, False)
    ],
    prevent_initial_call=True
)
@require_authorization
def on_click_load_gfa(n_clicks, checkbox_values, checkbox_ids, input_values, input_ids):
    # Get selected files
    selected_files = [c_id['index'] for c_val, c_id in zip(checkbox_values, checkbox_ids) if c_val]
    if not selected_files:
        return html.Div("❌ Please select at least one GFA file before importing.", style=error_style), [no_update for _ in checkbox_ids]
    invalid_files = [f for f in selected_files if not f.lower().endswith(".gfa")]
    if invalid_files:
        return html.Div(f"❌ Invalid file(s): {', '.join(invalid_files)}. Please select only .gfa files.",
                        style=error_style), [no_update for _ in checkbox_ids]

    #Get batch size from conf
    batch_size = get_db_load_gfa_batch_size()
    logger.info(f"Load GFA files {selected_files} with batch size {batch_size}")
    # Get the chromosome associated to file (if set)
    chromosome_dict = {i_id['index']: val for i_id, val in zip(input_ids, input_values)}
    list_chromosome_file = [chromosome_dict.get(f, "") for f in selected_files]
    if len(list_chromosome_file) > 1:
        for cf in list_chromosome_file:
            if cf is None or cf == "":
                return html.Div(
                    "❌ When multiple gfa are selected, it is required to set the chromosome for each of theses files (non null or empty value).",
                    style=error_style), [no_update for _ in checkbox_ids]

    for file_name, chromosome_file in zip(selected_files, list_chromosome_file):
        start_time = time.time()
        if chromosome_file == "":
            chromosome_file = None
        file_path = os.path.join(GFA_FOLDER, file_name)
        load_sequences(file_path, chromosome_file=chromosome_file, create=True)
        load_gfa_data_to_neo4j(file_path, chromosome_file=chromosome_file, batch_size=batch_size,
                                start_chromosome=None, create=True, haplotype=True, create_only_relations=False)
    logger.info(f"Graph from {file_name} loaded in {time.time() - start_time:.2f} s")
    logger.info("creating indexes")
    create_indexes(base=False, extend=True, genomes_index=True)
    (ret, msg_indexes) = wait_for_indexes()
    return html.Div(f"✅ GFA files loaded successfully: {', '.join(selected_files)}", style=success_style),  [ [] for _ in checkbox_ids ]



@app.callback(
    Output("gfa-message", "children", allow_duplicate=True),
    Output({'type': 'gfa-checkbox', 'index': ALL}, 'value', allow_duplicate=True),
    Input("btn-csv-import", "n_clicks"),
    State({'type': 'gfa-checkbox', 'index': ALL}, 'value'),
    State({'type': 'gfa-checkbox', 'index': ALL}, 'id'),
    State({'type': 'gfa-input', 'index': ALL}, 'value'),
    State({'type': 'gfa-input', 'index': ALL}, 'id'),
    running=[
        (Output("btn-csv-import", "disabled"), True, False)
    ],
    prevent_initial_call=True
)
@require_authorization
def on_click_csv_import(n_clicks, checkbox_values, checkbox_ids, input_values, input_ids):
    # Get selected files
    selected_files = [c_id['index'] for c_val, c_id in zip(checkbox_values, checkbox_ids) if c_val]
    if not selected_files:
        return html.Div("❌ Please select at least one GFA file before importing.", style=error_style), [no_update for _ in checkbox_ids]
    invalid_files = [f for f in selected_files if not f.lower().endswith(".gfa")]

    if invalid_files:
        return html.Div(f"❌ Invalid file(s): {', '.join(invalid_files)}. Please select only .gfa files.",
                        style=error_style), [no_update for _ in checkbox_ids]

    batch_size = get_db_load_gfa_batch_size()
    logger.info(f"Generate csv from GFA files {selected_files} with batch size {batch_size}")
    # Get the chromosome associated to file (if set)
    chromosome_dict = {i_id['index']: val for i_id, val in zip(input_ids, input_values)}
    list_chromosome_file = [chromosome_dict.get(f, "") for f in selected_files]

    if len(list_chromosome_file) > 1:
        for cf in list_chromosome_file:
            if cf is None or cf == "":
                return html.Div(
                    "❌ When multiple gfa are selected, it is required to set the chromosome for each of theses files (non null or empty value).",
                    style=error_style), [no_update for _ in checkbox_ids]

    for file_name, chromosome_file in zip(selected_files, list_chromosome_file):
        start_time = time.time()
        if chromosome_file != "":
            file_path = os.path.join(GFA_FOLDER, file_name)
            load_gfa_data_to_csv(file_path, import_dir="./data/import",
                                 chromosome_file=chromosome_file,
                                 chromosome_prefix=False,
                                 batch_size=batch_size,
                                 start_chromosome=None,
                                 haplotype=True)
        logger.info(f"CSV generation from {file_name} loaded in {time.time() - start_time:.2f} s")

    logger.info("✅ All GFA files loaded.")

    return html.Div(f"✅ GFA files loaded successfully: {', '.join(selected_files)}", style=success_style), [ [] for _ in checkbox_ids ]



@app.callback(
    Output("index_stats-message", "children", allow_duplicate=True),
    Input("btn-create-index", "n_clicks"),
    running=[
        (Output("btn-create-index", "disabled"), True, False)
    ],
    prevent_initial_call=True
)
@require_authorization
def on_click_create_index(n_clicks):
    logger.info("create indexes")
    create_indexes(base=True, extend=True, genomes_index=True)
    (ret, msg_indexes) = wait_for_indexes()
    logger.info("Indexes created")
    return html.Div(f"✅ Indexes creation command successfully done.", style=success_style)


@app.callback(
    Output("index_stats-message", "children", allow_duplicate=True),
    Input("btn-create-stats", "n_clicks"),
    running=[
        (Output("btn-create-stats", "disabled"), True, False)
    ],
    prevent_initial_call=True
)
@require_authorization
def on_click_create_stats(n_clicks):
    logger.info("create stats")
    create_stats_from_nodes()
    logger.info("Stats created")
    return html.Div(f"✅ Stats creation command successfully done.", style=success_style)



############# Annotations callbacks#################


@app.callback(
        Output('upload-annotations-output', 'children'),
        Input('upload-annotations-data', 'contents'),
        State('upload-annotations-data', 'filename'),
        State('upload-annotations-data', 'last_modified'),
        prevent_initial_call=True
    )
@require_authorization
def save_uploaded_files(list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        saved_files = []
        for content, filename in zip(list_of_contents, list_of_names):
            if not filename.endswith(('.gff', '.gff3', '.gtf')):
                continue
            try:
                data = content.encode("utf8").split(b";base64,")[1]
                file_path = os.path.join(ANNOTATIONS_FOLDER, filename)
                with open(file_path, "wb") as f:
                    f.write(base64.b64decode(data))
                saved_files.append(html.Li(f"File saved : {filename}"))
            except Exception as e:
                saved_files.append(html.Li(f"Error for saving file {filename} : {str(e)}"))

        return html.Ul(saved_files)
    return html.Div("No file.")

@app.callback(
    Output("annotation-message", "children", allow_duplicate=True),
    Input("dropdown-genome", "value"),
    prevent_initial_call=True
)
@require_authorization
def clear_genome_error_on_selection(genome):
    if genome:
        return ""
    dash.exceptions.PreventUpdate 


@app.callback(
    Output("annotation-message", "children"),
    Output({'type': 'annotation-checkbox', 'index': ALL}, 'value'),
    Input("btn-load-annotations-with-link", "n_clicks"),
    State({'type': 'annotation-checkbox', 'index': ALL}, 'value'),
    State({'type': 'annotation-checkbox', 'index': ALL}, 'id'),
    State({'type': 'annotation-dropdown', 'index': ALL}, 'value'),
    State({'type': 'annotation-dropdown', 'index': ALL}, 'id'),
    running=[
        (Output("btn-load-annotations-with-link", "disabled"), True, False)
    ],
    prevent_initial_call=True
)
@require_authorization
def on_click_load_annotations(n_clicks, checkbox_values, checkbox_ids, dropdown_values, dropdown_ids):
    triggered_id = ctx.triggered_id
    annotation_time = time.time()

    # Build list of selected files and map each file to its associated genome
    selected_files = []
    genome_mapping = {}
    for cb_val, cb_id, dd_val, dd_id in zip(checkbox_values, checkbox_ids, dropdown_values, dropdown_ids):
        file_name = cb_id['index']
        if cb_val:  # file is selected
            selected_files.append(file_name)
            genome_mapping[file_name] = dd_val

    # Error if no file selected
    if not selected_files:
        return html.Div("❌ Please select at least one annotation file before loading.", style=error_style), checkbox_values

    # Ensure each selected file has an associated genome
    for f in selected_files:
        genome = genome_mapping.get(f)
        if not genome:
            return html.Div(f"❌ Please select a genome for annotation file '{f}'.", style=error_style), checkbox_values

    # Check index state for each genome
    for genome in set(genome_mapping.values()):
        state_index = check_state_index("NodeIndex" + genome + "_position")
        if state_index is None:
            return html.Div(f"❌ Index for genome '{genome}' has not been created. Please create the index before loading annotations.", style=error_style), checkbox_values
        if int(state_index) != 100:
            return html.Div(f"❌ Index for genome '{genome}' is not completely created (state: {state_index}%). Please wait.", style=warning_style), checkbox_values

    # Validate file extensions
    invalid_files = [f for f in selected_files if not f.lower().endswith((".gff", ".gff3", ".gtf"))]
    if invalid_files:
        return html.Div(f"❌ Invalid file(s): {', '.join(invalid_files)}. Please select only .gff, .gff3, or .gtf files.", style=error_style), checkbox_values

    # Loop through each selected file and load annotations
    for f in selected_files:
        genome = genome_mapping[f]
        file_path = os.path.join(ANNOTATIONS_FOLDER, f)
        logger.info(f"Load annotations for file '{f}' with genome '{genome}'")
        load_annotations_neo4j(file_path, genome_ref=genome, single_chromosome=None)

    # Optionally link annotations if the "with-link" button was clicked
    if triggered_id == "btn-load-annotations-with-link":
        for genome in set(genome_mapping.values()):
            logger.info(f"Link annotations for genome '{genome}'")
            WARN = creer_relations_annotations_neo4j(genome)

    if WARN:
        return_message = html.Div(
        f"❌ Chromosome names mismatch between graph and annotation file. No annotation linked.",
        style=warning_style
    )
    else:
        return_message = html.Div(
            f"✅ Annotation(s) '{', '.join(selected_files)}' loaded with associated genomes.",
            style=success_style
        )

    # Return success message and reset all checkboxes (empty list)
    return return_message, [[] for _ in checkbox_values]


#
# @app.callback(
#     Output("annotation-message", "children"),
#     Output("annotations-files-selector", "value"),
#     Input("btn-load-annotations-with-link", "n_clicks"),
#     #Input("btn-load-only-annotations", "n_clicks"),
#     #Input("btn-link-annotations", "n_clicks"),
#     State("dropdown-genome", "value"),
#     State('annotations-files-selector', 'value'),
#     prevent_initial_call=True
# )
# @require_authorization
# def on_click_load_annotation(n_clicks_load_all, genome, annotation_file_names):
#     triggered_id = ctx.triggered_id
#     annotation_time = time.time()
#     logger.info("annotations file names : " + str(annotation_file_names))
#     logger.info("triggered id : " + str(triggered_id))
#     if triggered_id == "btn-load-annotations-with-link" or triggered_id == "btn-load-only-annotations":
#         if not annotation_file_names:
#             return html.Div("❌ Please select an annotation file before loading.", style=error_style), no_update
#         if not genome or genome == "":
#             return html.Div("❌ Please select a reference haplotype.", style=error_style), no_update
#         state_index = check_state_index("NodeIndex"+genome+"_position")
#         if state_index is None:
#             return html.Div(f"❌ Index {state_index} has not been created, please create index before loading annotations.", style=error_style), no_update
#         if int(state_index) != 100:
#             return html.Div(f"❌ Index {state_index} is not completly created (creation state : {state_index}%). Please wait until this index has been created.", style=warning_style), no_update
#
#         # Force to list
#         if isinstance(annotation_file_names, str):
#             annotation_file_names = [annotation_file_names]
#
#         # ✅ Check all files have .gtf / .gff3 / .gff extension
#         invalid_files = [f for f in annotation_file_names if not f.lower().endswith(".gff") and not f.lower().endswith(".gff3") and not f.lower().endswith(".gtf")]
#         if invalid_files:
#             return html.Div(f"❌ Invalid file(s): {', '.join(invalid_files)}. Please select only .gff, .gff3 or gtf files.", style=error_style), no_update
#         genomes_list = []
#         for f in annotation_file_names:
#             logger.info(f"Load annotations for file {f}")
#             state_index = check_state_index("NodeIndex"+genome+"_position")
#             if state_index is None:
#                 return html.Div(f"❌ Index {state_index} has not been created, please create index before loading annotations.", style=error_style), no_update
#             if int(state_index) != 100:
#                 return html.Div(f"❌ Index {state_index} is not completly created (creation state : {state_index}%). Please wait until this index has been created.", style=warning_style), no_update
#             file = os.path.join(ANNOTATIONS_FOLDER, f)
#             annotation_time = time.time()
#             load_annotations_neo4j(file, genome_ref = genome, single_chromosome = None)
#         logger.info("Annotations loaded in " + str(time.time()-annotation_time) + " s.")
#     if triggered_id == "btn-load-annotations-with-link" or triggered_id == "btn-link-annotations" :
#         logger.info("Link annotations")
#         annotation_relation_time = time.time()
#         creer_relations_annotations_neo4j(genome)
#         logger.info("Link annotations in " + str(time.time()-annotation_relation_time) + " s.")
#     if triggered_id == "btn-load-annotations-with-link" or triggered_id == "btn-load-only-annotations":
#         return html.Div(f"✅ Annotation '{annotation_file_names}' loaded for genome '{genome}'.", style=success_style), []
#     else:
#         return html.Div(f"✅ Annotation linked.", style=success_style),[]

############# Delete data callbacks#################

@app.callback(
    Output("delete-confirmation", "children"),
    Input("btn-delete", "n_clicks"),
    prevent_initial_call=True
)
@require_authorization
def delete_data_ask_confirmation(n_clicks):
    data_dir =  os.path.join(DATA_FOLDER, "databases/neo4j")
    transactions_dir = os.path.join(DATA_FOLDER,"transactions/neo4j")
    if n_clicks > 0:
        return html.Div([
            html.Div("⚠️ Confirm: this operation will delete all data and container."),
            html.Button("Confirm Delete", id="btn-confirm-delete", n_clicks=0, style={"marginTop": "5px", "color": "white", "backgroundColor": "red"})
        ])
    return ""

@app.callback(
    Output("delete-message", "children"),
    Output("delete-confirmation", "children", allow_duplicate=True),
    Output('db-management-page-store', 'data', allow_duplicate=True),
    Input("btn-confirm-delete", "n_clicks"),
    State('db-management-page-store', 'data'),
    running=[
        (Output("btn-delete", "disabled"), True, False),
        (Output("btn-confirm-delete", "disabled"), True, False)
    ],
    prevent_initial_call=True
)
@require_authorization
def confirm_delete_data(n_clicks, data):
    if not n_clicks:
        raise exceptions.PreventUpdate
    stop_container()
    logger.info(f"Deleting following data : {DATA_FOLDER} and {IMPORT_FOLDER} directory, {CONF_FILE} file.")
    try:
        if os.path.exists(DATA_FOLDER):
            shutil.rmtree(DATA_FOLDER)
            os.makedirs(DATA_FOLDER)
        if os.path.exists(IMPORT_FOLDER):
            shutil.rmtree(IMPORT_FOLDER)
            os.makedirs(IMPORT_FOLDER)
        if os.path.exists(CONF_FILE):
            #Reset container name in conf file
            keys_to_remove = ["container_name"]
            with open(CONF_FILE, "r") as f:
                conf = json.load(f)
                container_name = conf.get("container_name", "")
                if container_name is not None and container_name != "" :
                    remove_container(container_name)
            for key in keys_to_remove:
                conf.pop(key, None)
            with open(CONF_FILE, "w") as f:
                json.dump(conf, f, indent=4)
        if "container_name" in data:
            data.pop("container_name")
        return html.Div("✅ All data deleted successfully.", style=success_style), "", data
    except Exception as e:
        return html.Div(f"❌ Error while deleting data: {str(e)}", style=error_style), "", data


############# delete annotations callback ##################
@app.callback(
    Output("delete-annotations-confirmation", "children"),
    Input("btn-delete-annotations", "n_clicks"),
    prevent_initial_call=True
)
@require_authorization
def delete_annotations_ask_confirmation(n_clicks):
    if n_clicks > 0:
        return html.Div([
            html.Div("⚠️ Confirm: this operation will delete all annotations in database."),
            html.Button("Confirm Delete", id="btn-delete-annotations-confirmation", n_clicks=0, style={"marginTop": "5px", "color": "white", "backgroundColor": "red"})
        ])
    return ""

@app.callback(
    Output("delete-annotations-message", "children"),
    Output("delete-annotations-confirmation", "children", allow_duplicate=True),
    Input("btn-delete-annotations-confirmation", "n_clicks"),
    running=[
        (Output("btn-delete-annotations", "disabled"), True, False),
        (Output("btn-delete-annotations-confirmation", "disabled"), True, False)
    ],
    prevent_initial_call=True
)
@require_authorization
def confirm_delete_annotations(n_clicks):
    if not n_clicks:
        raise exceptions.PreventUpdate
    try:
        delete_annotations()
        return html.Div("✅ All annotations deleted successfully.", style=success_style), ""
    except Exception as e:
        return html.Div(f"❌ Error while deleting annotations: {str(e)}", style=error_style), ""



############# Container callbacks#################

@app.callback(
    Output('container-name-label', 'children'),
    Output('container-name-input', 'value'),
    Output('db-management-page-store', 'data'),
    Output("container-name-input", "disabled"),
    Output("btn-load-annotations-with-link", "disabled"),
    Output("btn-create-db", "style"),
    Output("btn-load-gfa", "style"),
    Output("btn-create-index", "disabled"),
    Output("btn-create-stats", "disabled"),
    Output("btn-dump-db", "disabled"),
    Input('db-management-page-store', 'data'),
    State('container-name-input', 'value'),
)
@require_authorization
def update_label(data, container_input):
    style_create_db = {"display": "inline-block"}
    style_add_gfa = {"display": "none"}
    if data is None :
        data = {}
    if "container_name" not in data :

        conf = load_config_from_json()
        if not conf or "container_name" not in conf or conf['container_name'] is None or conf['container_name'] == "":
            if container_input is not None and container_input != "":
                c_name = container_input
            else:
                c_name = "container_name"
            return f'No conf file found. Use "create new DB" procedure to generate it.', c_name, data, False, True, style_create_db, style_add_gfa, True, True, True
        else:
            style_create_db = {"display": "none"}
            style_add_gfa = {"display": "inline-block"}
            container_name = conf.get("container_name")
            data['container_name'] = get_container_name_no_prefix(container_name)
            return f'Container name : {container_name}', get_container_name_no_prefix(container_name), data, True, False, style_create_db, style_add_gfa, False, False, False
    else:
        style_create_db = {"display": "none"}
        style_add_gfa = {"display": "inline-block"}
        container_name = PREFIX_CONTAINER_NAME+data['container_name']
        return f"Container name : {container_name}", get_container_name_no_prefix(container_name), data, True, False, style_create_db, style_add_gfa, False, False, False
    
    
@app.callback(
    Output("create-db-confirmation", "children"),
    Input("btn-create-db", "n_clicks"),

    prevent_initial_call=True
)
@require_authorization
def create_db_ask_confirmation(n_clicks):
    if n_clicks > 0:
        return html.Div([
            html.Div("⚠️ Confirm: this operation will delete all data and container."),
            html.Button("Confirm creation of new DB", id="btn-confirm-create-db", n_clicks=0, style={"marginTop": "5px", "color": "white", "backgroundColor": "red"})
        ])
    return ""


@app.callback(
    Output("create-db-message", "children"),
    Output("create-db-confirmation", "children", allow_duplicate=True),
    Output('db-management-page-store', 'data', allow_duplicate=True),
    Output({'type': 'annotation-dropdown', 'index': ALL}, 'options'),
    Output({'type': 'gfa-checkbox', 'index': ALL}, 'value'),
    Output("container-name-input","value", allow_duplicate=True),
    Input("btn-confirm-create-db", "n_clicks"),
    State("container-name-input","value"),
    State("docker-image-dropdown","value"),
    State('db-management-page-store', 'data'),
    State('annotations-files-container', 'children'),
    State({'type': 'gfa-checkbox', 'index': ALL}, 'value'),
    State({'type': 'gfa-checkbox', 'index': ALL}, 'id'),
    State({'type': 'gfa-input', 'index': ALL}, 'value'),
    State({'type': 'gfa-input', 'index': ALL}, 'id'),
    running=[
        (Output("btn-create-db", "disabled"), True, False),
        (Output("btn-confirm-create-db", "disabled"), True, False)
    ],
    prevent_initial_call=True
)
@require_authorization
def confirm_create_db(n_clicks, container_name, docker_image, data, children, checkbox_values, checkbox_ids, input_values, input_ids):
    if data is None:
        data = {}

    n_dropdowns = len(children) if children else 0
    no_update_list = [no_update] * n_dropdowns

    if not n_clicks:
        raise exceptions.PreventUpdate

    # Check container name
    if not container_name:
        return html.Div("❌ No container name", style=error_style), "", data, no_update_list, checkbox_values, "container_name"
    else:
        container_name_prefixed = PREFIX_CONTAINER_NAME + container_name
    try:
        genomes_set = set()
        chromosomes_set = set()
        #Check if gfa files are selected
        selected_files = [c_id['index'] for c_val, c_id in zip(checkbox_values, checkbox_ids) if c_val]
        if selected_files:

            invalid_files = [f for f in selected_files if not f.lower().endswith(".gfa")]
            if invalid_files:
                return html.Div(f"❌ Invalid file(s): {', '.join(invalid_files)}. Please select only .gfa files.",
                                style=error_style), "", data, no_update_list, checkbox_values, container_name

            # Get batch size from conf
            batch_size = get_db_load_gfa_batch_size()
            logger.info(f"Generate csv from GFA files {selected_files} with batch size {batch_size}")
            # Get the chromosome associated to file (if set)
            chromosome_dict = {i_id['index']: val for i_id, val in zip(input_ids, input_values)}
            list_chromosome_file = [chromosome_dict.get(f, "") for f in selected_files]

            if len(list_chromosome_file) > 1:
                for cf in list_chromosome_file:
                    if cf is None or cf == "":
                        return html.Div(
                            "❌ When multiple gfa are selected, it is required to set the chromosome for each of theses files (non null or empty value).",
                            style=error_style), "", data, no_update_list, checkbox_values, container_name
            for file_name, chromosome_file in zip(selected_files, list_chromosome_file):
                start_time = time.time()
                if chromosome_file != "":
                    file_path = os.path.join(GFA_FOLDER, file_name)
                    genomes_analysed, chromosomes_analysed = load_gfa_data_to_csv(file_path, import_dir="./data/import",
                                         chromosome_file=chromosome_file,
                                         chromosome_prefix=False,
                                         batch_size=batch_size,
                                         start_chromosome=None,
                                         haplotype=True)
                    genomes_set = genomes_set | genomes_analysed
                    chromosomes_set = chromosomes_set | set(chromosomes_analysed)
                logger.info(f"CSV generation from {file_name} loaded in {time.time() - start_time:.2f} s")
        logger.info("All import files have been generated from gfa files in {time.time() - start_time:.2f} s, creating database.")
        creation_mode = create_db(container_name_prefixed, docker_image)
        # If creation by importing csv files it is necessary to create stats and indexes
        if creation_mode == "csv":
            stats = False
            if len(genomes_set) > 0 and len (chromosomes_set) > 0:
                logger.info("creating stats")
                create_stats(genomes_set, chromosomes_set)
                stats = True
            logger.info("creating base indexes")
            create_indexes(base=True, extend=True, genomes_index=False)
            if check_state_index("NodeIndexChromosome") is not None:
                t = 0
                while int(check_state_index("NodeIndexChromosome")) < 100 and t < MAX_TIME_INDEX:
                    time.sleep(10)
                    t += 10
                if not stats:
                    logger.info("creating stats from core node")
                    create_stats_from_nodes()
            logger.info("creating other indexes")
            create_indexes(base=False, extend=False, genomes_index=True)

        #Wait for all indexes are created
        (ret, msg_indexes) = wait_for_indexes()
        genomes = get_genomes() or []
        options_list = [
            [{"label": genome, "value": genome} for genome in genomes]
            for _ in range(n_dropdowns)
        ]
        data['container_name'] = container_name
        if ret == 0 :
            msg = "✅ DB and indexes successfully created"
        else:
            msg = f"❌ DB successfully created but error while creating index : {msg_indexes}"
        if len(selected_files) > 0:
            msg += f" with gfa files : {selected_files}"
        return html.Div(msg, style=success_style), "", data, options_list, [[] for _ in checkbox_values], container_name

    except Exception as e:
        logger.error(f"Error while creating database: {e}")
        return html.Div(f"❌ Error while creating database: {str(e)}", style=error_style), "", data, no_update_list, checkbox_values, container_name

@app.callback(
    Output("dump-message", "children", allow_duplicate=True),
    Input("btn-dump-db", "n_clicks"),
    State("docker-image-dropdown","value"),
    State('db-management-page-store', 'data'),
    running=[
        (Output("btn-dump-db", "disabled"), True, False)
    ],
    prevent_initial_call=True
)
@require_authorization
def dump_db_callback(n_clicks, docker_image, data):
    if not n_clicks:
        raise exceptions.PreventUpdate
    if os.path.exists(DUMP_FILE):
        os.remove(DUMP_FILE)

    if 'container_name' in data and data["container_name"] is not None and data["container_name"] != "":
        container_name = PREFIX_CONTAINER_NAME+data["container_name"]
    else:
        return html.Div("❌ No container name, you have to create the database before", style=error_style)
    try:
        dump_db(container_name, docker_image=DOCKER_IMAGE)
        return html.Div("✅ DB successfully dumped.", style=success_style)
    except Exception as e:
        return html.Div(f"❌ Error while dumping database: {str(e)}", style=error_style)