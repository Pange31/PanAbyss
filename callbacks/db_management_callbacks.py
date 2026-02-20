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


#Function to delete asynchron messages
def delete_messages(data):
    if "load_annotations_message" in data:
        data["load_annotations_message"] = ""
    if "db_creation_message" in data:
        data["db_creation_message"] = ""
    return data


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



#callback called on click on the load annotations button
@app.callback(
    Output("annotation-message", "children", allow_duplicate=True),
    Output("db-management-page-store", "data", allow_duplicate=True),
    Output("db-management-load-annotations-trigger", "data", allow_duplicate=True),
    Output("delete-annotations-message", "children", allow_duplicate=True),
    Input("btn-load-annotations-with-link", "n_clicks"),
    State({'type': 'annotation-checkbox', 'index': ALL}, 'value'),
    State({'type': 'annotation-checkbox', 'index': ALL}, 'id'),
    State({'type': 'annotation-dropdown', 'index': ALL}, 'value'),
    State({'type': 'annotation-dropdown', 'index': ALL}, 'id'),
    State('db-management-page-store', 'data'),
prevent_initial_call = True
)
@require_authorization
def on_click_load_annotations(n_clicks, annotation_checkbox_values, annotation_checkbox_ids, annotation_dropdown_values, annotation_dropdown_ids, data):
    triggered_id = ctx.triggered_id
    if not n_clicks:
        raise exceptions.PreventUpdate
    triggered_id = ctx.triggered_id
    if triggered_id not in ["btn-load-annotations-with-link"]:
        raise exceptions.PreventUpdate
    triggered_data = {"triggered_id": triggered_id}

    if data is None:
        data = {}

    data['annotation_checkbox_values'] = annotation_checkbox_values
    data['annotation_checkbox_ids'] = annotation_checkbox_ids
    data['annotation_dropdown_values'] = annotation_dropdown_values
    data['annotation_dropdown_ids'] = annotation_dropdown_ids
    data["load_annotations_status"] = "running"
    data["load_annotations_message"] = None
    return "", data, triggered_data, ""


#callback to load annotations
@app.callback(
    [Output('db-management-page-store', 'data', allow_duplicate=True),
     Output("global-notification", "data", allow_duplicate=True),
     Output("db-management-load-annotations-trigger", "data")],
    Input("db-management-load-annotations-trigger", "data"),
    State('db-management-page-store', 'data'),
    State({'type': 'annotation-checkbox', 'index': ALL}, 'value'),
    State({'type': 'annotation-checkbox', 'index': ALL}, 'id'),
    State({'type': 'annotation-dropdown', 'index': ALL}, 'value'),
    State({'type': 'annotation-dropdown', 'index': ALL}, 'id'),
    # running=[
    #     (Output("btn-load-annotations-with-link", "disabled"), True, False)
    # ],
    cancel=[Input("btn-cancel-load-annotations", "n_clicks")],
    background=True,
    suppress_callback_exceptions=True,
    prevent_initial_call=True
)
@require_authorization
def load_annotations_launch(trigger_data, data, checkbox_values, checkbox_ids, dropdown_values, dropdown_ids):
    triggered_id = trigger_data.get("triggered_id", None)
    if not triggered_id or not trigger_data :
        raise exceptions.PreventUpdate
    if data is None:
        data = {}
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
        toast_message = {
            "title": "Annotation loading error.",
            "message": "❌ Please select at least one annotation file before loading.",
            "type": "danger",
        }
        data["load_annotations_message"] = "❌ Please select at least one annotation file before loading."
        data["load_annotations_message_style"] = error_style
        data["load_annotations_status"] = "done"
        return data, toast_message, {}

    # Ensure each selected file has an associated genome
    for f in selected_files:
        genome = genome_mapping.get(f)
        if not genome:
            toast_message = {
                "title": "Annotation loading error.",
                "message": f"❌ Please select a genome for annotation file '{f}'.",
                "type": "danger",
            }
            data["load_annotations_message"] = f"❌ Please select a genome for annotation file '{f}'."
            data["load_annotations_message_style"] = error_style
            data["load_annotations_status"] = "done"
            return data, toast_message, {}

    # Check index state for each genome
    for genome in set(genome_mapping.values()):
        state_index = check_state_index("NodeIndex" + genome + "_position")
        if state_index is None:
            toast_message = {
                "title": "Annotation loading error.",
                "message": f"❌ Index for genome '{genome}' has not been created. Please create the index before loading annotations.",
                "type": "danger",
            }
            data["load_annotations_message"] = f"❌ Index for genome '{genome}' has not been created. Please create the index before loading annotations."
            data["load_annotations_message_style"] = error_style
            data["load_annotations_status"] = "done"
            return data, toast_message, {}
        if int(state_index) != 100:
            toast_message = {
                "title": "Annotation loading error.",
                "message": f"❌ Index for genome '{genome}' is not completely created (state: {state_index}%). Please wait.",
                "type": "danger",
            }
            data["load_annotations_message"] = f"❌ Index for genome '{genome}' is not completely created (state: {state_index}%). Please wait."
            data["load_annotations_message_style"] = error_style
            data["load_annotations_status"] = "done"
            return data, toast_message, {}

    # Validate file extensions
    invalid_files = [f for f in selected_files if not f.lower().endswith((".gff", ".gff3", ".gtf"))]
    if invalid_files:
        toast_message = {
            "title": "Annotation loading error.",
            "message": f"❌ Invalid file(s): {', '.join(invalid_files)}. Please select only .gff, .gff3, or .gtf files.",
            "type": "danger",
        }
        data["load_annotations_message"] = f"❌ Invalid file(s): {', '.join(invalid_files)}. Please select only .gff, .gff3, or .gtf files."
        data["load_annotations_message_style"] = error_style
        data["load_annotations_status"] = "done"
        return data, toast_message, {}

    # Loop through each selected file and load annotations
    for f in selected_files:
        genome = genome_mapping[f]
        file_path = os.path.join(ANNOTATIONS_FOLDER, f)
        logger.info(f"Load annotations for file '{f}' with genome '{genome}'")
        load_annotations_neo4j(file_path, genome_ref=genome, single_chromosome=None)

    # Optionally link annotations if the "with-link" button was clicked
    WARN = None
    #if triggered_id == "btn-load-annotations-with-link":
    for genome in set(genome_mapping.values()):
        logger.info(f"Link annotations for genome '{genome}'")
        WARN = creer_relations_annotations_neo4j(genome)

    if WARN:
        toast_message = {
            "title": "Annotation loading error.",
            "message": f"❌ Chromosome names mismatch between graph and annotation file. No annotation linked.",
            "type": "warning",
        }
        data["load_annotations_message"] = f"❌ Chromosome names mismatch between graph and annotation file. No annotation linked."
        data["load_annotations_message_style"] = warning_style
        data["load_annotations_status"] = "done"
        return data, toast_message, {}
    else:
        data['annotation_checkbox_values'] = [[] for _ in checkbox_values]
        toast_message = {
            "title": "Annotation loading.",
            "message": f"✅ Annotation(s) '{', '.join(selected_files)}' loaded with associated genomes.",
            "type": "success"
        }
        data["load_annotations_message"] = f"✅ Annotation(s) '{', '.join(selected_files)}' loaded with associated genomes."
        data["load_annotations_message_style"] = success_style
        data["load_annotations_status"] = "success"
        return data, toast_message, {}

    # Return success message and reset all checkboxes (empty list)
    data['annotation_checkbox_values'] = [[] for _ in checkbox_values]
    toast_message = {
        "title": "Annotation loading.",
        "message": f"✅ Annotation(s) '{', '.join(selected_files)}' loaded with associated genomes.",
        "type": "success"
    }
    data["load_annotations_message"] = f"✅ Annotation(s) '{', '.join(selected_files)}' loaded with associated genomes."
    data["load_annotations_message_style"] = success_style
    data["load_annotations_status"] = "success"
    return data, toast_message, {}


#Cancel annotation load button
@app.callback(
    Output('db-management-page-store', 'data', allow_duplicate=True),
    Input("btn-cancel-load-annotations", "n_clicks"),
    State('db-management-page-store', 'data'),
    prevent_initial_call=True,
)
def handle_cancel_annotation_load_click(n_clicks, data):
    if not n_clicks:
        return no_update

    data["load_annotations_message"] = (f"✅ Annotation loading cancelled.")
    data["load_annotations_message_style"] = success_style
    data["load_annotations_status"] = "done"
    return data


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
    Output("delete-message", "children", allow_duplicate=True),
    Output("delete-confirmation", "children", allow_duplicate=True),
    Output('db-management-page-store', 'data', allow_duplicate=True),
    Output("annotation-message", "children", allow_duplicate=True),
    Output("create-db-message", "children", allow_duplicate=True),
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
    data = delete_messages(data)
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
        if data is not None and "container_name" in data:
            data.pop("container_name")
        return html.Div("✅ All data deleted successfully.", style=success_style), "", data, "", ""
    except Exception as e:
        return html.Div(f"❌ Error while deleting data: {str(e)}", style=error_style), "", data, "", ""


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
    Output("delete-annotations-message", "children", allow_duplicate=True),
    Output("delete-annotations-confirmation", "children", allow_duplicate=True),
    Output("annotation-message", "children", allow_duplicate=True),
    Output('db-management-page-store', 'data', allow_duplicate=True),
    Input("btn-delete-annotations-confirmation", "n_clicks"),
    State('db-management-page-store', 'data'),
    running=[
        (Output("btn-delete-annotations", "disabled"), True, False),
        (Output("btn-delete-annotations-confirmation", "disabled"), True, False)
    ],
    prevent_initial_call=True
)
@require_authorization
def confirm_delete_annotations(n_clicks, data):
    data = delete_messages(data)
    if not n_clicks:
        raise exceptions.PreventUpdate
    try:
        delete_annotations()
        return html.Div("✅ All annotations deleted successfully.", style=success_style), "", "", data
    except Exception as e:
        return html.Div(f"❌ Error while deleting annotations: {str(e)}", style=error_style), "", "", data



############# Container callbacks#################


#Callback to create confirmation button after clicking on create db button
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


#Callback called when confirming database creation
@app.callback(
    Output("create-db-confirmation", "children", allow_duplicate=True),
    Output("db-management-page-store", "data", allow_duplicate=True),
    Output("db-management-create-db-trigger", "data", allow_duplicate=True),
    Output("delete-message", "children", allow_duplicate=True),
    Input("btn-confirm-create-db", "n_clicks"),
    State("container-name-input","value"),
    State("docker-image-dropdown","value"),
    State('db-management-page-store', 'data'),
    State('annotations-files-container', 'children'),
    State({'type': 'gfa-checkbox', 'index': ALL}, 'value'),
    State({'type': 'gfa-checkbox', 'index': ALL}, 'id'),
    State({'type': 'gfa-input', 'index': ALL}, 'value'),
    State({'type': 'gfa-input', 'index': ALL}, 'id'),
    prevent_initial_call=True
)
@require_authorization
def confirm_create_db(n_clicks, container_name, docker_image, data, children, checkbox_values, checkbox_ids, input_values, input_ids):
    if not n_clicks:
        raise exceptions.PreventUpdate
    triggered_id = ctx.triggered_id
    if triggered_id not in ["btn-confirm-create-db"]:
        raise exceptions.PreventUpdate
    triggered_data = {"triggered_id":triggered_id}

    # Check container name
    if not container_name:
        data["db_creation_message"] = "❌ No container name"
        data["db_creation_message_style"] = error_style
        return "",data, no_update, ""
    else:
        container_name_prefixed = PREFIX_CONTAINER_NAME + container_name

    if data is None:
        data = {}

    data['container_name'] = get_container_name_no_prefix(container_name)
    data['gfa_checkbox_values'] = checkbox_values
    data['gfa_checkbox_ids'] = checkbox_ids
    data['gfa_input_values'] = input_values
    data['gfa_input_ids'] = input_ids
    data['docker_image'] = docker_image
    data["create_db_status"]="running"
    data["db_creation_message"] = None
    return "", data,triggered_data, ""



@app.callback(
    [Output('db-management-page-store', 'data', allow_duplicate=True),
    Output("global-notification", "data", allow_duplicate=True),
    Output("db-management-create-db-trigger", "data")],
    Input("db-management-create-db-trigger", "data"),
    State('db-management-page-store', 'data'),
    State('annotations-files-container', 'children'),
    cancel=[Input("btn-cancel-create-db", "n_clicks")],
    background=True,
    suppress_callback_exceptions=True,
    prevent_initial_call=True
)
@require_authorization
def create_db_launch(trigger_data, data, children):
    triggered_id = trigger_data.get("triggered_id", None)
    if not triggered_id or not trigger_data:
        raise exceptions.PreventUpdate
    if data is None:
        data = {}
    n_dropdowns = len(children) if children else 0
    data['options_list'] = [no_update] * n_dropdowns


    container_name = data['container_name']
    checkbox_values = data['gfa_checkbox_values']
    checkbox_ids = data['gfa_checkbox_ids']
    input_values = data['gfa_input_values']
    input_ids = data['gfa_input_ids']
    docker_image = data['docker_image']

    # Check container name
    if not container_name:
        data["create_db_status"] = "done"
        data["db_creation_message"] = "❌ No container name"
        data["db_creation_message_style"] = error_style
        toast_message = {
            "title": "Database creation error.",
            "message": "❌ No container name",
            "type": "danger",
        }
        return data, toast_message, {}
    else:
        container_name_prefixed = PREFIX_CONTAINER_NAME + container_name
    try:
        genomes_set = set()
        chromosomes_stats = {}
        #Check if gfa files are selected
        selected_files = [c_id['index'] for c_val, c_id in zip(checkbox_values, checkbox_ids) if c_val]
        if selected_files:

            invalid_files = [f for f in selected_files if not f.lower().endswith(".gfa")]
            if invalid_files:
                data["create_db_status"] = "done"
                data["db_creation_message"] = f"❌ Invalid file(s): {', '.join(invalid_files)}. Please select only .gfa files."
                data["db_creation_message_style"] = error_style
                toast_message = {
                    "title": "Database creation error.",
                    "message": f"❌ Invalid file(s): {', '.join(invalid_files)}. Please select only .gfa files.",
                    "type": "danger",
                }
                return data, toast_message, {}

            # Get batch size from conf
            batch_size = get_db_load_gfa_batch_size()
            logger.info(f"Generate csv from GFA files {selected_files} with batch size {batch_size}")
            # Get the chromosome associated to file (if set)
            chromosome_dict = {i_id['index']: val for i_id, val in zip(input_ids, input_values)}
            list_chromosome_file = [chromosome_dict.get(f, "") for f in selected_files]
            if len(list_chromosome_file) > 1:
                for cf in list_chromosome_file:
                    if cf is None or cf == "":
                        data["create_db_status"] = "done"
                        data["db_creation_message"] = "❌ When multiple gfa are selected, it is required to set the chromosome for each of theses files (non null or empty value)."
                        data["db_creation_message_style"] = error_style
                        toast_message = {
                            "title": "Database creation error.",
                            "message": "❌ When multiple gfa are selected, it is required to set the chromosome for each of theses files (non null or empty value).",
                            "type": "danger",
                        }
                        return data, toast_message, {}
            for file_name, chromosome_file in zip(selected_files, list_chromosome_file):
                start_time = time.time()
                if chromosome_file != "":
                    chrom = chromosome_file
                else:
                    chrom = None
                file_path = os.path.join(GFA_FOLDER, file_name)
                genomes_analysed, chromosomes_analysed = load_gfa_data_to_csv(file_path, import_dir="./data/import",
                                     chromosome_file=chrom,
                                     chromosome_prefix=False,
                                     batch_size=batch_size,
                                     start_chromosome=None,
                                     haplotype=True)
                genomes_set = genomes_set | genomes_analysed
                chromosomes_stats = chromosomes_stats | chromosomes_analysed
                logger.info(f"CSV generation from {file_name} loaded in {time.time() - start_time:.2f} s")
        logger.info("All import files have been generated from gfa files in {time.time() - start_time:.2f} s, creating database.")
        creation_mode = create_db(container_name_prefixed, docker_image)
        # If creation by importing csv files it is necessary to create stats and indexes
        if creation_mode == "csv":
            stats = False
            if len(genomes_set) > 0 and len (chromosomes_stats) > 0:
                logger.info("creating stats")
                create_stats(genomes_set, chromosomes_stats)
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
        data["create_db_status"] = "success"
        data["db_creation_message"] = msg
        data["db_creation_message_style"] = success_style
        data['gfa_checkbox_values'] = [[] for _ in checkbox_values]
        data['options_list'] = options_list
        toast_message = {
            "title": "Database creation.",
            "message": msg,
            "type": "success",
        }
        return  data, toast_message, {}

    except Exception as e:
        logger.error(f"Error while creating database: {e}")
        data["create_db_status"] = "done"
        data["db_creation_message"] = f"❌ Error while creating database: {str(e)}"
        data["db_creation_message_style"] = error_style
        toast_message = {
            "title": "Database creation error.",
            "message": f"❌ Error while creating database: {str(e)}",
            "type": "danger",
        }
        return data, toast_message, {}


#Cancel database creation button
@app.callback(
    Output('db-management-page-store', 'data', allow_duplicate=True),
    Input("btn-cancel-create-db", "n_clicks"),
    State('db-management-page-store', 'data'),
    prevent_initial_call=True,
)
def handle_cancel_db_creation_click(n_clicks, data):
    if not n_clicks:
        return no_update
    logger.info(f"Deleting import data : {IMPORT_FOLDER} directory.")
    try:
        if os.path.exists(IMPORT_FOLDER):
            shutil.rmtree(IMPORT_FOLDER)
            os.makedirs(IMPORT_FOLDER)
        data["db_creation_message"] = (f"✅ DB creation cancelled. Warning : Depending on when the operation was stopped, "
                                       f"some data may have been created in the database. In this case, a full data reset may be required.")
        data["db_creation_message_style"] = success_style
    except Exception as e:
        data["db_creation_message"] = f"❌ Error while deleting data: {str(e)}"
        data["db_creation_message_style"] = error_style
    data["create_db_status"] = "done"
    return data

#Dump database callback
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


#Callback to refresh page after navigation or when data store has been changed (asynchron process)
@app.callback(
    Output('container-name-label', 'children'),
    Output('container-name-input', 'value'),
    #Output('db-management-page-store', 'data'),
    Output("container-name-input", "disabled"),
    Output("btn-create-db", "style"),
    Output("btn-load-gfa", "style"),
    Output("btn-create-index", "disabled"),
    Output("btn-create-stats", "disabled"),
    Output("btn-dump-db", "disabled"),
    Output("create-db-message", "children"),
    Output({'type': 'gfa-checkbox', 'index': ALL}, 'value'),
    Output("btn-create-db", "disabled"),
    Output("btn-load-gfa", "disabled"),
    Output("btn-cancel-create-db", "disabled"),
    Output("db-create-spinner", "style"),
    Output({'type': 'gfa-input', 'index': ALL}, 'value'),
    Output({'type': 'annotation-dropdown', 'index': ALL}, 'options'),
    Output({'type': 'annotation-checkbox', 'index': ALL}, 'value'),
    Output({'type': 'annotation-dropdown', 'index': ALL}, 'value'),
    Output("annotation-message", "children"),
    Output("btn-cancel-load-annotations", "disabled"),
    Output("db-load-annotation-spinner", "style"),
    Output("btn-load-annotations-with-link", "disabled"),
    Input('url', 'pathname'),
    Input('db-management-page-store', 'data'),
    State('container-name-input', 'value'),
    State({'type': 'gfa-checkbox', 'index': ALL}, 'value'),
    State({'type': 'annotation-checkbox', 'index': ALL}, 'value'),
    State({'type': 'annotation-dropdown', 'index': ALL}, 'value'),

    prevent_initial_call=False
)
@require_authorization
def update_page(pathname, data, container_input, checkbox_values, annotation_checkbox_values, annotations_ref_genome_values):
    if pathname != "/db_management":
        raise exceptions.PreventUpdate
    if data is None:
        data = {}

    #Part 1 : Management of database creation parts
    status_create_db = data.get("create_db_status", "done")
    create_db_message = ""
    if "db_creation_message" in data and data["db_creation_message"] is not None and data["db_creation_message"] != "":
        create_db_message = html.Div(data.get("db_creation_message", ""), style=data.get("db_creation_message_style", success_style))
    gfa_checkbox_values = data.get('gfa_checkbox_values',[[] for _ in checkbox_values])
    input_values = data.get('gfa_input_values', ["" for _ in checkbox_values])

    button_create_db = False
    button_load_gfa = False
    create_db_spinner_style = {"display": "none", "marginTop": "20px"}

    if status_create_db == "running":
        create_db_spinner_style = {"display": "block", "marginTop": "20px"}
        button_cancel = False
        button_create_db = True
        button_load_gfa = True
        if create_db_message == "" :
            create_db_message = html.Div("Creating database...", style=success_style)
    else:
        button_cancel = True


    #Part 2 : Management of annotations part
    annotation_checkbox_values = data.get('annotation_checkbox_values', [[] for _ in annotation_checkbox_values])
    annotations_ref_genome_values = data.get('annotation_dropdown_values', [[] for _ in annotations_ref_genome_values])
    annotation_message = ""
    if "load_annotations_message" in data and data["load_annotations_message"] is not None and data["load_annotations_message"] != "":
        annotation_message = html.Div(data.get("load_annotations_message", ""),
                                     style=data.get("load_annotations_message_style", success_style))

    status_load_annotations = data.get("load_annotations_status","done")
    create_load_annotation_spinner_style = {"display": "none", "marginTop": "20px"}
    button_cancel_annotation_disabled = True
    button_load_annotation_disabled = False
    if status_load_annotations == "running":
        create_load_annotation_spinner_style = {"display": "block", "marginTop": "20px"}
        button_cancel_annotation_disabled = False
        button_load_annotation_disabled = True
        if annotation_message == "":
            annotation_message = html.Div("Loading annotations...", style=success_style)
    style_create_db = {"display": "inline-block"}
    style_add_gfa = {"display": "none"}

    #Part 3 : getting conf file to check if the database already exists
    conf = load_config_from_json()
    style_create_db = {"display": "none"}
    style_add_gfa = {"display": "inline-block"}
    btn_create_index_disabled = False
    btn_create_state_disabled = False
    btn_dump_db_disabled = False
    container_name_input_disabled = True
    btn_load_annotations_disabled = False
    if not conf or "container_name" not in conf or conf['container_name'] is None or conf['container_name'] == "":
        style_create_db = {"display": "inline-block"}
        style_add_gfa = {"display": "none"}
        btn_create_index_disabled = True
        btn_create_state_disabled = True
        btn_dump_db_disabled = True
        container_name_input_disabled = False
        btn_load_annotations_disabled = True

        container_name = data.get("container_name", container_input if (container_input and container_input != "") else "container_name")
        msg_container = f'No conf file found. Use "create new DB" procedure to generate it.'
        annotations_option_list = []
    else:
        container_name = get_container_name_no_prefix(conf.get("container_name", "container_name"))
        msg_container = f"Container name: {container_name}"
    n_annotations_dropdown = len(annotations_ref_genome_values)
    if ("options_list" in data and isinstance(data["options_list"], list)
            and len(data["options_list"]) == n_annotations_dropdown):
        annotations_option_list = data["options_list"]
    else:
        genomes = get_genomes() or []
        annotations_option_list = [
            [{"label": g, "value": g} for g in genomes] if genomes else []
            for _ in range(n_annotations_dropdown)
        ]
        data["options_list"] = annotations_option_list


    return (msg_container, container_name,container_name_input_disabled,
                        style_create_db,style_add_gfa, btn_create_index_disabled,
                        btn_create_state_disabled, btn_dump_db_disabled,
                        create_db_message,  gfa_checkbox_values, button_create_db, button_load_gfa, button_cancel,
                        create_db_spinner_style, input_values, annotations_option_list, annotation_checkbox_values,
                        annotations_ref_genome_values,annotation_message,
                        button_cancel_annotation_disabled, create_load_annotation_spinner_style,
                        button_load_annotation_disabled)



# @app.callback(
#     Output('container-name-label', 'children'),
#     Output('container-name-input', 'value'),
#     #Output('db-management-page-store', 'data'),
#     Output("container-name-input", "disabled"),
#     Output("btn-load-annotations-with-link", "disabled"),
#     Output("btn-create-db", "style"),
#     Output("btn-load-gfa", "style"),
#     Output("btn-create-index", "disabled"),
#     Output("btn-create-stats", "disabled"),
#     Output("btn-dump-db", "disabled"),
#     Output("create-db-message", "children"),
#     Output({'type': 'gfa-checkbox', 'index': ALL}, 'value'),
#     Output("btn-create-db", "disabled"),
#     Output("btn-load-gfa", "disabled"),
#     Output("btn-cancel-create-db", "disabled"),
#     Output("db-create-spinner", "style"),
#     Output({'type': 'annotation-dropdown', 'index': ALL}, 'options'),
#     Output("annotation-message", "children"),
#     Input('url', 'pathname'),
#     Input('db-management-page-store', 'data'),
#     State('container-name-input', 'value'),
#     prevent_initial_call=False
# )
# @require_authorization
# def update_page(pathname, data, container_input, checkbox_values):
#     if pathname != "/db_management":
#         raise exceptions.PreventUpdate
#     if data is None:
#         data = {}
#     status_create_db = data.get("create_db_status", "done")
#     create_db_message = ""
#     if "db_creation_message" in data and data["db_creation_message"] is not None and data["db_creation_message"] != "":
#         create_db_message = html.Div(data.get("db_creation_message", ""), style=data.get("db_creation_message_style", success_style))
#     annotations_option_list = data.get('options_list',[])
#     gfa_checkbox_values = data.get('gfa_checkbox_values',[[] for _ in checkbox_values])
#
#
#     button_create_db = False
#     button_load_gfa = False
#     create_db_spinner_style = {"display": "none", "marginTop": "20px"}
#
#     if status_create_db == "running":
#         create_db_spinner_style = {"display": "block", "marginTop": "20px"}
#         button_cancel = False
#         button_create_db = True
#         button_load_gfa = True
#         if create_db_message == "" :
#             create_db_message = html.Div("Creating database...", style=success_style)
#     else:
#         button_cancel = True
#
#     style_create_db = {"display": "inline-block"}
#     style_add_gfa = {"display": "none"}
#
#     conf = load_config_from_json()
#     style_create_db = {"display": "none"}
#     style_add_gfa = {"display": "inline-block"}
#     btn_create_index_disabled = False
#     btn_create_state_disabled = False
#     btn_dump_db_disabled = False
#     container_name_input_disabled = True
#     btn_load_annotations_disabled = False
#     if not conf or "container_name" not in conf or conf['container_name'] is None or conf['container_name'] == "":
#         style_create_db = {"display": "inline-block"}
#         style_add_gfa = {"display": "none"}
#         btn_create_index_disabled = True
#         btn_create_state_disabled = True
#         btn_dump_db_disabled = True
#         container_name_input_disabled = False
#         btn_load_annotations_disabled = True
#
#         container_name = data.get("container_name", container_input if (container_input and container_input != "") else "container_name")
#         msg_container = f'No conf file found. Use "create new DB" procedure to generate it.'
#
#     else:
#         container_name = get_container_name_no_prefix(conf.get("container_name", "container_name"))
#         #data['container_name'] = get_container_name_no_prefix(container_name)
#         #container_name = PREFIX_CONTAINER_NAME + data['container_name']
#         #get_container_name_no_prefix(container_name)
#         msg_container = f"Container name: {container_name}"
#
#     return (f'No conf file found. Use "create new DB" procedure to generate it.', container_name, container_name_input_disabled,
#                         btn_load_annotations_disabled, style_create_db, style_add_gfa, btn_create_index_disabled,
#                         btn_create_state_disabled, btn_dump_db_disabled,
#                         create_db_message,  gfa_checkbox_values, button_create_db, button_load_gfa, button_cancel,
#                         create_db_spinner_style, annotations_option_list)

