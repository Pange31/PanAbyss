#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 21/08/2025

@author: fgraziani
"""

import dash
from dash import html, dcc, Input, Output, State, ctx
import dash_cytoscape as cyto

from neo4j_requests import *
from app import DB_VERSION

from neo4j_available_docker_images_conf import AVAILABLE_DOCKER_IMAGES



PREFIX_CONTAINER_NAME = "DB_"+ DB_VERSION + "_"
ANNOTATION_DIR = './data/annotations'
GFA_DIR = './data/gfa'


def list_annotation_files():
    """Return all .gff / .gtf / .gff3 files in the annotation directory."""
    return sorted([
        f for f in os.listdir(ANNOTATION_DIR)
        if f.lower().endswith(('.gff', '.gtf', '.gff3'))
    ])

def list_gfa_files():
    """Return all .gfa files in the annotation directory."""
    return sorted([
        f for f in os.listdir(GFA_DIR)
        if f.lower().endswith(('.gfa'))
    ])


# Layout de la page
def layout():
    genomes = get_genomes()
    if genomes is None :
        genomes = []
    return html.Div([
        html.H2("DB Management"),
        #Help section
        html.Details([
            html.Summary("ℹ️ Click here to display help"),
            html.Ul([
                html.Li("This page allows to manage database. Once the data have been created, this page is normally only used to reset the database. The creation of the database depends on the data :"),
                    html.Ul([
                        html.Li("Initial procedure :"),
                        html.Ul([
                            html.Li("Set the name of the neo4J container."),
                            html.Li("Select GFA files and set the associated chromosome in the following cases : "),
                            html.Ul([
                            html.Li("The file concern a unique chromosome."),
                            html.Li("Multiple gfa files are loaded, each file must refer ton only one chromosome."),
                            html.Li("No chromosome at all (set the value to 0 in this case)."),
                            html.Li("If there is only one GFA file with multiple chromosomes then there is nothing to set."),
                            ]),
                            html.Li("Click on Create new DB button."),
                        ]),
                        html.Li("Procedure to load further GFA files : "),
                        html.Ul([
                            html.Li(
                                "Once the database has been created, it is no longer possible to use the CSV generation process. In this case, you can use the **'Add data'** procedure by selecting the file and specifying the associated chromosome. However, this approach is **not recommended** for big gfa files, as it is slower than the initial loading process."),
                        ]),
                        html.Li("Dump procedure (for intermediate data volume): "),  
                        html.Ul([
                            html.Li("If a dump file has been generated (neo4j.dump file) just check if this file is into the /data/import directory, enter the neo4j container name and click on 'create new db' button (don't select gfa files)."),
                        ]),
                        html.Li("Load annotations files : "),  
                        html.Ul([
                            html.Li(
                                "To launch this operation the database must be created (else the button is disable)."),
                            html.Li("For this step it is necessary that indexes are created. It can take a few time after database creation."),
                            html.Li("For each annotation file to be loaded, select the file and specify the individual associated with that annotation file. Each annotation file must be linked to **only one** individual. Then, click the **'Load'** button."),
                            
                        ]),
                    ])
            ])
            ], style={"marginBottom": "20px"}),
    
        html.H3(id='container-name-label', children=""),


        html.Hr(),
        #############################" ---- Database creation Section ---- #####################################"

        html.H3("Create new Database"),
        html.Hr(),
        html.Div([
            html.H4("Important notes:"),
            html.Ul([
                html.Li("This operation allow to create a new database with gfa files."),
                html.Li("It is required to set a container name (name of the neo4j docker container)."),
                html.Li("If no gfa files selected, the it will use, dump (in priority) files or CSV files in ./data/import. If no data in /data/import then an empty database will be created."),
            ])
        ]),
        html.Br(),
        html.Div([
            html.Label("📦 Neo4j container name :  ", title='Set a name for your docker container.', ),
            html.Span(PREFIX_CONTAINER_NAME, style={
                'fontWeight': 'bold',
                'paddingRight': '5px'
            }),
            dcc.Input(
                id='container-name-input',
                type='text',
                debounce=True,
                style={'width': '400px'}
            )
        ], style={'marginBottom': '20px'}),
        html.Label("Select docker image to use :"),
        html.Div([
            dcc.Dropdown(
                id='docker-image-dropdown',
                options=[{'label': img, 'value': img} for img in AVAILABLE_DOCKER_IMAGES],
                value=AVAILABLE_DOCKER_IMAGES[0],
                clearable=False,
                style={'width': '400px'}
            ),
        ], style={'marginBottom': '20px'}),

        html.Br(),
        html.Div([
            html.H4("Instructions for loading a GFA files :"),
            html.Ul([
                html.Li("This step allow to create database for the first time."),
                html.Li("GFA files must be in the ./data/gfa directory. Small GFA can be directly uploaded by the drag/drop box."),
                html.Li("Once the database is created, it is necessary to reset all to recreate it."),
                html.Li("For multi GFA processing : each GFA must refer to a unique chromosome and the chromosome value must be set in the input box associated to each gfa file. If not set the value in the path / walk will be used."),
            ])
        ]),
        html.Br(),
        html.H4("GFA files loading (only for small gfa, big must be placed directly into the /data/gfa directory)"),
    
        dcc.Upload(
            id='upload-gfa-data',
            children=html.Div([
                'Dropdown or select ',
                html.A('GFA files')
            ]),
            style={
                'width': '50%',
                'height': '100px',
                'lineHeight': '100px',
                'borderWidth': '2px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px'
            },
            multiple=True
        ),
        dcc.Loading(
            id="loading-gfa-upload",
            #type="default",
            children=html.Div(id='upload-gfa-output'),
        ),

        html.Br(),

        html.Div([
            html.Div("GFA file", style={'fontWeight': 'bold', 'wordBreak': 'break-word'}),
            html.Div("Chromosome (if applicable)", title="If GFA concern only one chromosome, or if no chromosome, specify the chromosome value here (0 if no chromosome) or let this value unset if multiple chromosomes.", style={'fontWeight': 'bold'})
        ],
            style={
                'display': 'grid',
                'gridTemplateColumns': '350px 230px',
                'alignItems': 'center',
                'columnGap': '10px',
                'marginBottom': '8px'
            }),

        html.Div(
            id='gfa-files-container',
            children=[
                html.Div([
                    html.Div(
                        dcc.Checklist(
                            id={'type': 'gfa-checkbox', 'index': f},
                            options=[{'label': f, 'value': f}],
                            value=[],
                            labelStyle={
                                'display': 'inline-block',
                                'whiteSpace': 'normal',
                                'wordBreak': 'break-word',
                                'width': '100%'
                            }
                        )
                    ),

                    html.Div(
                        dcc.Input(
                            id={'type': 'gfa-input', 'index': f},
                            type='text',
                            placeholder='Enter chromosome (optional)',

                            style={'width': '100%'}
                        )
                    )
                ],
                    style={
                        'display': 'grid',
                        'gridTemplateColumns': '350px 230px',
                        'alignItems': 'center',
                        'columnGap': '10px',
                        'marginBottom': '6px'
                    })
                for f in list_gfa_files()
            ]
        ),


        html.Br(),

        html.Div(
            [
                html.P(
                    "⚠️ Warning: If your GFA file concerns only one chromosome, it is strongly recommended to enter the chromosome number.",
                    style={"margin": "4px 0"}
                ),
                html.P(
                    "If it does not correspond to any chromosome, enter 0.",
                    style={"margin": "4px 0"}
                ),
                html.P(
                    "If it contains multiple chromosomes (in this case only one file can be used), leave the field empty.",
                    style={"margin": "4px 0"}
                ),
            ],
            style={
                "border": "2px solid #f5c518",
                "backgroundColor": "#fff8e1",
                "padding": "10px",
                "borderRadius": "5px",
                "fontWeight": "bold",
                "color": "#8a6d3b",
                "fontSize": "14px",
                "marginBottom": "15px",
                "maxWidth": "50%"
            }
        ),

        html.Div([
            html.Button(
                "Create new DB",
                id="btn-create-db",
                n_clicks=0,
                style={"display": "inline-block", "marginRight": "15px"}
            ),
            html.Button(
                "Add data (DB must be created)",
                id="btn-load-gfa",
                n_clicks=0,
                style={"display": "none", "marginRight": "15px"}
            ),
            html.Button(
                "Cancel",
                id="btn-cancel-create-db",
                n_clicks=0,
                disabled=True,
                style={"display": "inline-block", "marginLeft": "15px"}
            )
        ], id="create-db-button-container"),
        html.Div(id="create-db-confirmation", style={"marginTop": "10px"}),
        html.Div(
            html.Div(className="custom-spinner"),
            id="db-create-spinner",
            style={"display": "none", "marginTop": "20px", "marginBottom": "20px"}
        ),
        # dcc.Loading(
        #     # type="circle",
        #     children=html.Div(id="create-db-message", style={"marginTop": "10px"})
        # ),
        html.Div(id="create-db-message"),
        dcc.Loading(
            # type="circle",
            children=html.Div(id="add-gfa-message", style={"marginTop": "10px"})
        ),


        ##################### ---- Annotation Upload Section ---- #################################"
        html.Br(),
        html.Hr(),
        html.H3("Annotations"),
        html.Hr(),

        html.Div([
            html.H4("Instructions for loading annotations:"),
            html.Ul([
                html.Li("Load annotation files – files must be located in the ./data/annotations directory."),
                html.Li("The database must already be loaded."),
                html.Li(
                    "Required indexes must be created and available in the database before loading annotations."), ])
        ]),
        html.Br(),
        html.H4(
            "Annotations files loading (only for small annotation files, big must be placed directly into the /data/annotations directory)"),

        dcc.Upload(
            id='upload-annotations-data',
            children=html.Div([
                'Dropdown or select ',
                html.A('GFF or GTF files')
            ]),
            style={
                'width': '50%',
                'height': '100px',
                'lineHeight': '100px',
                'borderWidth': '2px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px'
            },
            multiple=True
        ),
        dcc.Loading(
            id="loading-annotations-upload",
            # type="default",
            children=html.Div(id='upload-annotations-output'),
        ),

        html.Br(),
        html.Div([
            html.Div("Annotation file", style={'fontWeight': 'bold', 'wordBreak': 'break-word'}),
            html.Div("Associated genome", style={'fontWeight': 'bold'})
        ],
            style={
                'display': 'grid',
                'gridTemplateColumns': '350px 230px',
                'alignItems': 'center',
                'columnGap': '10px',
                'marginBottom': '8px'
            }),

        html.Div(
            id='annotations-files-container',
            children=[
                html.Div([

                    html.Div(
                        dcc.Checklist(
                            id={'type': 'annotation-checkbox', 'index': f},
                            options=[{'label': f, 'value': f}],
                            value=[],
                            labelStyle={
                                'display': 'inline-block',
                                'whiteSpace': 'normal',
                                'wordBreak': 'break-word',
                                'width': '100%'
                            }
                        )
                    ),

                    html.Div(
                        dcc.Dropdown(
                            id={'type': 'annotation-dropdown', 'index': f},
                            options=[{"label": genome, "value": genome} for genome in genomes],
                            placeholder="Select genome",
                            style={'width': '100%'}
                        )
                    )
                ],
                    style={
                        'display': 'grid',
                        'gridTemplateColumns': '350px 230px',
                        'alignItems': 'center',
                        'columnGap': '10px',
                        'marginBottom': '6px'
                    })
                for f in list_annotation_files()
            ]
        ),

        html.Br(),

        html.Div([
            html.Button("Load",
                        title="This will load annotations into database. Indexes must be created before (can take some time for big data).",
                        id="btn-load-annotations-with-link", n_clicks=0),
            html.Button("Cancel",
                        title="This will cancel load annotations operation.",
                        id="btn-cancel-load-annotations", disabled=True, n_clicks=0),

            html.Div(
                id="db-load-annotation-spinner",
                style={"display": "none", "marginTop": "20px", "marginBottom": "20px"},
                children=html.Div(className="custom-spinner"),
            ),

            html.Div(id="annotation-message"),
            # html.Button("Load only annotations", id="btn-load-only-annotations", n_clicks=0),
            # html.Button("Link annotations", id="btn-link-annotations", n_clicks=0),
        ]),

        ################ Database management #############################
        html.Br(),
        html.Hr(),
        html.H3("Management of DB"),
        html.Hr(),
        html.P("These operations enable specific management tasks:"),
        html.Ul([
            html.Li("Reset"),
            html.Li("CSV generation from GFAs"),
            html.Li("Creation of indexes and statistics in the database"),
            html.Li("Dump creation")
        ]),
        html.P("Normally, these operations do not need to be used since they are automatically launched when creating a new database."),
        html.Br(),
        html.Label("Reset all the database and configuration : this will delete data, configuration and neo4j container.",
                   style={'display': 'block', 'marginBottom': '8px'}),
        html.Button("Reset all", title="This will delete all data, conf and container.", id="btn-delete",
                    n_clicks=0),
        # Confirm deletion
        html.Div(id="delete-confirmation", style={"marginTop": "10px"}),
        dcc.Loading(
            # type="circle",
            children=html.Div(id="delete-message", style={"marginTop": "10px"})
        ),
        html.Br(),
        html.Label(
            "Reset all the annotations in the database.",
            style={'display': 'block', 'marginBottom': '8px'}),
        html.Button("Delete annotations", title="This will delete all annotations in database.", id="btn-delete-annotations",
                    n_clicks=0),
        # Confirm deletion
        html.Div(id="delete-annotations-confirmation", style={"marginTop": "10px"}),
        dcc.Loading(
            # type="circle",
            children=html.Div(id="delete-annotations-message", style={"marginTop": "10px"})
        ),
        html.Br(),
        html.Label(
            "Create csv import files : this will generate neo4j import csv files into /data/import directory. These files can be used to create a new DB.",
            style={'display': 'block', 'marginBottom': '8px'}),
        html.Button("Generate CSV Import file",
                    title="Generate csv files from gfa. These files can then be used to create the database. This process is automatic when creating a new database from gfa files.",
                    id="btn-csv-import", n_clicks=0, style={'marginRight': '10px'}),
        dcc.Loading(
            id="loading-gfa-msg",
            # type="circle",
            children=html.Div(id="gfa-message")
        ),
        html.Br(),
        html.Label("Create index or stats in database (if these steps have failed).",
                   style={'display': 'block', 'marginBottom': '8px'}),
        html.Div([
            html.Button("Create indexes",
                        title="This will generate / regenerate indexes. If the already exists there will be no action.",
                        id="btn-create-index", n_clicks=0,
                        style={'marginRight': '10px'}
            ),
            html.Button("Create stats", title="Recuperation procedure for stats node if not exist after loading data.",
                        id="btn-create-stats", n_clicks=0)
        ], style={'marginBottom': '10px'}),
        dcc.Loading(
            id="loading-index-stats-msg",
            # type="circle",
            children=html.Div(id="index_stats-message")
        ),
        html.Br(),
        html.Label("Dump database : this will generate a neo4j.dump in the /data/import directory. This file can be used to create a new DB.",
                   style={'display': 'block', 'marginBottom': '8px'}),

        html.Button("Dump DB", id="btn-dump-db", n_clicks=0),
        dcc.Loading(
            id="loading-dump-msg",
            # type="circle",
            children=html.Div(id="dump-message")
        ),
        html.Hr(style={"margin": "30px 0"}),

    

    ])



