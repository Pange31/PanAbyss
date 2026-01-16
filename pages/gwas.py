#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 21:55:56 2025

@author: fgraziani
"""

from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc

EXPORT_DIR = "./export/gwas/"

style_help = {
    "cursor": "pointer",
    "color": "black",
    "fontWeight": "bold",
    "marginLeft": "5px"
}



def layout():
    return html.Div([
        dcc.Store(id="gwas-page-store",storage_type="session"),
        html.H2("Shared Region Discovery"),

        html.Details([
            html.Summary("ℹ️ Click here to display help"),
            html.Ul([
                html.Li("Algorithm details : The principle is to determine the regions shared between a selection of haplotypes."
                        " There are two types of detection :"),
                    html.Ul([
                        html.Li("Shared nodes : Here, the objective is to detect the nodes shared by the selected haplotypes."
                                " Several parameters are involved:"),      
                            html.Ul([
                                html.Li("Min node size : a node will be detected only if it's size is superior to this value."),
                                html.Li("Max node size : a node will be detected only if it's size is inferior to this value. Set to zero or empty for no limitation."),
                                html.Li("Min (%) of selected haplotypes (= p): a node will be detected only if (p/100) * number of selected haplotypes are present on the node."
                                        "If set to zero it wil require at least one of the selected haplotypes."),
                                html.Li("Tolerance (%) (= t): a node with more than (t/100) * number of haplotypes present on this node will not be detected. If set to zero then nodes with a non selected haplotype will not be detected.")
                                
                            ]),
                            
                     html.Li("Deleted nodes : Here, the objective is to detect deletions shared by the selected haplotypes."
                             " This mode is activated only if 'include deletion' is checked."
                             " It detects nodes with the minimum of selected haplotypes and minimal size (see parameters) and at least one of the unselected haplotypes."
                             " If this node is following by a node with all the unselected haplotypes and only these haplotypes, and another following node with the selected haplotypes and almost on more haplotype, then a deletion node will be detected."
                             " In this case, the size of shared region will be incremented only by the first deleted node, the size of the other nodes will be ignored (but they will be displayed)."
                             " The following parameters are used :"),      
                         html.Ul([
                             html.Li("Min node size : a node will be detected only if it's size is superior to this value"),
                             html.Li("Max node size : a node will be detected only if it's size is inferior to this value. Set to zero or empty for no limitation."),
                             html.Li("Min (%) of selected haplotypes (= p): a node will be detected only if (p/100) * number of selected haplotypes are present on the node."
                                     "If set to zero it wil require at least one of the selected haplotypes."),
                             html.Li("Unselected haplotypes percentage (%) (= u): a node will be detected only if (u/100) * number of unselected haplotypes + (p/100) * number of selected haplotypes are present on the node."
                                     "If set to zero it wil require at least one of the unselected haplotypes."),
                             
                         ])
                            
                    ]),
                 html.Ul([
                     html.Li("general settings :"),
                         html.Ul([
                             html.Li("Size of region : this size is used to group nodes separated by less than this value (in bp)."),
                             html.Li("Limit search to chromosom : If a chromosom is selected, it will look for shared region only on this chromosom."),
                             html.Li("Reference haplotype : results will be displayed only for this haplotype, including annotations. If no one is selected then the first annotated haplotype will be displayed."
                                     "If the reference haplotype is not selected then the coordinates will be computed by selecting the first anchor node before / after the first haplotype. If no anchor is detected then the start will be set to 0 and the stop will be set to start."),
                             html.Li("Export to csv / export to csv with sequences : it will save the result into a csv file (without or with the sequences associated to each region). The file is located in the './export/gwas' directory."),
                             html.Li("Load csv : it allows to load the saved csv (it located in the '/gwas' directory or a selected file)."),
                             html.Li("First column : by clicking on the first columns it will display the region in the home page."),
                             html.Li("Region size column : by clicking on the size columns it will display the sequence associated to the region."),
                             ])
                     ])
            ])
        ], style={"marginBottom": "20px"}),


        # Area of genomes selection
        html.Div(id='genome-checkboxes'),
       html.H3("Select genomes : ", title="Select haplotypes for which you want to find shared regions."),

        dcc.Checklist(
            id='genome-list',
            options=[],
            value=[],
            labelStyle={'display': 'inline-block', 'marginRight':'10px'}
        ),
        
        html.Br(), 
        html.H3("Parameters : "),
        html.Div([
            html.Div([
                html.Label(
                    "Min node size to detect a shared region (integer) : ",
                    title="The nodes with a size below this value won't be detected by the process.",
                    style={'marginRight': '15px'}
                ),
                dcc.Input(id='gwas-min-node-size-int',style={'width': '80px', 'marginRight': '15px'},  type='number', step=1, value=10, debounce=True),
                html.Label(
                    "Max node size to detect a shared region (integer) : ",
                    title="The nodes with a size above this value won't be detected by the process. Set to 0 or empty for no max size.",
                    style={'marginRight': '15px'}
                ),
                dcc.Input(id='gwas-max-node-size-int',style={'width': '80px', 'marginRight': '15px'},  type='number', step=1, debounce=True),
                html.Label(
                    "Min (%) of selected haplotypes : ",
                     title="Min (%) of shared haplotypes = M. Number of selected haplotypes = N. To detect a shared node it must contains almost (M/100) x N of the selected haplotypes. If M = 0 then the minimum number of selected haplotypes will be 1.",
                     style={'marginRight': '15px'}),
                dcc.Input(id='gwas-min-percent_selected',style={'width': '80px', 'marginRight': '15px'}, type='number', step=1, value=100, debounce=True),
                html.Label(
                    "Tolerance (%) :  ",
                    title="Tolerance = T. Number of haplotypes on a node = n. To detect a shared node it must contains less than (T/100) x n of the non selected haplotypes. If T = 0 then detected nodes should contain only selected haplotypes.",
                    style={'marginRight': '15px'}
                ),
                dcc.Input(id='tolerance_percentage',style={'width': '80px', 'marginRight': '15px'}, type='number', step=1, value=0, debounce=True),
                html.Label(
                    "Max gap : ",
                    title="All the nodes detected will be grouped in larger regions. If two nodes are separated by less thant this value (in bp) they will be grouped in the same region.",
                     style={'marginRight': '15px'}
                ),
                dcc.Input(id='gwas-region-gap',style={'width': '80px', 'marginRight': '15px'}, type='number', step=1, value=10000, debounce=True),
                html.Div([
                    dcc.Checklist(
                        id='gwas-toggle-deletion',
                        options=[{'label':'', 'title':"If checked, the process will look for deletion node: it looks for nodes with minimal selected and unselected haplotypes followed by a node deleted for defined percentage of selected haplotype. The node size must be greater than min node size value.", 'value': 'show'}],
                        value=['show'],
                        style={'marginRight': '15px'}
                    ),
                    html.Label(
                        "Search for deletion (take more time). Unselected haplotypes percentage (%):  ",
                        style={"marginRight": "15px"},
                        title="If checked. Used to detect deleted nodes. Search for nodes with at least (T x number of unselected haplotypes / 100) unselected haplotypes + the defined percentage of selected haplotypes. For each node found, it looks for a deleted following node."
                    ),
                    dcc.Input(id='deletion-percentage',style={'width': '80px', 'marginRight': '15px'}, type='number', step=1, value=100, debounce=True),
                ], style={"display": "flex", "alignItems": "center", "marginRight": "20px"}),
            ], style={"display": "flex", "flexWrap": "wrap", "alignItems": "center", "marginRight": "20px"}
            ),
                dcc.Dropdown(
                    id='gwas_chromosomes_dropdown',
                    placeholder="Limit search to chromosome : ",
                    style={
                        "width": "250px",     
                        "minWidth": "150px",
                        "maxWidth": "100%",   
                        "flexShrink": 0
                    }
                ),
        ], style={"display": "flex", "flexWrap": "wrap", "flexDirection": "row", "alignItems": "center", "marginBottom": "20px"}),
        
        html.Br(), 
        html.Label(
            "Choose a reference haplotype :  ",
            title="Select the genome for which you want to view the results and obtain annotations. If no genome is selected, the result will be the first genome found with annotations. If there are no annotations, it will be the first genome found."
            ),
        html.Div(
                dcc.Dropdown(id='gwas_ref_genome_dropdown', placeholder="Reference haplotype : ", style={
                    "width": "250px",     
                    "minWidth": "150px",
                    "maxWidth": "100%",   
                    "flexShrink": 0
                })

        ),

    
        html.Button("Find shared regions", id='btn-find-shared', n_clicks=0, style={'margin': '15px 0'}),
        dcc.Loading(
            id="gwas_loading-spinner",
            #type="circle",  # 'default', 'circle', or 'dot'
            children=html.Div(id="load_spinner_zone")
        ),
    
        html.Div(id='shared-status', style={'marginBottom': '15px'}),
        html.Div(id='sequence-zone', style={"fontSize": "18px",
                                            "padding": "10px",
                                            "whiteSpace": "pre-wrap",
                                            "wordWrap": "break-word",
                                            "overflowWrap": "break-word"
                                            }),
        html.Label("Tip : Click on the size value to print the region sequence"), 
        html.Div([
            html.Button("💾 Export to CSV", title=f"Export the results in the table into a csv file in the {EXPORT_DIR} directory. Sequences won't be present.",id='save-csv-button', n_clicks=0),
            html.Button("💾 Export to CSV with sequences", title=f"Export the results in the table into a csv file in the {EXPORT_DIR} directory. Sequences will be present.", id='save-csv-with_seq-button', n_clicks=0),
            dcc.Download(id="download-csv"),
            html.Button("📂 Load csv", title="Load a csv file generated from this page.", id='load-csv-button', n_clicks=0),
            html.Div(id='save-feedback'),
            dcc.Upload(
                id='upload-csv',
                children=html.Div(['Glissez un fichier CSV ici ou cliquez pour sélectionner un fichier.']),
                style={
                    'display': 'none',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'padding': '10px',
                },
                multiple=False
            )
        ]),
        # Analyse array
        dash_table.DataTable(
            id='shared-region-table',
            columns=[
                {"name": "genome", "id": "genome"},
                {"name": "chromosome", "id": "chromosome"},
                {"name": "start", "id": "start"},
                {"name": "stop", "id": "stop"},
                {"name": "annotation before", "id": "annotation_before"},
                {"name": "annotations", "id": "annotations"},
                {"name": "annotation after", "id": "annotation_after"},
                {"name": "region size", "id": "region_size"},
                {"name": "shared size", "id": "shared_size"},
                {"name": "shared deleted nodes size", "id": "shared_deleted_size"},
                {"name": "Sequences", "id": "get_sequence", "presentation": "markdown"}
            ],
            data=[],
            sort_action='native',
            style_cell={
                'whiteSpace':'normal',
                'height':'auto',
                'textAlign':'left'},
            style_table={'overflowX': 'auto'},
            row_selectable='single',
            markdown_options={"html": True},
        ),
        

    html.Div(id='selected-region-output')
    ], style={'padding': '20px'})



