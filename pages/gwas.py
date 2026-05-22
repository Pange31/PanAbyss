#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 21:55:56 2025

@author: fgraziani
"""

from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc
from plotly.graph_objs.layout.ternary.aaxis import title

EXPORT_DIR = "./export/gwas/"

style_help = {
    "cursor": "pointer",
    "color": "black",
    "fontWeight": "bold",
    "marginLeft": "5px"
}

#Parameters style for the different gwas parameters
PARAM_STYLE = {
    "display": "flex",
    "flexDirection": "column",
    "gap": "6px",
    "minWidth": "130px",
    "padding": "10px",
    "border": "1px solid #e5e5e5",
    "borderRadius": "8px",
    "backgroundColor": "white",
    "boxShadow": "0 1px 3px rgba(0,0,0,0.05)",
}

def layout():
    return html.Div([
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
                            
                     html.Li("Deleted nodes: Here, the objective is to detect deletions shared by the selected haplotypes."
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
                     html.Li("General:"),
                         html.Ul([
                             html.Li("Size of region : this size is used to group nodes separated by less than this value (in bp)."),
                             html.Li("Limit search to chromosom : If a chromosom is selected, it will look for shared region only on this chromosom."),
                             html.Li("Reference haplotype : results will be displayed only for this haplotype, including annotations. If no one is selected then the first annotated haplotype will be displayed."
                                     "If the reference haplotype is not selected then the coordinates will be computed by selecting the first anchor node before / after the first haplotype. If no anchor is detected then the start will be set to 0 and the stop will be set to start."),
                             html.Li("Export to csv / export to csv with sequences : it will save the result into a csv file (without or with the sequences associated to each region). The file is located in the './export/gwas' directory."),
                             html.Li("Load csv : it allows to load the saved csv (it located in the '/gwas' directory or a selected file)."),
                             html.Li("First column : by clicking on the first columns it will display the region in the home page."),
                             html.Li("Sequences column : by clicking on this column it will display the sequence associated to the region."),
                             html.Li(
                                 "Use cache: if the same search has been computed then the results will be immediately retrieved."),
                             html.Li(
                                 "p-value: The p-value is calculated using a chi-square test on the contingency table of selected/non-selected individuals. The p-values are then aggregated on the region using the ACAT method, weighted by the sizes of the nodes.")
                             ]),

                     ])
            ])
        ], style={"marginBottom": "20px"}),


        # Area of genomes selection
        html.Div(id='genome-checkboxes'),
       html.H3("Select genomes : ", title="Select haplotypes for which you want to find shared regions."),

        # dcc.Checklist(
        #     id="genome-list",
        #     options=[],
        #     value=[],
        #     inline=True
        # ),
        dcc.Checklist(
            id="genome-list",
            options=[],
            value=[],
            inline=False,

            labelStyle={
                "display": "flex",
                "alignItems": "center",
                "gap": "5px",
                "fontSize": "14px",
                "whiteSpace": "nowrap"
            },

            style={
                "display": "grid",
                #"gridTemplateColumns": "repeat(auto-fit, minmax(120px, 1fr))",
                "gridTemplateColumns": "repeat(auto-fill, 250px)",
                "columnGap": "4px",
                "rowGap": "2px",
                "columnGap": "0px",
                "rowGap": "0px",
                "width": "100%"
            }
        ),

        html.Br(), 
        html.H3("Parameters : "),

        #Parameters section
        html.Div(
            [
                # =====================================================
                # Min node size
                # =====================================================
                html.Div(
                    [
                        html.Label(
                            "Min node size"
                        ),
                        html.Div(
                            dcc.Input(
                                id='gwas-min-node-size-int',
                                type='number',
                                step=1,
                                value=10,
                                debounce=False,
                                style={"width": "75px"}
                            ),
                            title="The nodes with a size below this value won't be detected by the process."
                        ),
                    ],
                    style=PARAM_STYLE
                ),

                # =====================================================
                # Max node size
                # =====================================================
                html.Div(
                    [
                        html.Label(
                            "Max node size",

                        ),
                        html.Div(
                            dcc.Input(
                                id='gwas-max-node-size-int',
                                type='number',
                                step=1,
                                debounce=False,
                                style={"width": "75px"}
                            ),
                            title="Set to 0 for no limit."
                        )
                    ],
                    style=PARAM_STYLE
                ),

                # =====================================================
                # Min selected %
                # =====================================================
                html.Div(
                    [
                        html.Label(
                            "Min selected (%)",

                        ),
                        html.Div(
                            dcc.Input(
                                id='gwas-min-percent_selected',
                                type='number',
                                step=1,
                                value=100,
                                debounce=False,
                                style={"width": "75px"}
                            ),
                            title="Minimum percentage of selected haplotypes."
                        )
                    ],
                    style=PARAM_STYLE
                ),

                # =====================================================
                # Tolerance
                # =====================================================
                html.Div(
                    [
                        html.Label(
                            "Tolerance (%)",

                        ),
                        html.Div(
                            dcc.Input(
                                id='tolerance_percentage',
                                type='number',
                                step=1,
                                value=0,
                                debounce=False,
                                style={"width": "75px"}
                            ),
                            title="Allowed percentage of non-selected haplotypes."
                        )
                    ],
                    style=PARAM_STYLE
                ),

                # =====================================================
                # Max gap
                # =====================================================
                html.Div(
                    [
                        html.Label(
                            "Max gap (bp)",

                        ),
                        html.Div(
                            dcc.Input(
                                id='gwas-region-gap',
                                type='number',
                                step=1,
                                value=10000,
                                debounce=False,
                                style={
                                    "width": "170px"
                                }
                            ),
                            title="Merge nodes that are separated by less than the GAP threshold."
                        )
                    ],
                    style=PARAM_STYLE
                ),

                # =====================================================
                # Cache
                # =====================================================
                html.Div(
                    [
                        html.Label("Cache"),
                        html.Div(
                            dcc.Checklist(
                                id='use-cache-checkbox',
                                options=[
                                    {
                                        'label': 'Use cache',
                                        'value': 'enabled'
                                    }
                                ],
                                value=['enabled']
                            ),
                            title="If enabled, previously computed results will be reused instead of rerunning the computation"
                        )
                    ],
                    style=PARAM_STYLE
                ),

                # =====================================================
                # Deletion search
                # =====================================================
                html.Div(
                    [
                        html.Label(
                            "Deletion search",
                            title="Search for deletion nodes."
                        ),
                        html.Div(
                            [
                                dcc.Checklist(
                                    id='gwas-toggle-deletion',
                                    options=[
                                        {
                                            'label': 'Enable',
                                            'value': 'show'
                                        }
                                    ],
                                    value=['show']
                                ),
                                html.Div(
                                    dcc.Input(
                                        id='deletion-percentage',
                                        type='number',
                                        step=1,
                                        value=100,
                                        debounce=False,
                                        placeholder="% unselected",
                                        style={"width": "85px"}
                                    ),
                                    title="A deletion is detected only when at least the specified percentage of non-selected individuals is present on it."
                                )
                            ],
                            style={
                                "display": "flex",
                                "alignItems": "center",
                                "gap": "10px",
                            }
                        ),

                    ],
                    style=PARAM_STYLE
                ),

                # =====================================================
                # Chromosome dropdown
                # =====================================================
                html.Div(
                    [

                        html.Label("Chromosome"),

                        dcc.Dropdown(
                            id='gwas_chromosomes_dropdown',
                            placeholder="All chromosomes",
                            style={
                                "width": "200px"
                            }
                        ),
                    ],
                    style=PARAM_STYLE
                ),
            ],
            #Global style for the whole parameters
            style={
                "display": "flex",
                "flexWrap": "wrap",
                "alignItems": "flex-start",
                "gap": "14px",
                "padding": "14px",
                "border": "1px solid #ddd",
                "borderRadius": "10px",
                "backgroundColor": "#fafafa",
                "maxWidth": "100%",
            }
        ),

        html.Br(),
        html.Label(
            "Choose a reference haplotype :  ",
            title="Select the genome for which you want to view the results and obtain annotations. If no genome is selected, the result will be the first genome found with annotations. If there are no annotations, it will be the first genome found."
            ),
        html.Div(
                dcc.Dropdown(id='gwas_ref_genome_dropdown',
                    options=[],
                    value=None,
                    placeholder="Reference haplotype : ",
                    clearable=False,
                    style={
                    "width": "250px",
                    "minWidth": "150px",
                    "maxWidth": "100%",
                    "flexShrink": 0
                })

        ),

    
        html.Button("Find shared regions", id='btn-find-shared', n_clicks=0, style={'margin': '15px 0', 'marginRight':'15px'}),
        html.Button("Cancel", id='btn-cancel-find-shared', disabled=True, n_clicks=0, style={'margin': '15px 0'}),
        html.Div(id='shared-status', style={'marginBottom': '15px'}),
        html.Div(
            id="gwas-progress-circle",
            style={"display": "none", "width": "100px", "height": "100px", "margin": "auto", "position": "relative",
                   "border-radius": "50%"},
            children=[
                html.Span(id="gwas-progress-text", children="0%",
                          style={"position": "absolute", "width": "100%", "text-align": "center",
                                 "line-height": "100px", "font-weight": "bold", "color": "#9d4edd"})
            ]
        ),
        html.Div(
            html.Div(className="custom-spinner"),
            id="spinner-container",
            style={"display": "none", "marginTop": "20px"}
        ),

        dcc.Interval(
            id='gwas-poll-interval',
            interval=5 * 1000,  # poll every 2 seconds
            n_intervals=0,
            disabled=True
        ),


        dcc.Graph(
            id="chromosome-graph",
            figure={},
            style={"display": "none"}
        ),


        html.Div(id='sequence-zone', style={"fontSize": "18px",
                                            "padding": "10px",
                                            "whiteSpace": "pre-wrap",
                                            "wordWrap": "break-word",
                                            "overflowWrap": "break-word"
                                            }),

        html.Div([

            html.Button(
                "💾 Export to CSV",
                id="save-csv-button",
                n_clicks=0,
                title=f"Export results into {EXPORT_DIR} (no sequences)",
                style={
                    "minWidth": "140px",
                    "height": "32px",
                    "whiteSpace": "nowrap",
                    "padding": "0 12px",
                    "fontSize": "12px",
                    "flexShrink": 0,
                }
            ),

            html.Button(
                "💾 Export with sequences",
                id="save-csv-with_seq-button",
                n_clicks=0,
                title=f"Export results into {EXPORT_DIR} (with sequences)",
                style={
                    "minWidth": "190px",
                    "height": "32px",
                    "whiteSpace": "nowrap",
                    "padding": "0 12px",
                    "fontSize": "12px",
                    "flexShrink": 0,
                }
            ),

            dcc.Download(id="download-csv"),
            dcc.Store(id="uploaded-file-store"),

            dcc.Upload(
                id="upload-csv",
                children=html.Div(
                    id="upload-csv-label",
                    children="📂 Drag & drop or click to select CSV",
                ),
                multiple=False,
                style={
                    "minWidth": "320px",
                    "width": "320px",
                    "height": "32px",

                    "display": "flex",
                    "alignItems": "center",
                    "justifyContent": "center",

                    "borderWidth": "1px",
                    "borderStyle": "dashed",
                    "borderRadius": "4px",

                    "fontSize": "11px",
                    "whiteSpace": "nowrap",
                    "overflow": "hidden",

                    "flexShrink": 0,
                    "boxSizing": "border-box",
                }
            ),

            html.Button(
                "Load CSV",
                id="run-csv-button",
                n_clicks=0,
                title="Run analysis on uploaded CSV",
                style={
                    "minWidth": "110px",
                    "height": "32px",
                    "whiteSpace": "nowrap",
                    "padding": "0 12px",
                    "fontSize": "12px",
                    "flexShrink": 0,
                }
            ),

            html.Div(id='save-feedback'),

        ], style={
            "display": "flex",
            "alignItems": "center",
            "gap": "8px",
            "flexWrap": "wrap",
            "width": "100%",
        }),


        # Analyse array
        dash_table.DataTable(
            id='shared-region-table',
            columns=[
                {"name": "Genome", "id": "genome"},
                {"name": "Chr", "id": "chromosome"},
                {"name": "Start", "id": "start"},
                {"name": "Stop", "id": "stop"},
                {"name": "Annotation before\ngene / distance", "id": "annotation_before", "presentation": "markdown"},
                {"name": "Annotations", "id": "annotations", "presentation": "markdown"},
                {"name": "Annotation after\ngene / distance", "id": "annotation_after", "presentation": "markdown"},
                {"name": "Region size", "id": "region_size", 'type': 'numeric'},
                {"name": "Shared size", "id": "shared_size", 'type': 'numeric'},
                {"name": "Shared deleted nodes size", "id": "shared_deleted_size", 'type': 'numeric'},
                {"name": "-log10(p-value)", "id": "pval", 'type': 'numeric'},
                {"name": "Score", "id": "score", 'type': 'numeric'},
                # {"name": "nodes number", "id": "nb_nodes_in_region"},
                # {"name": "start position mean", "id": "start_position_mean"},
                # {"name": "stop position mean", "id": "stop_position_mean"},
                #{"name": "Sequences", "id": "get_sequence", "presentation": "markdown", "selectable": True}
                {"name": "Sequences", "id": "get_sequence"}
            ],
            data=[],
            sort_action='native',
            filter_action='native',
            style_header={
                'whiteSpace': 'normal',
                'fontFamily': 'Inter, system-ui, sans-serif',
                'fontWeight': '600',
                'fontSize': '11px',

                'backgroundColor': '#f8fafc',
                'color': '#334155',

                'borderBottom': '1px solid #cbd5e1',
                'borderRight': '1px solid #e2e8f0',
                'borderLeft': '1px solid #e2e8f0',

                'padding': '6px 8px',
                'textAlign': 'center',
            },
            cell_selectable=True,
            style_cell={
                "whiteSpace": "pre-line",
                'userSelect': 'text',
                'textAlign': 'left',
                "height": "auto",
                "lineHeight": "15px",
                "textOverflow": "ellipsis",
                "minWidth": "0px",
                "maxWidth": "none",
                'fontSize': '11px',
                'whiteSpace': 'normal',
                "overflow": "hidden",
            },
            style_cell_conditional=[
                {
                    "if": {"column_id": "annotation_before"},
                    "minWidth": "120px",
                    "maxWidth": "120px",
                },
                {
                    "if": {"column_id": "annotation_after"},
                    "minWidth": "120px",
                    "maxWidth": "120px",
                },
            ],
            style_table={
                "overflowX": "auto",
                "width": "100%",
                "maxWidth": "100%",
                "border": "1px solid #e2e8f0",
                "borderRadius": "6px",
            },
            row_selectable='single',
            markdown_options={"html": True},
        ),
        

    html.Div(id='selected-region-output')
    ], style={'padding': '20px'})



