#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 13:52:04 2025

@author: fgraziani
"""

import dash_cytoscape as cyto
from dash import Dash, html,callback, dcc
import logging
from config import *
import dash_bootstrap_components as dbc

from pages.phylo_styles import stylesheet

logger = logging.getLogger("panabyss_logger")

PHYLO_BLOCK_TREE_RECOMPUTATION = get_phylo_conf()



# ============================================================
# Styles
# ============================================================

PARAM_CONTAINER_STYLE = {
    "display": "flex",
    "flexWrap": "wrap",
    "alignItems": "flex-start",
    "gap": "14px",
    "padding": "14px",
    "border": "1px solid #ddd",
    "borderRadius": "10px",
    "backgroundColor": "#fafafa",
    "width": "100%",
    "boxSizing": "border-box",
}

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
    "boxSizing": "border-box",
}

ACTION_BUTTON_STYLE = {
    "padding": "10px 18px",
    "fontSize": "15px",
    "minHeight": "42px",
}


def layout():

    return html.Div([

        html.H2("Phylogenetic tree"),

        # ============================================================
        # HELP SECTION
        # ============================================================

        html.Details([
            html.Summary("ℹ️ Click here to display help"),

            html.Ul([

                html.Li(
                    "This page allows to display 2 phylogenetic trees :"
                ),

                html.Ul([

                    html.Li(
                        "Load a newick file : this allows you to load a "
                        "file and display a reference tree, for example. "
                        "For that, juste drag / drop or select the newick file."
                    ),

                    html.Li(
                        "Plot global tree : this will compute a global tree "
                        "with a either a RAxML-NG (not recommended for big "
                        "pangenome) or a distance matrix and neighbor joining "
                        "(faster than raxml-ng but less accurate). "
                        "These two methods use a sample matrix of presence / "
                        "absence of a subset of sampled nodes. Nodes traversed "
                        "directly or in reverse are considered different nodes. "
                        "It is possible to select a chromosome to limit the "
                        "tree to this chromosome. If no chromosome selected "
                        "then the tree is computed on all chromosomes. "
                        "If the global tree has already been computed, it is "
                        "possible to load it directly by clicking on the "
                        "corresponding button (the button is hidden if no "
                        "existing tree)."
                    ),

                    html.Li(
                        "Plot tree of selected region : This allows you to "
                        "calculate a tree for the region currently being viewed "
                        "on the home page. It is therefore necessary to select "
                        "a region to view beforehand (on the home or gwas pages). "
                        "The tree is constructed based on a distance matrix. "
                        "This matrix is calculated using the Jaccard index, "
                        "taking into account the strand and repetition of each "
                        "node. If Weight by node size is checked, then nodes "
                        "size are used to weight the Jaccard index. The tree is "
                        "then calculated using the neighbor joining algorithm."
                    ),

                ])
            ])

        ], style={
            "marginBottom": "20px"
        }),


        # ============================================================
        # PARAMETERS — TWO COLUMNS
        # ============================================================

        html.Div([

            # ========================================================
            # LEFT COLUMN — GLOBAL TREE PARAMETERS
            # ========================================================

            html.Div([

                html.H3(
                    "Global tree parameters",
                    style={
                        "marginTop": "0",
                        "marginBottom": "8px"
                    }
                ),

                html.Div([

                    # ------------------------------------------------
                    # Global parameters
                    # ------------------------------------------------

                    html.Div([

                        # Chromosome
                        html.Div([
                            html.Label(
                                "Chromosome",
                                title=(
                                    "The regions / annotations will be "
                                    "related only to this chromosome."
                                )
                            ),

                            dcc.Dropdown(
                                id='phylogenetic_chromosomes_dropdown',
                                placeholder="All",
                                style={
                                    "width": "130px"
                                }
                            ),

                        ], style=PARAM_STYLE),


                        # Method
                        html.Div([
                            html.Label("Method"),

                            dcc.Dropdown(
                                id='method-dropdown',
                                options=[
                                    {
                                        'label': 'Neighbor Joining (fast method)',
                                        'value': 'nj'
                                    },
                                    {
                                        'label': 'RAxML-NG (Maximum Likelihood)',
                                        'value': 'raxml'
                                    }
                                ],
                                value='nj',
                                clearable=False,
                                style={
                                    "width": "100px"
                                }
                            ),

                        ], style=PARAM_STYLE),


                        # Cache
                        html.Div([
                            html.Label("Cache"),

                            html.Div(
                                dcc.Checklist(
                                    id="use-cache-checkbox",
                                    options=[
                                        {
                                            "label": "Use cache",
                                            "value": "use_cache"
                                        }
                                    ],
                                    value=["use_cache"]
                                ),

                                title=(
                                    "If same tree has already been computed "
                                    "retrieve the result."
                                )
                            )

                        ], style={
                            **PARAM_STYLE,
                            "display": (
                                "none"
                                if PHYLO_BLOCK_TREE_RECOMPUTATION
                                else "flex"
                            )
                        }),

                    ], style={
                        "display": "flex",
                        "flexWrap": "wrap",
                        "alignItems": "flex-start",
                        "gap": "8px",
                        "width": "100%",
                        "boxSizing": "border-box"
                    }),


                    # ------------------------------------------------
                    # Global action buttons
                    # ------------------------------------------------

                    html.Div([

                        html.Button(
                            "Plot global tree",
                            title=(
                                "This will compute the tree for the whole "
                                "pangenome or load last computed tree if exist..."
                            ),
                            n_clicks=0,
                            id="btn-plot-global-tree",
                            style={
                                **ACTION_BUTTON_STYLE,
                                "marginRight": "8px",
                                "marginBottom": "5px"
                            }
                        ),

                        html.Button(
                            "Cancel",
                            title=(
                                "This will cancel the global tree construction."
                            ),
                            disabled=True,
                            n_clicks=0,
                            id="btn-cancel-plot-global-tree",
                            style={
                                **ACTION_BUTTON_STYLE,
                                "marginRight": "8px",
                                "marginBottom": "5px"
                            }
                        ),

                        html.Button(
                            "Save tree",
                            title=(
                                "This will compute a new tree even it already exist."
                            ),
                            id="btn-save-global-tree",
                            n_clicks=0,
                            style={
                                **ACTION_BUTTON_STYLE,
                                "marginBottom": "5px"
                            }
                        ),

                        html.Div(
                            html.Div(className="custom-spinner"),
                            id="phylo-spinner-container",
                            style={
                                "display": "none",
                                "marginLeft": "10px"
                            }
                        )

                    ], style={
                        "display": "flex",
                        "flexWrap": "wrap",
                        "alignItems": "center",
                        "gap": "0",
                        "marginTop": "8px",
                        "width": "100%",
                        "boxSizing": "border-box"
                    })

                ], style={
                    **PARAM_CONTAINER_STYLE,
                    "width": "calc(100% - 10px)",
                    "minHeight": "165px",
                    "height": "auto",
                    "padding": "10px",
                    "gap": "6px",
                    "alignContent": "flex-start",
                    "boxSizing": "border-box",
                    "overflow": "visible"
                })

            ], style={
                "width": "calc(50% - 1px)",
                "minWidth": "0",
                "boxSizing": "border-box"
            }),


            # ========================================================
            # RIGHT COLUMN — LOCAL TREE PARAMETERS
            # ========================================================

            html.Div([

                html.H3(
                    "Local tree parameters",
                    style={
                        "marginTop": "0",
                        "marginBottom": "8px"
                    }
                ),

                html.Div([

                    # ------------------------------------------------
                    # Local parameter
                    # ------------------------------------------------

                    html.Div([

                        html.Label("Method"),

                        html.Div(
                            dcc.Checklist(
                                options=[
                                    {
                                        "label": "Weight by node size",
                                        "value": "weight_by_size"
                                    }
                                ],
                                value=[],
                                id="checkbox-weight-node-size"
                            ),

                            title=(
                                "If checked, the nodes size will be used "
                                "to weight the Jaccard index used to "
                                "compute distance matrix."
                            )
                        )

                    ], style=PARAM_STYLE),


                    # ------------------------------------------------
                    # Local action buttons
                    # ------------------------------------------------

                    html.Div([

                        html.Button(
                            "Plot local tree",
                            id="btn-plot-region",
                            n_clicks=0,
                            title=(
                                "Before using this button, data must be "
                                "displayed on home page. If data are displayed, "
                                "then this will computes the phylogenetic tree "
                                "of displayed haplotypes (see help for details)."
                            ),
                            style={
                                **ACTION_BUTTON_STYLE,
                                "marginRight": "8px",
                                "marginBottom": "5px"
                            }
                        ),

                        html.Button(
                            "Save tree",
                            title=(
                                "Tree will be saved into newick format."
                            ),
                            id="btn-save-tree",
                            style={
                                **ACTION_BUTTON_STYLE,
                                "marginBottom": "5px"
                            }
                        )

                    ], style={
                        "display": "flex",
                        "flexWrap": "wrap",
                        "alignItems": "center",
                        "gap": "0",
                        "marginTop": "5px",
                        "marginBottom": "5px",
                        "width": "100%",
                        "boxSizing": "border-box"
                    })

                ], style={
                    **PARAM_CONTAINER_STYLE,
                    "width": "calc(100% - 10px)",
                    "minHeight": "165px",
                    "height": "auto",
                    "padding": "10px",
                    "gap": "6px",
                    "alignContent": "flex-start",
                    "boxSizing": "border-box",
                    "overflow": "visible"
                })

            ], style={
                "width": "calc(50% - 1px)",
                "minWidth": "0",
                "boxSizing": "border-box"
            })

        ], style={
            "width": "100%",
            "display": "flex",
            "alignItems": "flex-start",
            "gap": "2px",
            "boxSizing": "border-box",
            "minWidth": "0"
        }),


        # ============================================================
        # OPTIONS
        # ============================================================

        html.Details(
            [

                html.Summary(
                    "Options",
                    style={
                        "cursor": "pointer",
                        "fontSize": "1.17em",
                        "fontWeight": "bold",
                        "marginTop": "0",
                        "marginBottom": "8px",
                        "padding": "0",
                        "userSelect": "none"
                    }
                ),

                html.Div(
                    [

                        # ================================================
                        # Load reference tree
                        # ================================================

                        html.Div(
                            [

                                html.Label(
                                    "Load a reference tree",
                                    style={
                                        "fontSize": "20px",
                                        "fontWeight": "500",
                                        "marginBottom": "6px",
                                        "display": "block"
                                    }
                                ),

                                dcc.Upload(
                                    id='upload-newick',
                                    children=html.Div(
                                        "Drag/drop or select a reference tree newick file",
                                        style={
                                            "whiteSpace": "nowrap",
                                            "overflow": "hidden",
                                            "textOverflow": "ellipsis",
                                            "padding": "0 8px",
                                            "boxSizing": "border-box"
                                        }
                                    ),
                                    style={
                                        "width": "400px",  # ← largeur du sélecteur
                                        "maxWidth": "100%",
                                        "height": "30px",
                                        "lineHeight": "28px",
                                        "borderWidth": "1px",
                                        "borderStyle": "dashed",
                                        "borderRadius": "5px",
                                        "textAlign": "center",
                                        "boxSizing": "border-box",
                                        "overflow": "hidden",
                                        "backgroundColor": "white"
                                    },
                                    multiple=False
                                ),

                                html.Div(
                                    id='upload-status',
                                    style={
                                        "marginTop": "4px",
                                        "lineHeight": "15px"
                                    }
                                )

                            ],
                            style={
                                "marginBottom": "15px"
                            }
                        ),

                        # ================================================
                        # Root tree
                        # ================================================

                        html.Div(
                            [
                                html.Label(
                                    "Root tree at:",
                                    style={
                                        "fontSize": "16px",
                                        "marginRight": "8px"
                                    }
                                ),

                                dcc.Dropdown(
                                    id="root-tree-dropdown",
                                    options=[
                                        {
                                            "label": "No rooting",
                                            "value": None
                                        },
                                        {
                                            "label": "Midpoint",
                                            "value": "midpoint"
                                        }
                                    ],
                                    value=None,
                                    clearable=False,
                                    style={
                                        "width": "300px"
                                    }
                                )
                            ],
                            style={
                                "display": "flex",
                                "alignItems": "center",
                                "marginBottom": "15px"
                            }
                        ),

                        # ================================================
                        # Colors
                        # ================================================

                        html.Div(
                            [

                                html.Label(
                                    "Colors",
                                    style={
                                        "fontSize": "20px",
                                        "fontWeight": "500",
                                        "marginBottom": "8px",
                                        "display": "block"
                                    }
                                ),

                                html.Div(
                                    id="individual-color-pickers",
                                    children=[],
                                    style={
                                        "display": "flex",
                                        "flexWrap": "wrap",
                                        "alignItems": "center",
                                        "gap": "5px 10px",
                                        "width": "100%"
                                    }
                                )

                            ],
                            style={
                                "width": "100%"
                            }
                        )

                    ],

                    style={
                        "width": "calc(100% - 10px)",
                        "padding": "14px",
                        "border": "1px solid #ddd",
                        "borderRadius": "10px",
                        "backgroundColor": "#fafafa",
                        "boxSizing": "border-box",
                        "marginTop": "5px"
                    }
                )

            ],

            open=False,

            style={
                "width": "100%",
                "marginTop": "2px",
                "marginBottom": "10px"
            }
        ),

        # ============================================================
        # POLLER
        # ============================================================

        dcc.Interval(
            id="phylo-job-poller",
            interval=5000,
            n_intervals=0,
            disabled=True
        ),


        # ============================================================
        # TREES — TWO COLUMNS
        # ============================================================

        html.Div([

            # ========================================================
            # LEFT — GLOBAL TREE
            # ========================================================

            html.Div([

                html.H3(
                    "Global tree",
                    style={
                        "marginTop": "2px",
                        "marginBottom": "8px"
                    }
                ),

                html.Div(
                    cyto.Cytoscape(
                        id='cytoscape-phylo',
                        elements=[],
                        stylesheet=stylesheet,
                        layout={
                            'name': 'preset'
                        },
                        style={
                            "width": "100%",
                            "height": "1000px"
                        },
                        zoomingEnabled=True,
                        userZoomingEnabled=True,
                        wheelSensitivity=0.1,
                        autounselectify=False,
                        boxSelectionEnabled=True,
                    ),
                    style={
                        'width': '100%',
                        'boxSizing': 'border-box',
                    }
                ),

            ], style={
                "width": "calc(50% - 1px)",
                "minWidth": "0",
                "boxSizing": "border-box"
            }),



            # ========================================================
            # RIGHT — LOCAL TREE
            # ========================================================

            html.Div([

                # ----------------------------------------------------
                # Technical components
                # ----------------------------------------------------

                dcc.Download(
                    id="download-tree"
                ),

                dcc.Loading(
                    id="loading-phylogenetic-msg",
                    children=html.Div(
                        id="phylogenetic-message"
                    )
                ),

                html.Div(
                    id='region-status',
                    style={
                        "margin": "6px 0"
                    }
                ),

                # ----------------------------------------------------
                # Local tree
                # ----------------------------------------------------

                html.H3(
                    "Local tree",
                    style={
                        "marginTop": "2px",
                        "marginBottom": "8px"
                    }
                ),
                html.Div(
                    cyto.Cytoscape(
                        id='cytoscape-phylo-region',
                        elements=[],
                        stylesheet=stylesheet,
                        layout={
                            'name': 'preset'
                        },
                        style={
                            "width": "100%",
                            "height": "1000px"
                        },
                        zoomingEnabled=True,
                        userZoomingEnabled=True,
                        wheelSensitivity=0.1,
                        autounselectify=False,
                        boxSelectionEnabled=True,
                    ),
                    style={
                        'width': '100%',
                        'boxSizing': 'border-box',
                    }
                ),


            ], style={
                "width": "calc(50% - 1px)",
                "minWidth": "0",
                "boxSizing": "border-box"
            })

        ], style={
            "width": "100%",
            "display": "flex",
            "alignItems": "flex-start",
            "gap": "13px",
            "boxSizing": "border-box",
            "minWidth": "0"
        })

    ],
    style={
        "width": "100%",
        "padding": "20px",
        "boxSizing": "border-box"
    })



