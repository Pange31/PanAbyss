import plotly.graph_objects as go


from neo4j_requests import *

import dash
from dash import Dash, dcc, html, Input, Output, State, callback_context, ctx, MATCH, ALL, no_update
import dash_svg as svg
from urllib.parse import parse_qs, urlparse
from dash.exceptions import PreventUpdate
import dash_cytoscape as cyto
import dash_bootstrap_components as dbc
import matplotlib.colors as mcolors

import pandas as pd
import numpy as np
import itertools
import json
import math

from app import *
from cache_manager import *

import os
from io import BytesIO
import base64
import logging

import urllib.parse


logger = logging.getLogger("panabyss_logger")

cyto.load_extra_layouts()

DEFAULT_SIZE_VALUE = 10
DEFAULT_SHARED_REGION_COLOR = "#008000"
DEFAULT_EXONS_COLOR = "#008000"
EXPORT_DIR = './export/graphs'

MAX_NODES_TO_VISUALIZE = get_max_nodes_to_visualize()
MAX_NODES_FROM_DB = get_max_nodes_from_db()
logger.debug(f"Conf max nodes to visualize: {MAX_NODES_TO_VISUALIZE} max nodes from db: {MAX_NODES_FROM_DB}")

#MAX_GAP is used to dash edges between nodes separated by more than this value
MAX_GAP = 50000

#MAX size of displayed nodes
S_MAX = 200
#MIN size of displayed nodes
S_MIN = 30

#MAX_EDGE WIDTH
MAX_EDGE_WIDTH = 4

def records_to_dataframe(nodes_data):
    rows = []
    for record in nodes_data:
        rows.append(nodes_data[record])
    return pd.DataFrame(rows)


def compute_stylesheet(color_number, nodes_names=False, node_shape_as_circle=False, nodes_size_scale=1):
    node_shape = 'circle' if node_shape_as_circle else 'round-rectangle'
    node_height = 'data(displayed_node_size)' if node_shape_as_circle else 18*nodes_size_scale
    if color_number > 1:
        stylesheet = [
            {
            'selector': 'node',
            'style': {
                'label': 'data(name)' if nodes_names else '',
                'backgroundColor':'data(color)',
                'text-opacity':1,
                'opacity':1,
                'shape': node_shape,
                'width':'data(displayed_node_size)',
                'height': node_height,
                #'height':'data(displayed_node_size)',
                'z-compound-depth': 'top'

                }
            },
            {
                'selector': '.main-node',
                'style': {'shape': node_shape}
            },
            {
                'selector': '.degenerate-node',
                'style': {'shape': 'ellipse'}
            },
            {
                'selector': 'edge',
                'style': {
                    'curve-style': 'unbundled-bezier',
                    'control-point-weights': [0.5],
                    'target-arrow-color': 'data(color)',
                    'target-arrow-shape': 'triangle',
                    'arrow-scale': 1,
                    'control-point-distances': [1],
                    'opacity':0.9,
                    'z-compound-depth': 'bottom'
                },
            },
            # {'selector': ':selected', 'style': {
            #     'backgroundColor': 'grey',
            #     'line-color': '#FF4136',
            #     'border-width': 3,
            #     'border-color': 'black'
            # }}
            {'selector': 'edge:selected', 'style': {
                'line-color': 'red',
                'background-color': 'red',
                'target-arrow-color': 'red',
                'width': 6,
                'opacity': 1,
                'z-index': 9999
            }},
            {'selector': 'edge:hover', 'style': {
                'line-color': 'red',
                'target-arrow-color': 'red',
                'width': 5,
                'opacity': 1,
            }},

            {'selector': 'node:selected', 'style': {
                'background-color': 'grey',
                'line-color': '#FF4136',
                'border-width': 3,
                'border-color': 'black'
            }},
            
            ]
        for i in range(1, color_number+1):
            sign = 1 if i % 2 == 0 else -1
            distance = sign * (20 + 10 * (i // 2))
            distance = 1 + 5*i
            stylesheet.append({
                'selector': f'.offset-{i}',
                'style': {'control-point-distances': [distance], 'opacity': 0.6, }
            })
    else:
        stylesheet = [
            {
                'selector': 'node',
                'style': {
                    'backgroundColor':'data(color)',
                    'label': 'data(name)' if nodes_names else '',
                    'text-opacity': 1,
                    'opacity':1,
                    'shape': node_shape,
                    'width':'data(displayed_node_size)',
                    'height': node_height,
                    #'height':'data(displayed_node_size)',
                    'z-compound-depth': 'top'
                }
            },
            {
                'selector': '.main-node',
                'style': {'shape': node_shape}
            },
            {
                'selector': '.degenerate-node',
                'style': {'shape': 'ellipse'}
            },
            {
                'selector': 'edge',
                'style': {
                    'line-color': '#A3C4BC',
                    'target-arrow-color': '#A3C4BC',
                    'target-arrow-shape': 'triangle',
                    'arrow-scale': 1,
                    'curve-style': 'straight',
                    'opacity':0.9,
                    'z-compound-depth': 'bottom'


                }
            },
            # {'selector': ':selected', 'style': {
            #     'background-color': 'grey',
            #     'line-color': '#FF4136',
            #     'border-width': 3,
            #     'border-color': 'black'
            # }}
            {'selector': 'edge:selected', 'style': {
                'line-color': 'red',
                'background-color': 'red',
                'target-arrow-color': 'red',
                'width': 6,
                'opacity': 1,
                'z-index': 9999
            }},
            {'selector': 'edge:hover', 'style': {
                'line-color': 'red',
                'target-arrow-color': 'red',
                'width': 5,
                'opacity': 1,
            }},

            {'selector': 'node:selected', 'style': {
                'background-color': 'grey',
                'line-color': '#FF4136',
                'border-width': 3,
                'border-color': 'black'
            }},
        ]
    return stylesheet


def flow_to_rgb(flow, node_style="default", exons_color=DEFAULT_EXONS_COLOR):
    if node_style=="exon":
        defined_color = hex_to_rgb_string(exons_color)
    else:
        r = int(255 * flow)
        g = int(0)
        b = int(255 * (1 - flow))
        defined_color = f'rgb({r},{g},{b})'
    return defined_color

def hex_to_rgb_string(hex_color):
    hex_color = hex_color.lstrip('#')
    if len(hex_color) != 6:
        raise ValueError(f"Invalid hex color: {hex_color}")
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return f'rgb({r}, {g}, {b})'


def get_color_palette(n):
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap("tab20")
    return [f'rgb({int(r*255)}, {int(g*255)}, {int(b*255)})' for r, g, b, _ in cmap(np.linspace(0, 1, n))]




#Merge the data of nodes to be removed into their predecessor node.
def merge_node_data(df, nodes_to_remove, predecessors):

    df = df.copy()

    # Fast lookup
    name_to_idx = dict(zip(df["name"], df.index))
    name_to_row = df.set_index("name").to_dict("index")

    for n in nodes_to_remove:

        if n not in predecessors:
            continue

        pred_list = list(predecessors[n])

        if not pred_list:
            continue

        pred_node = pred_list[0]

        visited_nodes = {n, pred_node}

        # Find first predecessor not removed
        while pred_node in nodes_to_remove and pred_node in predecessors:
            preds = list(predecessors.get(pred_node, []))
            if not preds:
                pred_node = None
                break

            pred_node = preds[0]

            if pred_node in visited_nodes and len(preds) > 1:
                pred_node = preds[1]

            if pred_node in visited_nodes:
                pred_node = None
                break

            visited_nodes.add(pred_node)

        if pred_node is None:
            continue

        if pred_node not in name_to_idx:
            continue

        if n not in name_to_row:
            continue

        idx_pred = name_to_idx[pred_node]

        # Node being merged (never modified)
        row_n = name_to_row[n]

        # SIZE
        size_n = row_n.get("size", 0) or 0
        size_pred = df.at[idx_pred, "size"]
        df.at[idx_pred, "size"] = size_pred + size_n


        # FEATURES
        f1 = row_n.get("features")
        f1 = f1 if isinstance(f1, list) else []

        f2 = df.at[idx_pred, "features"]
        f2 = f2 if isinstance(f2, list) else []

        df.at[idx_pred, "features"] = list(set(f1 + f2))


        # EXONS

        e1 = row_n.get("exons")
        e1 = e1 if isinstance(e1, list) else []

        e2 = df.at[idx_pred, "exons"]
        e2 = e2 if isinstance(e2, list) else []

        combined = e1 + e2

        seen_exons = set()
        merged_exons = []

        for exon in combined:

            if not isinstance(exon, dict):
                continue

            exon_id = exon.get("exon_id")

            if exon_id not in seen_exons:
                seen_exons.add(exon_id)
                merged_exons.append(exon)

        df.at[idx_pred, "exons"] = merged_exons


        # GENES

        g1 = row_n.get("genes_names")
        g1 = g1 if isinstance(g1, list) else []

        g2 = df.at[idx_pred, "genes_names"]
        g2 = g2 if isinstance(g2, list) else []

        df.at[idx_pred, "genes_names"] = list(set(g1 + g2))

        # GENOMES

        if "genomes" in df.columns:

            genomes_n = row_n.get("genomes")
            genomes_n = genomes_n if isinstance(genomes_n, list) else []

            genomes_pred = df.at[idx_pred, "genomes"]
            genomes_pred = genomes_pred if isinstance(genomes_pred, list) else []

            df.at[idx_pred, "genomes"] = list(
                set(genomes_n + genomes_pred)
            )

    # Remove compacted nodes
    df_compacted = df[~df["name"].isin(nodes_to_remove)].copy()

    # Normalize list columns
    list_cols = ["genomes", "features", "genes_names", "exons"]
    for col in list_cols:
        if col not in df_compacted.columns:
            continue

        df_compacted[col] = df_compacted[col].apply(
            lambda x: x if isinstance(x, list)
            else ([] if pd.isna(x) else [x])
        )

    return df_compacted



"""
This function compress a graph by removing
linear internal nodes, while allowing orientation inversions
across genomes.
"""
def graph_compression(df, flow_min=0, selected_genomes=None):

    genome_position_cols = [c for c in df.columns if c.endswith("_position")]
    if selected_genomes is not None:
        selected_position_cols = {
            f"{g}_position" for g in selected_genomes
        }
        genome_position_cols = [
            c for c in genome_position_cols
            if c in selected_position_cols
        ]

    flow_map = df.set_index("name")["flow"].to_dict()

    position_data = {
        col: df[["name", col]].dropna().sort_values(col)
        for col in genome_position_cols
    }

    nodes = df["name"].tolist()

    predecessors = {n: set() for n in nodes}
    successors = {n: set() for n in nodes}

    invalid_nodes = set()

    # STEP 1: BUILD GRAPH (UNCHANGED)
    for col, sub in position_data.items():

        ordered_nodes = sub["name"].tolist()

        for i in range(len(ordered_nodes) - 1):

            u = ordered_nodes[i]
            v = ordered_nodes[i + 1]

            if u not in invalid_nodes:
                successors[u].add(v)
                if len(successors[u]) > 2:
                    invalid_nodes.add(u)

            if v not in invalid_nodes:
                predecessors[v].add(u)
                if len(predecessors[v]) > 2:
                    invalid_nodes.add(v)

    # STEP 2: REMOVE INVALID NODES ONLY

    for n in invalid_nodes:
        predecessors.pop(n, None)
        successors.pop(n, None)

    valid_nodes = set(predecessors.keys()) | set(successors.keys())

    # STEP 3: DETECTION
    nodes_to_remove = set()

    for node in valid_nodes:

        preds = predecessors.get(node, set())
        succs = successors.get(node, set())

        curr_flow = flow_map.get(node, 0)

        if curr_flow < flow_min:
            continue

        if not preds or not succs:
            continue

        if len(preds) > 2 or len(succs) > 2:
            continue

        # CASE 1: chain
        if len(preds) == 1 and len(succs) == 1:
            pred = next(iter(preds))
            # IMPORTANT FIX: safe access
            if len(successors.get(pred, set())) == 1:
                nodes_to_remove.add(node)

        # CASE 2: symmetric
        elif len(preds) == 2 and len(succs) == 2:
            neighbors = preds | succs
            if len(neighbors) == 2:

                n1, n2 = tuple(neighbors)

                n1_neighbors = successors.get(n1, set()) | predecessors.get(n1, set())
                n2_neighbors = successors.get(n2, set()) | predecessors.get(n2, set())

                if len(n1_neighbors) == 2 and len(n2_neighbors) == 2:
                    nodes_to_remove.add(node)

    # STEP 4: CLEAN ONLY FOR MERGE
    df_compacted = merge_node_data(df, nodes_to_remove, predecessors)

    return df_compacted



def compute_graph_elements(data, ref_genome, selected_genomes, size_min, all_genomes, all_chromosomes,
                           specifics_genomes=None, color_genomes=[], x_max=1000, y_max=1000, labels=True,
                           min_shared_genome=100, tolerance=0, color_shared_regions=DEFAULT_SHARED_REGION_COLOR,
                           exons=False, exons_color=DEFAULT_EXONS_COLOR, colored_edges_size=5,
                           compression=False, min_flow_compression_value=0, max_nodes_to_visualize=MAX_NODES_TO_VISUALIZE,
                           nodes_size_scale=1, genes_color=None):
    logger.debug(f"Compute elements with ref genome {ref_genome} node min size : {size_min}")

    legend_nodes_size_dict = {
        "size_min": "1",
        "size_max": ""
    }
    if data != None and len(data) > 0:
        if ref_genome is not None and ref_genome != "":
            position_field = ref_genome+"_position"
        else:
            position_field = "mean_pos"
        df = records_to_dataframe(data)
        n_rows = len(df)
        logger.debug(f"Nb rows before filter : {n_rows}")
        df = df[df["size"] >= size_min].copy()
        n_rows = len(df)
        if df.empty:
            return [], 0, legend_nodes_size_dict
        #If compression is set to true then it will compress the graph
        #Nodes with a single predecessor and a single successor, and for which the predecessor
        #has only one outgoing node, are removed
        if compression:
            df = graph_compression(df, min_flow_compression_value, selected_genomes)
        n_rows = len(df)
        logger.debug(f"Nb rows after filter and compression : {n_rows}")
        if n_rows > max_nodes_to_visualize:
            return {}, n_rows, legend_nodes_size_dict

        selected_set = set(selected_genomes)

        df = df[df["genomes"].apply(
            lambda g: bool(set(g) & selected_set)
            if isinstance(g, (list, set, tuple))
            else False
        )].copy()
        #logger.debug(f"Compute elements - dataframe")
        def mean_position(row):
            positions = [row.get(f"{g}_node")
                         for g in row["genomes"] if f"{g}_node" in row]
            positions = [p for p in positions if p is not None]
            return np.mean(positions) if positions else 0

        size_max = df['size'].max()
        size_min = df['size'].min()

        df["displayed_node_size"] = (S_MIN + (df["size"]-size_min)*(S_MAX-S_MIN)/(size_max-size_min))*nodes_size_scale

        legend_nodes_size_dict = {
                "size_min": str(size_min),
                "size_max": str(size_max)
        }


        df["mean_pos"] = df.apply(mean_position, axis=1)
        x_min, x_max_data = df["mean_pos"].min(), df["mean_pos"].max()
        df["x"] = ((df["mean_pos"] - x_min) / max(1,(x_max_data - x_min))) * x_max



        X_SPAN = 40000
        GAP = 40

        min_x = df["mean_pos"].min()
        max_x = df["mean_pos"].max()

        df["x_mean"] = ((df["mean_pos"] - min_x)/ (max(1,max_x - min_x))) * X_SPAN - X_SPAN / 2
        df = df.sort_values("x_mean").reset_index(drop=True)
        #df["x"] = 0.0
        for n in range(len(df)):
            if n == 0:
                df.at[n, "x"] = df.at[n, "x_mean"]
            else:
                prev_x = df.at[n - 1, "x"]
                prev_radius = df.at[n - 1, "displayed_node_size"] / 2
                curr_radius = df.at[n, "displayed_node_size"] / 2
                min_distance = prev_radius + curr_radius + GAP
                delta_x_mean = df.at[n, "x_mean"] - df.at[n - 1, "x_mean"]
                df.at[n, "x"] = prev_x + max(delta_x_mean, min_distance)

        Y_MAX = 2500  # max y coordinate
        Y_MIN = 150  # min y gap between main axe
        FLOW_CENTER = 0.95
        ALPHA = 2.5

        df["y"] = 0.0

        mask = df["flow"] < FLOW_CENTER
        df.loc[mask, "y"] = Y_MIN + (Y_MAX - Y_MIN) * (df.loc[mask, "flow"] / FLOW_CENTER) ** ALPHA

        #df.loc[df.index % 2 == 0, "y"] *= -1

        #This line is due to a pb with cytoscape : if no y with negative and positive values then
        #graph disappear when zooming
        #the first node y value is set negative to avoid this problem
        df.at[df.index[0], "y"] = -Y_MAX
        #logger.debug(df[["x", "y"]].describe())


        df["genome_key"] = df["genomes"].apply(lambda g: "".join(sorted(g)))
        genome_keys = sorted(
            df["genome_key"].drop_duplicates(), key=lambda x: (-len(x), x))
        #y_positions = {k: 1 for k in enumerate(genome_keys)}
        #df["y"] = df["genome_key"].map(y_positions)

        color_map = {k: c for k, c in zip(
            genome_keys, get_color_palette(len(genome_keys)))}

        nodes = []

        size_max_noeud = 10
        #logger.debug(f"Compute elements - begin loop")

        for _, row in df.iterrows():
            node_style = "default"
            classes = "main-node" if row['ref_node'] == row['name'] else "degenerate-node"
            if exons :
                if "features" in row and "exon" in row['features']:
                    node_style="exon"
            node_color = flow_to_rgb(row['flow'],node_style,exons_color)

            # if genes_color:
            #     genes = row.get("genes_names") or []
            #     for gene in genes:
            #         if gene.lower() in genes_color:
            #             node_color = genes_color[gene.lower()]
            #             break

            data_nodes = {
                'data': {
                    'id': row.get('name'),
                    'name': row.get('name'),
                    'displayed_node_size': row.get('displayed_node_size'),
                    'ref_node': row.get('ref_node'),
                    'size': row.get('size'),
                    'flow': row.get('flow'),
                    'genomes': row.get('genomes'),
                    'chromosome': row.get('chromosome'),
                    'sequence': row.get('sequence'),
                    'genes_names': row.get('genes_names'),
                    'features': row.get('features'),
                    'exons': row.get('exons'),
                    'color': node_color,
                    'position': row.get(position_field)
                },
                'position': {'x': row['x'], 'y': row['y']},
                'classes':classes
            }

            nodes.append(data_nodes)



        edges = []
        edges_dict = {}

        for genome in selected_genomes:
            #logger.debug(f"Compute elements - genome {genome}")

            mask = df["genomes"].apply(lambda g: genome in g)
            nodes_g = df.loc[mask, :]

            col = f"{genome}_position"
            if col not in nodes_g.columns:
                continue

            nodes_g = nodes_g.loc[nodes_g[col].notnull()]
            if nodes_g.empty:
                continue

            nodes_g = nodes_g.sort_values(by=col, ascending=True, ignore_index=True)

            #compute source and targets
            sources = nodes_g["name"].iloc[:-1].values
            targets = nodes_g["name"].iloc[1:].values
            positions = nodes_g[col].values
            sizes = nodes_g["size"].values
            annotations = nodes_g["genes_names"].values

            for i, (source, target) in enumerate(zip(sources, targets)):
                edge_key = tuple(sorted([source, target]))
                if edge_key not in edges_dict:
                    edges_dict[edge_key] = {
                        "direction": {f"{source}->{target}": 1},
                        "direction_genome": {genome: f"{source}->{target}"},
                        "flow": 1,
                        "annotations": set(annotations[i] + annotations[i + 1]),
                        "genomes": [genome],
                        "SV": False
                    }
                else:
                    e = edges_dict[edge_key]
                    e["direction"][f"{source}->{target}"] = e["direction"].get(f"{source}->{target}", 0) + 1
                    e["direction_genome"][genome] = f"{source}->{target}"
                    e["flow"] += 1
                    if genome not in e["genomes"]:
                        e["genomes"].append(genome)
                    e["annotations"].update(annotations[i] + annotations[i + 1])

                if genome == ref_genome and "direction_genome_ref" not in edges_dict[edge_key]:
                    edges_dict[edge_key]["direction_genome_ref"] = f"{source}->{target}"

                if (positions[i + 1] - positions[i]) > MAX_GAP and MAX_GAP > sizes[i]:
                    edges_dict[edge_key]["SV"] = True

        #logger.debug(f"Compute elements - specific color")
        # Specific colors for a selected set of genomes
        colored_genomes = {}
        for g, c in zip(all_genomes, color_genomes):
            if c != '#000000':
                colored_genomes[g] = c
        for (s, t), dic in edges_dict.items():
            link_color = 'gray'
            flow = dic["flow"]
            virtual_flow = flow
            if specifics_genomes is not None and len(specifics_genomes) > 0:
                min_required_shared = max(
                    ceil(min_shared_genome * len(specifics_genomes) / 100), 1)
                max_allowed_extra = ceil(tolerance * len(dic["genomes"]) / 100)
                set_specifics_genomes = set(specifics_genomes)
                set_genomes = set(dic["genomes"])

                intersect = set_specifics_genomes.intersection(set_genomes)
                if len(intersect) >= min_required_shared and len(dic["genomes"]) - len(intersect) <= max_allowed_extra:
                    link_color = color_shared_regions
                    virtual_flow = len(all_genomes)
            label = ""
            label_color = "black"
            if labels :
                if not genes_color or len(genes_color) == 0:
                    first_label = True
                    for a in dic["annotations"]:
                        if first_label:
                            label += str(a)
                            first_label = False
                        else:
                            label += ", " + str(a)
                else :
                    first_label = True
                    for a in dic["annotations"]:
                        if a.lower() in genes_color:
                            if first_label:
                                label += str(a)
                            else :
                                label += ", " + str(a)
                            label_color = genes_color[a.lower()]
                            break

            if dic["SV"]:
                line_style = "dashed"
            else:
                line_style = "solid"
            if "direction_genome_ref" in dic :
                source = dic["direction_genome_ref"].split("->")[0]
                target = dic["direction_genome_ref"].split("->")[1]
            else:
                max_dir = 0
                for d in dic["direction"]:
                    if dic["direction"][d] > max_dir:
                        source = d.split("->")[0]
                        target = d.split("->")[1]
                        max_dir = dic["direction"][d]

            edges.append({
                'data': {
                    'source': source,
                    'target': target,
                    'flow': flow,
                    'genomes': dic["genomes"]
                },
                'style': {
                    'line-color': link_color,
                    'line-style':line_style,
                    'target-arrow-color': link_color,
                    'label': label,
                    'color': label_color,
                    'text-rotation': 'autorotate',
                    'text-margin-y': -20,
                    'width': (virtual_flow+int(0.2*len(all_genomes)))/len(all_genomes)*MAX_EDGE_WIDTH
                }
            })
            if len(colored_genomes) > 0:
                i = 0
                for g in list(colored_genomes.keys()):
                    if g in dic["genomes"]:
                        source = dic["direction_genome"][g].split("->")[0]
                        target = dic["direction_genome"][g].split("->")[1]
                        sign = 1 if i % 2 == 0 else -1
                        distance = sign * (20 + 10 * (i // 2))
                        if dic["SV"]:
                            line_style = "dashed"
                        else:
                            line_style = "solid"
                        edges.append({
                            'data': {
                                'source': source,
                                'target': target,
                                'flow': 1,
                                'genomes': [g]
                            },
                            'classes': f'offset-{i}',
                            'style': {
                                'line-color': colored_genomes[g],
                                'line-style' : line_style,
                                'target-arrow-color': colored_genomes[g],
                                'width': colored_edges_size

                            }
                        })
                        i += 1
        logger.debug("nb nodes : " + str(len(nodes)) +
              " - nb edges : " + str(len(edges)))
        return nodes + edges, len(nodes), legend_nodes_size_dict
    else:
        return [], 0, legend_nodes_size_dict


legend = html.Div(
    style={
        "bottom": "20px",
        "left": "20px",
        "background": "rgba(255,255,255,0.95)",
        "border": "1px solid #ccc",
        "borderRadius": "6px",
        "padding": "10px",
        "width": "400px",
        "fontSize": "14px"
    },
    children=[

        html.Div("Legend", style={"fontWeight": "bold", "marginBottom": "8px"}),

        # === NODE SIZE ===

        #New version of node size legend
        svg.Svg(width="320", height="70", children=[

            # ===== MIN =====
            svg.Rect(
                x="10",
                y="20",
                width=str(S_MIN),
                height="18",
                fill="red",
                stroke="black",
                strokeWidth="1",
                rx="4",
                ry="4"
            ),
            svg.Text(
                id="legend-size-min",
                children=str(S_MIN),
                x=str(10 + S_MIN / 2),
                y="29",
                fontSize="11",
                fill="black",
                fontWeight="bold",
                textAnchor="middle",
                dominantBaseline="middle"
            ),

            # ===== ARROW =====
            svg.Text(
                "→",
                x=str(20 + S_MIN + 10),
                y="29",
                fontSize="18",
                fill="black",
                textAnchor="middle",
                dominantBaseline="middle"
            ),

            # ===== MAX =====
            svg.Rect(
                x=str(50 + S_MIN),
                y="20",
                width=str(S_MAX),
                height="18",
                fill="red",
                stroke="black",
                strokeWidth="1",
                rx="4",
                ry="4"
            ),
            svg.Text(
                id="legend-size-max",
                children=str(S_MAX),
                x=str(50 + S_MIN + S_MAX / 2),
                y="29",
                fontWeight="bold",
                fontSize="11",
                fill="black",
                textAnchor="middle",
                dominantBaseline="middle"
            ),

            # ===== LABEL =====
            svg.Text(
                "Node size: min → max (sequence length)",
                x="10",
                y="60",
                fontSize="14",
                fill="black"
            )
        ]),


        # === NODE COLOR (BLUE → RED) ===

        #new version with bar
        svg.Svg(width="300", height="60", children=[

            svg.Rect(x="10", y="10", width="45", height="20", fill="blue"),
            svg.Rect(x="55", y="10", width="45", height="20", fill="#2f2f9f"),
            svg.Rect(x="100", y="10", width="45", height="20", fill="#7f007f"),
            svg.Rect(x="145", y="10", width="45", height="20", fill="#bf003f"),
            svg.Rect(x="190", y="10", width="45", height="20", fill="red"),

            svg.Rect(
                x="10", y="10",
                width="225",
                height="20",
                fill="none",
                stroke="black",
                strokeWidth="0.5"
            ),

            svg.Text(
                "Color: few individuals → all individuals",
                x="10",
                y="45",
                fontSize="14",
                dominantBaseline="middle"
            )
        ]),


        # === REPEATED NODE ===
        svg.Svg(width="300", height="70", children=[
            svg.Ellipse(cx="40", cy="25", rx="28", ry="12", fill="red", stroke="black"),
            svg.Text("Repeated node", x="80", y="30", fontSize="14")
        ]),


        # === EDGE WIDTH ===
        svg.Svg(width="360", height="30", children=[
            svg.Line(
                x1="10", y1="15", x2="60", y2="15",
                stroke="#333",
                strokeWidth="5"
            ),
            svg.Text(
                "Edge width ∝ number of individuals",
                x="80", y="20",
                fontSize="14"
            )
        ]),

        # === DASHED EDGE ===
        svg.Svg(width="360", height="30", children=[
            svg.Line(
                x1="10", y1="15", x2="60", y2="15",
                stroke="#333",
                strokeWidth="3",
                strokeDasharray="6,4"
            ),
            svg.Text(
                "Dashed edge: long genomic distance",
                x="80", y="20", fontSize="14"
            )
        ]),
    ]
)



def layout(data=None, initial_size_limit=10):
    all_genomes = get_genomes()
    all_genomes.sort()
    all_chromosomes = get_chromosomes()
    features = get_annotations_features()
    if all_genomes :
        max_label_len = max(len(s) for s in all_genomes)
    else :
        max_label_len = 30
    min_item_width = min(max_label_len * 8 + 15, 350)
    if data != None:
        elements, nodes_count, legend_nodes_size_dict = compute_graph_elements(
            data, all_genomes, initial_size_limit, all_genomes, all_chromosomes, [], [])
    else:
        elements = []
    size_max = 500

    return html.Div([
        dcc.Store(id='update_graph_command_storage', storage_type='memory'),
        html.Div([
            html.H2("PANABYSS"),

            html.Details([
                html.Summary("ℹ️ Click here to display help"),
                html.P("PanAbyss is an application that allows you to view and manipulate pan-genomic data. "
                       "This page allows you to view portions of the pan-genome:"),
                html.Ul([
                    html.Li(
                        "Haplotype selection: choose the haplotype that will be used for annotations and coordinates."),
                    html.Li(
                        "Select a chromosome: searches will only be performed on this chromosome."),
                    html.Li("Choose from the following options:"),
                    html.Ul([
                        html.Li("Select the start and end of the region on the selected chromosome. If there is no defined end, then the searched region will be from the defined start to the end of the pangenome. The region should not be too large, otherwise the display will take too long or will not be possible."),
                        html.Li("Search by annotation with a gene name or ID.")
                    ]),
                    html.Li("Select the haplotypes to be viewed: it is possible to exclude some haplotypes. "
                            "In this case, nodes containing only these haplotypes will not be displayed."),
                    html.Li("You can download the current graph by clicking on 'as jpg', 'as png' or 'as svg' button. "
                            "Graphs png and jpg are saved into ./export/graphs directory, while svg graphs are downloaded.")
                ]),
                html.P("Display settings:"),
                html.Ul([
                    html.Li(
                        "The minimum node size slider allows you to hide nodes smaller than this value."),
                    html.Li(
                        "Select layout: allows you to change the representation algorithm if it is not suitable."),
                    html.Li(
                        "Hide labels: allows you to hide labels (they can take up too much space)."),
                    html.Li(
                        "Show exons: allows you to color the node linked to an exon annotation (color is defined by user)."),
                    html.Li(
                        "Colored edges size: set the edge size for colored haplotypes."),
                    html.Li(
                        "Graph compression: when a min node size is set greater than 0, this will compact the linear parts of the graph (i.e. nodes connected to exactly 2 nodes will be compacted. It is possible to set the minimal percentage of individuals for compressing a node (for instance to compress only the core nodes)."),
                    html.Li("Search shared paths : if the box is checked, the display is modified to allow the selection of haplotypes for which you want to view shared links. This display can be configured:"),
                    html.Ul([
                        html.Li(
                            "Selection of haplotypes for which common links are sought"),
                        html.Li("Selection of the minimum percentage of selected haplotypes that must be present on the link "
                                "(for example, if the value is set to 50 and 10 haplotypes are selected, at least 5 of the haplotypes must be on the link). "
                                "If the value is zero, then at least one of the selected haplotypes will be required."),
                        html.Li("Tolerance: links containing fewer than [(tolerance / 100) × number of haplotypes passing through this node] "
                                "unselected haplotypes will be reported."),
                        html.Li("Link color: choice of color for reported links")
                    ]),
                    html.Li(
                        "Search shared paths : if the box is unchecked, then it is possible to select a specific color for a each haplotype."),
                    html.Li(
                        "If one of this option is modified, click on the Update graph button."),
                ]),
                html.P("Display description :"),
                html.Ul([
                    html.Li(
                        "By clicking on a node or a link, data of this node will be displayed under the 'update graph' button."),
                    html.Li(
                        "If annotations exists in the visualized region they will be concatenated under the node / link description area."),
                    html.Li("Graph description :"),
                    html.Ul([
                        html.Li(
                            "Node shape : a node is drawn as a round rectangle (can be configured in about page), unless if it's a repeated node in which case it will be displayed as an ellipse."),
                        html.Li("Node color : The color of the nodes ranges from blue to red. The bluer the color, the less common the node (node associated with only one or a small number of haplotypes). Conversely, red nodes are those of the core genome shared by all haplotypes."),
                        html.Li(
                            "Node size : the size of a node is proportionnal to the size of the sequence associated."),
                        html.Li(
                            "Link size : the size of a link is proportional to the number of haplotypes passing through that link."),
                        html.Li(
                            "Link shape : the dashed links represent structural variations for which nodes are likely missing for the associated individuals."),
                        html.Li(
                            "Link direction : The direction is defined as follows: if the reference individual traverses the link, "
                            "their direction is used; otherwise, the direction most commonly observed among the individuals passing through the link is chosen. "
                            "For color-coded individuals, the specific direction for each is indicated by their assigned color."),
                    ]),
                    html.Li("Zoom :"),
                    html.Ul([
                        html.Li(
                            "To zoom in : first select the nodes you want to zoom in by holding down the left mouse button and Shift key. Then push the 'zoom on selection' button."),
                        html.Li(
                            "Zoom out 2000 bp : zoom out the region by starting 1000 bp before the current region and 1000 bp after."),
                        html.Li(
                            "To retrieve the initial region just push the 'Reset zoom' button."),
                    ]),
                ]),
                html.P("Query params:"),
                html.Ul([
                    html.Li([
                        "It is possible to directly visualize a region by specifying query parameters in the request."]),

                    html.Li([
                    "For example, to visualize the region from 10,330,185 to 10,529,460 on chromosome 21 of the haplotype CHM13_0, "
                    "you can use the following query parameter: ",
                    html.Code("?haplotype=CHM13_0&chromosome=21&start=10330185&end=10529460")
                        ]),
                    html.Li(["Here is the list of available parameters:"]),
                    html.Ul([
                        html.Li([
                            "Mandatory parameters:",
                            html.Ul([
                                html.Li("haplotype"),
                                html.Li("chromosome")
                            ])
                        ]),
                        html.Li([
                            "Search parameters:",
                            html.Ul([
                                html.Li("start / end"),
                                html.Li("geneName"),
                                html.Li("geneId")
                            ])
                        ])
                    ])

                ]),
            ], style={"marginBottom": "20px"}),


        ]),
        # Upper block : settings
        html.Div([


            # Left block
            html.Div([
                html.Div([
                    html.Div([
                        html.Label("Reference haplotype", title="Select an haplotype to search / display annotation and to define genomic coordinates to search region.",
                                   style={'display': 'block', 'marginBottom': '5px'}),
                        dcc.Dropdown(
                            id='genomes-dropdown',
                            options=[{'label': genome, 'value': genome}
                                     for genome in all_genomes],
                            # value=all_genomes[0],
                            clearable=False,
                            style={'width': '300px'}
                        )
                    ], style={'marginRight': '30px'}),
                    html.Div([
                        html.Label("Chromosome", title="The regions / annotations will be related only to this chromosome.",
                                   style={'display': 'block', 'marginBottom': '5px'}),
                        dcc.Dropdown(
                            id='chromosomes-dropdown',
                            options=[{'label': chrom, 'value': chrom}
                                     for chrom in all_chromosomes],
                            clearable=False,
                            placeholder="Choose a chromosome",
                            # value="1",
                            style={'width': '200px'}
                        )
                    ])
                ], style={'display': 'flex', 'padding': '20px', 'border': '1px solid #ccc', 'minWidth': '300px', 'boxSizing': 'border-box'}),


                html.Div([
                    html.Div([
                        #Div search and display
                        html.Div([
                                #div search options
                            html.H5("Search region", style={'textAlign': 'left', 'marginBottom': '15px'}),
                            html.Div([
                                # html.Label("Start : ", title="Start on the selected haplotype / chromosome."),
                                # dcc.Input(id='start-input', type='number', style={'width': '100px', 'marginRight': '10px'}),
                                # html.Label("End : ", title="End on the selected haplotype / chromosome."),
                                # dcc.Input(id='end-input', type='number', style={'width': '100px', 'marginRight': '20px'})

                                html.Div([
                                    html.Div([
                                        html.Label("Start : ", title="Start on the selected haplotype / chromosome."),
                                        dcc.Input(id='start-input', type='text', inputMode='numeric', pattern='[0-9]*', style={'width': '140px','minWidth': '120px'})
                                    ], style={'display': 'flex', 'alignItems': 'center', 'gap': '5px'}),

                                    html.Div([
                                        html.Label("End : ", title="End on the selected haplotype / chromosome."),
                                        dcc.Input(id='end-input', type='text', inputMode='numeric', pattern='[0-9]*', style={'width': '140px','minWidth': '120px'})
                                    ], style={'display': 'flex', 'alignItems': 'center', 'gap': '5px'})

                                ], style={
                                    'display': 'flex',
                                    'flexWrap': 'wrap',
                                    'gap': '10px'
                                })


                            ], style={'marginBottom': '10px'}),


                            html.Div(
                                style={
                                    'display': 'flex',
                                    'flexDirection': 'column',
                                    'gap': '6px',
                                    'marginBottom': '20px',
                                    'width': '100%'
                                },

                                children=[

                                    # Label
                                    html.Label(
                                        "Annotations :",
                                        title="Search by annotations."
                                    ),

                                    html.Div(
                                        style={
                                            'display': 'flex',
                                            'flexDirection': 'row',
                                            'alignItems': 'center',
                                            'gap': '12px',
                                            'flexWrap': 'wrap'
                                        },

                                        children=[

                                            dcc.Dropdown(
                                                id='features-dropdown',

                                                options=[
                                                    {'label': f, 'value': f}
                                                    for f in features
                                                ],
                                                # options=[
                                                #             {'label': f"{feature}_name", 'value': f"{feature}_name"}
                                                #             for feature in features
                                                #         ] + [
                                                #             {'label': f"{feature}_id", 'value': f"{feature}_id"}
                                                #             for feature in features
                                                #         ],
                                                value='gene' if 'gene' in features else None,
                                                clearable=False,
                                                placeholder="Feature",

                                                style={
                                                    'width': 'max-content',
                                                    'minWidth': '160px',
                                                    'maxWidth': '240px'
                                                }
                                            ),

                                            dcc.Input(
                                                id='feature-input',
                                                type='text',
                                                placeholder='Name or id to search',
                                                debounce=True,

                                                style={
                                                    'width': '200px',
                                                    'minWidth': '140px',
                                                    'maxWidth': '260px',

                                                    'height': '38px',
                                                    'padding': '0 10px',
                                                    'boxSizing': 'border-box'
                                                }
                                            )
                                        ]
                                    )
                                ]
                            ),



                            html.Button('Search', id='search-button',
                                        n_clicks=0, style={'marginTop': '10px'}),
                            dcc.Loading(
                                id="loading-search-msg",
                                #type="circle",
                                children=html.Div(id="search-message")
                            ),
                        ],style={
                            'width': '45%',
                            'padding': '10px',
                            'boxSizing': 'border-box'
                        }),
                        html.Div([
                            #Div displayed region
                            html.H4("Displayed region", style={'marginBottom': '10px'}),
                            html.Div([
                                html.Div(id='displayed-region-container', children=[
                                    html.Div("No region selected.", id='no-region-msg',
                                             style={'fontStyle': 'italic', 'color': '#777'})
                                ])
                            ], style={
                                'border': '1px solid #ccc',
                                'borderRadius': '5px',
                                'padding': '10px',
                                'backgroundColor': '#f9f9f9',
                                'minHeight': '100px'
                            }),
                        ],
                            style={
                                'width': '45%',
                                'padding': '10px',
                                'borderLeft': '2px solid #ddd',
                                'boxSizing': 'border-box'
                            }),
                    ],style={
                        'display': 'flex',
                        'flexDirection': 'row',
                        'justifyContent': 'space-between',
                        'alignItems': 'flex-start',
                        'width': '100%',
                    }),
                    html.Div([
                        html.Div([
                            html.Label(
                                "Haplotypes to visualize :",
                                title="Nodes containing only unselected haplotypes won't be displayed.",
                                style={'marginBottom': '5px'}
                            ),

                            html.Button(
                                "Select all",
                                id="select_all_genomes",
                                n_clicks=0,
                                style={"marginLeft": "10px", "marginRight": "5px"}
                            ),

                            html.Button(
                                "Unselect all",
                                id="unselect_all_genomes",
                                n_clicks=0
                            ),
                        ], style={"display": "flex", "alignItems": "center"}),

                        dcc.Checklist(
                            id="genome_selector",
                            options=[{"label": g, "value": g} for g in all_genomes],
                            value=all_genomes,
                            labelStyle={
                                "margin": "0px",
                                "padding": "0px",
                                "lineHeight": "1",
                                "fontSize": "14px",
                                "display": "flex",
                                "alignItems": "center",
                                "gap": "4px",
                                "whiteSpace": "nowrap",
                            },
                            style={
                                "display": "grid",
                                "gridTemplateColumns": f"repeat(auto-fit, minmax({min_item_width}px, max-content))",
                                "columnGap": "12px",
                                "rowGap": "2px",
                                "width": "100%",
                            }
                        )
                    ],
                        style={'marginBottom': '20px'}),



                ],
                    style={'marginBottom': '20px'}
                )
            ], style={'flex': '1', 'padding': '20px', 'border': '1px solid #ccc'}),

            # Right block
            html.Div([

                html.Div([
                    html.Label("Minimal size of nodes :",
                               title="Nodes with size under this value won't be displayed."),

                    dcc.Slider(
                        id='size_slider',
                        min=1,
                        max=size_max,
                        step=1,
                        marks={i: str(i) for i in range(0, size_max + 1, int(size_max/10))},
                        value=DEFAULT_SIZE_VALUE,
                        updatemode="mouseup",
                        tooltip={"placement": "bottom",
                                 "always_visible": False},
                    ),

                    html.Div(id='size_stats', style={'marginTop': '10px'})

                ]),
                html.Div(id='nb-noeuds', style={'margin': '10px'}),

                html.Div([

                    html.Div([
                        html.Div([
                            html.H4(
                                "Layout",
                                title="The layout computes the graph...",
                                style={'margin': '0'}
                            ),
                            dcc.Dropdown(
                                id='layout-dropdown',
                                options=[
                                    {'label': 'fcose', 'value': 'fcose'},
                                    {'label': 'dagre', 'value': 'dagre'},
                                    {'label': 'preset', 'value': 'preset'}
                                ],
                                value='fcose',
                                clearable=False,
                                style={'minWidth': '120px'}
                            ),
                        ], style={
                            'display': 'flex',
                            'alignItems': 'center',
                            'gap': '8px',
                            'whiteSpace': 'nowrap'
                        }),

                        html.Div([
                            dcc.Checklist(
                                id='graph-compression',
                                options=[
                                    {
                                        'label': 'Compress nodes with at least ',
                                        'value': 'graph_compression'
                                    }
                                ],
                                value=[],
                                style={'display': 'inline-block'}
                            ),

                            dcc.Input(
                                id='min-flow',
                                type='number',
                                min=0,
                                max=100,
                                step=1,
                                value=0,
                                style={
                                    'width': '70px',
                                    'marginLeft': '5px'
                                }
                            ),

                            html.Span(" % of individuals. "),
                        ]),


                    ], style={
                        'display': 'flex',
                        'alignItems': 'center',
                        'flexWrap': 'wrap',
                        'gap': '12px',
                        'marginBottom': '10px'
                    }),


                    #Size sliders
                    html.Div([
                        # Colored edge size slider
                        html.Div([
                            html.Label("Colored edges size", style={'marginBottom': '0'}),
                            dcc.Slider(
                                id='colored-edge-size-slider',
                                min=1,
                                max=20,
                                step=1,
                                value=5,
                                marks={i: str(i) for i in range(0, 25, 5)}
                            )
                        ], style={
                            'display': 'flex',
                            'alignItems': 'center',
                            'gap': '15px',
                            'width': '300px',
                            'paddingRight': '20px',
                            'marginBottom': '10px'
                        }),

                        # Node scale size slider
                        html.Div([
                            html.Label("Nodes size scale", style={'marginBottom': '0'}),
                            dcc.Slider(
                                id='node-size-scale-slider',
                                min=0.1,
                                max=5,
                                step=0.1,
                                value=1,
                                marks={
                                    0.1: '0.1',
                                    1: '1',
                                    2: '2',
                                    3: '3',
                                    4: '4',
                                    5: '5'
                                }
                            )
                        ], style={
                            'display': 'flex',
                            'alignItems': 'center',
                            'gap': '15px',
                            'width': '300px',
                            'paddingRight': '20px',
                            'marginBottom': '10px'
                        }),
                    ],
                        style={
                            'display': 'flex',
                            'flexWrap': 'wrap',
                            'gap': '20px',
                            'alignItems': 'center'
                        }),

                    # === Second line: checkboxes and color picker ===
                    html.Div([
                        dcc.Checklist(
                            options=[{
                                'label': 'Hide annotations',
                                'value': 'hide',
                                'title': 'Uncheck if annotations takes too much space on the graph.'
                            }],
                            id='show-labels',
                            style={'marginRight': '30px'}
                        ),
                        dcc.Checklist(
                            options=[{
                                'label': 'Nodes names',
                                'value': 'nodes_names',
                                'title': 'Check to display nodes names.'
                            }],
                            id='nodes-names',
                            style={'marginRight': '30px'},
                            value=[]
                        ),

                        dcc.Checklist(
                            options=[{
                                'label': 'Show exons',
                                'value': 'exons',
                                'title': 'Check to display exons.'
                            }],
                            id='show-exons',
                            style={'marginRight': '10px'},
                            value=[]
                        ),

                        dbc.Input(
                            id='exon-color-picker',
                            type='color',
                            value=DEFAULT_EXONS_COLOR,
                            style={'width': '25px', 'minWidth':'25px', 'height': '25px', 'marginRight': '10px'}
                        ),





                    ], style={'display': 'flex', 'alignItems': 'center', 'gap': '8px'}),
                ], style={'display': 'flex', 'flexDirection': 'column', 'marginBottom':'10px'}),


                html.Div(
                    style={
                        'display': 'flex',
                        'alignItems': 'center',
                        'gap': '30px',  # espace entre checklist et le mini-div
                        'marginBottom': '20px'
                    },
                    children=[
                        dcc.Checklist(
                            options=[{'label': 'Shared paths', 'value': 'shared'}],
                            id='shared-mode'
                        )


                    ]
                ),
            
                html.Div([

                    html.Div([
                        html.Label(
                            "Vizualize shared paths :", title="Detect the links through which the selected haplotypes pass - the following input fields allow you to refine the search. ", style={'marginBottom': '20px'}),

                        dcc.Checklist(
                            id="specific-genome_selector",
                            options=[{"label": g, "value": g} for g in all_genomes],
                            inline=False,

                            labelStyle={
                                "display": "flex",
                                "alignItems": "center",
                                "gap": "4px",
                                "fontSize": "14px",
                                "whiteSpace": "nowrap",
                                "lineHeight": "1"
                            },

                            style={
                                "display": "grid",
                                "gridTemplateColumns": f"repeat(auto-fit, minmax({min_item_width}px, max-content))",
                                "columnGap": "2px",
                                "rowGap": "0px",
                                "width": "100%",
                                "lineHeight": "1",
                                "marginBottom": "5px"
                            }
                        ),
                        html.Label("Min (%) of shared haplotypes : ", title="Min (%) of shared haplotypes = M. Number of selected haplotypes = N. To detect a shared node it must contains almost (M/100) x N of the selected haplotypes. If M = 0 then the minimum number of selected haplotypes will be 1."),
                        dcc.Input(id='min_shared_genomes-input', type='number',
                                  value=100, style={'width': '100px', 'marginRight': '10px', "marginBottom": "5px"}),
                        html.Label("Tolerance (%) : ", title="Tolerance = T. Number of haplotypes on a node = n. To detect a shared node it must contains less than (T/100) x n of the non selected haplotypes. If T = 0 then detected nodes should contain only selected haplotypes."),
                        dcc.Input(id='tolerance-input', type='number', value=0,
                                  style={'width': '100px', 'marginRight': '20px'}),

                        html.Label(
                            "Link color", title="This color will be used to color links between detected shared nodes."),

                        dbc.Input(
                            id='shared-region-color-picker',
                            type='color',
                            # value=DEFAULT_SHARED_REGION_COLOR,
                            style={'width': '25px',
                                   'height': '25px', 'marginLeft': '10px'}
                        ),

                    ], id='shared-checklist-container'),

                    html.Div([
                        html.Div([
                            dbc.Input(
                                id={'type': 'color-picker', 'index': i},
                                type='color',
                                value="#000000",
                                style={
                                    'width': '22px',
                                    'height': '22px',
                                    'minWidth': '22px',
                                    'padding': '0',
                                    'border': 'none',
                                    'borderRadius': '5px',
                                    'overflow': 'hidden',
                                    'cursor': 'pointer'
                                }
                            ),
                            html.Label(
                                s,
                                style={
                                    "marginLeft": "5px",
                                    "fontSize": "12px",
                                    "whiteSpace": "nowrap"
                                }
                            )
                        ],
                            style={
                                'display': 'flex',
                                'alignItems': 'center',
                                'width': f'{min_item_width + 15}px',
                                'padding': '1px 2px',
                                "marginBottom": "5px",
                                "overflow": "visible"
                            })
                        for i, s in enumerate(all_genomes)
                    ],

                        style={
                            'display': 'flex',
                            'flexWrap': 'wrap',
                            'gap': '2px 6px'
                        }, id='color-picker-container')

                ], id='sample-controls'),
                html.Button("Update graph", id="update-btn",
                            n_clicks=0, style={'marginTop': '10px'}),

                html.Div(html.H4(id='node-info', style={'margin': '10px'})),
                html.Div(html.Label("Annotations in the region:", title="Compiles all annotations for the displayed nodes.", style={
                         'marginBottom': '5px'})),
                html.Div(html.H4(id='annotations-info',
                         style={'margin': '10px'})),
                html.Div(id="gene-color-picker-container")
            ], style={'flex': '1', 'padding': '20px', 'border': '1px solid #ccc', 'marginLeft': '20px', 'minWidth': '300px',
                      'boxSizing': 'border-box', 'display': 'flex', 'flexDirection': 'column'})
        ], style={
            'display': 'flex',
            'flexDirection': 'row',
            'flexWrap': 'wrap',
            'justifyContent': 'space-between'
        }),
        html.Div([
            html.Div([
                # Button on top of the graph
                html.Button(
                    "🔍 Zoom on selection",
                    id='btn-zoom',
                    title='Before using this button, nodes must be selected by holding left mouse button and Shift key.',
                    style={
                        'whiteSpace': 'nowrap'
                    }
                ),
                html.Button(
                    "🔄 Reset Zoom",
                    id='btn-reset-zoom',
                    style={
                        'whiteSpace': 'nowrap'
                    }
                ),
                html.Button(
                    "🔄 Zoom out 2000 bp",
                    id='btn-zoom-out',
                    style={
                        'whiteSpace': 'nowrap'
                    }
                ),
                html.Button(
                    "Toggle Legend",
                    id='btn-toggle-legend',
                    style={
                        'whiteSpace': 'nowrap'
                    }
                ),
                # Download graph div
                html.Div([
                    html.Label(
                        'Download graph:',
                        title="PNG and JPG files will be saved into ./export/graphs directory while SVG files are downloaded.",
                        style={"marginLeft": "10px"}),
                    html.Button("as jpg", id="btn-save-jpg",
                                style={"marginLeft": "10px"}),
                    html.Button("as png", id="btn-save-png",
                                style={"marginLeft": "10px"}),
                    html.Button("as svg", id="btn-save-svg",
                                style={"marginLeft": "10px"}),
                    dcc.Download(id="download-graph"),
                    # html.Button("as svg", id="btn-save-svg", style={"marginLeft":"10px"}),
                    html.Div(id="dummy-output")
                ]),
            ], style={'marginBottom': '10px', 'display': 'flex', 'gap': '10px'}),

            html.Div(
                legend,
                id='legend-div',
                style={'marginBottom': '10px','display': 'none'}
            )
        ]),

        # Graph block
        # cyto.Cytoscape(
        #     id='graph',
        #     style={'width': '100%', 'height': '1000px'},
        #     zoomingEnabled=True,
        #     userZoomingEnabled=True,
        #     userPanningEnabled=True,
        #     wheelSensitivity=0.1,
        #     boxSelectionEnabled=True,
        #
        # )

        html.Div(
            [
                cyto.Cytoscape(
                    id='graph',
                    style={
                        'width': '100%',
                        'height': '1000px',
                    },
                    zoomingEnabled=True,
                    userZoomingEnabled=True,
                    userPanningEnabled=True,
                    wheelSensitivity=0.1,
                    boxSelectionEnabled=True,
                )
            ],
            style={
                'border': '1px solid #e5e7eb',
                'borderRadius': '12px',
                'padding': '10px',
                'backgroundColor': '#f8f9fa',
                'boxShadow': '0 2px 6px rgba(0,0,0,0.04)',
            }
        )


        ])

#callback to show / hide legend
@app.callback(
    Output('legend-div', 'style'),
    Input('btn-toggle-legend', 'n_clicks'),
    State('legend-div', 'style')
)
def toggle_legend(n_clicks, current_style):
    if n_clicks is None:
        return current_style
    new_display = 'none' if current_style.get('display', 'block') != 'none' else 'block'
    current_style['display'] = new_display
    return current_style


#Callback to selected all / unselect all genomes
@app.callback(
    Output("genome_selector", "value"),
    Input("select_all_genomes", "n_clicks"),
    Input("unselect_all_genomes", "n_clicks"),
    State("genome_selector", "options"),
    prevent_initial_call=True
)
def update_genome_selection(n_select, n_unselect, options):
    ctx = dash.callback_context

    if not ctx.triggered:
        return dash.no_update

    button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    all_values = [opt["value"] for opt in options]

    if button_id == "select_all_genomes":
        return all_values

    elif button_id == "unselect_all_genomes":
        return []

    return dash.no_update

# Callback to get nodes or link info when clicking on it
@app.callback(
    Output('node-info', 'children'),
    Input('graph', 'tapNodeData'),
    Input('graph', 'tapEdgeData')
)
def display_element_data(node_data, edge_data):
    triggered_id = ctx.triggered_id

    # -------------------------
    # EDGE
    # -------------------------
    if triggered_id == 'graph' and "prop_id" in ctx.triggered[0] and ctx.triggered[0]['prop_id'] == 'graph.tapEdgeData' and edge_data:
        return html.Div([
            html.Div(f"Selected link : {edge_data.get('source')} → {edge_data.get('target')}"),
            html.Div(f"• Flow : {edge_data.get('flow')}"),
            html.Div(f"• Haplotypes : {', '.join(edge_data.get('genomes', []))}")
        ])

    # -------------------------
    # NODE
    # -------------------------
    elif triggered_id == 'graph' and "prop_id" in ctx.triggered[0] and ctx.triggered[0]['prop_id'] == 'graph.tapNodeData' and node_data:
        exon_spans = []
        for exon in sorted(
                node_data.get("exons", []),
                key=lambda e: e.get("exon_id") or ""
        ):
            exon_id = exon.get("exon_id")
            if exon_id is None:
                continue

            tooltip_lines = []
            # transcripts
            transcripts = exon.get("transcript_ids", [])
            if transcripts:
                tooltip_lines.append("Transcripts: " + ", ".join(transcripts))

            # coordinates
            tooltip_lines.append(f"{exon.get('start')} - {exon.get('end')}")

            exon_spans.append(
                html.Span(
                    exon_id,
                    title="\n".join(tooltip_lines),
                    style={
                        "marginRight": "6px",
                        "textDecoration": "underline",
                        "cursor": "pointer",
                        "fontSize": "14px"
                    }
                )

            )
        return html.Div([
            html.B(
                f"Selected node : {node_data.get('label', node_data.get('name'))}"
                f" - Ref node : {node_data.get('ref_node')}"
            ),
            html.Div([
                html.Span(f"• Size : {node_data.get('size')} | "),
                html.Span(f"Position : {node_data.get('position')} | "),
                html.Span(f"Flow : {node_data.get('flow')}")
            ]),
            html.Div([
                html.Span(f"• Haplotypes : {', '.join(node_data.get('genomes', []))}")
            ]),
            html.Div([
                html.Span(f"• Genes : {', '.join(node_data.get('genes_names', []))}")
            ]),
            html.Div([
                html.Span(f"• Features : {', '.join(node_data.get('features', []))}")
            ]),
            #Exon with tooltip
            # html.Div([
            #     html.Span("• Exons : "),
            #     html.Span(exon_spans)
            # ]),
            html.Div([
                html.Span("• Exons : "),
                html.Div(
                    exon_spans,
                    style={
                        "display": "flex",
                        "flexWrap": "wrap",
                        "gap": "6px",
                        "marginTop": "4px"
                    }
                )
            ]),
            html.Br(),
            html.Div([
                html.B("Sequence (first 1000 bp only):")
            ]),
            html.Pre(
                node_data.get("sequence"),
                style={
                    "whiteSpace": "pre-wrap",
                    "overflowWrap": "break-word"
                }
            )
        ])

    return "Click on a node or link to display data."


#Function to construct the region information
def get_displayed_div(start, end, feature_name, feature_value):
    no_display = True
    def info_line(label, value):
        if value is not None and value != "":
            return html.Div([
                html.Span(f"{label}: ", style={'font-weight': 'bold', 'marginRight': '5px'}),
                html.Span(str(value))
            ], style={'margin-bottom': '5px'})
    info_rows = []

    if start is not None :
        start_txt = str(start)
        no_display = False
    else:
        start_txt = str(0)
    if end is not None:
        end_txt = str(end)
        no_display = False
    else:
        end_txt = "-"
    info_rows.append(
        html.Div([
            html.Span("Region: ", style={'font-weight': 'bold', 'marginRight': '5px'}),
            html.Span(f"{start_txt} — {end_txt}")
        ], style={'margin-bottom': '5px'})
    )
    if feature_name and feature_name != "" and feature_value and feature_value != "":
        info_rows.append(info_line(feature_name, feature_value))
        no_display = False
    if no_display :
        displayed_div = html.Div("No region selected.", style={'fontStyle': 'italic', 'color': '#777'})
    else:
        displayed_div = info_rows
    return displayed_div


#Callback to set size of nodes in legend
@app.callback(
    Output("legend-size-min", "children"),
    Output("legend-size-max", "children"),
    Input("home-page-legend-store", "data")
)
def update_legend_size(legend_data):
    if not legend_data:
        return "0", "0", "0"

    # sécurité si structure inattendue
    try:
        min_val = legend_data.get("size_min", "1")
        max_val = legend_data.get("size_max", "")
    except AttributeError:
        return "0", ""

    return (
        str(min_val),
        str(max_val)
    )

#This function build annotations to display for the global region
def build_annotations(nodes_data, genes_color=None):
    genes_set = set()
    genes_to_transcripts = {}

    region_transcripts = {}

    for node_data in nodes_data.values():

        # get genes names
        if "genes_names" in node_data:
            for a in node_data["genes_names"]:
                genes_set.add(a.upper())

        # get transcripts
        for transcript in node_data.get("transcripts", []):

            transcript_id = transcript.get("transcript_id")
            if transcript_id is None:
                continue
            transcript_id = transcript_id.upper()

            if transcript_id not in region_transcripts:
                transcript_gene = transcript.get("gene_name")
                if transcript_gene is not None:
                    transcript_gene = transcript_gene.upper()

                region_transcripts[transcript_id] = {
                    "start": transcript["start"],
                    "end": transcript["end"],
                    "exons": [],
                    "gene_name": transcript_gene
                }

                # mapping gene -> transcripts (tooltip)
                if transcript_gene is not None:
                    if transcript_gene not in genes_to_transcripts:
                        genes_to_transcripts[transcript_gene] = []
                    genes_to_transcripts[transcript_gene].append(transcript_id)

        # get exons
        for exon in node_data.get("exons", []):

            exon_id = exon.get("exon_id")
            if exon_id is None:
                continue
            exon_id = exon_id.upper()
            transcript_ids = exon.get("transcript_ids", [])

            for transcript_id in transcript_ids:

                if transcript_id is None:
                    continue

                if transcript_id not in region_transcripts:
                    continue

                region_transcripts[transcript_id]["exons"].append({
                    "exon_id": exon_id,
                    "start": exon["start"],
                    "end": exon["end"]
                })

    # construct ordered exon list (toujours utile pour futur usage)
    for transcript_id, transcript_data in region_transcripts.items():

        exons = transcript_data.get("exons", [])
        unique_exons = {}
        for exon in exons:
            exon_id = exon.get("exon_id")
            if exon_id is None:
                continue

            if exon_id not in unique_exons:
                unique_exons[exon_id] = exon
            else:
                existing = unique_exons[exon_id]
                existing["start"] = min(existing["start"], exon["start"])
                existing["end"] = max(existing["end"], exon["end"])

        exons = list(unique_exons.values())

        # sort
        exons.sort(key=lambda x: (x["start"], x["end"], x["exon_id"]))

        exon_ids = [e["exon_id"] for e in exons]
        transcript_data["exon_chain"] = " - ".join(exon_ids)
        transcript_data["exons"] = exon_ids

    genes_list = sorted(list(genes_set))

    #Gene with tooltip
    genes_html = []

    for gene in genes_list:
        transcripts_for_gene = genes_to_transcripts.get(gene, [])
        tooltip_lines = ["transcript_id (start - end) exon_id"]
        for tid in transcripts_for_gene:
            t = region_transcripts.get(tid)

            if t is None:
                continue

            tooltip_lines.append(
                f"{tid} ({t['start']}-{t['end']}) {t.get('exon_chain', '')}"
            )
        tooltip_text = "\n".join(tooltip_lines) if tooltip_lines else "No transcripts"

    #     genes_html.append(
    #         html.Span(
    #             gene,
    #             title=tooltip_text,
    #             style={
    #                 "textDecoration": "underline",
    #                 "cursor": "pointer",
    #                 "whiteSpace": "nowrap"
    #             }
    #         )
    #     )

        genes_html.append(
            html.Div(
                [
                    dbc.Input(
                        id={'type': 'gene-color', 'gene': gene},
                        type='color',
                        value=genes_color.get(gene.lower(), "#000000") if genes_color else "#000000",
                        style={
                            'width': '22px',
                            'height': '22px',
                            'minWidth': '22px',
                            'padding': '0',
                            'border': 'none',
                            'borderRadius': '5px',
                            'overflow': 'hidden',
                            'cursor': 'pointer'
                        }
                    ),
                    html.Span(
                        gene,
                        title=tooltip_text,
                        style={
                            "marginLeft": "5px",
                            "fontSize": "12px",
                            "whiteSpace": "nowrap",
                            "textDecoration": "underline",
                            "cursor": "pointer"
                        }
                    )
                ],
                style={
                    'display': 'flex',
                    'alignItems': 'center',
                    'padding': '1px 2px',
                    "marginBottom": "5px"
                }
            )
        )

    # annotations_html = html.Div([
    #     html.B("Genes: "),
    #     html.Div(
    #         genes_html,
    #         style={
    #             "display": "flex",
    #             "flexWrap": "wrap",
    #             "gap": "6px"
    #         }
    #     )
    # ])
    annotations_html = html.Div([
        html.B("Genes: "),
        html.Div(
            genes_html,
            style={
                'display': 'flex',
                'flexWrap': 'wrap',
                'gap': '2px 6px'
            },
            id='gene-color-picker-container'
        )
    ])

    return annotations_html


# Main callback to update graph when changing size, or selecting genomes, etc.
@app.callback(
    Output("graph", "elements"),
    Output("nb-noeuds", 'children'),
    #Output('shared_storage_nodes', 'data', allow_duplicate=True),
    Output('search-message', 'children'),
    Output('annotations-info', 'children'),
    Output('graph', 'stylesheet'),
    Output('graph', 'layout', allow_duplicate=True),
    Output('home-page-store', 'data', allow_duplicate=True),
    Output('graph', 'selectedNodeData'),
    Output('graph', 'selectedEdgeData'),
    Output('zoom_shared_storage_nodes', 'data', allow_duplicate=True),
    Output('start-input', 'value', allow_duplicate=True),
    Output('end-input', 'value', allow_duplicate=True),
    Output('features-dropdown', 'value', allow_duplicate=True),
    Output('feature-input', 'value', allow_duplicate=True),
    Output('displayed-region-container', 'children'),
    Output("phylogenetic-page-store", "data", allow_duplicate=True),
    Output('sequences-page-store', 'data', allow_duplicate=True),
    Output('home-page-legend-store', 'data'),
    State('genome_selector', 'value'),
    State('shared-mode', 'value'),
    State('specific-genome_selector', 'value'),
    State({'type': 'color-picker', 'index': ALL}, 'value'),
    State('show-labels', 'value'),
    Input('update-btn', 'n_clicks'),
    Input('btn-zoom', 'n_clicks'),
    Input('btn-zoom-out', 'n_clicks'),
    Input('btn-reset-zoom', 'n_clicks'),
    State('graph', 'selectedNodeData'),
    State('size_slider', 'value'),
    State('home-page-store', 'data'),
    Input('search-button', 'n_clicks'),
    Input('update_graph_command_storage', 'data'),
    State('start-input', 'value'),
    State('end-input', 'value'),
    State('features-dropdown', 'value'),
    State('feature-input', 'value'),
    State('genomes-dropdown', 'value'),
    State('chromosomes-dropdown', 'value'),
    State('shared_storage', 'data'),
    State('shared_storage_nodes', 'data'),
    State('min_shared_genomes-input', 'value'),
    State('tolerance-input', 'value'),
    State('shared-region-color-picker', 'value'),
    State('zoom_shared_storage_nodes', 'data'),
    State('show-exons', 'value'),
    State('exon-color-picker', 'value'),
    State('layout-dropdown', 'value'),
    State("phylogenetic-page-store", "data"),
    State('sequences-page-store', 'data'),
    State('colored-edge-size-slider', 'value'),
    State('node-size-scale-slider', 'value'),
    State('graph-compression', 'value'),
    State('min-flow', 'value'),
    State('nodes-names', 'value'),
    State('global_parameters', 'data'),
    State({"type": "gene-color", "gene": ALL}, "value"),
    State({"type": "gene-color", "gene": ALL}, "id"),
    prevent_initial_call=True
)
def update_graph(selected_genomes, shared_mode, specifics_genomes, color_genomes, show_labels, 
                 update_n_clicks, zoom_clicks, zoom_out_clicks, reset_zoom_bouton_clicks, 
                 selected_nodes_data, size_slider, home_data_storage, n_clicks, update_graph_command_storage, start, end,
                 feature_name, feature_value, genome, chromosome, data_storage, data_storage_nodes,
                 min_shared_genome, tolerance, shared_regions_link_color, zoom_shared_storage, 
                 show_exons, exons_color, layout_choice, phylo_data, sequences_data, colored_edges_size, nodes_size_scale,
                 graph_compression_value, min_flow_compression, nodes_names_value, global_parameters,
                 genes_color_values, genes_color_ids):
    if genome is not None and chromosome is not None:
        if "nodes_cache_id" not in data_storage_nodes:
            raise PreventUpdate
        nodes_cache_id = data_storage_nodes["nodes_cache_id"]
        cached = get_session_cache(nodes_cache_id)
        ctx = dash.callback_context
        return_metadata = {"return_code":"", "flow":None, "nodes_number":0, "removed_genomes":None}
        message = ""
        nodes = {}
        genes_color = {}
        start_value = None
        end_value = None
        new_request = False
        #alternatif genome in case of zoom on a region that doesn't contain the reference genome
        alt_genome = ""
        if home_data_storage is None:
            home_data_storage = {}
        home_data_storage["genome_zoom"] = None
        home_data_storage["zoom"] = False
        triggered_id = ctx.triggered_id
        if (triggered_id == "search-button" and n_clicks > 0
                and ((start is None or start == "") or (end is None or end == ""))
                and (feature_value is None or feature_value == "")):
            return (no_update, no_update, f"❌ You must set start / end or feature value.", no_update, no_update,
                    no_update, no_update, no_update, no_update, no_update,
                    no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update)

        max_nodes_to_visualize = MAX_NODES_TO_VISUALIZE
        if global_parameters and "max_nodes_to_visualize" in global_parameters:
            max_nodes_to_visualize = global_parameters["max_nodes_to_visualize"]
        max_nodes_from_db = MAX_NODES_FROM_DB
        if global_parameters and "max_nodes_from_db" in global_parameters:
            max_nodes_from_db = global_parameters["max_nodes_from_db"]
        node_shape_as_circle = False
        if global_parameters and "circle" in global_parameters and global_parameters["circle"] == True:
            node_shape_as_circle = True
        if home_data_storage is None:
            home_data_storage = {}
        size_slider_val = DEFAULT_SIZE_VALUE
        if size_slider is None :
            if home_data_storage is not None and 'slider_value' in home_data_storage:
                size_slider_val = home_data_storage["slider_value"]
        else:
            size_slider_val = size_slider
            home_data_storage["slider_value"] = size_slider

        #Checks if min node size has been decreased : if so it is required to get data from database
        if size_slider_val is not None and "current_size" in home_data_storage and home_data_storage["current_size"] > size_slider_val and cached["min_node_size"] > size_slider_val:
            logger.debug(f"Min node size has been set to {size_slider_val} and is lower than old value {home_data_storage['current_size']} - nodes will be updated from database.")
            new_request = True
            if "start" in home_data_storage :
                start_value = home_data_storage["start"]
            if "end" in home_data_storage :
                end_value = home_data_storage["end"]
            if feature_name in home_data_storage and feature_value in home_data_storage :
                feature_name = home_data_storage["feature_name"]
                feature_value = home_data_storage["feature_value"]
        else:
            if "current_size" in  home_data_storage :
                logger.debug(f"Min node size {size_slider_val} - old value {home_data_storage['current_size']}.")
            else:
                logger.debug(f"min node size : {size_slider_val}")

        if "search_return_metadata" in home_data_storage and home_data_storage["search_return_metadata"] is not None :
            return_metadata = home_data_storage["search_return_metadata"]
            home_data_storage["search_return_metadata"] = None
        # save the parameters into store
        if genome is not None:
            home_data_storage["selected_genome"] = genome
        if "selected_genome" in home_data_storage:
            genome = home_data_storage["selected_genome"]
        if chromosome is not None:
            home_data_storage["selected_chromosome"] = chromosome
        if shared_regions_link_color is not None:
            home_data_storage["shared_regions_link_color"] = shared_regions_link_color

        if start is not None:
            home_data_storage["start"] = start
            start_value = start
            home_data_storage["feature_name"] = ""
            home_data_storage["feature_value"] = ""

        if end is not None:
            home_data_storage["end"] = end
            end_value = end
            home_data_storage["feature_name"] = ""
            home_data_storage["feature_value"] = ""

        if feature_name is not None and feature_name != "" and feature_value is not None and feature_value != "":
            home_data_storage["feature_name"] = feature_name
            home_data_storage["feature_value"] = feature_value

        if min_shared_genome is None:
            min_shared_genome = 100
        if tolerance is None:
            tolerance = 0
        if color_genomes is not None:
            home_data_storage["color_genomes"] = color_genomes
        if specifics_genomes is not None:
            home_data_storage["specifics_genomes"] = specifics_genomes

        for value, id_dict in zip(genes_color_values, genes_color_ids):
            if not isinstance(id_dict, dict):
                continue
            gene = id_dict["gene"]
            if not gene:
                continue
            if value and value.lower() != "#000000":
                genes_color[gene.lower()] = value

        compression = 'graph_compression' in graph_compression_value
        min_flow_compression_value = 0
        if compression:
            try:
                min_flow_compression_value = float(min_flow_compression)/100
            except (TypeError, ValueError):
                min_flow_compression_value = 0
            if not (0 <= min_flow_compression_value <= 1):
                min_flow_compression_value = 0
        else:
            min_flow_compression_value = 0

        nodes_names = 'nodes_names' in nodes_names_value
        # zoom on selected nodes
        zoom_shared_storage_out = zoom_shared_storage or {}
        zoom_ranges = {}
        if triggered_id == "btn-zoom":
            #In case of zoom => get the selected nodes to prepare a new request
            selected_nodes_name = set()
            if selected_nodes_data is not None and len(selected_nodes_data) > 0:
                selected_nodes_name = set([node['name'] for node in selected_nodes_data])
            else:
                raise PreventUpdate
            #Check if it is the first zoom to store it

            if "nodes" in cached:
                nodes = cached["nodes"]
                if zoom_shared_storage_out or len(zoom_shared_storage_out) == 0:
                    #First zoom => store the old data to retrieve them when reset zoom
                    cached["zoom"] = cached["nodes"]
                zoom_shared_storage_out["start"] = home_data_storage["start"]
                zoom_shared_storage_out["end"] = home_data_storage["end"]

                position_field = genome + "_position"
                selected_positions =set()
                for k,node in nodes.items():
                    if node["name"] in selected_nodes_name and position_field in node :
                        selected_positions.add(node[position_field])
                    if alt_genome is None or alt_genome == "":
                        if node["name"] in selected_nodes_name:
                            alt_genome = node["genomes"][0]
                if len(selected_positions) == 0 and alt_genome != "":
                    #Try to switch to another reference genome
                    genome = alt_genome
                    position_field = genome + "_position"
                    for n in nodes:
                        node = nodes[n]
                        if node["name"] in selected_nodes_name and position_field in node :
                            selected_positions.add(node[position_field])
                if len(selected_positions) > 0:
                    start_value = min(selected_positions)
                    home_data_storage["start"] = start_value
                    end_value = max(selected_positions)
                    home_data_storage["end"] = end_value
                    logger.debug(f"Zoom - start : {start_value} - end : {end_value}")
                else:
                    logger.debug(f"No position found in the selected nodes for the reference genome {genome}")

        if triggered_id == "btn-zoom-out":
            if "start" in home_data_storage and home_data_storage["start"] is not None \
                and "end" in home_data_storage and home_data_storage["end"] is not None:
                start_value = max(0,int(home_data_storage["start"]) - 1000)
                end_value = int(home_data_storage["end"]) + 1000
                home_data_storage["start"] = start_value
                home_data_storage["end"] = end_value
                home_data_storage["min_node_size"] = size_slider_val
            else:
                raise PreventUpdate
        if triggered_id == "btn-reset-zoom":
            if len(zoom_shared_storage_out) > 0:
                logger.debug(f"reset zoom to {zoom_shared_storage_out['start']} - {zoom_shared_storage_out['end']}")
                start_value = zoom_shared_storage_out["start"]
                end_value = zoom_shared_storage_out["end"]
                home_data_storage["start"] = start_value
                home_data_storage["end"] = end_value
            else:
                logger.debug(f"No zoom, can't reset zoom")
                raise PreventUpdate



        logger.debug("update graph : " + str(ctx.triggered[0]['prop_id']))
        stylesheet = []
        if shared_mode and 'shared' in shared_mode:
            specifics_genomes_list = specifics_genomes
            color_genomes_list = []
        else:
            specifics_genomes_list = []
            color_genomes_list = color_genomes
        labels = True
        if show_labels and 'hide' in show_labels:
            labels = False
        exons=False
        if show_exons and 'exons' in show_exons:
            exons = True
        all_genomes = data_storage["genomes"]
        all_chromosomes = data_storage["chromosomes"]
        #Checks if it is required to request database
        if (triggered_id== "search-button" and n_clicks > 0)  \
            or triggered_id in ["btn-zoom", "btn-reset-zoom", "btn-zoom-out"] \
            or (triggered_id == "update_graph_command_storage" and update_graph_command_storage is not None) \
            or new_request:
            new_data = {}
            if triggered_id == "search-button":
                #New search => reset genes color
                genes_color = {}
            #Delete local phylo graph if exists
            if start_value and end_value:
                if (triggered_id == "btn-reset-zoom"
                    and "zoom" in cached
                    and len(cached["zoom"]) > 0
                    and home_data_storage.get("min_node_size",0) == size_slider_val):
                    logger.debug(
                        f"Retrieve {len(cached['zoom'])} nodes before zoom.")
                else:
                    logger.debug(f"Getting data from database from {start_value} to {end_value} on chr {chromosome} for genome {genome}")
            else:
                if feature_name and feature_value:
                    logger.debug(
                        f"Getting data from database for {feature_name} {feature_value} for genome {genome}")
            if phylo_data is not None and "newick_region" in phylo_data:
                phylo_data["newick_region"] = None
            if sequences_data is not None :
                sequences_data = {}
            if start_value is not None:
                use_anchor = True
                if triggered_id == "btn-zoom":
                    use_anchor = False
                    home_data_storage["zoom"] = True
                    home_data_storage["genome_zoom"] = alt_genome
                if (triggered_id == "btn-reset-zoom"
                    and "zoom" in cached
                    and len(cached["zoom"]) > 0
                    and home_data_storage.get("min_node_size",0) == size_slider_val):
                    new_data = cached["zoom"]
                    return_metadata["return_code"] = "OK"
                    return_metadata['nodes_number'] = len(new_data)
                else:
                    new_data, return_metadata = get_nodes_by_region(
                            genome, chromosome=chromosome, start=start_value, end=end_value, use_anchor=use_anchor, min_node_size=size_slider_val,
                            max_nodes_number=max_nodes_from_db, selected_genomes=selected_genomes)

                #data_storage_nodes = new_data
                logger.debug("len new_data : " + str(len(new_data)))
            else:
                if (feature_name is not None and feature_name != "" and feature_value is not None and feature_value != "") and chromosome is not None:
                        new_data,return_metadata = get_nodes_by_feature(
                            genome, chromosome=chromosome, feature= feature_name, value=feature_value, min_node_size=size_slider_val,
                            max_nodes_number=max_nodes_from_db, selected_genomes=selected_genomes, use_anchor=True)
                else:
                    new_data, return_metadata = get_nodes_by_region(
                        genome, chromosome=chromosome, start=0, end=end, min_node_size=size_slider_val,
                        max_nodes_number=max_nodes_from_db, selected_genomes=selected_genomes)



            if triggered_id in ["btn-reset-zoom", "search-button"] and nodes_cache_id in zoom_shared_storage_out:
                cached["zoom"] = {}
            cached["min_node_size"] = size_slider_val
            cached["nodes"] = new_data
            nodes_cache.set(nodes_cache_id, cached, expire=8 * 3600)

            # Get the start / end value when graph is updated
            genome_position = genome + "_position"
            nodes_with_position = [node for node in new_data.values() if genome_position in node]
            if len(nodes_with_position) > 0:
                min_node = min(nodes_with_position, key=lambda x: x[genome_position])
                max_node = max(nodes_with_position, key=lambda x: x[genome_position])
                max_node_size = max_node.get("size", None)
                start_value = min_node.get(genome_position, None)
                end_value = max_node.get(genome_position) + max_node_size
                home_data_storage["start"] = start_value
                home_data_storage["end"] = end_value
                logger.debug(f"start value : {start_value} - end value : {end_value}")
            #data_storage_nodes = new_data

            nodes = new_data

            elements, nodes_count, legend_nodes_size_dict = compute_graph_elements(new_data, genome, selected_genomes, size_slider_val, all_genomes,
                                              all_chromosomes, specifics_genomes_list,
                                              color_genomes_list, labels=labels, min_shared_genome=min_shared_genome,
                                              tolerance=tolerance, color_shared_regions=shared_regions_link_color,
                                              exons=exons, exons_color=exons_color, colored_edges_size=colored_edges_size,
                                              compression = compression, min_flow_compression_value = min_flow_compression_value,
                                              max_nodes_to_visualize=max_nodes_to_visualize, nodes_size_scale=nodes_size_scale,
                                              genes_color=genes_color)
            home_data_storage["current_size"] = size_slider_val
            home_data_storage["min_node_size"] = size_slider_val
            if triggered_id == "search-button":
                zoom_shared_storage_out = {}
                message = html.Div("❌ Error.", style=warning_style)
            if len(elements) == 0 and nodes_count == 0:
                start_value = None
                end_value = None
                home_data_storage["start"] = start_value
                home_data_storage["end"] = end_value

            if len(elements) == 0 and new_data and nodes_count > 0:
                return_metadata["return_code"] = "WIDE"

            if new_data is not None and return_metadata["return_code"] == "OK":
                message = html.Div(f"✅ Region has been successfully found, number of node {return_metadata['nodes_number']}.", style=success_style)
            elif new_data is not None and return_metadata["return_code"] == "ZOOM":
                message = html.Div(
                    f"⚠️ Region is too large, data has been filtered and only nodes with flow > {return_metadata['flow']} are displayed. Nodes number : {return_metadata['nodes_number']}.",
                    style=warning_style)
            elif new_data is not None and return_metadata["return_code"] == "FILTER" and 'removed_genomes' in return_metadata:
                message = html.Div(
                    f"⚠️ Region is too large, some individuals have been removed from search : {return_metadata['removed_genomes']}. These individuals won't be complete. Total nodes number : {return_metadata['nodes_number']}.",
                    style=warning_style)
            elif new_data is not None and return_metadata["return_code"] == "FILTER":
                message = html.Div(
                    f"⚠️ Region is too large, some individuals have been removed from search and the graph is partail. Total nodes number : {return_metadata['nodes_number']}.",
                    style=warning_style)
            elif new_data is not None and return_metadata["return_code"] == "PARTIAL":
                message = html.Div(
                    f"✅ Region has been successfully found (without core anchor), number of node {return_metadata['nodes_number']}.",
                    style=success_style)
            elif return_metadata["return_code"] == "OK":
                message = html.Div(
                    f"✅ Region has been successfully found, number of node {return_metadata['nodes_number']}.",
                    style=success_style)
            elif return_metadata["return_code"] == "WIDE":
                message = html.Div("⚠️ Region is too wide and cannot be displayed.", style=warning_style)
            elif return_metadata["return_code"] == "NO_DATA":
                if feature_name is not None and feature_name != "" and feature_value is not None and feature_value != "":
                    message = html.Div(f"❌ No nodes associated to {feature_name} {feature_value} found.", style=error_style)
                else:
                    message = html.Div("❌ Region not found.", style=error_style)
            else:
                message = html.Div("❌ Error.", style=error_style)

        else:
            #No needs to get data from database
            start_value = home_data_storage.get("start",None)
            end_value = home_data_storage.get("end",None)
            feature_name = home_data_storage.get("feature_name", "")
            feature_value = home_data_storage.get("feature_value", "")
            nodes = cached.get("nodes", {})
            elements, nodes_count, legend_nodes_size_dict = compute_graph_elements(nodes, genome, selected_genomes, size_slider_val, all_genomes,
                                              all_chromosomes, specifics_genomes_list,
                                              color_genomes_list, labels=labels, min_shared_genome=min_shared_genome,
                                              tolerance=tolerance, color_shared_regions=shared_regions_link_color,
                                              exons=exons, exons_color=exons_color, colored_edges_size=colored_edges_size,
                                              compression = compression, min_flow_compression_value = min_flow_compression_value,
                                              max_nodes_to_visualize=max_nodes_to_visualize, nodes_size_scale=nodes_size_scale,
                                              genes_color=genes_color)

            if len(elements) == 0 and nodes_count > 0:
                message = html.Div("⚠️ Region is too wide and cannot be displayed.", style=warning_style)

        defined_color = 0
        if color_genomes is not None:
            for c in color_genomes:
                if c != "#000000":
                    defined_color += 1
        stylesheet = compute_stylesheet(defined_color, nodes_names, node_shape_as_circle,nodes_size_scale=nodes_size_scale)
        count = len(elements)

        annotations = ""
        if nodes != None and len(nodes) > 0:
            annotations = build_annotations(nodes, genes_color)
        #default layout is fcose
        layout = {
            'name': 'fcose',
            'maxIterations': 100000,
            'maxSimulationTime': 5000,
            'quality': "proof",
            'fit': True
        }
        if layout_choice:
            match (layout_choice):
                case "dagre":
                    layout = {'name': 'dagre',
                              'rankDir': "RL",
                              'nodeDimensionsIncludeLabels': True,
                              'fit': True
                              }
                case "preset":
                    layout = {'name': 'preset',
                              'fit': True
                              }
        #displayed region construction:
        displayed_div = get_displayed_div(start_value, end_value, feature_name, feature_value)
        return (elements,f"{nodes_count} displayed nodes", message, annotations, stylesheet,
                layout, home_data_storage, [], [], zoom_shared_storage_out,
                None, None, feature_name, "", displayed_div, phylo_data, sequences_data, legend_nodes_size_dict)
    else:
        return ([], "", f"❌ No data loaded, first load a gfa into DB management page.", "", no_update,
                no_update, no_update, [], [], no_update,
                None, None, feature_name, "", "", no_update, no_update, no_update)


# color picker
@app.callback(
    Output('color-picker-container', 'style'),
    Output('shared-checklist-container', 'style'),
    Input('shared-mode', 'value')
)
def toggle_inputs(shared_mode):
    if shared_mode and 'shared' in shared_mode:
        return {'display': 'none'}, {'display': 'block'}
    else:
        return {'display': 'flex', 'flexWrap': 'wrap'}, {'display': 'none'}


@app.callback(
    Output('query-params-store', 'data'),
    Output('url', 'search'),
    Input('url', 'search'),
    prevent_initial_call=True
)
def process_query_params(search):
    if not search:
        raise PreventUpdate
    logger.debug(f"Query params {search}")
    params = parse_qs(urlparse(search).query)
    url_hap = params.get('haplotype', [None])[0]
    url_chromosome = (
        params.get('chr', [None])[0]
        or params.get('chromosome', [None])[0]
    )
    url_feature_name = params.get('featureName', [None])[0]
    url_feature_value = params.get('featureValue', [None])[0]
    url_start = params.get('start', [None])[0]
    url_end = params.get('end', [None])[0]

    valid_query = (
        url_hap is not None
        and url_chromosome is not None
        and (
            (url_start is not None and url_end is not None)
            or
            (url_feature_name is not None and url_feature_value is not None)
        )
    )
    if not valid_query:
        raise PreventUpdate
    query_data = {
        "url_hap": url_hap,
        "url_chromosome": url_chromosome,
        "url_feature_name": url_feature_name,
        "url_feature_value": url_feature_value,
        "url_start": url_start,
        "url_end": url_end,
        "force_refresh": time.time()
    }
    logger.debug(query_data)
    # clear URL after processing
    return query_data, ""


@app.callback(
    Output('size_slider', 'value'),
    Output('chromosomes-dropdown', 'value'),
    Output('genomes-dropdown', 'value'),
    Output('start-input', 'value'),
    Output('end-input', 'value'),
    Output('features-dropdown', 'value'),
    Output('feature-input', 'value'),
    Output('shared-region-color-picker', 'value'),
    Output('specific-genome_selector', 'value'),
    Output('update_graph_command_storage', 'data'),
    Input('url', 'pathname'),
    Input('home-page-store', 'data'),
    Input('query-params-store', 'data'),
    State('shared_storage', 'data'),
    State('genomes-dropdown', 'options'),
    State('chromosomes-dropdown', 'options'),
    State('specific-genome_selector', 'value'),
    State('demo_store', 'data'),
)
def update_parameters_on_page_load(pathname,data,
                                   query_params,shared_data,
                                   options_genomes,options_chromosomes,
                                   specifics_genomes, demo_data):

    if pathname != "/" and pathname != "/demo":
        raise PreventUpdate
    slider_value = DEFAULT_SIZE_VALUE

    selected_genome = None
    selected_chromosome = None

    start_input = None
    end_input = None

    features = shared_data.get("features", None)

    if features is not None and "gene" in features:
        feature_name = "gene"
    else:
        feature_name = None

    feature_value = ""
    shared_regions_link_color = DEFAULT_SHARED_REGION_COLOR
    update_graph_command_storage = dash.no_update
    if specifics_genomes is not None:
        selected_shared_genomes = specifics_genomes
    else:
        selected_shared_genomes = []

    if data is None:
        data = {}

    #Saved values
    if "slider_value" in data and data["slider_value"] is not None:
        slider_value = data["slider_value"]
    if "shared_regions_link_color" in data:
        shared_regions_link_color = data["shared_regions_link_color"]
    if "specifics_genomes" in data:
        selected_shared_genomes = data["specifics_genomes"]

    #Query params
    no_query_params = True

    if query_params is not None:
        no_query_params = False
        selected_genome = query_params.get("url_hap")
        selected_chromosome = query_params.get("url_chromosome")
        url_feature_name = query_params.get("url_feature_name")
        url_feature_value = query_params.get("url_feature_value")
        url_start = query_params.get("url_start")
        url_end = query_params.get("url_end")
        if url_feature_name is not None and url_feature_value:
            feature_name = url_feature_name
            feature_value = url_feature_value

        else:
            start_input = int(url_start)
            end_input = int(url_end)

        update_graph_command_storage = query_params
    elif (pathname == "/demo" and demo_data and len(demo_data) > 0 and "url_start" in demo_data
          and "url_end" in demo_data and "url_chromosome" in demo_data and "url_hap" in demo_data):
        no_query_params = False
        selected_genome = demo_data.get("url_hap", options_genomes[0]["value"])
        selected_chromosome = demo_data.get("url_chromosome", options_chromosomes[0]["value"])
        start_input = demo_data.get("url_start")
        end_input = demo_data.get("url_end")
        update_graph_command_storage = demo_data

    if no_query_params:
        if "selected_genome" in data:
            selected_genome = data["selected_genome"]
        else:
            if options_genomes:
                selected_genome = options_genomes[0]["value"]

        if "selected_chromosome" in data:
            selected_chromosome = data["selected_chromosome"]
        else:
            if options_chromosomes:
                selected_chromosome = options_chromosomes[0]["value"]
        start_input = None
        end_input = None
        if features is not None and "gene" in features:
            feature_name = "gene"
        else:
            feature_name = None
        feature_value = ""

    return (slider_value, selected_chromosome, selected_genome, start_input,
        end_input, feature_name, feature_value, shared_regions_link_color,
        selected_shared_genomes, update_graph_command_storage)



######## download graph callbacks ###############
@app.callback(
    Output('dummy-output', 'children'),
    Output("download-graph", "data", allow_duplicate=True),
    Input('graph', 'imageData'),
    State('chromosomes-dropdown', 'value'),
    State('genomes-dropdown', 'value'),
    State('home-page-store', 'data'),
    prevent_initial_call=True
)
def save_image_to_file(image_data, chromosome, genome, data):
    start = data.get("start", "")
    end =  data.get("end", "")
    if not image_data:
        raise PreventUpdate
    # check directory exists
    os.makedirs(EXPORT_DIR, exist_ok=True)


    #header, base64_data = image_data.split(',', 1)
    header, content = image_data.split(',', 1)
    # Get image format from header
    if 'image/png' in header:
        ext = 'png'
        image_bytes = base64.b64decode(content)
    elif 'image/jpeg' in header:
        ext = 'jpg'
        image_bytes = base64.b64decode(content)
    elif 'image/svg+xml' in header:
        ext = 'svg'
        if ";base64" in header:
            # SVG base64
            image_bytes = base64.b64decode(content)
        else:
            # SVG texte / URL-encoded
            svg_text = urllib.parse.unquote(content)
            image_bytes = svg_text.encode("utf-8")
    else:
        raise ValueError("Unsupported image format")


    file_name = "graph_"+str(genome)+"_chr_"+str(chromosome) + \
        "_start_"+str(start)+"_end_"+str(end)+"."+ext
    #image_bytes = base64.b64decode(base64_data)
    if not SERVER_MODE:
        save_path = os.path.join(os.getcwd(), EXPORT_DIR, file_name)
        with open(save_path, 'wb') as f:
            f.write(image_bytes)
    
        return f"Image downloaded in {save_path}", None
    else:
        return "File downloaded.", dcc.send_bytes(lambda io: io.write(image_bytes), file_name)


@app.callback(
    Output("graph", "generateImage"),
    Input("btn-save-jpg", "n_clicks"),
    Input("btn-save-png", "n_clicks"),
    Input("btn-save-svg", "n_clicks")
)
def trigger_image_save(n_clicks_jpg, n_clicks_png, n_clicks_svg):
    if not ctx.triggered_id:
        raise PreventUpdate

    # Get image format
    fmt = ctx.triggered_id.split('-')[-1]

    if fmt == 'svg':
        return {'type': fmt, 'action': 'download'}
    else:
        return {'type': fmt, 'action': 'store'}

