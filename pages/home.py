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

#MAX_GAP is used to dash edges between nodes separated by more than this value
MAX_GAP = 50000


def records_to_dataframe(nodes_data):
    rows = []
    for record in nodes_data:
        rows.append(nodes_data[record])
    return pd.DataFrame(rows)


def compute_stylesheet(color_number):
    if color_number > 1:
        stylesheet = [
            {
            'selector': 'node',
            'style': {
                'label': 'data(label)',
                'backgroundColor':'data(color)',
                'text-opacity':1,
                'opacity':1,
                'width':'data(displayed_node_size)',
                'height':'data(displayed_node_size)',
                'z-compound-depth': 'top'

                }
            },
            {
                'selector': '.main-node',
                'style': {'shape': 'circle'}
            },
            {
                'selector': '.degenerate-node',
                'style': {'shape': 'square'}
            },
            {
                'selector': 'edge',
                'style': {
                    'curve-style': 'unbundled-bezier',
                    'control-point-weights': [0.5],
                    'target-arrow-color': 'data(color)',
                    'target-arrow-shape': 'triangle',
                    'arrow-scale': 0.5,
                    'control-point-distances': [1],
                    'opacity':0.9,
                    'z-compound-depth': 'bottom'
                },
            },
            {'selector': ':selected', 'style': {
                'backgroundColor': 'grey',
                'line-color': '#FF4136',
                'border-width': 3,
                'border-color': 'black'
            }}  
            
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
                    'label': 'data(label)',
                    'text-opacity': 1,
                    'opacity':1,
                    'width':'data(displayed_node_size)',
                    'height':'data(displayed_node_size)',
                    'z-compound-depth': 'top'
                }
            },
            {
                'selector': '.main-node',
                'style': {'shape': 'circle'}
            },
            {
                'selector': '.degenerate-node',
                'style': {'shape': 'square'}
            },
            {
                'selector': 'edge',
                'style': {
                    'line-color': '#A3C4BC',
                    'target-arrow-color': '#A3C4BC',
                    'target-arrow-shape': 'triangle',
                    'arrow-scale': 0.5,
                    'curve-style': 'straight',
                    'opacity':0.9,
                    'z-compound-depth': 'bottom'


                }
            },
            {'selector': ':selected', 'style': {
                'background-color': 'grey',
                'line-color': '#FF4136',
                'border-width': 3,
                'border-color': 'black'
            }}
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
    for n in nodes_to_remove:
        if n not in predecessors:
            continue  # No predecessor, cannot merge

        pred_node = list(predecessors[n])[0]
        visited_nodes = set([n, pred_node])
        #Get the first predecessor not removed
        while pred_node in nodes_to_remove and pred_node in predecessors:
            pred_node = list(predecessors[pred_node])[0]
            if (pred_node in visited_nodes and len(list(predecessors[pred_node])) > 1):
                pred_node = list(predecessors[pred_node])[1]
            if (pred_node in visited_nodes) :
                pred_node = None
                break
            visited_nodes.add(pred_node)

            if pred_node is None:
                break

        if pred_node is None:
            continue

        idx = df.loc[df['name'] == pred_node].index[0]
        # Update size
        df.loc[df['name'] == pred_node, 'size'] += df.loc[df['name'] == n, 'size'].values[0]

        # add features
        f1 = df.loc[df['name'] == n, 'features'].iloc[0]
        f2 = df.loc[df['name'] == pred_node, 'features'].iloc[0]

        features = list(set((f1 if isinstance(f1, list) else []) + (f2 if isinstance(f2, list) else [])))
        df.at[idx, 'features'] = features

        # Concatenate annotations
        a1 = df.loc[df['name'] == n, 'annotations'].iloc[0]
        a2 = df.loc[df['name'] == pred_node, 'annotations'].iloc[0]
        annotations = list(set((a1 if isinstance(a1, list) else []) + (a2 if isinstance(a2, list) else [])))
        df.at[idx, 'annotations'] = annotations

    # Remove the compacted nodes
    df_compacted = df[~df['name'].isin(nodes_to_remove)].copy()

    return df_compacted


"""
This function compress a graph by removing
linear internal nodes, while allowing orientation inversions
across genomes.
"""
def graph_compression(df):

    genome_position_cols = [c for c in df.columns if c.endswith("_position")]

    # Predecessors and successors (directed graph)
    predecessors = {n: set() for n in df["name"]}
    successors = {n: set() for n in df["name"]}

    # Nodes with more than one successor / predecessor (or two in case of reverse nodes)
    invalid_nodes = set()

    # ------------------------------------------------------------------
    # Step 1: Build the directed graph genome by genome
    # ------------------------------------------------------------------
    for genome_position in genome_position_cols:
        # Keep only nodes existing in this genome
        sub = df[["name", genome_position]].dropna(subset=[genome_position])
        sub = sub.sort_values(by=genome_position)
        ordered_nodes = sub["name"].tolist()

        # Connect consecutive nodes
        for u, v in zip(ordered_nodes[:-1], ordered_nodes[1:]):
            if u not in invalid_nodes:
                successors[u].add(v)
                if len(successors[u]) > 2:
                    invalid_nodes.add(u)

            if v not in invalid_nodes:
                predecessors[v].add(u)
                if len(predecessors[v]) > 2:
                    invalid_nodes.add(v)

    # ------------------------------------------------------------------
    # Step 2: Remove invalid nodes from the graph
    # ------------------------------------------------------------------
    for n in invalid_nodes:
        predecessors.pop(n, None)
        successors.pop(n, None)

    # ------------------------------------------------------------------
    # Step 3: Identify removable nodes
    # ------------------------------------------------------------------
    nodes_to_remove = set()
    all_nodes = set(predecessors.keys()) | set(successors.keys())

    for node in all_nodes:
        preds = predecessors.get(node, set())
        succs = successors.get(node, set())

        #checks the flow of previous node, it must be equal to current node
        #this is due to different start coordinates on the different haplotypes
        curr_flow = df.loc[df['name'] == node, 'flow'].iloc[0]
        for pred in preds:
            pred_flow = df.loc[df['name'] == pred, 'flow'].iloc[0]
            if pred_flow != curr_flow:
                continue
            continue
        # Node must be locally linear
        if len(preds) > 2 or len(succs) > 2 or len(preds) == 0 or len(succs) == 0:
            continue

        #Node is compressed if there is only one predecessor and one successor
        #and if the predecessor has only one successor
        if len(preds) == 1 and len(succs) == 1:
            if len(successors.get(list(preds)[0],[])) == 1:
                nodes_to_remove.add(node)
        else:
            #Checks if node is traversed in both direction
            #Node is compressed if it has only two neighbors
            #and their neighbors have only two neighbors too
            if len(preds) == 2 and len(succs) == 2:
                neighbors = preds | succs
                if len(neighbors) == 2:
                    pred_neighbors = successors.get(list(neighbors)[0],set()) | predecessors.get(list(neighbors)[0],set())
                    succ_neighbors = successors.get(list(neighbors)[1], set()) | predecessors.get(list(neighbors)[1], set())
                    if len(pred_neighbors) == 2 and len(succ_neighbors) == 2:
                        nodes_to_remove.add(node)

    # ------------------------------------------------------------------
    # Step 4: Remove nodes from the dataframe
    # ------------------------------------------------------------------
    df_compacted = merge_node_data(df, nodes_to_remove, predecessors)

    return df_compacted

def compute_graph_elements(data, ref_genome, selected_genomes, size_min, all_genomes, all_chromosomes,
                           specifics_genomes=None, color_genomes=[], x_max=1000, y_max=1000, labels=True,
                           min_shared_genome=100, tolerance=0, color_shared_regions=DEFAULT_SHARED_REGION_COLOR,
                           exons=False, exons_color=DEFAULT_EXONS_COLOR, colored_edges_size=5,
                           compression=False):
    logger.debug(f"Compute elements with ref genome {ref_genome}")

    if data != None and len(data) > 0:
        if ref_genome is not None and ref_genome != "":
            position_field = ref_genome+"_position"
        else:
            position_field = "mean_pos"
        df = records_to_dataframe(data)
        df = df[df["size"] >= size_min].copy()
        #If compression is set to true then it will compress the graph
        #Nodes with a single predecessor and a single successor, and for which the predecessor
        #has only one outgoing node, are removed
        if compression:
            df = graph_compression(df)


        df = df[df["genomes"].apply(lambda g: any(
            x in selected_genomes for x in g))].copy()
        #logger.debug(f"Compute elements - dataframe")
        def mean_position(row):
            positions = [row.get(f"{g}_node")
                         for g in row["genomes"] if f"{g}_node" in row]
            positions = [p for p in positions if p is not None]
            return np.mean(positions) if positions else 0

        size_max = df['size'].max()
        size_min = df['size'].min()
        s_max = 200
        s_min = 30
        df["displayed_node_size"] = s_min + (df["size"]-size_min)*(s_max-s_min)/(size_max-size_min)


        df["mean_pos"] = df.apply(mean_position, axis=1)
        x_min, x_max_data = df["mean_pos"].min(), df["mean_pos"].max()
        df["x"] = ((df["mean_pos"] - x_min) /
                   (x_max_data - x_min + 1e-6)) * x_max



        X_SPAN = 40000
        GAP = 40

        min_x = df["mean_pos"].min()
        max_x = df["mean_pos"].max()

        df["x_mean"] = ((df["mean_pos"] - min_x)/ (max_x - min_x)) * X_SPAN - X_SPAN / 2

        df = df.sort_values("x_mean").reset_index(drop=True)
        df["x"] = 0.0

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

        df.loc[mask, "y"] = (
                Y_MIN
                + (Y_MAX - Y_MIN)
                * (df.loc[mask, "flow"] / FLOW_CENTER) ** ALPHA
        )


        #df.loc[df.index % 2 == 0, "y"] *= -1

        #This line is due to a pb with cytoscape : if no y with negative and positive values then
        #grpah disappear when zooming
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
                    'annotations': row.get('annotations'),
                    'features': row.get('features'),
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

            # Préfiltrage plus rapide
            mask = df["genomes"].apply(lambda g: genome in g)
            nodes_g = df.loc[mask, :]

            col = f"{genome}_position"
            if col not in nodes_g.columns:
                continue

            nodes_g = nodes_g.loc[nodes_g[col].notnull()]
            if nodes_g.empty:
                continue

            # Tri sans copy
            nodes_g = nodes_g.sort_values(by=col, ascending=True, ignore_index=True)

            # Calcul des sources et cibles via décalage
            sources = nodes_g["name"].iloc[:-1].values
            targets = nodes_g["name"].iloc[1:].values
            positions = nodes_g[col].values
            sizes = nodes_g["size"].values
            annotations = nodes_g["annotations"].values

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
            if labels:

                first_label = True

                for a in dic["annotations"]:
                    if first_label:
                        label += str(a)
                        first_label = False
                    else:
                        label += ", " + str(a)
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
                    'text-rotation': 'autorotate',
                    'width': (virtual_flow+int(0.2*len(all_genomes)))/len(all_genomes)*10
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
        return nodes + edges, len(nodes)
    else:
        return [], 0


legend = html.Div(
    style={
        "bottom": "20px",
        "left": "20px",
        "background": "rgba(255,255,255,0.95)",
        "border": "1px solid #ccc",
        "borderRadius": "6px",
        "padding": "10px",
        "width": "320px",
        "fontSize": "12px"
    },
    children=[

        html.Div("Legend", style={"fontWeight": "bold", "marginBottom": "8px"}),

        # === NODE SIZE ===
        svg.Svg(width="300", height="40", children=[
            svg.Circle(cx="35", cy="20", r="6", fill="blue"),
            svg.Circle(cx="35", cy="20", r="10", fill="blue", opacity="0.5"),
            svg.Circle(cx="35", cy="20", r="14", fill="blue", opacity="0.4"),

            svg.Text("Node size ∝ sequence length", x="70", y="25")
        ]),

        # === NODE COLOR (BLUE → RED) ===
        svg.Svg(width="300", height="30", children=[
            svg.Circle(cx="20", cy="15", r="8", fill="blue"),
            svg.Circle(cx="40", cy="15", r="8", fill="#7f007f"),
            svg.Circle(cx="60", cy="15", r="8", fill="red"),
            svg.Text(
                "Color: few individuals → all individuals",
                x="80", y="20"
            )
        ]),

        # === REPEATED NODE ===
        svg.Svg(width="300", height="30", children=[
            svg.Rect(x="25", y="5", width="20", height="20", fill="red"),
            svg.Text("Repeated node", x="80", y="20")
        ]),

        # === EDGE WIDTH ===
        svg.Svg(width="300", height="30", children=[
            svg.Line(
                x1="10", y1="15", x2="60", y2="15",
                stroke="#333",
                strokeWidth="5"
            ),
            svg.Text(
                "Edge width ∝ number of individuals",
                x="80", y="20"
            )
        ]),

        # === DASHED EDGE ===
        svg.Svg(width="300", height="30", children=[
            svg.Line(
                x1="10", y1="15", x2="60", y2="15",
                stroke="#333",
                strokeWidth="3",
                strokeDasharray="6,4"
            ),
            svg.Text(
                "Dashed edge: long genomic distance",
                x="80", y="20"
            )
        ]),
    ]
)


def layout(data=None, initial_size_limit=10):
    all_genomes = get_genomes()
    all_genomes.sort()
    all_chromosomes = get_chromosomes()
    features = get_annotations_features()
    if data != None:
        elements, nodes_count = compute_graph_elements(
            data, all_genomes, initial_size_limit, all_genomes, all_chromosomes, [], [])
    else:
        elements = []
    size_max = 500

    return html.Div([
        dcc.Store(id='zoom_shared_storage_nodes', storage_type='memory'),
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
                        "Graph compression: when a min node size is set greater than 0, this will compact the linear parts of the graph (i.e. nodes connected to exactly 2 nodes will be compacted."),
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
                            "Node shape : a node is drawn as a circle, unless if it's a repeated node in which case it will be displayed as a square."),
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
                                html.Label("Start : ", title="Start on the selected haplotype / chromosome."),
                                dcc.Input(id='start-input', type='number', style={'width': '100px', 'marginRight': '10px'}),
                                html.Label("End : ", title="End on the selected haplotype / chromosome."),
                                dcc.Input(id='end-input', type='number', style={'width': '100px', 'marginRight': '20px'})
                            ], style={'marginBottom': '10px'}),
                            html.Div(
                                style={
                                    'display': 'flex',
                                    'alignItems': 'center',
                                    'gap': '15px',
                                    'marginBottom': '20px'
                                },
                                children=[

                                    html.Div(
                                        style={'width': '300px'},
                                        children=[
                                            dcc.Dropdown(
                                                id='features-dropdown',
                                                options=[
                                                            {'label': f"{feature}_name", 'value': f"{feature}_name"}
                                                            for feature in features
                                                        ] + [
                                                            {'label': f"{feature}_id", 'value': f"{feature}_id"}
                                                            for feature in features
                                                        ],
                                                value='gene_name' if 'gene' in features else None,
                                                clearable=False,
                                                placeholder="Choose a feature to search"
                                            )
                                        ]
                                    ),

                                    dcc.Input(
                                        id='feature-input',
                                        type='text',
                                        placeholder='Value to search',
                                        debounce=True,
                                        style={
                                            'height': '38px',
                                            'lineHeight': '38px',
                                            'padding': '0 10px',
                                            'boxSizing': 'border-box'
                                        }
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
                        html.Label("Haplotypes to visualize :", title="Nodes containing only unselected haplotypes won't be displayed.", style={
                                   'marginBottom': '5px'}),
                        dcc.Checklist(
                            id="genome_selector",
                            options=[{"label": g, "value": g}
                                     for g in all_genomes],
                            value=all_genomes,
                            inline=True
                        )

                    ],
                        style={'marginBottom': '20px'}),
                    # Download graph div
                    html.Div([
                        html.Label(
                            'Download graph:', title="PNG and JPG files will be saved into ./export/graphs directory while SVG files are downloaded.", style={"marginLeft": "10px"}),
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
                    # dcc.Loading(id="loading-spinner", type="circle",
                    #             children=html.Div(id="output-zone"))
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
                        min=0,
                        max=size_max,
                        step=1,
                        marks={i: str(i) for i in range(0, size_max + 1, int(size_max/10))},
                        value=DEFAULT_SIZE_VALUE,
                        tooltip={"placement": "bottom",
                                 "always_visible": False},
                    ),

                    html.Div(id='size_stats', style={'marginTop': '10px'})

                ]),
                html.Div(id='nb-noeuds', style={'margin': '10px'}),
                html.Div([
                    # === First line: Min node size + Layout + Dropdown ===
                    html.Div([
                        html.Div(
                            id='size-output',
                            children='Min node size : 10',
                            style={'marginRight': '20px'}
                        ),
                        html.H4(
                            "Layout",
                            title="The layout computes the graph. If the graph is not readable, it is possible to modify the algorithm. Dagre layout will display more linear graphs and fcose will display more compact graphs.",
                            style={'marginRight': '10px'}
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
                            style={'width': '120px', 'display': 'inline-block', 'marginRight':'10px'}
                        ),
                        html.Div(
                            style={'display': 'flex', 'alignItems': 'center', 'gap': '5px'},  # petit gap ici
                            children=[
                                html.Label("Colored edges size"),
                                html.Div(
                                    style={'width': '280px'},
                                    children=[
                                        dcc.Slider(
                                            id='colored-edge-size-slider',
                                            min=1,
                                            max=20,
                                            step=1,
                                            value=5,
                                            marks={i: str(i) for i in range(0, 25, 5)}
                                        )
                                    ]
                                )
                            ]
                        ),
                    ], style={'display': 'flex', 'alignItems': 'center', 'marginBottom': '10px', 'gap': '8px'}),

                    # === Second line: checkboxes and color picker ===
                    html.Div([
                        dcc.Checklist(
                            options=[{
                                'label': 'Hide labels',
                                'value': 'hide',
                                'title': 'Uncheck if labels takes too much space on the graph.'
                            }],
                            id='show-labels',
                            style={'marginRight': '30px'}
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
                            style={'width': '25px', 'height': '25px', 'marginRight': '10px'}
                        ),

                        dcc.Checklist(
                            options=[{
                                'label': 'Graph compression',
                                'value': 'graph_compression',  # <-- obligatoire pour Checklist
                                'title': 'Check to compact the graph nodes.'
                            }],
                            id='graph-compression',
                            style={'marginRight': '10px'},
                            value=[]
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
                            options=[{"label": g, "value": g}
                                     for g in all_genomes],
                            # value=[],
                            inline=True
                        ),
                        html.Label("Min (%) of shared haplotypes : ", title="Min (%) of shared haplotypes = M. Number of selected haplotypes = N. To detect a shared node it must contains almost (M/100) x N of the selected haplotypes. If M = 0 then the minimum number of selected haplotypes will be 1."),
                        dcc.Input(id='min_shared_genomes-input', type='number',
                                  value=0, style={'width': '100px', 'marginRight': '10px'}),
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
                            html.Label(s),
                            dbc.Input(
                                id={'type': 'color-picker', 'index': i},
                                type='color',
                                value="#000000",
                                style={'width': '25px',
                                       'height': '25px', 'marginLeft': '10px'}
                            )

                        ], style={'marginRight': '10px'}) for i, s in enumerate(all_genomes)
                    ], style={'display': 'flex',
                              'flexWrap': 'wrap',
                              'gap': '5px'}, id='color-picker-container')

                ], id='sample-controls'),
                html.Button("Update graph", id="update-btn",
                            n_clicks=0, style={'marginTop': '10px'}),

                html.Div(html.H4(id='node-info', style={'margin': '10px'})),
                html.Div(html.Label("Annotations :", title="Compiles all annotations for the displayed nodes.", style={
                         'marginBottom': '5px'})),
                html.Div(html.H4(id='annotations-info',
                         style={'margin': '10px'}))
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
                    title='Before using this button, nodes must be selected by holding left mouse button and Shift key.'
                ),
                html.Button(
                    "🔄 Reset Zoom",
                    id='btn-reset-zoom',
                ),
                html.Button(
                    "🔄 Zoom out 2000 bp",
                    id='btn-zoom-out',
                )
            ], style={'marginBottom': '10px', 'display': 'flex', 'gap': '10px'}),

            html.Div(
                legend,
                style={'marginBottom': '10px'}
            )
        ]),

        # Graph block
        cyto.Cytoscape(
            id='graph',
            style={'width': '100%', 'height': '1000px'},
            zoomingEnabled=True,
            userZoomingEnabled=True,
            userPanningEnabled=True,
            wheelSensitivity=0.1,
            boxSelectionEnabled=True,
            
        )
        ])


# Callback to get nodes or link info when clicking on it
@app.callback(
    Output('node-info', 'children'),
    Input('graph', 'tapNodeData'),
    Input('graph', 'tapEdgeData')
)
def display_element_data(node_data, edge_data):
    triggered_id = ctx.triggered_id

    if triggered_id == 'graph' and ctx.triggered[0]['prop_id'] == 'graph.tapEdgeData' and edge_data:
        return (
            f"Selected link : {edge_data.get('source')} → {edge_data.get('target')}\n"
            f"• Flow : {edge_data.get('flow')}\n"
            f"• Haplotypes : {', '.join(edge_data.get('genomes', []))}"
        )
    elif triggered_id == 'graph' and ctx.triggered[0]['prop_id'] == 'graph.tapNodeData' and node_data:
        return (
            f"Selected node : {node_data.get('label', node_data.get('name'))}\n"
            f"• Size : {node_data.get('size')}\n"
            f"• Position : {node_data.get('position')}\n"
            f"• Flow : {node_data.get('flow')}\n"
            f"• Ref node : {node_data.get('ref_node')}\n"
            f"• Haplotypes : {', '.join(node_data.get('genomes', []))}"
            f"• Annotations : {', '.join(node_data.get('annotations', []))}"
            f"• Features : {', '.join(node_data.get('features', []))}"
            f"• Sequence (first 1000 bp only) : {node_data.get('sequence')}\n"
        )
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


# Main callback to update graph when changing size, or selecting genomes, etc.
@app.callback(
    Output("graph", "elements"),
    Output("nb-noeuds", 'children'),
    Output('shared_storage_nodes', 'data', allow_duplicate=True),
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
    State('graph-compression', 'value'),
    prevent_initial_call=True
)
def update_graph(selected_genomes, shared_mode, specifics_genomes, color_genomes, show_labels, 
                 update_n_clicks, zoom_clicks, zoom_out_clicks, reset_zoom_bouton_clicks, 
                 selected_nodes_data, size_slider, home_data_storage, n_clicks, update_graph_command_storage, start, end,
                 feature_name, feature_value, genome, chromosome, data_storage, data_storage_nodes,
                 min_shared_genome, tolerance, shared_regions_link_color, zoom_shared_storage, 
                 show_exons, exons_color, layout_choice, phylo_data, sequences_data, colored_edges_size,
                 graph_compression_value):
    if genome is not None and chromosome is not None:
        ctx = dash.callback_context
        return_metadata = {"return_code":"", "flow":None, "nodes_number":0, "removed_genomes":None}
        message = ""
        start_value = None
        end_value = None
        new_request = False
        triggered_id = ctx.triggered_id
        logger.debug(f"{triggered_id} update")
        if home_data_storage is None:
            home_data_storage = {}
        if size_slider is None :
            if home_data_storage is not None and 'slider_value' in home_data_storage:
                size_slider_val = home_data_storage['slider_value']
            else:
                size_slider_val = DEFAULT_SIZE_VALUE
        else:
            size_slider_val = size_slider
        #Checks if min node size has been decreased : if so it is required to get data from database
        if size_slider_val is not None and "current_size" in home_data_storage and home_data_storage["current_size"] > size_slider_val:
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
        compression = 'graph_compression' in graph_compression_value
        # zoom on selected nodes
        zoom_shared_storage_out = zoom_shared_storage or {}
        if triggered_id == "btn-zoom":
            if selected_nodes_data is not None and len(selected_nodes_data) > 0:
                selected_nodes_name = set([node['name'] for node in selected_nodes_data])
            if len(zoom_shared_storage_out) ==0:
                zoom_shared_storage_out["start"] = home_data_storage["start"]
                zoom_shared_storage_out["end"] = home_data_storage["end"]
            position_field = genome + "_position"
            selected_positions =set()
            for n in data_storage_nodes:
                node = data_storage_nodes[n]
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
                start_value = max(0,home_data_storage["start"] - 1000)
                end_value = home_data_storage["end"] + 1000
                home_data_storage["start"] = start_value
                home_data_storage["end"] = end_value
        if triggered_id == "btn-reset-zoom":
            if len(zoom_shared_storage_out) > 0:
                logger.debug(f"reset zoom to {zoom_shared_storage_out['start']} - {zoom_shared_storage_out['end']}")
                start_value = zoom_shared_storage_out["start"]
                end_value = zoom_shared_storage_out["end"]
                zoom_shared_storage_out = {}
                home_data_storage["start"] = start_value
                home_data_storage["end"] = end_value
            else:
                if ("feature_name" in home_data_storage and home_data_storage["feature_name"] is not None and home_data_storage["feature_name"] != ""
                    and "feature_value" in home_data_storage and home_data_storage["feature_value"] is not None and home_data_storage["feature_value"] != ""):
                    feature_name = home_data_storage["feature_name"]
                    feature_value = home_data_storage["feature_value"]
                    logger.debug(f"No zoom, display {feature_name} {feature_value}")
                elif "start" in home_data_storage and "end" in home_data_storage:
                    start_value = home_data_storage["start"]
                    end_value = home_data_storage["end"]
                    logger.debug(f"No zoom, display region {start_value} - {end_value}")


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
            #Delete local phylo graph if exists
            logger.debug("Getting data from database")
            if phylo_data is not None and "newick_region" in phylo_data:
                phylo_data["newick_region"] = None
            if sequences_data is not None :
                sequences_data = {}
            if start_value is not None:
                use_anchor = True
                if triggered_id == "btn-zoom":
                    use_anchor = False
                new_data, return_metadata = get_nodes_by_region(
                        genome, chromosome=chromosome, start=start_value, end=end_value, use_anchor=use_anchor, min_node_size=size_slider_val)
                #data_storage_nodes = new_data
                logger.debug("len new_data : " + str(len(new_data)))
            else:
                if (feature_name is not None and feature_name != "" and feature_value is not None and feature_value != "") and chromosome is not None:
                        new_data,return_metadata = get_nodes_by_feature(
                            genome, chromosome=chromosome, feature= feature_name, value=feature_value, min_node_size=size_slider_val)
                else:
                    new_data, return_metadata = get_nodes_by_region(
                        genome, chromosome=chromosome, start=0, end=end, min_node_size=size_slider_val)

            # Get the start / end value when graph is updated
            genome_position = genome + "_position"
            nodes_with_position = [node for node in new_data.values() if genome_position in node]
            if len(nodes_with_position) > 1:
                print(f"new data length {len(new_data)}")
                min_node = min(nodes_with_position, key=lambda x: x[genome_position])
                max_node = max(nodes_with_position, key=lambda x: x[genome_position])

                max_node_size = max_node.get("size", None)
                start_value = min_node.get(genome_position, None)
                end_value = max_node.get(genome_position) + max_node_size
                home_data_storage["start"] = start_value
                home_data_storage["end"] = end_value
                logger.debug(f"start value : {start_value} - end value : {end_value}")
            data_storage_nodes = new_data
            elements, nodes_count = compute_graph_elements(new_data, genome, selected_genomes, size_slider_val, all_genomes,
                                              all_chromosomes, specifics_genomes_list,
                                              color_genomes_list, labels=labels, min_shared_genome=min_shared_genome,
                                              tolerance=tolerance, color_shared_regions=shared_regions_link_color,
                                              exons=exons, exons_color=exons_color, colored_edges_size=colored_edges_size,
                                              compression = compression)
            home_data_storage["current_size"] = size_slider_val
            if triggered_id == "search-button":
                zoom_shared_storage_out = {}
                message = html.Div("❌ Error.", style=warning_style)
            if len(elements) == 0:
                start_value = None
                end_value = None
                home_data_storage["start"] = start_value
                home_data_storage["end"] = end_value


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
                    f"⚠️ No core genome anchor found near searched region, the result is partial and some nodes may be absent for some individuals. Total nodes number : {return_metadata['nodes_number']}.",
                    style=warning_style)
            elif return_metadata["return_code"] == "OK":
                message = html.Div(
                    "✅ Region has been successfully found.",
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
            logger.debug(f"min node size : {size_slider_val}")
            elements, nodes_count = compute_graph_elements(data_storage_nodes, genome, selected_genomes, size_slider_val, all_genomes,
                                              all_chromosomes, specifics_genomes_list,
                                              color_genomes_list, labels=labels, min_shared_genome=min_shared_genome,
                                              tolerance=tolerance, color_shared_regions=shared_regions_link_color,
                                              exons=exons, exons_color=exons_color, colored_edges_size=colored_edges_size,
                                              compression = compression)

        defined_color = 0
        if color_genomes is not None:
            for c in color_genomes:
                if c != "#000000":
                    defined_color += 1
        stylesheet = compute_stylesheet(defined_color)
        count = len(elements)
        annotations = ""
        set_annot = set()
        if data_storage_nodes != None:
            for n in data_storage_nodes:
                if "annotations" in data_storage_nodes[n]:
                    for a in data_storage_nodes[n]["annotations"]:
                        set_annot.add(a)
        for a in set_annot:
            annotations += str(a) + "\n"

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
        return (elements,f"{nodes_count} displayed nodes", data_storage_nodes, message, annotations, stylesheet,
                layout, home_data_storage, [], [], zoom_shared_storage_out,
                None, None, feature_name, "", displayed_div, phylo_data, sequences_data)
    else:
        return ([], "", no_update, f"❌ No data loaded, first load a gfa into DB management page.", "", no_update,
                no_update, no_update, [], [], no_update,
                None, None, feature_name, "", "", no_update, no_update)


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
    Output('home-page-store', 'data', allow_duplicate=True),
    Output("size-output", 'children'),
    Input('size_slider', 'value'),
    State('home-page-store', 'data'),
    prevent_initial_call=True
)
def save_slider_value(size_slider_val, data):
    if data is None:
        data = {}
    if size_slider_val is None:
        data['slider_value'] = DEFAULT_SIZE_VALUE
    else:
        data['slider_value'] = size_slider_val
    return data, f"Min node size  : {data['slider_value']}"

# Restore value after navigation
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
    Input('url', 'search'),
    Input('home-page-store', 'data'),
    State('shared_storage', 'data'),
    State('genomes-dropdown', 'options'),
    State('chromosomes-dropdown', 'options'),
    State('specific-genome_selector', 'value')
)
def update_parameters_on_page_load(pathname, search, data, shared_data, options_genomes, options_chromosomes, specifics_genomes):
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
    if "slider_value" in data and data["slider_value"] is not None:
        slider_value = data["slider_value"]
    if "shared_regions_link_color" in data:
        shared_regions_link_color = data["shared_regions_link_color"]
    if "specifics_genomes" in data:
        selected_shared_genomes = data["specifics_genomes"]
    #Get query param if setted
    #Query params exemple : ?haplotype=Korso_0&chromosome=1&featureName=gene_name&featureValue=BolK_1g00590
    #?haplotype=Korso_0&chromosome=1&start=100000&end=1090000
    no_query_params = True
    if search:
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
        if (url_hap is not None and url_chromosome is not None
                and ((url_start is not None and url_end is not None)
                     or (url_feature_name is not None and url_feature_value is not None))):
            no_query_params = False
            selected_genome = url_hap
            selected_chromosome = url_chromosome
            if url_feature_name is not None and url_feature_value:
                feature_name = url_feature_name
                feature_value = url_feature_value
            else:
                start_input = int(url_start)
                end_input = int(url_end)
            update_graph_command_storage = {"force_refresh": time.time(), "url_hap":url_hap, "url_chromosome":url_chromosome, "url_feature_name":url_feature_name,
                                            "url_feature_value":url_feature_value, "url_start":url_start, "url_end":url_end}
            logger.debug(update_graph_command_storage)
    if no_query_params :
        logger.debug("No query params")
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
    return slider_value, selected_chromosome, selected_genome, start_input, end_input, feature_name, feature_value, shared_regions_link_color, selected_shared_genomes, update_graph_command_storage


# Algorithm cytoscape choice
# @app.callback(
#     Output('graph', 'layout'),
#     Input('layout-dropdown', 'value')
# )
# def toggle_layout(layout_choice):
#     if layout_choice:
#         match(layout_choice):
#             case "dagre":
#                 return {
#                     'name': 'dagre',
#                     'rankDir': "RL",
#                     'nodeDimensionsIncludeLabels': True
#                 }
#             case "preset":
#                 return {
#                     'name': 'dagre',
#                     'rankDir': "RL",
#                     'nodeDimensionsIncludeLabels': True
#                 }
#
#     #Default layout is fcose
#     return {
#         'name': 'fcose',
#         'maxIterations': 100000,
#         'maxSimulationTime': 5000,
#         # 'nodeRepulsion': 10000,
#         # 'gravity': 0.1,
#         # 'gravityRangeCompound': 1.5,
#         # 'idealEdgeLength': 100,
#         # 'componentSpacing': 100,
#         # 'nodeDimensionsIncludeLabels': True,
#         # 'edgeElasticity': 0.1,
#         # 'nestingFactor': 0.8,
#         # 'tile': True,
#         'quality': "proof",
#         'fit': True
#     }


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
    print(f"start : {start}")
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

