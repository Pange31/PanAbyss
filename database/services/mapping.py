"""
Created on Tue Aug 19 16:05:16 2026

@author: fgraziani

This file contains the function to map a sequence on the pangenome graph
"""

from database.services.neo4j_kmer_requests import *
from database.services.neo4j_requests import get_nodes_from_sequence_id
import random
import time
import matplotlib.pyplot as plt
import plotly.graph_objects as go


"""
This function takes the nodes detected with their length, kmer counts and coverage and:
    1)It selects the top seeds (min min_seeds seeds and maximum max_seeds) by coverage
      If the min coverage associated to these seeds is 100% then it can return more than max_seeds 
      (all the 100% coverage seeds will be returned)
    2) It selects maximum max_selected_nodes nodes with the best coverage
Returns :
    seeds, selected_nodes: 
        2 dictionaries with node id as key and count, size and coverage as value
    
"""
def select_seeds(
        nodes_data,
        min_seeds: int = 10,
        max_seeds: int = 100,
        max_selected_nodes: int = 10000,
):
    # --------------------------------------------------------
    # Empty input.
    # --------------------------------------------------------

    n_nodes = nodes_data.size

    if n_nodes == 0:
        return {}, {}

    coverage = nodes_data["coverage"]

    # ========================================================
    # Select seeds.
    # ========================================================

    # We first try progressively larger top percentages until
    # we have at least 10 candidates.
    #
    # The percentage refers to the number of nodes, ranked by
    # decreasing coverage.
    # ========================================================

    min_seeds = min(min_seeds, max_seeds)

    # Sort by decreasing coverage.
    order = np.argsort(coverage)[::-1]

    # --------------------------------------------------------
    # Find the initial percentile.
    #
    # Start with the top 10%, then progressively increase
    # until we have at least min_seeds nodes.
    # --------------------------------------------------------

    if n_nodes <= min_seeds:
        seed_indices = np.arange(n_nodes)

    else:
        seed_indices = None

        for percentage in range(10, 101, 10):
            n_top = int(
                np.ceil(n_nodes * percentage / 100)
            )

            n_top = max(
                min_seeds,
                min(n_top, n_nodes),
            )

            threshold = coverage[
                order[n_top - 1]
            ]

            candidate_indices = np.flatnonzero(
                coverage >= threshold
            )

            if candidate_indices.size >= min_seeds:
                seed_indices = candidate_indices
                break

        # This is only a safeguard; 100% should always
        # produce at least min_seeds if enough nodes exist.
        if seed_indices is None:
            seed_indices = np.arange(n_nodes)

    # --------------------------------------------------------
    # Restrict to max_seeds.
    # --------------------------------------------------------

    if seed_indices.size > max_seeds:

        max_coverage = coverage[order[0]]

        # Special case:
        # if more than max_seeds nodes have 100% coverage,
        # keep all of them.
        max_coverage_indices = np.flatnonzero(
            coverage == max_coverage
        )

        if (
            max_coverage >= 100.0
            and max_coverage_indices.size > max_seeds
        ):
            seed_indices = max_coverage_indices

        else:
            # max_seeds-th best coverage.
            threshold = coverage[
                order[max_seeds - 1]
            ]

            # Keep all nodes tied at the threshold.
            seed_indices = np.flatnonzero(
                coverage >= threshold
            )

    selected_seed_nodes = nodes_data[seed_indices]

    seeds = {}

    for node in selected_seed_nodes:
        node_id = int(node["node_id"])

        seeds[node_id] = {
            "count": int(node["count"]),
            "size": int(node["size"]),
            "coverage": float(node["coverage"]),
        }

    logger.info(
        f"Number of seeds: {len(seeds)}"
    )

    # ========================================================
    # Select filtered nodes.
    # ========================================================

    if n_nodes <= max_selected_nodes:
        filtered_indices = np.arange(n_nodes)

    else:
        # Keep strictly max_selected_nodes best nodes.
        best_indices = np.argpartition(
            coverage,
            -max_selected_nodes,
        )[-max_selected_nodes:]

        filtered_indices = best_indices

    filtered_nodes_data = nodes_data[filtered_indices]

    selected_nodes = {}

    for node in filtered_nodes_data:
        node_id = int(node["node_id"])

        selected_nodes[node_id] = {
            "count": int(node["count"]),
            "size": int(node["size"]),
            "coverage": float(node["coverage"]),
        }

    logger.info(
        f"Number of selected nodes: "
        f"{len(selected_nodes)}"
    )

    return seeds, selected_nodes




"""
Return nodes supported by at least
map_kmer_threshold distinct k-mers.
"""
def get_seeds(
        seq: str,
        k: int,
        canonical=True,
        max_seeds: int = 100,
        min_coverage = 20,
        min_node_size = None,
        max_selected_nodes=10000,
):
    if not min_node_size:
        min_node_size = k+4
    start_time = time.time()

    if k <= 0:
        raise ValueError(
            f"k must be > 0, got {k}"
        )

    if not 1 <= max_seeds :
        raise ValueError(f"max_seeds must be greater than 1 ")

    if len(seq) < k:
        return {}

    # --------------------------------------------------------
    # Encode the complete sequence.
    # --------------------------------------------------------

    kmer_codes = encode_kmers(
        seq,
        k=k,
        canonical=canonical,
    )

    logger.info(
        f"kmers: {len(kmer_codes)}. "
        f"Time to encode: "
        f"{time.time() - start_time:.2f}s."
    )

    if kmer_codes.size == 0:
        return {}

    # --------------------------------------------------------
    # Remove duplicate encoded k-mers.
    # --------------------------------------------------------

    kmer_codes = np.unique(
        kmer_codes
    )

    logger.info(
        f"kmers unique: {len(kmer_codes)}"
    )

    # --------------------------------------------------------
    # Look up k-mers and count supporting nodes.
    # --------------------------------------------------------

    load_kmer_index(
        k=k,
        canonical=canonical,
    )


    t0 = time.time()

    node_kmer_counts = (
        get_encoded_kmer_node_counts_batch(
            kmer_codes,
            k=k,
            canonical=canonical,
        )
    )

    logger.info(
        f"Time to lookup and count nodes: "
        f"{time.time() - t0:.2f}s."
    )

    if not node_kmer_counts:
        return {}

    # --------------------------------------------------------
    # Batch lookup of node sizes.
    # --------------------------------------------------------

    t0 = time.time()

    node_ids = np.fromiter(
        node_kmer_counts.keys(),
        dtype=np.uint64,
    )

    counts = np.fromiter(
        node_kmer_counts.values(),
        dtype=np.uint64,
    )

    node_sizes = get_node_sizes_batch(
        k=k,
        node_ids=node_ids
    )

    logger.info(
        f"Time to lookup node sizes: "
        f"{time.time() - t0:.2f}s."
    )

    # --------------------------------------------------------
    # Calculate k-mer coverage for each node.
    # --------------------------------------------------------

    valid = (
        (node_sizes > 0)
        & (node_sizes >= min_node_size)
    )

    valid_node_ids = node_ids[valid]
    valid_counts = counts[valid]
    valid_sizes = node_sizes[valid]

    max_kmer_counts = (
        valid_sizes - k + 1
    )

    coverage = (
        100.0
        * valid_counts
        / max_kmer_counts
    )


    mask = coverage > min_coverage

    nodes_data_dtype = np.dtype([
        ("node_id", "<u8"),
        ("count", "<u8"),
        ("size", "<u8"),
        ("coverage", "<f8"),
    ])

    nodes = np.empty(
        np.count_nonzero(mask),
        dtype=nodes_data_dtype,
    )

    nodes["node_id"] = valid_node_ids[mask]
    nodes["count"] = valid_counts[mask]
    nodes["size"] = valid_sizes[mask]
    nodes["coverage"] = coverage[mask]

    logger.info(
        f"Nodes passing min_node_size: "
        f"{len(valid_node_ids)}"
    )

    logger.info(
        f"Nodes passing min_coverage: "
        f"{len(nodes)}"
    )



    # --------------------------------------------------------
    # Select nodes
    # --------------------------------------------------------

    seeds, filtered_nodes = select_seeds(nodes,max_seeds=max_seeds, max_selected_nodes=max_selected_nodes)


    return seeds, filtered_nodes, nodes


"""
Find dense genomic regions using a fixed sliding window.

A region is retained only if it contains at least one seed.

Regions are sorted by:
    1. number of seeds (descending)
    2. number of nodes (descending)

Parameters
----------
nodes : dict
    Nodes indexed by sequence_id.

seeds : dict
    Selected seed nodes indexed by sequence_id.

window_size : int
    Size of the sliding window in bp.

min_nodes : int
    Minimum number of nodes required in a window.

Returns
-------
regions : list[dict]
    Dense regions containing at least one seed.
"""
def find_dense_regions(
        nodes,
        seeds,
        sequence_size=90000,
        min_nodes=5,
        sequence_type="DNA",
        dna_window_factor=2.0,
        rna_window_factor=10.0
):

    sequence_type = sequence_type.upper()

    if sequence_type == "DNA":
        window_size = int(
            sequence_size * dna_window_factor
        )

    elif sequence_type == "RNA":
        window_size = int(
            sequence_size * rna_window_factor
        )

    else:
        raise ValueError(
            f"Unknown sequence_type: {sequence_type}. "
            f"Expected 'DNA' or 'RNA'."
        )
    # --------------------------------------------------------
    # Retrieve Neo4j nodes.
    # --------------------------------------------------------

    neo4j_nodes = get_nodes_from_sequence_id(
        list(nodes.keys())
    )

    if not neo4j_nodes:
        return []

    # --------------------------------------------------------
    # Build set of seed sequence IDs.
    # --------------------------------------------------------

    seed_ids = set(seeds.keys())

    # --------------------------------------------------------
    # Group nodes by chromosome.
    # --------------------------------------------------------

    chromosomes = {}

    for node in neo4j_nodes:

        chromosome = node["node"]["chromosome"]

        chromosomes.setdefault(
            chromosome,
            [],
        ).append(node)

    regions = []

    # ========================================================
    # Process each chromosome independently.
    # ========================================================

    for chromosome, chromosome_nodes in chromosomes.items():

        # ----------------------------------------------------
        # Sort nodes by genomic position.
        # ----------------------------------------------------

        chromosome_nodes.sort(
            key=lambda x: x["node"]["position_mean"]
        )

        positions = np.asarray(
            [
                node["node"]["position_mean"]
                for node in chromosome_nodes
            ],
            dtype=np.int64,
        )

        chromosome_min = positions[0]
        chromosome_max = positions[-1]

        # ----------------------------------------------------
        # Sliding window.
        # ----------------------------------------------------

        step = window_size // 2

        window_start = chromosome_min

        while window_start <= chromosome_max:

            window_end = (
                window_start + window_size
            )

            # ------------------------------------------------
            # Find nodes inside the window.
            # ------------------------------------------------

            left = np.searchsorted(
                positions,
                window_start,
                side="left",
            )

            right = np.searchsorted(
                positions,
                window_end,
                side="right",
            )

            n_window_nodes = right - left

            # ------------------------------------------------
            # Apply minimum node count.
            # ------------------------------------------------

            if n_window_nodes >= min_nodes:

                window_nodes = chromosome_nodes[
                    left:right
                ]

                # --------------------------------------------
                # Identify seeds in this window.
                # --------------------------------------------

                window_seed_nodes = [
                    node
                    for node in window_nodes
                    if node["sequence_id"] in seed_ids
                ]

                n_seeds = len(
                    window_seed_nodes
                )

                # --------------------------------------------
                # Keep only regions containing seeds.
                # --------------------------------------------

                if n_seeds > 0:

                    total_size = sum(
                        node["node"]["size"]
                        for node in window_nodes
                    )

                    genome_positions = {}

                    for node in window_nodes:

                        node_data = node["node"]

                        for genome in node_data.get("genomes", []):

                            position_key = f"{genome}_position"

                            position = node_data.get(position_key)

                            if position is None:
                                continue

                            if genome not in genome_positions:
                                genome_positions[genome] = {
                                    "min": position,
                                    "max": position,
                                }

                            else:
                                genome_positions[genome]["min"] = min(
                                    genome_positions[genome]["min"],
                                    position,
                                )

                                genome_positions[genome]["max"] = max(
                                    genome_positions[genome]["max"],
                                    position,
                                )

                    regions.append({
                        "chromosome": chromosome,
                        "position_min": int(window_start),
                        "position_max": int(window_end),

                        "n_seeds": int(n_seeds),
                        "n_nodes": int(n_window_nodes),

                        "total_size": int(total_size),

                        "genome_positions": genome_positions,

                        "seeds": [
                            node["node"]["name"]
                            for node in window_seed_nodes
                        ],

                        "nodes": [
                            node["node"]["name"]
                            for node in window_nodes
                        ],
                    })

            window_start += step

    # --------------------------------------------------------
    # Sort regions:
    #
    #   1. number of seeds
    #   2. number of nodes
    # --------------------------------------------------------

    regions.sort(
        key=lambda x: (
            x["n_seeds"],
            x["n_nodes"],
        ),
        reverse=True,
    )

    return regions


def find_regions(
        seq: str,
        k: int=13,
        canonical=True,
        max_seeds: int = 100,
        min_coverage = 20,
        min_node_size = None,
        max_selected_nodes=10000,

        ):


    start_time = time.time()
    size_seq = len(seq)
    seeds, filtered_nodes, nodes = get_seeds(seq, k=k, canonical=canonical, max_seeds=max_seeds,
                                             min_coverage=min_coverage, min_node_size=min_node_size,
                                             max_selected_nodes=max_selected_nodes)
    regions = find_dense_regions(filtered_nodes, seeds, sequence_size=size_seq, sequence_type="DNA", min_nodes= 5)

    logger.info(f"{len(regions)} regions found in {time.time()-start_time} s.")
    return regions




def mutate_adn(sequence, taux, seed=None):
    """
    Transforme une séquence ADN par mutations, délétions, insertions
    et inversions.

    Paramètres
    ----------
    sequence : str
        Séquence ADN contenant A, C, G et T.
    taux : float
        Pourcentage de bases à modifier, entre 0 et 100.
        Par exemple 10 signifie environ 10 %.
    seed : int, optionnel
        Graine pour rendre le résultat reproductible.

    Retour
    ------
    str
        Séquence ADN transformée.
    """

    if seed is not None:
        random.seed(seed)

    sequence = sequence.upper()

    # Vérification de la séquence
    if not all(base in "ACGT" for base in sequence):
        raise ValueError("La séquence doit uniquement contenir A, C, G et T.")

    if not 0 <= taux <= 100:
        raise ValueError("Le taux doit être compris entre 0 et 100.")

    if len(sequence) == 0:
        return sequence

    # Nombre approximatif de positions à traiter
    nombre_operations = round(len(sequence) * taux / 100)

    bases = "ACGT"

    for _ in range(nombre_operations):

        # On choisit aléatoirement une opération
        operation = random.choice([
            "mutation",
            "deletion",
            "insertion",
            "inversion"
        ])

        if len(sequence) == 0:
            break

        # ---------------------------------------------------------
        # MUTATION : A -> C, par exemple
        # ---------------------------------------------------------
        if operation == "mutation":

            position = random.randrange(len(sequence))
            ancienne_base = sequence[position]

            nouvelles_bases = [
                b for b in bases
                if b != ancienne_base
            ]

            nouvelle_base = random.choice(nouvelles_bases)

            sequence = (
                sequence[:position]
                + nouvelle_base
                + sequence[position + 1:]
            )

        # ---------------------------------------------------------
        # DELETION : suppression d'une base
        # ---------------------------------------------------------
        elif operation == "deletion":

            position = random.randrange(len(sequence))

            sequence = (
                sequence[:position]
                + sequence[position + 1:]
            )

        # ---------------------------------------------------------
        # INSERTION : ajout d'une base aléatoire
        # ---------------------------------------------------------
        elif operation == "insertion":

            position = random.randrange(len(sequence) + 1)
            nouvelle_base = random.choice(bases)

            sequence = (
                sequence[:position]
                + nouvelle_base
                + sequence[position:]
            )

        # ---------------------------------------------------------
        # INVERSION : inversion d'un segment
        # ---------------------------------------------------------
        elif operation == "inversion":

            if len(sequence) >= 2:

                debut = random.randrange(len(sequence))
                fin = random.randrange(debut + 1, len(sequence) + 1)

                segment = sequence[debut:fin]
                segment_inverse = segment[::-1]

                sequence = (
                    sequence[:debut]
                    + segment_inverse
                    + sequence[fin:]
                )

    return sequence

