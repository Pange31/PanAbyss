import logging


import numpy as np
from types import SimpleNamespace

from database.services.mapping.index_kmer_requests import encode_kmers

from database.services.mapping.Seed import Seed

logger = logging.getLogger("panabyss_logger")


"""
Test the synteny calculation between two DNA sequences
without requiring a Neo4j connection.

Parameters
----------
query_sequence : str
    Query DNA sequence.

node_sequence : str
    DNA sequence simulating a Neo4j node.
    It can be shorter or longer than the query.

k : int
    K-mer size.

canonical : bool
    Whether canonical k-mers are used.

max_kmer_distance : int
    Maximum allowed distance between consecutive diagonals.

max_query_occurrences : int
    Maximum number of query occurrences allowed for a k-mer.

timings : bool
    If True, return the execution timings.

Returns
-------
seed : Seed
    Temporary Seed object containing the synteny score,
    anchors and other intermediate results.

timing_results : dict or None
    Execution times for encoding, anchors and chains.
"""
def compute_synteny(
    query_sequence: str,
    node_sequence: str,
    k: int = 13,
    canonical: bool = True,
    max_kmer_distance: int = 20,
    max_query_occurrences: int = 50,
    timings: bool = False,
):


    # --------------------------------------------------
    # 1. Encode query k-mers
    # --------------------------------------------------

    query_codes, query_orientations = encode_kmers(
        query_sequence,
        k=k,
        canonical=canonical,
        orientations=True,
    )

    query_codes = np.ascontiguousarray(
        query_codes,
        dtype=np.uint32,
    )

    query_orientations = np.ascontiguousarray(
        query_orientations,
        dtype=np.int8,
    )

    query_positions = np.arange(
        query_codes.size,
        dtype=np.int64,
    )

    # Sort query k-mers by their codes
    order = np.argsort(query_codes, kind="stable")

    query = SimpleNamespace(
        sequence=query_sequence,
        k=k,
        canonical=canonical,
        has_kmers=query_codes.size > 0,
        sorted_kmer_codes=np.ascontiguousarray(
            query_codes[order],
            dtype=np.uint32,
        ),
        sorted_positions=np.ascontiguousarray(
            query_positions[order],
            dtype=np.int64,
        ),
        sorted_orientations=np.ascontiguousarray(
            query_orientations[order],
            dtype=np.int8,
        ),
    )

    # --------------------------------------------------
    # 2. Simulate a Seed / Neo4j sequence node
    # --------------------------------------------------

    seed = Seed(
        sequence_id=0,
        count=1,
        size=len(node_sequence),
        coverage=0.0,
        query_coverage=0.0,
        sequence=node_sequence,
    )

    # --------------------------------------------------
    # 3. Calculate synteny
    # --------------------------------------------------

    timing_results = {} if timings else None

    synteny = seed.calculate_synteny(
        query=query,
        max_kmer_distance=max_kmer_distance,
        max_query_occurrences=max_query_occurrences,
        timings=timing_results,
    )

    return seed, timing_results


"""
This function generate classical cases to test synteny:
    - test = "manual" => in this case the query and node sequences are directly used
    - test = "full_truncated" => synteny is fully conserved but there are absent kmers
        query = k1 + k2 + k3 + k4 + k5 + k6 + k7 + k8 + k9 + k10 + k11 + k12 + k13 + k14 + k15
        node =  k1 + ka + k2 + kb + k5 + kc + kd + ke + k9 + k10 + k11 + k12 + kf +  kg +  k15
        with ka-g kmer of size between 5 - 14 and absent from query sequence
    - test = "full_truncated_inversion" => synteny is fully conserved but there are absent kmers and an inversion
        query = k1 + k2 + k3 + k4 + k5 + k6 + k7 + k8 + k9 + k10 + k11 + k12 + k13 + k14 + k15
        node =  k1 + ka + k2 + kb + k5 + kc + kd + ke + (k9 + k10 + k11 + k12)inv + kf +  kg +  k15
        with ka-g kmer of size between 5 - 14 and absent from query sequence
    - test = "fragmented" => synteny is not fully conserved, there are 3 distinct groups
        query = k1 + k2 + k3 + k4 + k5 + k6 + k7 + k8 + k9 + k10 + k11 + k12 + k13 + k14 + k15
        node =  k14 + k15 + ka + k2 + k3 + k4 + kb + k1 
        with ka-b sequences of size > 20 bp
"""
def test_synteny(test="manual", query=None, node=None):

    anchors = [
        "AAGCCCAATAAAC",  # A1
        "GCTCTGACTGGCC",  # A2
        "TACGGTACCTGAA",  # A3
        "CGATGCTAGCTAC",  # A4
        "TTACCGGATCGTA",  # A5
        "ACGTAGCTTACCG",  # A6
        "GGCATACCTAGCT",  # A7
        "CTAGGCTAACGTA",  # A8
    ]

    query_spacers = [
        "AA", "CCCCC", "GG", "TTTTTTTT",
        "AC", "GTGTGT", "CAAA"
    ]

    node_spacers = [
        "T", "GAG", "CTCTCTCT", "A",
        "TCGTA", "AG", "ATATATATAT"
    ]

    # Reverse-complement a DNA sequence.
    def reverse_complement(sequence):
        return sequence.translate(
            str.maketrans("ACGT", "TGCA")
        )[::-1]

    match test:

        case "full_truncated":

            # Build the query sequence.
            query = anchors[0]

            for spacer, anchor in zip(query_spacers, anchors[1:]):
                query += spacer + anchor

            # Build the node sequence in the original anchor order.
            node = anchors[0]

            for spacer, anchor in zip(node_spacers, anchors[1:]):
                node += spacer + anchor

            node_anchor_sequences = anchors.copy()

        case "full_truncated_inversion":

            # Modify the spacer after A4 to avoid a boundary k-mer match.
            inversion_query_spacers = query_spacers.copy()
            inversion_query_spacers[3] = "GTTTTTTT"

            # Build the query sequence.
            query = anchors[0]

            for spacer, anchor in zip(inversion_query_spacers, anchors[1:]):
                query += spacer + anchor

            # Build the first four anchors in their original orientation.
            node = (
                anchors[0]
                + node_spacers[0] + anchors[1]
                + node_spacers[1] + anchors[2]
                + node_spacers[2] + anchors[3]
            )

            # Reverse-complement the complete region from A5 to A8.
            node += (
                reverse_complement(anchors[7])
                + reverse_complement(node_spacers[6])
                + reverse_complement(anchors[6])
                + reverse_complement(node_spacers[5])
                + reverse_complement(anchors[5])
                + reverse_complement(node_spacers[4])
                + reverse_complement(anchors[4])
                + reverse_complement(node_spacers[3])
            )

            # Define the expected anchor sequences in the node.
            node_anchor_sequences = [
                anchors[0],
                anchors[1],
                anchors[2],
                anchors[3],
                reverse_complement(anchors[4]),
                reverse_complement(anchors[5]),
                reverse_complement(anchors[6]),
                reverse_complement(anchors[7]),
            ]

        case "fragmented":

            # Build the query sequence with all eight anchors in their original order.
            query = anchors[0]

            for spacer, anchor in zip(query_spacers, anchors[1:]):
                query += spacer + anchor

            # Build the node sequence with three fragmented groups:
            #
            # Group 1: A7 -> A8
            # Group 2: A2 -> A3 -> A4
            # Group 3: A1
            #
            # Groups are separated by more than 20 bp.

            long_spacer = "N" * 30

            node = (
                    anchors[6]
                    + node_spacers[6]
                    + anchors[7]
                    + long_spacer
                    + anchors[1]
                    + node_spacers[1]
                    + anchors[2]
                    + node_spacers[2]
                    + anchors[3]
                    + long_spacer
                    + anchors[0]
            )

            # Keep all eight anchors to preserve the correspondence
            # with the query. Missing anchors will have position -1.
            node_anchor_sequences = anchors.copy()
        case "manual":

            # Use the manually provided sequences.
            if query is None or node is None:
                raise ValueError(
                    "Both query and node must be provided in manual mode."
                )

            node_anchor_sequences = anchors

        case _:
            raise ValueError(f"Unknown test case: {test}")

    # Calculate anchor positions in both sequences.
    query_positions = [
        query.find(anchor) for anchor in anchors
    ]

    node_positions = [
        node.find(anchor) for anchor in node_anchor_sequences
    ]

    # Display distances between consecutive anchors.
    print("\nDistances between anchors:")

    for i in range(1, len(anchors)):

        if query_positions[i] == -1 or query_positions[i - 1] == -1:
            query_distance = "N/A"
        else:
            query_distance = query_positions[i] - query_positions[i - 1]

        if node_positions[i] == -1 or node_positions[i - 1] == -1:
            node_distance = "N/A"
        else:
            node_distance = node_positions[i] - node_positions[i - 1]

        print(
            f"Anchor {i}: "
            f"Query = {query_distance} bp | "
            f"Node = {node_distance} bp"
        )

    # Compute synteny.
    seed, timings = compute_synteny(
        query_sequence=query,
        node_sequence=node,
        k=13,
        canonical=True,
        max_kmer_distance=20,
        timings=True,
    )

    print("\nDetected anchors:")

    for i in range(len(seed.anchor_query)):
        q_pos = seed.anchor_query[i]
        n_pos = seed.anchor_node[i]

        print(
            f"Anchor {i + 1}: "
            f"Query position = {q_pos} | "
            f"Node position = {n_pos} | "
            f"Orientation = {seed.anchor_orientation[i]} | "
            f"Score = {seed.anchor_scores[i]}"
        )

    print("\n========== SYNTENY TEST ==========")
    print("Test:", test)
    print("Query:", query)
    print("Node: ", node)
    print("Query length:", len(query))
    print("Node length:", len(node))
    print("Synteny:", seed.synteny)
    print("Number of anchors:", len(seed.anchor_query))
    print("Anchor scores:", seed.anchor_scores)
    print("Timings:", timings)
    print("==================================\n")

    return query, node, seed