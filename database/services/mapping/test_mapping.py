import logging


import numpy as np
import random

import os
import csv
import tkinter as tk
from tkinter import filedialog

from types import SimpleNamespace

from database.services.mapping.index_kmer_requests import encode_kmers

from database.services.mapping.Seed import Seed

from database.services.mapping.mapping import find_regions

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

"""
Transform a DNA sequence through mutations, deletions, insertions,
and inversions.

Parameters
----------
sequence : str
    DNA sequence containing A, C, G, and T.
rate : float
    Percentage of bases to modify, between 0 and 100.
    For example, 10 means approximately 10%.
seed : int, optional
    Random seed to make the result reproducible.

Returns
-------
str
    Transformed DNA sequence.
"""
def mutate_dna(sequence, mutation_rate, seed=None):

    if seed is not None:
        random.seed(seed)

    sequence = sequence.upper()

    # Validate the sequence
    if not all(base in "ACGT" for base in sequence):
        raise ValueError("The sequence must contain only A, C, G, and T.")

    if not 0 <= mutation_rate <= 100:
        raise ValueError("The rate must be between 0 and 100.")

    if len(sequence) == 0:
        return sequence

    # Approximate number of operations to perform
    num_operations = round(len(sequence) * mutation_rate / 100)

    bases = "ACGT"

    for _ in range(num_operations):

        # Randomly select an operation
        operation = random.choice([
            "mutation",
            "deletion",
            "insertion",
            "inversion"
        ])

        if len(sequence) == 0:
            break

        # ---------------------------------------------------------
        # MUTATION: A -> C, for example
        # ---------------------------------------------------------
        if operation == "mutation":

            position = random.randrange(len(sequence))
            old_base = sequence[position]

            possible_bases = [
                base for base in bases
                if base != old_base
            ]

            new_base = random.choice(possible_bases)

            sequence = (
                    sequence[:position]
                    + new_base
                    + sequence[position + 1:]
            )

        # ---------------------------------------------------------
        # DELETION: remove a base
        # ---------------------------------------------------------
        elif operation == "deletion":

            position = random.randrange(len(sequence))

            sequence = (
                    sequence[:position]
                    + sequence[position + 1:]
            )

        # ---------------------------------------------------------
        # INSERTION: add a random base
        # ---------------------------------------------------------
        elif operation == "insertion":

            position = random.randrange(len(sequence) + 1)
            new_base = random.choice(bases)

            sequence = (
                    sequence[:position]
                    + new_base
                    + sequence[position:]
            )

        # ---------------------------------------------------------
        # INVERSION: reverse a segment
        # ---------------------------------------------------------
        elif operation == "inversion":

            if len(sequence) >= 2:
                start = random.randrange(len(sequence))
                end = random.randrange(start + 1, len(sequence) + 1)

                segment = sequence[start:end]
                reversed_segment = segment[::-1]

                sequence = (
                        sequence[:start]
                        + reversed_segment
                        + sequence[end:]
                )

    return sequence



"""
This function will charge a selected csv file containing entries like : gene_name;sequence
Then it will search for gene into the gene names of this csv file
Then if mutation_rate is not None it will mutate the sequence with the desired mutation rate
Then it will search this sequence into the pangenome graph
"""
def search_gene(gene, mutation_rate=None):

    # --------------------------------------------------
    # 1. Validate mutation rate
    # --------------------------------------------------

    if mutation_rate is not None:
        if not isinstance(mutation_rate, (int, float)):
            raise TypeError("mutation_rate must be a number or None.")

        if not 0 <= mutation_rate <= 100:
            raise ValueError("mutation_rate must be between 0 and 100.")

    # --------------------------------------------------
    # 2. Select CSV file
    # --------------------------------------------------

    root = tk.Tk()
    root.withdraw()

    filepath = filedialog.askopenfilename(
        initialdir=os.getcwd(),
        title="Select gene CSV file",
        filetypes=[
            ("CSV files", "*.csv"),
            ("All files", "*.*"),
        ],
    )

    root.destroy()

    if not filepath:
        raise FileNotFoundError("No CSV file was selected.")


    # --------------------------------------------------
    # 3. Load available genes and search for the requested gene
    # --------------------------------------------------

    available_genes = []
    query = None
    # Increase the maximum CSV field size.
    csv.field_size_limit(10 ** 8)
    with open(filepath, "r", encoding="utf-8-sig", newline="") as csvfile:

        reader = csv.DictReader(csvfile, delimiter=";")

        if reader.fieldnames is None:
            raise ValueError("The CSV file is empty.")

        columns = {
            column.strip().lower(): column
            for column in reader.fieldnames
        }

        if "sequence" not in columns:
            raise ValueError("The CSV file must contain a 'sequence' column.")

        gene_column = columns.get("gene") or columns.get("gene_name")

        if gene_column is None:
            raise ValueError(
                "The CSV file must contain a 'gene' or 'nom_gene' column."
            )

        sequence_column = columns["sequence"]

        for row in reader:

            row_gene = row[gene_column].strip()
            available_genes.append(row_gene)

            if row_gene.lower() == gene.strip().lower():
                query = row[sequence_column].strip()

    # --------------------------------------------------
    # 4. Display available genes
    # --------------------------------------------------

    print("\n========== AVAILABLE GENES ==========")

    columns_per_row = 4

    for i in range(0, len(available_genes), columns_per_row):
        row_genes = available_genes[i:i + columns_per_row]

        print(" | ".join(f"{name:<20}" for name in row_genes))

    print(f"\nTotal genes: {len(available_genes)}")
    print("=====================================\n")

    # --------------------------------------------------
    # 5. Check requested gene
    # --------------------------------------------------

    if query is None:
        raise ValueError(f"Gene '{gene}' was not found in the CSV file.")


    # --------------------------------------------------
    # 6. Apply mutations if requested
    # --------------------------------------------------

    if mutation_rate is not None:

        mutate_sequence = mutate_dna(
            query,
            mutation_rate,
        )

    else:
        mutate_sequence = query

    # --------------------------------------------------
    # 7. Find regions
    # --------------------------------------------------

    regions, selected_seeds, all_seeds = find_regions(mutate_sequence)

    return regions, selected_seeds, all_seeds