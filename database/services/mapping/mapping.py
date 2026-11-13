"""
Created on 12/11/2026

@author: fgraziani

This file implements main functions to map a query sequence on the pangenome graph
Main algorithm steps :

1. Retrieves kmer of the query in the graph kmer indexes
2. Filter the detected nodes to have a minimum query coverage or a minimum node coverage
3. Selects seeds (with a maximum seeds parameters set by default to 10000):
    3.a Select the best query coverage nodes (only for nodes with sufficient query coverage)
    3.b Select the nodes with higher node coverage to complet the max_seeds seeds
4. For each seed compute a synteny score (with a LIS algorithm)
5. Cluster seeds by genomic distance for each genome
6. Merge clustered if they share a minimum of common nodes (set by default to 20%) => regions
7. Compute a global synteny score on each region and general statistics

Notes
----------
- Small kmer are used (13 by default and max 15) to detect small nodes and to be resilient to mutations
"""
import logging
import time

import numpy as np

from database.services.mapping.index_kmer_requests import get_encoded_kmer_node_counts_batch, get_node_sizes_batch, encode_kmers
from database.services.neo4j_requests import get_nodes_from_sequence_id, get_sequence_from_ids

from database.services.mapping.Query import Query
from database.services.mapping.Seed import Seed, calculate_chains_numba

logger = logging.getLogger("panabyss_logger")


"""
Select seed indices based on query and node coverage.

Nodes passing at least one coverage threshold are retained.
Nodes meeting the query coverage threshold have priority and
are ranked by decreasing query coverage, with node coverage
as a secondary criterion.

Remaining slots are filled with other retained nodes ranked
by decreasing node coverage.

Parameters
----------
node_coverage : np.ndarray
    Node coverage percentages.
query_coverage : np.ndarray
    Query coverage percentages.
min_node_coverage : float
    Minimum node coverage required to retain a node.
min_query_coverage : float
    Minimum query coverage required to retain a node and
    qualify it for priority selection.
max_selected_seeds : int
    Maximum number of selected nodes.

Returns
-------
candidate_indices : np.ndarray
    Indices of nodes passing at least one coverage threshold.
selected_indices : np.ndarray
    Indices of selected nodes relative to candidate_indices,
    ordered by selection priority.
"""
def select_seed_indices(
    node_coverage,
    query_coverage,
    min_node_coverage,
    min_query_coverage,
    max_selected_seeds,
):

    if max_selected_seeds <= 0 or len(node_coverage) == 0:
        return (
            np.empty(0, dtype=np.int64),
            np.empty(0, dtype=np.int64),
        )

    # --------------------------------------------------------
    # Filter candidates.
    # --------------------------------------------------------

    coverage_mask = (
        (node_coverage >= min_node_coverage)
        | (query_coverage >= min_query_coverage)
    )

    candidate_indices = np.flatnonzero(coverage_mask)

    if len(candidate_indices) == 0:
        return (
            candidate_indices,
            np.empty(0, dtype=np.int64),
        )

    candidate_node_coverage = node_coverage[candidate_indices]
    candidate_query_coverage = query_coverage[candidate_indices]

    # --------------------------------------------------------
    # Priority 1: query coverage.
    # --------------------------------------------------------

    priority_mask = candidate_query_coverage >= min_query_coverage

    priority_indices = np.flatnonzero(priority_mask)

    #Sort by query coverage (only nodes with query coverage >= min_query_coverage) and then by node_coverage
    priority_order = np.lexsort((
        -candidate_node_coverage[priority_indices],
        -candidate_query_coverage[priority_indices],
    ))

    priority_indices = priority_indices[priority_order]
    priority_indices = priority_indices[:max_selected_seeds]

    # --------------------------------------------------------
    # Priority 2: node coverage.
    # --------------------------------------------------------

    remaining_slots = max_selected_seeds - len(priority_indices)

    if remaining_slots <= 0:
        return candidate_indices, priority_indices

    remaining_indices = np.flatnonzero(~priority_mask)

    #sort by node coverage and order by decreasing order
    remaining_order = np.argsort(candidate_node_coverage[remaining_indices])[::-1]

    remaining_indices = remaining_indices[remaining_order]
    remaining_indices = remaining_indices[:remaining_slots]

    selected_indices = np.concatenate((
        priority_indices,
        remaining_indices,
    ))

    return candidate_indices, selected_indices


"""
Find candidate sequences in the pangenome graph that match a query sequence.

The function:
    1. Builds a Query object from the input sequence.
    2. Retrieves graph nodes matching the query k-mers.
    3. Filters nodes by sequence size.
    4. Computes node and query coverage for each candidate.
    5. Retains nodes meeting at least one coverage threshold.
    6. Selects up to max_selected_seeds candidates, prioritizing nodes
       meeting the minimum query coverage threshold.
    7. Retrieves the selected nodes' k-mers and DNA sequences.
    8. Creates Seed objects and calculates synteny.
    9. Retrieves all physical Neo4j Node occurrences.

Parameters
----------
seq : str
    Input DNA sequence.

k : int
    K-mer size.

canonical : bool, default=True
    Whether canonical k-mer encoding is used.

min_node_coverage : float, default=20
    Minimum node coverage percentage required to retain a candidate.

min_query_coverage : float, default=20
    Minimum query coverage percentage required to retain a candidate.
    Candidates meeting this threshold are given selection priority.

min_node_size : int or None, default=None
    Minimum node sequence size.
    If None, the minimum size is set to k + 4.

max_selected_seeds : int, default=10000
    Maximum number of candidates selected for Seed construction
    and synteny calculation.

max_kmer_distance : int, default=20
    Maximum difference between consecutive anchor diagonals
    when clustering anchors into collinear chains.

max_query_occurrences : int, default=50
    Maximum number of occurrences of a query k-mer allowed
    during anchor construction.

timings : bool, default=False
    Whether to display execution timings.

Returns
-------
selected_seeds : dict
    Dictionary mapping Sequence node IDs to Seed objects.

all_seeds : np.ndarray
    Structured array containing all candidates passing the node size
    and coverage filters, before limiting the number of selected seeds.

query : Query
    Prepared query object.
"""
def get_seeds(
    seq: str,
    k: int,
    canonical: bool = True,
    min_node_coverage: float = 20,
    min_query_coverage: float = 20,
    min_node_size: int = None,
    max_selected_seeds: int = 10000,
    max_kmer_distance: int = 20,
    max_query_occurrences: int = 50,
    timings: bool = False,
):

    t_total = time.perf_counter()

    # --------------------------------------------------------
    # Validate input parameters.
    # --------------------------------------------------------

    if not isinstance(seq, str):
        raise TypeError("seq must be a string.")

    if k <= 0:
        raise ValueError("k must be > 0.")

    if max_selected_seeds < 0:
        raise ValueError("max_selected_seeds must be >= 0.")

    if min_node_size is None or min_node_size < k:
        min_node_size = k + 4

    empty_dtype = np.dtype([
        ("node_id", "<u8"),
        ("count", "<u8"),
        ("size", "<u8"),
        ("coverage", "<f8"),
        ("query_coverage", "<f8"),
    ])

    empty_seeds = np.empty(0, dtype=empty_dtype)

    if len(seq) < k:
        logger.warning("Query sequence is shorter than k.")
        return {}, empty_seeds, None

    # --------------------------------------------------------
    # Prepare query.
    # --------------------------------------------------------

    t0 = time.perf_counter()

    query = Query(
        sequence=seq,
        k=k,
        canonical=canonical,
    )

    t_query = time.perf_counter() - t0

    if not query.has_kmers:
        logger.warning("No valid k-mers found in query.")
        return {}, empty_seeds, query

    # Use unique query k-mers for index lookup.
    unique_query_kmer_codes = query.unique_kmer_codes

    # --------------------------------------------------------
    # Retrieve matching graph nodes.
    # --------------------------------------------------------

    t0 = time.perf_counter()

    node_ids, node_counts, node_kmer_index = (
        get_encoded_kmer_node_counts_batch(
            unique_query_kmer_codes,
            k=k,
            canonical=canonical,
            return_node_kmer_index=True,
            timings=timings
        )
    )

    t_index = time.perf_counter() - t0

    if len(node_ids) == 0:
        logger.info("No graph nodes match the query k-mers.")
        return {}, empty_seeds, query

    # --------------------------------------------------------
    # Retrieve node sizes.
    # --------------------------------------------------------

    t0 = time.perf_counter()

    node_sizes = get_node_sizes_batch(k=k,node_ids=node_ids)

    t_sizes = time.perf_counter() - t0

    # --------------------------------------------------------
    # Filter nodes by size.
    # --------------------------------------------------------

    t0 = time.perf_counter()

    valid_mask = (node_sizes > 0) & (node_sizes >= min_node_size)


    valid_node_ids = node_ids[valid_mask]
    valid_counts = node_counts[valid_mask]
    valid_sizes = node_sizes[valid_mask]

    if len(valid_node_ids) == 0:
        logger.info("No nodes passed the minimum size filter.")
        return {}, empty_seeds, query


    # --------------------------------------------------------
    # Calculate coverage.
    # --------------------------------------------------------

    max_kmer_counts = valid_sizes - k + 1

    node_coverage = (
        100.0 * valid_counts / max_kmer_counts
    )

    query_coverage = (
        100.0 * valid_counts / len(unique_query_kmer_codes)
    )

    # --------------------------------------------------------
    # Filter and select candidates.
    # --------------------------------------------------------

    t0 = time.perf_counter()

    candidate_indices, selected_indices = select_seed_indices(
        node_coverage=node_coverage,
        query_coverage=query_coverage,
        min_node_coverage=min_node_coverage,
        min_query_coverage=min_query_coverage,
        max_selected_seeds=max_selected_seeds,
    )

    t_filtering = time.perf_counter() - t0

    if len(candidate_indices) == 0:
        logger.info("No nodes passed the coverage filters.")
        return {}, empty_seeds, query

    # --------------------------------------------------------
    # Build candidate array.
    # --------------------------------------------------------

    all_seeds = np.empty(
        len(candidate_indices),
        dtype=empty_dtype,
    )

    all_seeds["node_id"] = valid_node_ids[candidate_indices]
    all_seeds["count"] = valid_counts[candidate_indices]
    all_seeds["size"] = valid_sizes[candidate_indices]
    all_seeds["coverage"] = node_coverage[candidate_indices]
    all_seeds["query_coverage"] = query_coverage[candidate_indices]

    # --------------------------------------------------------
    # Retrieve selected candidates.
    # --------------------------------------------------------

    t0 = time.perf_counter()

    selected_data = all_seeds[selected_indices]

    # --------------------------------------------------------
    # Prepare the query-matching k-mers for selected seeds.
    # --------------------------------------------------------

    unique_nodes, starts, ends, sorted_kmers = node_kmer_index

    selected_node_ids = selected_data["node_id"]

    # Find the position of each selected node in the sorted index.
    indexes = np.searchsorted(unique_nodes, selected_node_ids)

    # Identify positions that are within the bounds of unique_nodes.
    valid = indexes < len(unique_nodes)

    # For positions within bounds, verify that the node ID actually matches.
    # searchsorted returns the insertion position, which does not guarantee
    # that the searched node exists in the index.
    valid_positions = np.flatnonzero(valid)

    valid[valid_positions] = (
            unique_nodes[indexes[valid_positions]]
            == selected_node_ids[valid_positions]
    )

    selected_kmers = {}

    for node_id, index, is_valid in zip(selected_node_ids, indexes, valid):

        node_id = int(node_id)

        if is_valid:
            selected_kmers[node_id] = np.ascontiguousarray(
                sorted_kmers[starts[index]:ends[index]],
                dtype=np.uint32,
            )
        else:
            selected_kmers[node_id] = np.empty(0, dtype=np.uint32)

    selected_ids = selected_node_ids.tolist()

    t_selection = time.perf_counter() - t0

    if not selected_ids:
        logger.info("No seeds selected.")
        return {}, all_seeds, query

    # --------------------------------------------------------
    # Retrieve selected sequences.
    # --------------------------------------------------------

    t0 = time.perf_counter()

    node_sequences = get_sequence_from_ids(selected_ids)

    t_sequences = time.perf_counter() - t0

    if node_sequences is None:
        raise RuntimeError(
            "Unable to retrieve selected sequences from Neo4j."
        )

    # --------------------------------------------------------
    # Create Seed objects.
    # --------------------------------------------------------

    selected_seeds = {}

    missing_sequences = 0

    for node in selected_data:

        sequence_id = int(node["node_id"])

        node_sequence = node_sequences.get(sequence_id)

        # Support string keys if the Neo4j service returns them.
        if node_sequence is None:
            node_sequence = node_sequences.get(str(sequence_id))

        if node_sequence is None:
            missing_sequences += 1
            continue

        seed = Seed(
            sequence_id=sequence_id,
            count=int(node["count"]),
            size=int(node["size"]),
            coverage=float(node["coverage"]),
            query_coverage=float(node["query_coverage"]),
            sequence=node_sequence,
            query_coverage_count=int(node["count"]),
            kmer_codes=selected_kmers[sequence_id],
        )

        selected_seeds[sequence_id] = seed

    if missing_sequences:
        logger.warning(
            "Could not retrieve %d selected sequences.",
            missing_sequences,
        )

    if not selected_seeds:
        logger.warning("No selected sequences could be retrieved.")
        return {}, all_seeds, query

    # --------------------------------------------------------
    # Calculate synteny.
    # --------------------------------------------------------

    t0 = time.perf_counter()

    synteny_timings = {
        "encoding": 0.0,
        "anchors": 0.0,
        "chains": 0.0,
    }

    for seed in selected_seeds.values():

        seed.calculate_synteny(
            query=query,
            max_kmer_distance=max_kmer_distance,
            max_query_occurrences=max_query_occurrences,
            timings=synteny_timings if timings else None,
        )

    t_synteny = time.perf_counter() - t0

    # --------------------------------------------------------
    # Retrieve physical Neo4j Node occurrences.
    # --------------------------------------------------------

    t0 = time.perf_counter()

    neo4j_nodes = get_nodes_from_sequence_id(list(selected_seeds.keys()))
    t_neo4j = time.perf_counter() - t0
    t0 = time.perf_counter()
    if neo4j_nodes:

        for result in neo4j_nodes:

            sequence_id = int(result["sequence_id"])

            seed = selected_seeds.get(sequence_id)

            if seed is None:
                continue

            # Keep all properties of the physical Node occurrence.
            node_properties = dict(result["node"])

            seed.add_neo4j_node(node_properties)

    t_occurrences = time.perf_counter() - t0

    # --------------------------------------------------------
    # Log results.
    # --------------------------------------------------------

    logger.info(
        "Seed detection completed: %d candidates, %d selected.",
        len(all_seeds),
        len(selected_seeds),
    )

    if timings:

        logger.debug("\n============= GET SEEDS TIMINGS =============")
        logger.debug(f"Query k-mers:          {query.n_kmers:,}")
        logger.debug(f"Unique query k-mers:   {query.n_unique_kmers:,}")
        logger.debug(f"Candidate nodes:       {len(node_ids):,}")
        logger.debug(f"Filtered nodes:        {len(all_seeds):,}")
        logger.debug(f"Selected seeds:        {len(selected_seeds):,}")
        logger.debug("---------------------------------------------")
        logger.debug(f"Neo4j nodes retrieval: {t_neo4j:.3f} s")
        logger.debug(f"Query preparation:     {t_query:.3f} s")
        logger.debug(f"Index retrieval:       {t_index:.3f} s")
        logger.debug(f"Node sizes:            {t_sizes:.3f} s")
        logger.debug(f"Filtering:             {t_filtering:.3f} s")
        logger.debug(f"Seed selection:        {t_selection:.3f} s")
        logger.debug(f"Sequence retrieval:    {t_sequences:.3f} s")
        logger.debug(f"Synteny calculation:   {t_synteny:.3f} s")
        logger.debug(f"Occurrence retrieval:  {t_occurrences:.3f} s")

        if timings:
            logger.debug("---------------------------------------------")
            logger.debug(
                f"  K-mer encoding:      "
                f"{synteny_timings['encoding']:.3f} s"
            )
            logger.debug(
                f"  Anchor construction: "
                f"{synteny_timings['anchors']:.3f} s"
            )
            logger.debug(
                f"  Chain calculation:   "
                f"{synteny_timings['chains']:.3f} s"
            )

        logger.debug("---------------------------------------------")
        logger.debug(
            f"Total:                 "
            f"{time.perf_counter() - t_total:.3f} s"
        )
        logger.debug("=============================================\n")

    return selected_seeds, all_seeds, query


"""
Disjoint-set data structure used to build genomic clusters with Union-Find algorithm.

Each physical Neo4j Node occurrence is represented by
its unique occurrence name.
"""
# class UnionFind:
#
#
#     def __init__(self, elements):
#
#         self.parent = {
#             element: element
#             for element in elements
#         }
#
#         self.rank = {
#             element: 0
#             for element in elements
#         }
#
#     def find(self, element):
#
#         if self.parent[element] != element:
#             self.parent[element] = self.find(
#                 self.parent[element]
#             )
#
#         return self.parent[element]
#
#     def union(self, element_a, element_b):
#
#         root_a = self.find(element_a)
#         root_b = self.find(element_b)
#
#         if root_a == root_b:
#             return
#
#         if self.rank[root_a] < self.rank[root_b]:
#             root_a, root_b = root_b, root_a
#
#         self.parent[root_b] = root_a
#
#         if self.rank[root_a] == self.rank[root_b]:
#             self.rank[root_a] += 1




"""
Merge genome clusters based on their shared occurrence IDs.

Two clusters are merged when the number of shared occurrences,
divided by the size of the larger cluster, reaches the specified
overlap threshold.
Uses a Union-Find data structure to merge connected clusters.

Clusters are merged transitively: if cluster A overlaps with B,
and B overlaps with C, all three clusters are merged, even if
A and C do not directly overlap.

Parameters
----------
    genome_clusters (list[set]): Collection of genome clusters,
        where each cluster contains occurrence IDs.
    min_overlap (float, optional): Minimum overlap ratio required
        to merge two clusters. Calculated as the number of shared
        occurrences divided by the size of the larger cluster.
        Defaults to 0.10.

Returns
----------
    list[set]: Merged clusters, with each cluster containing the
        union of occurrence IDs from all merged clusters.
"""

def merge_clusters_by_overlap(genome_clusters, min_overlap=0.10):

    n_clusters = len(genome_clusters)

    if n_clusters <= 1:
        return list(genome_clusters)

    cluster_sizes = [len(cluster) for cluster in genome_clusters]

    # Build inverted index: occurrence ID -> cluster IDs.
    occurrence_to_clusters = {}

    for cluster_id, cluster in enumerate(genome_clusters):

        for occurrence_id in cluster:

            if occurrence_id not in occurrence_to_clusters:
                occurrence_to_clusters[occurrence_id] = []

            occurrence_to_clusters[occurrence_id].append(cluster_id)

    # Count shared occurrences between cluster pairs.
    shared_counts = {}

    for cluster_ids in occurrence_to_clusters.values():

        n = len(cluster_ids)

        for i in range(n):

            cluster_a = cluster_ids[i]

            for j in range(i + 1, n):

                cluster_b = cluster_ids[j]
                pair = (cluster_a, cluster_b)

                shared_counts[pair] = (
                    shared_counts.get(pair, 0) + 1
                )

    # Union-Find using lists.
    parent = list(range(n_clusters))
    rank = [0] * n_clusters

    def find(x):

        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]

        return x

    def union(a, b):

        root_a = find(a)
        root_b = find(b)

        if root_a == root_b:
            return

        if rank[root_a] < rank[root_b]:
            parent[root_a] = root_b

        elif rank[root_a] > rank[root_b]:
            parent[root_b] = root_a

        else:
            parent[root_b] = root_a
            rank[root_a] += 1

    # Merge sufficiently overlapping clusters.
    for (cluster_a, cluster_b), shared_count in shared_counts.items():

        denominator = max(cluster_sizes[cluster_a],cluster_sizes[cluster_b])

        if denominator == 0:
            continue

        overlap = shared_count / denominator

        if overlap >= min_overlap:
            union(cluster_a, cluster_b)

    # Reconstruct merged clusters.
    merged_clusters = {}

    for cluster_id, cluster in enumerate(genome_clusters):

        root = find(cluster_id)

        if root not in merged_clusters:
            merged_clusters[root] = set()

        merged_clusters[root].update(cluster)

    return list(merged_clusters.values())

"""
Calculate a global synteny score for a cluster of Seeds.

The procedure:
    1. Collect anchors from all Seeds.
    2. Keep the best Seed for each query position.
    3. Order Seeds according to their median query position.
    4. Build virtual concatenated Seed coordinates.
    5. Calculate synteny using the existing chain algorithm.

Parameters
----------
seeds : list[Seed]
    Unique Seeds belonging to the cluster.

max_kmer_distance : int
    Maximum distance used for diagonal grouping.

Returns
-------
float
    Global cluster synteny score.
"""
def calculate_cluster_synteny(
    seeds,
    total_query_kmer_count,
    max_kmer_distance=20,
):

    if not seeds:
        return 0.0

    # --------------------------------------------------------
    # Collect anchors from all Seeds.
    # --------------------------------------------------------

    candidates = []

    for seed in seeds:

        if seed.anchor_query.size == 0:
            continue

        for i in range(seed.anchor_query.size):

            candidates.append({
                "query": int(seed.anchor_query[i]),
                "node": int(seed.anchor_node[i]),
                "orientation": int(seed.anchor_orientation[i]),
                "seed_id": seed.sequence_id,
                "anchor_score": int(seed.anchor_scores[i]),
            })

    if not candidates:
        return 0.0

    # --------------------------------------------------------
    # Keep the best anchor for each query position.
    # --------------------------------------------------------

    best_anchors = {}

    for anchor in candidates:

        query_position = anchor["query"]
        previous = best_anchors.get(query_position)

        if previous is None:
            best_anchors[query_position] = anchor

        elif anchor["anchor_score"] > previous["anchor_score"]:
            best_anchors[query_position] = anchor

    selected_anchors = sorted(
        best_anchors.values(),
        key=lambda anchor: anchor["query"],
    )

    if len(selected_anchors) < 2:
        return 0.0

    # --------------------------------------------------------
    # Reconstruct virtual node coordinates.
    #
    # Within the same Seed:
    #     Preserve the original node distance.
    #
    # Between different Seeds:
    #     Use the distance observed in the query.
    #
    # Coordinates are reconstructed independently
    # for each orientation.
    # --------------------------------------------------------

    n = len(selected_anchors)

    anchor_query = np.array(
        [anchor["query"] for anchor in selected_anchors],
        dtype=np.int64,
    )

    anchor_node = np.zeros(n, dtype=np.int64)

    anchor_orientation = np.array(
        [anchor["orientation"] for anchor in selected_anchors],
        dtype=np.int8,
    )

    for orientation in (1, -1):

        orientation_indices = [
            i
            for i in range(n)
            if anchor_orientation[i] == orientation
        ]

        if not orientation_indices:
            continue

        # First anchor establishes the virtual coordinate origin.
        first_index = orientation_indices[0]
        anchor_node[first_index] = 0

        for j in range(1, len(orientation_indices)):

            previous_index = orientation_indices[j - 1]
            current_index = orientation_indices[j]

            previous = selected_anchors[previous_index]
            current = selected_anchors[current_index]

            query_distance = current["query"] - previous["query"]

            if current["seed_id"] == previous["seed_id"]:

                # Same Seed: preserve the original node distance.
                node_distance = current["node"] - previous["node"]

            else:
                # Different Seeds: reconstruct the distance
                # using the query coordinates.
                node_distance = (
                    query_distance
                    if orientation == 1
                    else -query_distance
                )

            anchor_node[current_index] = anchor_node[previous_index] + node_distance

    # --------------------------------------------------------
    # Calculate global synteny.
    # --------------------------------------------------------
    synteny, _ = calculate_chains_numba(
        np.ascontiguousarray(anchor_query),
        np.ascontiguousarray(anchor_node),
        np.ascontiguousarray(anchor_orientation),
        max_kmer_distance,
        total_query_kmer_count
    )

    return float(synteny)
#
# def calculate_cluster_synteny_old(
#     seeds,
#     max_kmer_distance=20,
# ):
#
#
#     if not seeds:
#         return 0.0
#
#     # --------------------------------------------------------
#     # Collect anchors from all Seeds.
#     # --------------------------------------------------------
#
#     candidates = []
#
#     for seed in seeds:
#
#         if seed.anchor_query.size == 0:
#             continue
#
#         for i in range(seed.anchor_query.size):
#
#             candidates.append({
#                 "query": int(seed.anchor_query[i]),
#                 "node": int(seed.anchor_node[i]),
#                 "orientation": int(seed.anchor_orientation[i]),
#                 "seed_id": seed.sequence_id,
#                 "anchor_score": seed.anchor_scores[i]
#             })
#
#     if not candidates:
#         return 0.0
#
#     # --------------------------------------------------------
#     # Keep the best anchor for each query position.
#     # --------------------------------------------------------
#
#     best_anchors = {}
#
#     for anchor in candidates:
#         query_position = anchor["query"]
#         previous = best_anchors.get(query_position)
#
#         if previous is None:
#             best_anchors[query_position] = anchor
#
#         elif anchor["anchor_score"] > previous["anchor_score"]:
#             best_anchors[query_position] = anchor
#
#     selected_anchors = list(best_anchors.values())
#
#     if len(selected_anchors) < 2:
#         return 0.0
#
#     # --------------------------------------------------------
#     # Group anchors by Seed.
#     # --------------------------------------------------------
#
#     seed_anchors = {}
#
#     for anchor in selected_anchors:
#         seed_id = anchor["seed_id"]
#         seed_anchors.setdefault(seed_id,[],).append(anchor)
#
#     # --------------------------------------------------------
#     # Determine the order of Seeds.
#     #
#     # Seeds are ordered according to the median position
#     # of their retained anchors in the Query.
#     # --------------------------------------------------------
#
#     seed_order = sorted(
#         seed_anchors,
#         key=lambda seed_id: np.median([
#             anchor["query"]
#             for anchor in seed_anchors[seed_id]
#         ]),
#     )
#
#     # --------------------------------------------------------
#     # Build virtual concatenated Seed coordinates.
#     #
#     # Each Seed occupies a separate coordinate interval.
#     # --------------------------------------------------------
#
#     seed_lookup = {
#         seed.sequence_id: seed
#         for seed in seeds
#     }
#
#     anchor_query = []
#     anchor_node = []
#     anchor_orientation = []
#
#     offset = 0
#
#     for seed_id in seed_order:
#
#         anchors = seed_anchors[seed_id]
#
#         seed = seed_lookup[seed_id]
#
#         # Determine the virtual interval size.
#         seed_size = (
#             int(np.max(seed.anchor_node)) + 1
#             if seed.anchor_node.size > 0
#             else 0
#         )
#
#         for anchor in anchors:
#             anchor_query.append(anchor["query"])
#             anchor_node.append(offset + anchor["node"])
#             anchor_orientation.append(anchor["orientation"])
#
#         offset += seed_size
#
#     # --------------------------------------------------------
#     # Convert to contiguous NumPy arrays.
#     # --------------------------------------------------------
#
#     anchor_query = np.ascontiguousarray(
#         anchor_query,
#         dtype=np.int64,
#     )
#
#     anchor_node = np.ascontiguousarray(
#         anchor_node,
#         dtype=np.int64,
#     )
#
#     anchor_orientation = np.ascontiguousarray(
#         anchor_orientation,
#         dtype=np.int8,
#     )
#
#     # --------------------------------------------------------
#     # Calculate global synteny.
#     # --------------------------------------------------------
#     synteny, _ = calculate_chains_numba(
#         anchor_query,
#         anchor_node,
#         anchor_orientation,
#         max_kmer_distance,
#     )
#
#     return float(synteny)


"""
Cluster physical Neo4j node occurrences into genomic regions and
calculate their associated synteny and k-mer statistics.

Physical node occurrences from the selected Seeds are first grouped
by chromosome. Within each chromosome, overlapping or nearby genomic
intervals are clustered independently for each genome.

Clusters sharing physical node occurrences are then merged to obtain
consolidated genomic regions. Occurrences without genomic position
information are preserved as singleton clusters.

For each resulting region, the function calculates:
    - The number of unique Seeds and physical node occurrences.
    - The genomic position ranges across genomes.
    - The total physical size of the region.
    - The sum of individual Seed synteny scores.
    - The global cluster synteny score, calculated from the anchors
      of all Seeds using calculate_cluster_synteny().
    - The number of unique and shared k-mers.
    - The percentage of query k-mers covered by the region.

Regions are filtered according to their minimum number of physical
node occurrences and sorted by global synteny score, number of
unique Seeds and total physical size.

Parameters
----------
selected_seeds : dict
    Dictionary mapping sequence IDs to their corresponding Seed
    objects. Each Seed must provide its Neo4j node occurrences,
    synteny score and matched k-mer codes.

query : Query
    Query object providing the total number of unique query k-mers
    through its n_unique_kmers attribute.

max_distance : int, default (None) = 10% of query length
    Maximum allowed genomic distance between consecutive intervals
    when building clusters independently for each genome.

min_nodes : int, default=1
    Minimum number of physical Neo4j node occurrences required
    for a region to be retained.

min_nodes : int, default=0.1
    Minimum synteny score for a region to be retained.

max_kmer_distance : int, default=100
    Maximum allowed diagonal distance passed to
    calculate_cluster_synteny() when computing the global
    synteny score of each region.

min_overlap : float, default=0.10
    Minimum overlap ratio required when merging clusters that
    share physical node occurrences. Must be between 0 and 1.

Returns
-------
list of dict
    Regions sorted by decreasing global synteny score. Each region
    contains chromosome information, genomic position ranges,
    Seed and node counts, total size, synteny scores, k-mer
    statistics and physical node occurrence names.

"""

def cluster_selected_seeds(
    selected_seeds,
    query,
    max_distance=None,
    min_nodes=1,
    min_synteny_score = 0.1,
    max_kmer_distance=100,
    min_overlap=0.10,
    print_chr_stats=False,
):

    if max_distance is None or max_distance < 0:
        if max_distance is None:
            min_d = min(len(query.sequence), 1000)
            max_distance = max(min_d, 0.1 * len(query.sequence))
        else:
            raise ValueError("max_distance must be >= 0.")

    if min_nodes < 1:
        raise ValueError("min_nodes must be >= 1.")

    if not 0.0 <= min_synteny_score <= 1.0:
        raise ValueError("min_synteny_score must be between 0 and 1")

    if not 0 <= min_overlap <= 1:
        raise ValueError("min_overlap must be between 0 and 1.")

    if not selected_seeds:
        return []

    # --------------------------------------------------------
    # Group physical occurrences by chromosome.
    # --------------------------------------------------------

    chromosomes = {}
    seen_occurrences = set()

    for sequence_id, seed in selected_seeds.items():

        for node_data in seed.neo4j_nodes:

            occurrence_id = node_data.get("name")

            if occurrence_id is None:
                raise ValueError(
                    "A physical Neo4j Node occurrence has no 'name'."
                )

            # Avoid duplicate physical occurrences.
            if occurrence_id in seen_occurrences:
                continue

            seen_occurrences.add(occurrence_id)

            node_data = dict(node_data)
            node_data["sequence_id"] = int(sequence_id)

            chromosome = node_data.get("chromosome")

            if chromosome is None:
                continue

            if chromosome not in chromosomes:
                chromosomes[chromosome] = {}

            chromosomes[chromosome][occurrence_id] = node_data

    regions = []
    chromosome_stats = {}

    # ========================================================
    # Process each chromosome independently.
    # ========================================================

    for chromosome, nodes_by_id in chromosomes.items():

        chromosome_stats[chromosome] = {
            "n_regions": 0,
            "regions": [],
        }

        occurrence_ids = list(nodes_by_id.keys())

        # ----------------------------------------------------
        # Build genomic position indexes.
        # ----------------------------------------------------

        genome_index = {}
        node_positions = {}

        for occurrence_id, node_data in nodes_by_id.items():

            node_positions[occurrence_id] = {}

            genomes = node_data.get("genomes", [])

            if isinstance(genomes, str):
                genomes = [genomes]

            node_size = int(node_data.get("size", 0))

            for genome in genomes:
                position = node_data.get(f"{genome}_position")

                if position is None:
                    continue

                position = int(position)

                if genome in node_positions[occurrence_id]:
                    continue

                if genome not in genome_index:
                    genome_index[genome] = []

                genome_index[genome].append((
                    position,
                    position + node_size,
                    occurrence_id,
                ))

                node_positions[occurrence_id][genome] = position

        # ====================================================
        # 1. Build interval clusters independently per genome.
        # ====================================================

        genome_clusters = []

        for genome, entries in genome_index.items():

            entries.sort(key=lambda item: item[0])

            current_component = set()
            component_max_end = None

            for start, end, occurrence_id in entries:

                if component_max_end is None:

                    current_component.add(occurrence_id)
                    component_max_end = end

                    continue

                distance = start - component_max_end

                if distance <= max_distance:

                    current_component.add(occurrence_id)
                    if end > component_max_end:
                        component_max_end = end

                else:

                    genome_clusters.append(current_component)

                    current_component = {occurrence_id}
                    component_max_end = end

            if current_component:
                genome_clusters.append(current_component)

        # ====================================================
        # 2. Merge clusters sharing physical occurrences.
        # ====================================================

        merged_clusters = merge_clusters_by_overlap(
            genome_clusters,
            min_overlap=min_overlap,
        )

        # ====================================================
        # 3. Build region dictionaries.
        # ====================================================

        for cluster_set in merged_clusters:

            # Preserve original occurrence order.
            cluster_ids = [
                occurrence_id
                for occurrence_id in occurrence_ids
                if occurrence_id in cluster_set
            ]

            if len(cluster_ids) < min_nodes:
                continue

            cluster_nodes = [
                nodes_by_id[occurrence_id]
                for occurrence_id in cluster_ids
            ]

            # ------------------------------------------------
            # Retrieve unique Sequence IDs.
            # ------------------------------------------------

            sequence_ids = list(dict.fromkeys(
                nodes_by_id[occurrence_id]["sequence_id"]
                for occurrence_id in cluster_ids
            ))

            # ------------------------------------------------
            # Calculate genome-specific position ranges.
            # ------------------------------------------------

            genome_positions = {}

            for occurrence_id in cluster_ids:

                for genome, position in node_positions[
                    occurrence_id
                ].items():

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

            # ------------------------------------------------
            # Calculate region score.
            # ------------------------------------------------

            region_sum_score = sum(
                selected_seeds[sequence_id].synteny
                for sequence_id in sequence_ids
                if sequence_id in selected_seeds
            )

            # ------------------------------------------------
            # Calculate global cluster synteny.
            # ------------------------------------------------

            cluster_seeds = [
                selected_seeds[sequence_id]
                for sequence_id in sequence_ids
                if sequence_id in selected_seeds
            ]
            total_query_kmer_count = len(query.sequence) - query.k + 1
            region_score = calculate_cluster_synteny(
                seeds=cluster_seeds,
                total_query_kmer_count = total_query_kmer_count,
                max_kmer_distance=max_kmer_distance,
            )

            #Filter low synteny score
            if region_score < min_synteny_score:
                continue

            # ------------------------------------------------
            # Calculate region position boundaries.
            # ------------------------------------------------

            positions_mean = [
                nodes_by_id[occurrence_id].get("position_mean")
                for occurrence_id in cluster_ids
                if nodes_by_id[occurrence_id].get(
                    "position_mean"
                ) is not None
            ]

            region_position_min = (
                min(positions_mean)
                if positions_mean
                else None
            )

            region_position_max = (
                max(positions_mean)
                if positions_mean
                else None
            )

            # ------------------------------------------------
            # Calculate total physical size.
            # ------------------------------------------------

            # total_size = sum(
            #     int(node.get("size", 0))
            #     for node in cluster_nodes
            # )

            genome_sizes = [
                positions["max"] - positions["min"]
                for positions in genome_positions.values()
            ]

            total_size = max(genome_sizes) if genome_sizes else 0

            # ------------------------------------------------
            # Calculate region k-mer metrics.
            # ------------------------------------------------

            sequence_kmers = []

            for sequence_id in sequence_ids:

                seed = selected_seeds.get(sequence_id)

                if seed is None:
                    continue

                if seed.matched_kmer_codes.size > 0:

                    sequence_kmers.append(
                        np.unique(
                            np.asarray(
                                seed.matched_kmer_codes,
                                dtype=np.uint32,
                            )
                        )
                    )

            if sequence_kmers:

                all_region_kmers = np.concatenate(sequence_kmers)

                unique_region_kmers, kmer_sequence_counts = (
                    np.unique(
                        all_region_kmers,
                        return_counts=True,
                    )
                )

                n_kmers = len(unique_region_kmers)
                n_shared_kmers = int(np.count_nonzero(kmer_sequence_counts > 1))

            else:

                unique_region_kmers = np.empty(0,dtype=np.uint32,)

                n_kmers = 0
                n_shared_kmers = 0

            # ------------------------------------------------
            # Calculate query k-mer coverage.
            # ------------------------------------------------

            n_query_kmers = query.n_unique_kmers

            query_kmer_coverage = (
                100.0 * n_kmers / n_query_kmers
                if n_query_kmers > 0
                else 0.0
            )

            # ------------------------------------------------
            # Build occurrence names.
            # ------------------------------------------------

            occurrence_names = [
                nodes_by_id[occurrence_id].get(
                    "name",
                    str(occurrence_id),
                )
                for occurrence_id in cluster_ids
            ]

            # ------------------------------------------------
            # Store region.
            # ------------------------------------------------

            region = {
                "chromosome": chromosome,

                "position_min": region_position_min,
                "position_max": region_position_max,

                # Unique DNA sequences.
                "n_seeds": len(sequence_ids),

                # Physical Neo4j occurrences.
                "n_nodes": len(cluster_ids),

                "total_size": int(total_size),
                "score": float(region_score),
                "sum_score": float(region_sum_score),
                "genome_positions": genome_positions,

                "nodes": occurrence_names,

                # Region k-mer information.
                "n_kmers": int(n_kmers),
                "query_kmer_coverage": float(query_kmer_coverage),
                "n_shared_kmers": int(n_shared_kmers),
                "kmer_codes": unique_region_kmers,
            }

            regions.append(region)

            chromosome_stats[chromosome]["regions"].append({
                "score": float(region_score),
                "n_nodes": len(cluster_ids),
                "total_size": int(total_size),
                "n_kmers": int(n_kmers),
                "query_kmer_coverage": float(
                    query_kmer_coverage
                ),
                "n_shared_kmers": int(n_shared_kmers),
                "position_min": region_position_min,
                "position_max": region_position_max,
            })

            chromosome_stats[chromosome]["n_regions"] += 1

    # --------------------------------------------------------
    # Sort regions by score, number of seeds and total size.
    # --------------------------------------------------------

    regions.sort(
        key=lambda region: (
            region["score"],
            region["n_seeds"],
            region["total_size"],
        ),
        reverse=True,
    )

    # --------------------------------------------------------
    # Print chromosome statistics.
    # --------------------------------------------------------

    logger.debug("\n========== CHROMOSOME REGION STATISTICS ==========")

    for chromosome, stats in chromosome_stats.items():

        logger.debug(f"\nChromosome: {chromosome}")
        logger.debug(f"Number of regions: {stats['n_regions']}")

        sorted_stats = sorted(
            stats["regions"],
            key=lambda region: region["score"],
            reverse=True,
        )

        for index, region in enumerate(
            sorted_stats,
            start=1,
        ):

            logger.debug(
                f"  Region {index}: "
                f"score={region['score']:.3f}, "
                f"nodes={region['n_nodes']}, "
                f"n_kmers={region['n_kmers']}, "
                f"total size={region['total_size']}, "
                f"query_kmer_coverage="
                f"{region['query_kmer_coverage']:.2f}%, "
                f"shared_kmers={region['n_shared_kmers']}, "
                f"start={region['position_min']}, "
                f"end={region['position_max']}"
            )

    logger.debug("\n=================================================")

    return regions



"""
Main entry point for query mapping and genomic region detection.

Identifies candidate DNA sequences matching the query, evaluates
their coverage and synteny, selects the best candidate seeds, and
clusters their physical Neo4j occurrences into genomic regions
based on their genomic positions.

The clustering process groups nearby occurrences independently
for each genome and merges clusters sharing physical occurrences
across genomes.

Parameters
----------
seq : str
    Input DNA query sequence.

k : int, default=13
    K-mer length used for sequence matching.

canonical : bool, default=True
    Whether to use canonical k-mers, treating reverse-complement
    k-mers as equivalent.

min_coverage : float, default=20
    Minimum query and node k-mer coverage percentage required
    for candidate selection.

min_node_size : int or None, default=None
    Minimum physical node size required for candidate selection.
    If None, no additional minimum node size is imposed.

max_selected_seeds : int, default=10000
    Maximum number of candidate seeds retained for further
    analysis.

max_distance : int, default (None) = 10% of query length
    Maximum genomic distance, in base pairs, between consecutive
    intervals for them to belong to the same genomic cluster.

min_nodes : int, default=1
    Minimum number of distinct physical Neo4j occurrences
    required for a genomic region to be retained.
    
min_synteny_score : int, default=0.1 (min 0, max 1)
    Minimum of synteny score
    required for a genomic region to be retained.
    
max_kmer_distance : int, default=20
    Maximum k-mer distance allowed when evaluating synteny
    between candidate seeds.

max_query_occurrences : int, default=50
    Maximum number of occurrences of a query k-mer considered
    during candidate detection.

print_stats : bool, default=True
    Whether to display timing information and chromosome-level
    region statistics.

Returns
-------
regions : list of dict
    Detected genomic regions, sorted by decreasing synteny score.
    Each dictionary contains:

    - chromosome : Chromosome identifier.
    - position_min : Minimum genomic position across all genomes.
    - position_max : Maximum genomic position across all genomes.
    - n_seeds : Number of unique DNA sequences represented.
    - n_nodes : Number of distinct physical Neo4j occurrences.
    - total_size : Maximum genomic span (in base pairs) among
      the genomes represented in the region.
    - score : Global cluster synteny score.
    - sum_score : Sum of individual seed synteny scores.
    - genome_positions : Minimum and maximum genomic positions
      for each genome.
    - nodes : Names of the physical Neo4j occurrences.
    - n_kmers : Number of unique k-mers represented in the region.
    - query_kmer_coverage : Percentage of unique query k-mers
      represented in the region.
    - n_shared_kmers : Number of k-mers shared by multiple
      candidate sequences.
    - kmer_codes : Unique k-mer codes represented in the region.

selected_seeds : dict
    Selected Seed objects indexed by Sequence node ID.
    These seeds are used for genomic clustering and region
    characterization.

all_seeds : np.ndarray
    All candidate seeds passing the initial filtering steps
    performed during query mapping.

query : Query or None
    Prepared query object containing the query k-mer information.
    None if no valid query could be constructed.

Notes
-----
If no valid query or selected seeds are found, the function
returns an empty region list together with the available
candidate information.

The genomic region size (total_size) is defined as the maximum
difference between the minimum and maximum genomic positions
among all genomes represented in the region.
"""
def find_regions(
    seq: str,
    k: int = 13,
    canonical: bool = True,
    min_coverage: float = 20,
    min_node_size: int = None,
    max_selected_seeds: int = 10000,
    max_distance: int = None,
    min_nodes: int = 1,
    min_synteny_score=0.1,
    max_kmer_distance: int = 20,
    max_query_occurrences: int = 50,
    print_stats: bool = True,
):


    t_total = time.perf_counter()

    # --------------------------------------------------------
    # Detect candidate seeds.
    # --------------------------------------------------------

    selected_seeds, all_seeds, query = get_seeds(
        seq=seq,
        k=k,
        canonical=canonical,
        min_node_coverage=min_coverage,
        min_query_coverage=min_coverage,
        min_node_size=min_node_size,
        max_selected_seeds=max_selected_seeds,
        max_kmer_distance=max_kmer_distance,
        max_query_occurrences=max_query_occurrences,
        timings=print_stats,
    )

    if query is None or not selected_seeds:
        logger.info("No selected seeds found.")
        return [], selected_seeds, all_seeds, query

    # --------------------------------------------------------
    # Cluster physical occurrences into genomic regions.
    # --------------------------------------------------------

    t0 = time.perf_counter()

    regions = cluster_selected_seeds(
        selected_seeds=selected_seeds,
        query=query,
        max_distance=max_distance,
        min_nodes=min_nodes,
        min_synteny_score=min_synteny_score,
        max_kmer_distance=max_kmer_distance,
        print_chr_stats=print_stats
    )

    t_clustering = time.perf_counter() - t0

    logger.info(
        "%d regions found from %d selected seeds in %.3f seconds.",
        len(regions),
        len(selected_seeds),
        time.perf_counter() - t_total,
    )

    if print_stats:
        logger.info(f"Region clustering: {t_clustering:.3f} s")
        logger.info(
            f"Total mapping:     "
            f"{time.perf_counter() - t_total:.3f} s"
        )

    return regions, selected_seeds, all_seeds



