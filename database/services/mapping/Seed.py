"""
Created on 12/11/2026

@author: fgraziani
"""
import numpy as np
from numba import njit

from database.services.mapping.index_kmer_requests import encode_kmers


"""
Build anchor correspondences between query elements and node elements.

For each node, the function searches for matching query elements based
on their codes. Each valid match generates an anchor linking a query
position to a node index, together with their relative orientation.

Only matches whose number of occurrences in the query does not exceed
`max_query_occurrences` are retained. This filtering helps avoid
ambiguous matches caused by codes that appear too frequently in the
query.

The function uses a two-pass strategy:
    1. Count the exact number of valid anchors to allocate output arrays
       with the required size.
    2. Populate the preallocated arrays with query positions, node
       indices, and relative orientations.

The implementation is optimized for performance using Numba's JIT
compilation (`@njit`) and avoids dynamic array resizing.

Parameters
----------
query_codes : numpy.ndarray
    Sorted array of codes associated with query elements.
    Multiple occurrences of the same code are allowed.

query_positions : numpy.ndarray
    Array containing the original positions or indices of query
    elements. Its length must match `query_codes`.

query_orientations : numpy.ndarray
    Array containing the orientation values associated with query
    elements. Used only when `has_query_orientations` is True.
    Its length must match `query_codes` when orientations are enabled.

has_query_orientations : bool
    Indicates whether query orientation information is available.

    If True, the relative orientation is computed as the product of
    the query orientation and the corresponding node orientation.

    If False, all generated anchors are assigned a default orientation
    of +1.

node_codes : numpy.ndarray
    Array containing the codes associated with node elements.
    Each element is compared against `query_codes` to identify matches.

node_orientations : numpy.ndarray
    Array containing the orientation values associated with node
    elements. Its length must match `node_codes`.

max_query_occurrences : int
    Maximum allowed number of occurrences of a given code in the query.
    Codes occurring more frequently than this threshold are ignored
    when generating anchors.

Returns
-------
anchor_query : numpy.ndarray
    One-dimensional array of dtype int64 containing the query positions
    associated with each generated anchor.

anchor_node : numpy.ndarray
    One-dimensional array of dtype int64 containing the indices of the
    corresponding nodes.

anchor_orientation : numpy.ndarray
    One-dimensional array of dtype int8 containing the relative
    orientation of each anchor.

    Possible values:
        +1 : Same orientation.
        -1 : Opposite orientation.

Notes
-----
- The three returned arrays always have the same length.
- Each anchor represents a correspondence between one query element
  and one node element sharing the same code.
- A node may generate multiple anchors if its code appears multiple
  times in the query.
- Nodes with no matching query code are ignored.
- Nodes whose matching code occurs more than `max_query_occurrences`
  times in the query are ignored.
- The function assumes that `query_codes` is sorted in ascending order,
  as required by `numpy.searchsorted`.
- If no valid anchors are found, three empty arrays with the expected
  data types are returned.

Complexity
----------
Let:
    N = number of nodes
    Q = number of query elements
    A = number of generated anchors

Time complexity:
    O(N * log(Q) + A)

The algorithm performs two binary searches per node during each of
the two passes, resulting in O(N * log(Q) + A) overall complexity.

Space complexity:
    O(A)

The output arrays are preallocated using the exact number of valid
anchors, avoiding unnecessary memory allocations or resizing.

Implementation Details
----------------------
The relative orientation is calculated as:

    relative_orientation = query_orientation * node_orientation

The resulting value is normalized to either +1 or -1:

    orientation = -1 if relative_orientation < 0 else +1

When query orientations are unavailable, the relative orientation
defaults to +1.
"""
@njit(cache=True)
def build_anchors_numba(
    query_codes,
    query_positions,
    query_orientations,
    has_query_orientations,
    node_codes,
    node_orientations,
    max_query_occurrences,
):
    n_node = node_codes.size
    n_anchors = 0

    # First pass: calculate the exact number of anchors.
    for i in range(n_node):
        code = node_codes[i]
        left = np.searchsorted(query_codes, code, side="left")
        right = np.searchsorted(query_codes, code, side="right")
        occurrences = right - left

        if occurrences == 0 or occurrences > max_query_occurrences:
            continue

        n_anchors += occurrences

    if n_anchors == 0:
        return (
            np.empty(0, dtype=np.int64),
            np.empty(0, dtype=np.int64),
            np.empty(0, dtype=np.int8),
        )

    anchor_query = np.empty(n_anchors, dtype=np.int64)
    anchor_node = np.empty(n_anchors, dtype=np.int64)
    anchor_orientation = np.empty(n_anchors, dtype=np.int8)

    # Second pass: populate the anchor arrays.
    index = 0

    for i in range(n_node):
        code = node_codes[i]
        left = np.searchsorted(query_codes, code, side="left")
        right = np.searchsorted(query_codes, code, side="right")
        occurrences = right - left

        if occurrences == 0 or occurrences > max_query_occurrences:
            continue

        node_orientation = node_orientations[i]

        for j in range(left, right):
            if has_query_orientations:
                relative_orientation = (
                    query_orientations[j] * node_orientation
                )
            else:
                relative_orientation = 1

            orientation = -1 if relative_orientation < 0 else 1

            anchor_query[index] = query_positions[j]
            anchor_node[index] = i
            anchor_orientation[index] = orientation

            index += 1

    return anchor_query, anchor_node, anchor_orientation


"""
Chaining algorithm on a node
Calculate a chain consistency score from a set of anchor correspondences.

The function identifies coherent chains of anchors by grouping them
according to their diagonal positions and orientations, then computing
the Longest Increasing Subsequence (LIS) within each group.

Anchors are grouped dynamically
based on the distance between consecutive diagonal values.

Two consecutive anchors belong to the same cluster if the absolute
difference between their diagonal values is less than or equal to
`max_kmer_distance`.

The function uses a two-stage strategy:

    1. Diagonal clustering:
       Anchors are separated by orientation, sorted by their actual
       diagonal values, and grouped according to the maximum allowed
       distance between consecutive diagonals.

    2. Chain detection:
       Within each cluster, anchors are sorted by node position and
       the Longest Increasing Subsequence (LIS) of their query
       positions is computed.

Reverse-oriented anchors are handled by negating their query positions
before LIS computation, allowing both orientations to use the same
increasing subsequence algorithm.

The final score is calculated as the square root of the sum of squared
LIS lengths divided by the total number of input anchors.

Parameters
----------
anchor_query : numpy.ndarray
    One-dimensional array of dtype int64 containing the query positions
    associated with each anchor.

anchor_node : numpy.ndarray
    One-dimensional array of dtype int64 containing the node positions
    associated with each anchor.

anchor_orientation : numpy.ndarray
    One-dimensional array of dtype int8 containing the relative
    orientation of each anchor.

    Possible values:
        +1 : Forward orientation.
        -1 : Reverse orientation.

max_kmer_distance : int
    Maximum allowed difference between the diagonal values of two
    consecutive anchors belonging to the same cluster.

Returns
-------
score : float
    Chain consistency score calculated as:

        score = sqrt(sum(LIS_length ** 2) / n_anchors)

    where:
        - LIS_length is the length of the longest increasing
          subsequence identified within a cluster.
        - n_anchors is the total number of input anchors.

    Only clusters containing at least two anchors are considered
    when accumulating the score.

    Returns 0.0 if the input contains no anchors.

Notes
-----
- All input arrays must have the same length.
- Anchor orientations are expected to be either +1 or -1.
- Forward and reverse orientations are processed independently.
- Forward-oriented anchors use the following diagonal:

      diagonal = node_position - query_position

- Reverse-oriented anchors use the following diagonal:

      diagonal = node_position + query_position

Algorithm
---------
The function consists of four main steps:

1. Orientation separation
   Anchors are processed independently for forward (+1) and reverse
   (-1) orientations.

2. Diagonal calculation and sorting
   Each anchor is assigned a diagonal value according to its
   orientation. Anchors are then sorted by their actual diagonal
   values.

3. Dynamic diagonal clustering
   Consecutive anchors are grouped when the difference between their
   diagonal values is less than or equal to `max_kmer_distance`.

   A new cluster is created whenever this condition is not satisfied.

4. Longest Increasing Subsequence and score aggregation
   Within each cluster, anchors are sorted by node position.
   The LIS of their corresponding query positions is computed using
   a binary search approach.

"""
@njit(cache=True)
def calculate_chains_numba(
    anchor_query,
    anchor_node,
    anchor_orientation,
    max_kmer_distance,
):
    n = anchor_query.size

    if n == 0:
        return 0.0, np.empty(0, dtype=np.int64)

    diagonals = np.empty(n, dtype=np.int64)
    anchor_indices = np.empty(n, dtype=np.int64)

    total_squared_size = 0.0

    # Store the local synteny score of each anchor.
    anchor_scores = np.zeros(n, dtype=np.int64)

    # =========================================================
    # STEP 1: Process each orientation independently.
    # =========================================================
    for orientation_value in (1, -1):

        # Collect anchors and calculate their diagonal values.
        orientation_count = 0

        for anchor_index in range(n):

            if anchor_orientation[anchor_index] != orientation_value:
                continue

            query_position = anchor_query[anchor_index]
            node_position = anchor_node[anchor_index]

            if orientation_value == 1:
                diagonal = node_position - query_position
            else:
                diagonal = node_position + query_position

            diagonals[orientation_count] = diagonal
            anchor_indices[orientation_count] = anchor_index

            orientation_count += 1

        if orientation_count < 2:
            continue

        # Sort anchors by their diagonal values.
        sorted_positions = np.argsort(diagonals[:orientation_count])

        # =====================================================
        # STEP 2: Identify clusters and compute their LIS.
        # =====================================================
        cluster_start = 0

        # The last iteration ensures that the final cluster
        # is also processed.
        for i in range(1, orientation_count + 1):

            if (
                i == orientation_count
                or diagonals[sorted_positions[i]]
                - diagonals[sorted_positions[i - 1]]
                > max_kmer_distance
            ):

                cluster_end = i
                cluster_size = cluster_end - cluster_start

                # =================================================
                # STEP 3: Compute the LIS for the current cluster.
                # =================================================
                if cluster_size >= 2:

                    cluster_anchor_indices = anchor_indices[
                        sorted_positions[cluster_start:cluster_end]
                    ]

                    # Sort anchors by their node positions.
                    node_order = np.argsort(
                        anchor_node[cluster_anchor_indices]
                    )

                    tails = np.empty(cluster_size, dtype=np.int64)

                    # Store the anchor index associated with
                    # each LIS length.
                    tails_indices = np.empty(
                        cluster_size, dtype=np.int64
                    )

                    # Store the predecessor of each anchor.
                    predecessors = np.full(
                        cluster_size, -1, dtype=np.int64
                    )

                    lis_length = 0

                    # -------------------------------------------------
                    # Compute the LIS and store its predecessors.
                    # -------------------------------------------------
                    for j in range(cluster_size):

                        anchor_index = cluster_anchor_indices[
                            node_order[j]
                        ]

                        query_position = anchor_query[anchor_index]

                        # Reverse the query coordinate for reverse
                        # orientation.
                        if orientation_value == -1:
                            query_position = -query_position

                        # Binary search in the tails array.
                        search_start = 0
                        search_end = lis_length

                        while search_start < search_end:

                            search_middle = (
                                search_start + search_end
                            ) // 2

                            if tails[search_middle] < query_position:
                                search_start = search_middle + 1
                            else:
                                search_end = search_middle

                        # Store the predecessor of the current anchor.
                        if search_start > 0:
                            predecessors[j] = tails_indices[
                                search_start - 1
                            ]

                        # Update the LIS.
                        tails[search_start] = query_position
                        tails_indices[search_start] = j

                        if search_start == lis_length:
                            lis_length += 1

                    # -------------------------------------------------
                    # Calculate the global synteny score.
                    # -------------------------------------------------
                    if lis_length >= 2:

                        total_squared_size += (
                            lis_length * lis_length
                        )

                        # -------------------------------------------------
                        # Reconstruct the selected LIS.
                        # -------------------------------------------------
                        current_index = tails_indices[
                            lis_length - 1
                        ]

                        while current_index != -1:

                            anchor_index = cluster_anchor_indices[
                                node_order[current_index]
                            ]

                            # All anchors belonging to this LIS
                            # receive the same local score.
                            anchor_scores[anchor_index] = lis_length

                            current_index = predecessors[
                                current_index
                            ]

                # Start the next cluster.
                cluster_start = i

    # =========================================================
    # STEP 4: Calculate the final global chain consistency score.
    # =========================================================
    synteny = np.sqrt(total_squared_size / n)

    return synteny, anchor_scores


"""
Represents a unique candidate sequence (Neo4j Sequence node).

A sequence can have multiple physical occurrences in the graph.
These occurrences are stored in neo4j_nodes, each preserving
its own name and genomic location properties.
"""
class Seed:
    def __init__(
        self,
        sequence_id: int,
        count: int,
        size: int,
        coverage: float,
        query_coverage: float,
        sequence: str = None,
        query_coverage_count: int = None,
        kmer_codes=None,
    ):
        # Sequence node attributes.
        self.sequence_id = int(sequence_id)
        self.count = int(count)
        self.size = int(size)
        self.coverage = float(coverage)
        self.query_coverage = float(query_coverage)

        # Sequence data.
        self.sequence = sequence
        self.synteny = 0.0
        self.anchor_scores = None

        # K-mer codes and orientations.
        self.kmer_codes = (
            np.ascontiguousarray(kmer_codes, dtype=np.uint32)
            if kmer_codes is not None
            else np.empty(0, dtype=np.uint32)
        )
        self.matched_kmer_codes = self.kmer_codes.copy()

        self.kmer_orientations = np.empty(
            0,
            dtype=np.int8,
        )

        # Anchors used for cluster-level synteny.
        self.anchor_query = np.empty(0, dtype=np.int64)
        self.anchor_node = np.empty(0, dtype=np.int64)
        self.anchor_orientation = np.empty(0, dtype=np.int8)

        # Physical Neo4j Node occurrences associated with this sequence.
        self.neo4j_nodes = []
        self.node_names = []

        # Additional query coverage information.
        self.query_coverage_count = query_coverage_count

    # ========================================================
    # Properties
    # ========================================================

    @property
    def node_id(self):
        """
        Return the Sequence node identifier.

        This identifier is different from the name of a physical
        Node occurrence.
        """
        return self.sequence_id

    @property
    def name(self):
        """Return the names of all physical Node occurrences."""
        return self.node_names

    # ========================================================
    # Neo4j node management
    # ========================================================

    def add_neo4j_node(self, properties):
        """
        Add a physical Node occurrence while preserving all
        its Neo4j properties.
        """
        data = dict(properties)

        node_name = data.get("name")

        if node_name is not None and node_name not in self.node_names:
            self.node_names.append(node_name)

        self.neo4j_nodes.append(data)

    # ========================================================
    # Synteny calculation
    # ========================================================
    """
    Calculate the synteny score between this Seed and a Query.

    The calculation consists of three main steps:
    1. Encode the Seed sequence into k-mers.
    2. Build anchors between matching query and Seed k-mers.
    3. Calculate collinear chains from these anchors.

    Optional timing information can be collected for each step.
    """
    def calculate_synteny(
        self,
        query,
        max_kmer_distance: int = 20,
        max_query_occurrences: int = 50,
        timings: dict = None,
    ) -> float:

        import time

        if max_kmer_distance < 0:
            raise ValueError("max_kmer_distance must be >= 0")

        if max_query_occurrences < 1:
            raise ValueError("max_query_occurrences must be >= 1")

        if not query.has_kmers or self.sequence is None:
            self.synteny = 0.0
            return self.synteny

        # ----------------------------------------------------
        # Encode Seed k-mers.
        # ----------------------------------------------------

        t0 = time.perf_counter()

        node_codes, node_orientations = encode_kmers(
            self.sequence,
            k=query.k,
            canonical=query.canonical,
            orientations=True,
        )

        node_codes = np.ascontiguousarray(
            node_codes,
            dtype=np.uint32,
        )

        node_orientations = np.ascontiguousarray(
            node_orientations,
            dtype=np.int8,
        )

        # Store the k-mers and orientations associated with this Seed.
        self.kmer_codes = node_codes
        self.kmer_orientations = node_orientations

        if timings is not None:
            timings["encoding"] = timings.get("encoding", 0.0) + (
                time.perf_counter() - t0
            )

        if node_codes.size == 0:
            self.synteny = 0.0
            return self.synteny

        # ----------------------------------------------------
        # Build anchors.
        # ----------------------------------------------------

        t0 = time.perf_counter()

        anchors = build_anchors_numba(
            query.sorted_kmer_codes,
            query.sorted_positions,
            query.sorted_orientations,
            True,
            node_codes,
            node_orientations,
            max_query_occurrences,
        )

        if timings is not None:
            timings["anchors"] = timings.get("anchors", 0.0) + (
                time.perf_counter() - t0
            )

        anchor_query, anchor_node, anchor_orientation = anchors
        # Store anchors for subsequent cluster-level calculations.
        self.anchor_query = anchor_query
        self.anchor_node = anchor_node
        self.anchor_orientation = anchor_orientation

        if anchor_query.size == 0:
            self.synteny = 0.0
            return self.synteny

        # ----------------------------------------------------
        # Calculate collinear chains.
        # ----------------------------------------------------

        t0 = time.perf_counter()

        self.synteny, self.anchor_scores = calculate_chains_numba(
            self.anchor_query,
            self.anchor_node,
            self.anchor_orientation,
            max_kmer_distance
        )

        if timings is not None:
            timings["chains"] = timings.get("chains", 0.0) + (
                time.perf_counter() - t0
            )

        return self.synteny