"""
Created on Tue Aug 19 16:05:16 2026

@author: fgraziani

This file contains the function to access the kmer index
"""

import json
from pathlib import Path
import numpy as np
from utils.kmer import encode_kmers
import matplotlib.pyplot as plt
import logging

logger = logging.getLogger("panabyss_logger")


# Configuration
DEFAULT_INDEX_PATH = "data/indexes"
INDEX_CATALOG_FILENAME = "kmer.indexes.json"
MAX_KMER_SIZE = 15
DEFAULT_KMER_SIZE = 13


# Global index cache
KMER_INDEX = None
POSTINGS_INDEX = None
DISCARDED_INDEX = None
NODE_SIZE_INDEX = None

# Metadata of the currently loaded index.
INDEX_K = None
INDEX_CANONICAL = None

# Must be identical to the writer.
INDEX_DTYPE = np.dtype([
    ("offset", "<u8"),
    ("count", "<u2"),
])

NODE_SIZE_DTYPE = np.dtype([
    ("node_id", "<u8"),
    ("size", "<u8"),
])


"""
Unload an index (used after index creation)
"""
def unload_index():

    global KMER_INDEX
    global POSTINGS_INDEX
    global DISCARDED_INDEX
    global INDEX_K
    global INDEX_CANONICAL
    global NODE_SIZE_INDEX

    KMER_INDEX = None
    POSTINGS_INDEX = None
    DISCARDED_INDEX = None
    NODE_SIZE_INDEX = None

    INDEX_K = None
    INDEX_CANONICAL = None


"""
Load the PanAbyss k-mer index catalog.

If the catalog does not exist yet, return an empty catalog.
"""
def load_index_catalog(output_dir=DEFAULT_INDEX_PATH):
    catalog_path = Path(output_dir) / INDEX_CATALOG_FILENAME

    if not catalog_path.exists():
        return {
            "application": "PanAbyss",
            "indexes": [],
        }

    with open(
            catalog_path,
            "r",
            encoding="utf-8",
    ) as f:

        catalog = json.load(f)

    if catalog.get("application") != "PanAbyss":
        raise ValueError(
            f"Invalid index catalog: "
            f"expected application='PanAbyss', "
            f"got {catalog.get('application')!r}"
        )

    if "indexes" not in catalog:
        catalog["indexes"] = []

    return catalog


"""
Find the index metadata corresponding to
the requested k-mer size and canonical mode.
"""
def find_index_metadata(
        k,
        canonical,
        output_dir=DEFAULT_INDEX_PATH,
):
    catalog = load_index_catalog(output_dir)

    for entry in catalog["indexes"]:

        if (
                entry["kmer_size"] == k
                and entry["canonical"] == canonical
        ):
            return entry

    raise FileNotFoundError(
        f"No k-mer index found for "
        f"k={k}, canonical={canonical}"
    )


"""
Load the k-mer index corresponding to (k, canonical).

The loaded resources are stored in global variables:

    KMER_INDEX
    POSTINGS_INDEX
    DISCARDED_INDEX

If the requested index is already loaded, nothing is done.
"""

def load_kmer_index(
        k,
        canonical,
        output_dir=DEFAULT_INDEX_PATH,
):
    global KMER_INDEX
    global POSTINGS_INDEX
    global DISCARDED_INDEX
    global INDEX_K
    global INDEX_CANONICAL

    # --------------------------------------------------------
    # Already loaded
    # --------------------------------------------------------

    if (
            KMER_INDEX is not None
            and INDEX_K == k
            and INDEX_CANONICAL == canonical
    ):
        return

    # --------------------------------------------------------
    # Find the files in the catalog.
    # --------------------------------------------------------

    metadata = find_index_metadata(
        k=k,
        canonical=canonical,
        output_dir=output_dir,
    )

    output_dir = Path(output_dir)

    index_path = output_dir / metadata["kmer_index_filename"]

    postings_path = output_dir / metadata["kmer_postings_filename"]

    discarded_path = output_dir / metadata["kmer_discarded_filename"]

    # --------------------------------------------------------
    # Check files.
    # --------------------------------------------------------

    for path in (
            index_path,
            postings_path,
            discarded_path,
    ):
        if not path.exists():
            raise FileNotFoundError(
                f"Index file not found: {path}"
            )

    # Number of possible k-mers.
    num_kmers = 4 ** k

    # Load main index.
    #
    # The index is memory-mapped and therefore does not
    # require loading the complete file into RAM.
    KMER_INDEX = np.memmap(
        index_path,
        dtype=INDEX_DTYPE,
        mode="r",
        shape=(num_kmers,),
    )

    # Load discarded k-mers.
    DISCARDED_INDEX = np.fromfile(
        discarded_path,
        dtype=np.uint32,
    )

    # Open postings file.
    # # access with seek().
    # POSTINGS_INDEX = open(
    #     postings_path,
    #     "rb",
    # )
    POSTINGS_INDEX = np.memmap(
        postings_path,
        dtype=np.uint64,
        mode="r",
    )

    INDEX_K = k
    INDEX_CANONICAL = canonical



def load_node_size_index(
        k: int,
        output_dir=DEFAULT_INDEX_PATH,
):
    global NODE_SIZE_INDEX

    if NODE_SIZE_INDEX is not None and INDEX_K == k:
        return

    output_dir = Path(output_dir)

    catalog = load_index_catalog(output_dir)

    # Find the first index entry matching k.
    matching_index = next(
        (
            entry
            for entry in catalog.get("indexes", [])
            if entry.get("kmer_size") == k
            and entry.get("node_size_index_filename")
        ),
        None,
    )

    if matching_index is None:
        raise FileNotFoundError(
            f"No node size index registered for k={k} "
            f"in {output_dir / INDEX_CATALOG_FILENAME}"
        )

    filename = matching_index["node_size_index_filename"]
    path = output_dir / filename

    if not path.exists():
        raise FileNotFoundError(
            f"Node size index registered in catalog but "
            f"file does not exist: {path}"
        )

    NODE_SIZE_INDEX = np.memmap(
        path,
        dtype=NODE_SIZE_DTYPE,
        mode="r",
    )



"""
Return the sequence size associated with a node ID.

Returns
-------
int
    Sequence size.

None
    If the node is not present in the index.
"""
def get_node_size(k: int, node_id: int):
    load_node_size_index()
    if NODE_SIZE_INDEX.size == 0:
        return None

    # Binary search on node_id.
    node_ids = NODE_SIZE_INDEX["node_id"]

    position = np.searchsorted(
        node_ids,
        np.uint64(node_id),
    )

    if position >= NODE_SIZE_INDEX.size:
        return None

    if node_ids[position] != node_id:
        return None

    return int(
        NODE_SIZE_INDEX["size"][position]
    )

"""
Look up a k-mer in the currently selected PanAbyss index.

Returns:
    (TAG, list[int])

    TAG = "ERROR" / "FOUND" / "NOT_FOUND" / "DISCARDED"

    If TAG = "FOUND" then the list is not empty and
        contains the node IDs associated with the k-mer
        else list is empty.
"""
def get_kmer_nodes(
        kmer,
        canonical=True,
        output_dir=DEFAULT_INDEX_PATH,
):
    global KMER_INDEX
    global POSTINGS_INDEX
    global DISCARDED_INDEX

    TAG = "NOT_FOUND"
    nodes_list = []

    # Validate k-mer size.
    k = len(kmer)

    if not 1 <= k <= MAX_KMER_SIZE:
        raise ValueError(
            f"k must be between 1 and "
            f"{MAX_KMER_SIZE}, got {k}"
        )

    # Load index if necessary.
    load_kmer_index(
        k=k,
        canonical=canonical,
        output_dir=output_dir,
    )

    # Encode k-mer (same encoder than for index creation)
    encoded, orientations = encode_kmers(
        kmer.encode("ascii"),
        k=k,
        canonical=canonical,
    )

    if encoded.size != 1:
        TAG = "ERROR"
        return TAG, []

    kmer_code = int(
        encoded[0]
    )

    # Read index entry.
    entry = KMER_INDEX[kmer_code]

    offset = int(
        entry["offset"]
    )

    count = int(
        entry["count"]
    )

    # --------------------------------------------------------
    # count == 0 means either:
    #
    #   - k-mer does not exist
    #   - k-mer was discarded
    #
    # We therefore check the discarded index.
    # --------------------------------------------------------

    if count == 0:
        TAG = "NOT_FOUND"
        position = np.searchsorted(
            DISCARDED_INDEX,
            np.uint32(kmer_code),
        )

        if (
                position < DISCARDED_INDEX.size
                and DISCARDED_INDEX[position]
                == kmer_code
        ):
            TAG = "DISCARDED"

        return TAG, []

    TAG = "FOUND"
    # Read postings.

    nodes = POSTINGS_INDEX[
        offset: offset + count
    ]

    return TAG, nodes.tolist()



"""
Look up multiple encoded k-mers in the currently selected
PanAbyss index.

Parameters
----------
kmer_codes : np.ndarray
    One-dimensional array of encoded k-mers with dtype
    uint32.

k : int
    K-mer size.

canonical : bool
    Whether canonical k-mer encoding is used.

output_dir : str or Path
    Directory containing the k-mer index.

Returns
-------
list[tuple[str, list[int]]]
    One result per input k-mer.

    Each result is:

        ("FOUND", [node_ids])
        ("NOT_FOUND", [])
        ("DISCARDED", [])
"""
def get_encoded_kmer_nodes_batch(
        kmer_codes,
        k,
        canonical=True,
        output_dir=DEFAULT_INDEX_PATH,
):


    global KMER_INDEX
    global POSTINGS_INDEX
    global DISCARDED_INDEX

    # --------------------------------------------------------
    # Validate k.
    # --------------------------------------------------------

    if not 1 <= k <= MAX_KMER_SIZE:
        raise ValueError(
            f"k must be between 1 and "
            f"{MAX_KMER_SIZE}, got {k}"
        )

    # --------------------------------------------------------
    # Convert input to uint32.
    #
    # np.asarray() does not copy if kmer_codes already has
    # the correct dtype.
    # --------------------------------------------------------

    kmer_codes = np.asarray(
        kmer_codes,
        dtype=np.uint32,
    )

    if kmer_codes.ndim != 1:
        raise ValueError(
            "kmer_codes must be one-dimensional."
        )

    if kmer_codes.size == 0:
        return []

    # --------------------------------------------------------
    # Load the corresponding index.
    # --------------------------------------------------------

    load_kmer_index(
        k=k,
        canonical=canonical,
        output_dir=output_dir,
    )

    # --------------------------------------------------------
    # Batch lookup in the direct index.
    #
    # This is the main optimization:
    #
    #     KMER_INDEX[kmer_codes]
    #
    # performs all index lookups using NumPy advanced
    # indexing instead of one Python lookup per k-mer.
    # --------------------------------------------------------

    entries = KMER_INDEX[
        kmer_codes
    ]

    counts = entries["count"]
    offsets = entries["offset"]

    # --------------------------------------------------------
    # Prepare results.
    # --------------------------------------------------------

    results = [
        None
        for _ in range(kmer_codes.size)
    ]

    # --------------------------------------------------------
    # Process FOUND k-mers.
    # --------------------------------------------------------

    found_positions = np.flatnonzero(
        counts > 0
    )

    for position in found_positions:

        position = int(position)

        offset = int(
            offsets[position]
        )

        count = int(
            counts[position]
        )

        # ----------------------------------------------------
        # Read postings.
        # ----------------------------------------------------


        nodes = POSTINGS_INDEX[
            offset: offset + count
        ]

        results[position] = (
            "FOUND",
            nodes.tolist(),
        )

    # --------------------------------------------------------
    # Process k-mers with count == 0.
    #
    # A zero count means either:
    #
    #   - k-mer was not found
    #   - k-mer was discarded because it exceeded
    #     max_postings during index construction
    # --------------------------------------------------------

    zero_positions = np.flatnonzero(
        counts == 0
    )

    if zero_positions.size > 0:

        zero_codes = kmer_codes[
            zero_positions
        ]

        # ----------------------------------------------------
        # DISCARDED_INDEX is sorted, so use binary search.
        # ----------------------------------------------------

        discarded_positions = np.searchsorted(
            DISCARDED_INDEX,
            zero_codes,
        )

        valid_positions = (
            discarded_positions
            < DISCARDED_INDEX.size
        )

        # ----------------------------------------------------
        # Default all zero-count k-mers to NOT_FOUND.
        # ----------------------------------------------------

        for position in zero_positions:

            results[int(position)] = (
                "NOT_FOUND",
                [],
            )

        # ----------------------------------------------------
        # Replace NOT_FOUND with DISCARDED where appropriate.
        # ----------------------------------------------------

        if np.any(valid_positions):

            candidate_indices = np.flatnonzero(
                valid_positions
            )

            candidate_discarded_positions = (
                discarded_positions[
                    valid_positions
                ]
            )

            candidate_codes = zero_codes[
                valid_positions
            ]

            matches = (
                DISCARDED_INDEX[
                    candidate_discarded_positions
                ]
                == candidate_codes
            )

            for candidate_index in np.flatnonzero(
                    matches
            ):

                original_position = int(
                    zero_positions[
                        candidate_indices[
                            candidate_index
                        ]
                    ]
                )

                results[original_position] = (
                    "DISCARDED",
                    [],
                )

    return results

"""
Look up multiple encoded k-mers and directly count
the number of distinct queried k-mers supporting
each node.
Parameters
----------
kmer_codes : np-array
    Kmer codes of the query sequence

k : int
    k-mer size

canonical : bool
    True to use canonical index, False otherwise.

output_dir : str
    Path to the index.
    
return_node_kmer_index : bool
    If true returns the k-mer found associated to each node, otherwise returns only nodes and kmer found count
    
Returns
-------
If return_node_kmer_index = True:
    Tuple
        node_kmer_counts -> (node_names, starts, ends, kmers)
        Tuple: (node names array, node starts array, node ends array, array of kmers)
        To find the kmers associated to a given node, for example node 42 in node names array
            => node_kmer_counts[3][starts[42]:ends[42]]
Else:   
    dict
        node_id -> number of supporting k-mers
"""
def get_encoded_kmer_node_counts_batch(
    kmer_codes,
    k,
    canonical=True,
    output_dir=DEFAULT_INDEX_PATH,
    return_node_kmer_index=False,
    timings=False
):
    import time

    global KMER_INDEX
    global POSTINGS_INDEX

    t_total = time.perf_counter()

    kmer_codes = np.asarray(kmer_codes, dtype=np.uint32)

    if kmer_codes.ndim != 1:
        raise ValueError("kmer_codes must be one-dimensional.")

    if kmer_codes.size == 0:
        if return_node_kmer_index:
            return (
                np.empty(0, dtype=np.uint64),
                np.empty(0, dtype=np.uint64),
                None,
            )
        return {}

    # ------------------------------------------------------------
    # Load index
    # ------------------------------------------------------------

    t0 = time.perf_counter()

    load_kmer_index(
        k=k,
        canonical=canonical,
        output_dir=output_dir,
    )

    t_load = time.perf_counter() - t0

    # ------------------------------------------------------------
    # Retrieve k-mer index entries
    # ------------------------------------------------------------

    t0 = time.perf_counter()

    entries = KMER_INDEX[kmer_codes]

    counts = entries["count"]
    offsets = entries["offset"]

    found_mask = counts > 0

    if not np.any(found_mask):
        print("[TIMING] Index loading:       {:.3f} s".format(t_load))
        print("[TIMING] No matching k-mers.")

        if return_node_kmer_index:
            return (
                np.empty(0, dtype=np.uint64),
                np.empty(0, dtype=np.uint64),
                None,
            )

        return {}

    offsets = offsets[found_mask].astype(np.uint64, copy=False)
    counts = counts[found_mask].astype(np.uint64, copy=False)
    found_kmers = kmer_codes[found_mask]

    total_postings = int(np.sum(counts))

    t_lookup = time.perf_counter() - t0

    # ------------------------------------------------------------
    # Allocate arrays
    # ------------------------------------------------------------

    t0 = time.perf_counter()

    all_nodes = np.empty(total_postings, dtype=np.uint64)

    if return_node_kmer_index:
        all_kmers = np.empty(total_postings, dtype=np.uint32)

    t_allocation = time.perf_counter() - t0

    # ------------------------------------------------------------
    # Read postings
    # ------------------------------------------------------------

    t0 = time.perf_counter()

    destination = 0

    order = np.argsort(offsets)

    offsets_sorted = offsets[order]
    counts_sorted = counts[order]
    kmers_sorted = found_kmers[order]

    for offset, count, kmer_code in zip(
        offsets_sorted,
        counts_sorted,
        kmers_sorted,
    ):
        offset = int(offset)
        count = int(count)

        nodes = POSTINGS_INDEX[offset:offset + count]

        all_nodes[destination:destination + count] = nodes

        if return_node_kmer_index:
            all_kmers[destination:destination + count] = kmer_code

        destination += count

    t_postings = time.perf_counter() - t0

    # ------------------------------------------------------------
    # Sort nodes and associated k-mers
    # ------------------------------------------------------------

    t0 = time.perf_counter()

    order = np.argsort(all_nodes)

    sorted_nodes = all_nodes[order]

    if return_node_kmer_index:
        sorted_kmers = all_kmers[order]

    t_sort = time.perf_counter() - t0

    # ------------------------------------------------------------
    # Group postings by node
    # ------------------------------------------------------------

    t0 = time.perf_counter()

    starts = np.concatenate((
        np.array([0], dtype=np.int64),
        np.flatnonzero(
            sorted_nodes[1:] != sorted_nodes[:-1]
        ) + 1,
    ))

    ends = np.concatenate((
        starts[1:],
        np.array([len(sorted_nodes)], dtype=np.int64),
    ))

    unique_nodes = sorted_nodes[starts]
    node_counts = ends - starts

    t_grouping = time.perf_counter() - t0

    # ------------------------------------------------------------
    # Build compact node-to-k-mer index
    # ------------------------------------------------------------

    node_kmer_index = None

    if return_node_kmer_index:
        node_kmer_index = (
            unique_nodes,
            starts,
            ends,
            sorted_kmers,
        )

    else:
        node_kmer_counts = dict(
            zip(
                unique_nodes.tolist(),
                node_counts.tolist(),
            )
        )

    # ------------------------------------------------------------
    # Print timings
    # ------------------------------------------------------------
    if timings:
        logger.debug("\n========== K-MER INDEX TIMINGS ==========")
        logger.debug(f"Input k-mers:          {len(kmer_codes):,}")
        logger.debug(f"Found k-mers:          {len(found_kmers):,}")
        logger.debug(f"Total postings:        {total_postings:,}")
        logger.debug(f"Unique nodes:          {len(unique_nodes):,}")
        logger.debug("-----------------------------------------")
        logger.debug(f"Index loading:         {t_load:.3f} s")
        logger.debug(f"Index lookup:          {t_lookup:.3f} s")
        logger.debug(f"Array allocation:      {t_allocation:.3f} s")
        logger.debug(f"Postings extraction:   {t_postings:.3f} s")
        logger.debug(f"Sorting:               {t_sort:.3f} s")
        logger.debug(f"Grouping:              {t_grouping:.3f} s")
        logger.debug("-----------------------------------------")
        logger.debug(
            f"Total function time:   "
            f"{time.perf_counter() - t_total:.3f} s"
        )
        logger.debug("=========================================\n")

    if return_node_kmer_index:
        return unique_nodes, node_counts, node_kmer_index

    return node_kmer_counts


"""
Plot the distribution of k-mer posting counts
from the currently loaded KMER_INDEX.

K-mer count ranges:
    1 - 5
    6 - 10
    11 - 50
    51 - 100
    101 - 200
    201 - 500
    501 - 1000
"""
def plot_kmer_count_histogram():

    ################################## Prepare data ######################
    if KMER_INDEX is None:
        raise RuntimeError(
            "K-mer index has not been initialized."
        )


    # Extract counts.
    counts = KMER_INDEX["count"]

    # Convert to int64 for safe comparisons.
    counts = np.asarray(
        counts,
        dtype=np.int64,
    )

    # 0 count will get all the discarded kmer => filter 0 count
    counts = counts[counts > 0]

    # Define histo ranges.
    ranges = [
        ("<= 5", counts <= 5),
        ("6 - 10", (counts > 5) & (counts <= 10)),
        ("11 - 50", (counts > 10) & (counts <= 50)),
        ("51 - 100", (counts > 50) & (counts <= 100)),
        ("101 - 200", (counts > 100) & (counts <= 200)),
        ("201 - 500", (counts > 200) & (counts <= 500)),
        ("501 - 1000", (counts > 500) & (counts <= 1000)),
    ]

    labels = []
    values = []

    for label, mask in ranges:
        labels.append(label)
        values.append(
            int(np.count_nonzero(mask))
        )

    ################################## PLOT ######################

    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(
        labels,
        values,
        color="steelblue",
    )

    ax.set_xlabel("Number of postings per k-mer")

    ax.set_ylabel("Number of k-mers")

    ax.set_title(
        f"K-mer posting count distribution "
        f"(k={INDEX_K}, "
        f"canonical={INDEX_CANONICAL})"
    )

    ax.grid(axis="y", alpha=0.3)


    # Display values above bars.
    for bar, value in zip(
            bars,
            values,
    ):
        ax.text(
            bar.get_x()
            + bar.get_width() / 2,
            bar.get_height(),
            f"{value:,}",
            ha="center",
            va="bottom",
        )

    plt.tight_layout()
    plt.show()

"""
Return node sizes for an array of node IDs.

Parameters
----------
k : int
    kmer size
node_ids : array-like
    Node IDs to look up.

Returns
-------
np.ndarray
    Node sizes.
    Returns 0 for node IDs that are not present.
"""
def get_node_sizes_batch(k:int, node_ids):

    global NODE_SIZE_INDEX

    load_node_size_index(k=k)

    node_ids = np.asarray(
        node_ids,
        dtype=np.uint64,
    )

    if node_ids.size == 0:
        return np.empty(
            0,
            dtype=np.uint64,
        )

    indexed_ids = NODE_SIZE_INDEX["node_id"]

    positions = np.searchsorted(
        indexed_ids,
        node_ids,
    )

    sizes = np.zeros(
        node_ids.size,
        dtype=np.uint64,
    )

    valid = positions < indexed_ids.size

    if not np.any(valid):
        return sizes

    valid_indices = np.flatnonzero(valid)
    valid_positions = positions[valid]

    matches = (
        indexed_ids[valid_positions]
        == node_ids[valid]
    )

    sizes[
        valid_indices[matches]
    ] = NODE_SIZE_INDEX["size"][
        valid_positions[matches]
    ]

    return sizes