"""
Created on Tue Aug 18 18:05:16 2026

@author: fgraziani

This file contains the function to index a pangenome from neo4j database
"""

import json

from neo4j.exceptions import Neo4jError
from database.driver.neo4j_driver import get_driver
import logging

logger = logging.getLogger("panabyss_logger")

from numba import njit

import time
import shutil
from pathlib import Path

import numpy as np
from utils.base_utils import get_current_version
from utils.kmer import encode_kmers
from database.services.neo4j_kmer_requests import load_index_catalog



##########################################
# Configuration
##########################################

MAX_KMER_SIZE = 15
MAX_POSTINGS = 1000
DEFAULT_KMER_SIZE = 13

DEFAULT_INDEX_PATH = "data/indexes"
INDEX_CATALOG_FILENAME = "kmer.indexes.json"

# Final index:
# offset : byte offset in kmer.postings
# count  : number of node IDs
INDEX_DTYPE = np.dtype([
    ("offset", "<u8"),
    ("count", "<u2"),
])

#Size index:
# node_id
# size
NODE_SIZE_DTYPE = np.dtype([
    ("node_id", "<u8"),
    ("size", "<u8"),
])

"""
Add or update a k-mer index entry in the PanAbyss
index catalog.

An existing entry with the same (k, canonical)
combination is replaced.
"""


def register_kmer_index(
        output_dir,
        k: int,
        canonical: bool,
        filenames
):
    output_dir = Path(output_dir)

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    catalog_path = Path(output_dir) / INDEX_CATALOG_FILENAME
    catalog = load_index_catalog(output_dir)
    entry = {
        "application_version":
            get_current_version(),
        "kmer_size":
            k,
        "canonical":
            bool(canonical),
        **filenames,
    }


    # Remove an existing entry for the same k/canonical combination.
    catalog["indexes"] = [
        existing
        for existing in catalog["indexes"]
        if not (
                existing.get("kmer_size") == k
                and existing.get("canonical") == bool(canonical)
        )
    ]


    # Add the new entry.
    catalog["indexes"].append(entry)


    # Sort catalog for readability.
    catalog["indexes"].sort(
        key=lambda x: (
            x["kmer_size"],
            not x["canonical"],
        )
    )

    # Write atomically.
    temporary_path = catalog_path.with_suffix(".tmp")

    with open(
            temporary_path,
            "w",
            encoding="utf-8",
    ) as f:
        json.dump(
            catalog,
            f,
            indent=4,
            ensure_ascii=False,
        )

        f.write("\n")

    temporary_path.replace(
        catalog_path
    )

    return catalog_path



"""
Stream Sequence nodes from Neo4j.
fetch_size is the buffer size to transfert the neo4j results

Expected properties:
    n.id       : integer
    n.sequence : DNA string
"""
def stream_sequences(fetch_size: int = 20000, k: int = DEFAULT_KMER_SIZE):
    driver = get_driver()

    if driver is None:
        return

    query = """
    MATCH (n:Sequence)
    WHERE id(n) IS NOT NULL
      AND n.sequence IS NOT NULL
      AND size(n.sequence) >= $k
    RETURN id(n) AS node_id, n.sequence AS sequence
    """

    with driver.session(
            default_access_mode="READ",
            fetch_size=fetch_size,
    ) as session:

        result = session.run(query, k=k)

        # IMPORTANT:
        # iterate directly over result instead of list(result)
        # so Neo4j streams records instead of materializing everything.
        for record in result:
            # yield record["node_id"], record["sequence"]

            sequence = record["sequence"]

            if isinstance(sequence, str):
                sequence = sequence.encode("ascii")

            yield record["node_id"], sequence


"""
Build the direct k-mer index using range-based temporary partitions.

PASS 1:
    - Read each sequence exactly once.
    - Encode canonical k-mers.
    - Remove duplicate k-mers within each node.
    - Route sorted k-mers directly to fixed k-mer ranges.
    - Write (kmer, node_id) pairs to temporary binary files.

PASS 2: (No second scan of Neo4j nodes)
    - Load one partition at a time.
    - Sort pairs by (kmer, node_id).
    - Build the final index and postings.

Partitions use fixed k-mer ranges. They are intentionally not
balanced by occurrence count because balancing would require an
additional global counting/repartitioning pass.

Final files:
    kmer.index
    kmer.postings
    kmer.discarded
"""


def build_chunked_direct_index(
        output_dir=DEFAULT_INDEX_PATH,
        k=DEFAULT_KMER_SIZE,
        num_partitions=32,
        max_postings=MAX_POSTINGS,
        fetch_size=20000,
        pair_buffer_size=2_000_000,
        canonical=True
):
    # =========================================================
    # Validation
    # =========================================================

    if not 1 <= k <= MAX_KMER_SIZE:
        raise ValueError(
            f"k must be between 1 and {MAX_KMER_SIZE}"
        )

    if num_partitions < 1:
        raise ValueError(
            "num_partitions must be >= 1"
        )

    if max_postings > 65535:
        raise ValueError(
            "max_postings must be <= 65535"
        )

    if pair_buffer_size <= 0:
        raise ValueError(
            "pair_buffer_size must be > 0"
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )
    node_size_index_path = output_dir / "nodes.size.index"

    build_node_size_index = not node_size_index_path.exists()

    num_kmers = 4 ** k


    # Paths
    temp_dir = output_dir / f".kmer_partitions_k{k}"

    if temp_dir.exists():
        shutil.rmtree(temp_dir)

    temp_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    mode = (
        "canonical"
        if canonical
        else "noncanonical"
    )

    index_path = output_dir / f"kmer.{mode}.{k}.index"
    postings_path = output_dir / f"kmer.{mode}.{k}.postings"
    discarded_path = output_dir / f"kmer.{mode}.{k}.discarded"

    filenames = {"kmer_index_filename": f"kmer.{mode}.{k}.index",
                 "kmer_discarded_filename": f"kmer.{mode}.{k}.discarded",
                 "kmer_postings_filename": f"kmer.{mode}.{k}.postings"
                 }

    # =========================================================
    # Pair representation
    #
    # kmer    : uint32
    # node_id : uint64
    #
    # 12 bytes per pair.
    # =========================================================

    pair_dtype = np.dtype([
        ("kmer", "<u4"),
        ("node_id", "<u8"),
    ])

    # =========================================================
    # PASS 1
    # =========================================================

    logger.info("=== PASS 1 : Encoding and writing chunks ===")
    logger.info(f"K-mers          : {num_kmers:,}")
    logger.info(f"Partitions      : {num_partitions}")

    logger.info(f"Pair buffer     : {pair_buffer_size:,}")

    # ---------------------------------------------------------
    # Create partition files.
    #
    # Partition p contains the k-mer range:
    #
    #     [p * num_kmers / P,
    #      (p + 1) * num_kmers / P)
    #
    # Because kmers returned by np.unique() are sorted, each
    # sequence can be split into contiguous slices.
    # ---------------------------------------------------------

    partition_paths = [
        temp_dir / f"partition_{p:03d}.bin"
        for p in range(num_partitions)
    ]

    partition_files = [
        open(
            path,
            "wb",
            buffering=16 * 1024 * 1024,
        )
        for path in partition_paths
    ]

    node_size_file = None
    if build_node_size_index:
        logger.info(f"Node size index does not exist: {node_size_index_path}")

        node_size_file = open(
            node_size_index_path,
            "wb",
            buffering=16 * 1024 * 1024,
        )


    # ---------------------------------------------------------
    # Reusable pair buffers.
    # ---------------------------------------------------------

    partition_buffers = [
        np.empty(
            pair_buffer_size,
            dtype=pair_dtype,
        )
        for _ in range(num_partitions)
    ]

    partition_sizes = np.zeros(
        num_partitions,
        dtype=np.int64,
    )

    # =========================================================
    # PASS 1 loop
    # =========================================================

    pass1_start = time.perf_counter()

    sequence_count = 0
    pair_count = 0

    # ---------------------------------------------------------
    # Precompute partition boundaries in k-mer space.
    #
    # This avoids performing the multiplication/division
    # for every k-mer.
    # ---------------------------------------------------------

    partition_boundaries = (
            np.arange(
                num_partitions + 1,
                dtype=np.uint64,
            )
            * np.uint64(num_kmers)
            // np.uint64(num_partitions)
    )

    try:

        for node_id, sequence in stream_sequences(
                fetch_size=fetch_size,
                k=k,
        ):

            sequence_count += 1

            # -------------------------------------------------
            # Build node size index if necessary.
            # -------------------------------------------------
            if build_node_size_index:
                node_size_record = np.array(
                    (
                        np.uint64(node_id),
                        np.uint64(len(sequence)),
                    ),
                    dtype=NODE_SIZE_DTYPE,
                )

                node_size_file.write(
                    node_size_record.tobytes()
                )

            # -------------------------------------------------
            # Encode canonical k-mers.
            # -------------------------------------------------

            kmers = encode_kmers(
                sequence,
                k,
                canonical=canonical
            )

            if kmers.size == 0:
                continue

            # -------------------------------------------------
            # Remove duplicate k-mers from this node.
            #
            # np.unique() returns sorted values.
            # -------------------------------------------------

            kmers = np.unique(kmers)

            n_kmers = kmers.size

            if n_kmers == 0:
                continue

            pair_count += n_kmers

            # -------------------------------------------------
            # Find the first k-mer belonging to each partition.
            #
            # Since kmers is sorted, searchsorted gives us
            # contiguous slices directly.
            #
            # This avoids:
            #
            #   - partition_ids
            #   - flatnonzero
            #   - concatenate
            #   - boolean masks
            # -------------------------------------------------

            boundaries = np.searchsorted(
                kmers,
                partition_boundaries[1:-1],
                side="left",
            )

            previous = 0

            # -------------------------------------------------
            # Write each non-empty contiguous slice.
            # -------------------------------------------------

            for p in range(num_partitions):

                if p < num_partitions - 1:
                    end = int(boundaries[p])
                else:
                    end = n_kmers

                if end <= previous:
                    continue

                selected_kmers = kmers[
                    previous:end
                ]

                n = end - previous

                current_size = int(
                    partition_sizes[p]
                )

                # -------------------------------------------------
                # Flush the current buffer if necessary.
                # -------------------------------------------------

                if current_size + n > pair_buffer_size:

                    if current_size > 0:
                        partition_files[p].write(
                            partition_buffers[p][
                                :current_size
                            ].tobytes()
                        )

                        current_size = 0
                        partition_sizes[p] = 0

                # -------------------------------------------------
                # Very large sequence segment.
                # -------------------------------------------------

                if n >= pair_buffer_size:

                    direct_pairs = np.empty(
                        n,
                        dtype=pair_dtype,
                    )

                    direct_pairs["kmer"] = selected_kmers
                    direct_pairs["node_id"] = np.uint64(
                        node_id
                    )

                    partition_files[p].write(
                        direct_pairs.tobytes()
                    )

                    del direct_pairs

                else:

                    # -------------------------------------------------
                    # Copy k-mers and node ID into the reusable buffer.
                    # -------------------------------------------------

                    buffer = partition_buffers[p]

                    buffer[
                        current_size:
                        current_size + n
                    ]["kmer"] = selected_kmers

                    buffer[
                        current_size:
                        current_size + n
                    ]["node_id"] = np.uint64(
                        node_id
                    )

                    partition_sizes[p] = (
                            current_size + n
                    )

                previous = end

            # -------------------------------------------------
            # Progress.
            # -------------------------------------------------

            if sequence_count % 500000 == 0:
                elapsed = (
                        time.perf_counter()
                        - pass1_start
                )

                logger.info(
                    f"[PASS1 "
                    f"{sequence_count:,}] "
                    f"elapsed={elapsed:.2f}s | "
                    f"pairs={pair_count:,}"
                )

    finally:
        # -----------------------------------------------------
        # Flush all remaining buffers.
        # -----------------------------------------------------
        for p in range(num_partitions):

            size = int(
                partition_sizes[p]
            )

            if size > 0:
                partition_files[p].write(
                    partition_buffers[p][
                        :size
                    ].tobytes()
                )

            partition_files[p].close()

        # -----------------------------------------------------
        # Finalize node size index.
        # -----------------------------------------------------
        if node_size_file is not None:
            node_size_file.close()

    pass1_time = (
            time.perf_counter()
            - pass1_start
    )

    del partition_buffers
    del partition_sizes
    del partition_boundaries

    logger.info("PASS 1 complete.")

    logger.info(f"Sequences : {sequence_count:,}")
    logger.info(f"Pairs     : {pair_count:,}")
    logger.info(f"Time      : {pass1_time:.2f}s")

    # =========================================================
    # Partition statistics
    # =========================================================

    logger.info("=== Partition sizes ===")

    partition_pair_counts = np.zeros(
        num_partitions,
        dtype=np.uint64,
    )

    for p, path in enumerate(
            partition_paths
    ):
        size = path.stat().st_size

        count = (
                size // pair_dtype.itemsize
        )

        partition_pair_counts[p] = count

        logger.info(
            f"Partition {p:02d}: "
            f"{count:,} pairs | "
            f"{size / 1024 ** 3:.2f} GB"
        )

    del partition_pair_counts

    # =========================================================
    # PASS 2
    # =========================================================

    logger.info("=== PASS 2 : Sorting chunks and building index ===")

    pass2_start = time.perf_counter()

    # ---------------------------------------------------------
    # Final index.
    # ---------------------------------------------------------

    index = np.memmap(
        index_path,
        dtype=INDEX_DTYPE,
        mode="w+",
        shape=(num_kmers,),
    )

    index["offset"] = 0
    index["count"] = 0

    # ---------------------------------------------------------
    # Discarded k-mers.
    # ---------------------------------------------------------

    discarded = np.zeros(
        num_kmers,
        dtype=np.bool_,
    )

    # ---------------------------------------------------------
    # Postings output buffer.
    # ---------------------------------------------------------

    postings_block_size = 8_000_000

    postings_buffer = np.empty(
        postings_block_size,
        dtype=np.uint64,
    )

    postings_buffer_size = 0

    postings_file = open(
        postings_path,
        "wb",
        buffering=16 * 1024 * 1024,
    )

    total_postings = 0
    discarded_count = 0

    def flush_postings_buffer():
        nonlocal postings_buffer_size

        if postings_buffer_size == 0:
            return

        postings_file.write(
            postings_buffer[
                :postings_buffer_size
            ].tobytes()
        )

        postings_buffer_size = 0

    # =========================================================
    # Process partitions
    # =========================================================

    try:

        for p, path in enumerate(
                partition_paths
        ):

            partition_start_time = (
                time.perf_counter()
            )

            logger.info(f"Processing partition {p + 1}/{num_partitions}...")

            file_size = path.stat().st_size

            if file_size == 0:
                logger.info("  Empty")

                continue

            # -------------------------------------------------
            # Load complete partition.
            # -------------------------------------------------

            t0 = time.perf_counter()

            pairs = np.fromfile(
                path,
                dtype=pair_dtype,
            )

            load_time = (
                    time.perf_counter()
                    - t0
            )

            logger.info(f"  Loaded {len(pairs):,} pairs in {load_time:.2f}s")

            if pairs.size == 0:
                del pairs
                continue

            # =================================================
            # Sort by (kmer, node_id)
            # =================================================

            t0 = time.perf_counter()

            order = np.lexsort(
                (
                    pairs["node_id"],
                    pairs["kmer"],
                )
            )
            #Change to these lines to sort only on kmer key
            # order = np.argsort(
            #     pairs["kmer"],
            #     kind="quicksort",
            # )

            pairs = pairs[order]

            del order

            sort_time = (
                    time.perf_counter()
                    - t0
            )

            logger.info(f"  Sorted in {sort_time:.2f}s")

            # =================================================
            # Process sorted pairs
            # =================================================

            t0 = time.perf_counter()

            kmers = pairs["kmer"]
            nodes = pairs["node_id"]

            n_pairs = pairs.size

            # -------------------------------------------------
            # Locate k-mer group boundaries.
            # -------------------------------------------------

            if n_pairs == 1:

                group_starts = np.array(
                    [0],
                    dtype=np.int64,
                )

            else:

                group_starts = np.flatnonzero(
                    kmers[1:] != kmers[:-1]
                ).astype(
                    np.int64,
                    copy=False,
                ) + 1

                group_starts = np.concatenate(
                    (
                        np.array(
                            [0],
                            dtype=np.int64,
                        ),
                        group_starts,
                    )
                )

            n_groups = group_starts.size

            # -------------------------------------------------
            # Process each k-mer group.
            # -------------------------------------------------

            for group_index in range(
                    n_groups
            ):

                start = int(
                    group_starts[group_index]
                )

                if group_index + 1 < n_groups:
                    end = int(
                        group_starts[
                            group_index + 1
                            ]
                    )
                else:
                    end = n_pairs

                kmer = int(
                    kmers[start]
                )

                count = end - start

                # -------------------------------------------------
                # Discard overrepresented k-mers.
                # -------------------------------------------------

                if count > max_postings:
                    discarded[kmer] = True
                    discarded_count += 1

                    continue

                # -------------------------------------------------
                # Record index entry.
                # -------------------------------------------------

                index[kmer]["offset"] = (
                    total_postings
                )

                index[kmer]["count"] = (
                    count
                )

                # -------------------------------------------------
                # Copy nodes to the postings output buffer.
                # -------------------------------------------------

                source_pos = start
                remaining = count

                while remaining > 0:

                    available = (
                            postings_block_size
                            - postings_buffer_size
                    )

                    take = min(
                        remaining,
                        available,
                    )

                    postings_buffer[
                        postings_buffer_size:
                        postings_buffer_size + take
                    ] = nodes[
                        source_pos:
                        source_pos + take
                    ]

                    postings_buffer_size += take
                    source_pos += take
                    remaining -= take
                    total_postings += take

                    if (
                            postings_buffer_size
                            == postings_block_size
                    ):
                        flush_postings_buffer()

            process_time = (
                    time.perf_counter()
                    - t0
            )

            logger.info(f"  Processed in {process_time:.2f}s")

            del pairs
            del kmers
            del nodes
            del group_starts

            logger.info(f"  Partition total: {time.perf_counter() - partition_start_time:.2f}s")

    finally:

        flush_postings_buffer()
        postings_file.close()

    del postings_buffer

    # =========================================================
    # Finalize index
    # =========================================================

    index.flush()
    del index

    # =========================================================
    # Write discarded k-mers
    # =========================================================

    logger.info("Writing discarded k-mers...")

    t0 = time.perf_counter()

    discarded_kmers = np.flatnonzero(
        discarded
    ).astype(
        np.uint32,
        copy=False,
    )

    discarded_kmers.tofile(
        discarded_path
    )

    discarded_write_time = (
            time.perf_counter()
            - t0
    )

    discarded_count = discarded_kmers.size

    del discarded
    del discarded_kmers

    # =========================================================
    # Cleanup
    # =========================================================

    shutil.rmtree(
        temp_dir,
        ignore_errors=True,
    )

    # =========================================================
    # Timing
    # =========================================================

    pass2_time = (
            time.perf_counter()
            - pass2_start
    )

    total_time = (
            time.perf_counter()
            - pass1_start
    )

    # =========================================================
    # Final report
    # =========================================================

    logger.info("=== Direct index complete ===")
    logger.info(f"Sequences          : {sequence_count:,}")
    logger.info(f"Pairs              : {pair_count:,}")
    logger.info(f"Retained postings  : {total_postings:,}")
    logger.info(f"Discarded k-mers   : {discarded_count:,}")

    logger.info("=== Timing ===")

    logger.info(f"PASS 1             : {pass1_time:.2f}s")
    logger.info(f"PASS 2             : {pass2_time:.2f}s")
    logger.info(f"Discarded write    : {discarded_write_time:.2f}s")

    logger.info(f"TOTAL              : {total_time:.2f}s")
    logger.info(f"Index size         : {index_path.stat().st_size / 1024 ** 3:.2f} GB")
    logger.info(f"Postings size      : {postings_path.stat().st_size / 1024 ** 3:.2f} GB")

    logger.info(f"Discarded size     : {discarded_path.stat().st_size / 1024 ** 2:.2f} MB")

    # =========================================================
    # Update index catalog
    # =========================================================

    register_kmer_index(
        output_dir=output_dir,
        k=k,
        canonical=canonical,
        filenames=filenames
    )

    return (
        index_path,
        postings_path,
        discarded_path,
    )

