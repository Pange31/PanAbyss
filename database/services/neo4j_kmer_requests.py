import json
from pathlib import Path
import numpy as np
from utils.kmer import encode_kmer_numba
import matplotlib.pyplot as plt

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

# Metadata of the currently loaded index.
INDEX_K = None
INDEX_CANONICAL = None

# Must be identical to the writer.
INDEX_DTYPE = np.dtype([
    ("offset", "<u8"),
    ("count", "<u2"),
])

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
    catalog = load_index_catalog(
        output_dir
    )

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
    # access with seek().
    POSTINGS_INDEX = open(
        postings_path,
        "rb",
    )

    INDEX_K = k
    INDEX_CANONICAL = canonical


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
    encoded = encode_kmer_numba(
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
    POSTINGS_INDEX.seek(
        offset * np.dtype(np.uint64).itemsize
    )

    nodes = np.fromfile(
        POSTINGS_INDEX,
        dtype=np.uint64,
        count=count,
    )

    return TAG, nodes.tolist()


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