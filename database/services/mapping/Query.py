"""
Created on 12/11/2026

@author: fgraziani
"""
import numpy as np

from database.services.mapping.index_kmer_requests import encode_kmers

"""
Represents the query sequence and prepares the indexes
required for synteny calculations with different Seeds.

K-mers are encoded only once during initialization.
The sorted indexes allow efficient retrieval of k-mer
positions and orientations.
"""
class Query:
    def __init__(self, sequence, k, canonical=True):

        if not isinstance(sequence, str):
            raise TypeError("Query sequence must be a string.")

        if not isinstance(k, int) or k <= 0:
            raise ValueError("k must be a positive integer.")

        self.sequence = sequence
        self.k = k
        self.canonical = canonical

        # Encode query k-mers.
        self.kmer_codes, self.kmer_orientations = encode_kmers(
            sequence,
            k=k,
            canonical=canonical,
            orientations=True,
        )

        # Convert outputs to NumPy arrays.
        self.kmer_codes = np.asarray(
            self.kmer_codes,
            dtype=np.uint32,
        )

        self.kmer_orientations = np.asarray(
            self.kmer_orientations,
            dtype=np.int8,
        )

        # Generate the original positions of query k-mers.
        self.positions = np.arange(
            len(self.kmer_codes),
            dtype=np.int32,
        )

        # Convert outputs to NumPy arrays.
        self.kmer_codes = np.asarray(
            self.kmer_codes,
            dtype=np.uint64,
        )

        self.kmer_orientations = np.asarray(
            self.kmer_orientations,
            dtype=np.int8,
        )

        self.positions = np.asarray(
            self.positions,
            dtype=np.int32,
        )

        # ----------------------------------------------------
        # Prepare sorted k-mer indexes.
        # ----------------------------------------------------

        if len(self.kmer_codes) > 0:

            order = np.argsort(
                self.kmer_codes,
                kind="stable",
            )

            self.sorted_kmer_codes = np.ascontiguousarray(
                self.kmer_codes[order],
                dtype=np.uint64,
            )

            self.sorted_positions = np.ascontiguousarray(
                self.positions[order],
                dtype=np.int32,
            )

            self.sorted_orientations = np.ascontiguousarray(
                self.kmer_orientations[order],
                dtype=np.int8,
            )

            self.unique_kmer_codes = np.unique(
                self.sorted_kmer_codes
            )

        else:

            self.sorted_kmer_codes = np.empty(
                0,
                dtype=np.uint64,
            )

            self.sorted_positions = np.empty(
                0,
                dtype=np.int32,
            )

            self.sorted_orientations = np.empty(
                0,
                dtype=np.int8,
            )

            self.unique_kmer_codes = np.empty(
                0,
                dtype=np.uint64,
            )

    # ========================================================
    # Properties
    # ========================================================

    @property
    def has_kmers(self):
        """Return whether the query contains any k-mers."""

        return len(self.kmer_codes) > 0

    @property
    def n_kmers(self):
        """Return the total number of k-mers in the query."""

        return len(self.kmer_codes)

    @property
    def n_unique_kmers(self):
        """Return the number of distinct k-mers in the query."""

        return len(self.unique_kmer_codes)

    @property
    def synteny_index(self):
        """
        Return the arrays required for synteny calculations.

        K-mers are sorted by code, with their corresponding
        positions and orientations.
        """

        return (
            self.sorted_kmer_codes,
            self.sorted_positions,
            self.sorted_orientations,
        )
