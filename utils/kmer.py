from numba import njit
import numpy as np

"""
Encode a DNA sequence into uint32 k-mers.

If canonical=True:
    return min(forward, reverse-complement).

If canonical=False:
    return the k-mer in the original sequence orientation.

Invalid bases split the sequence into independent segments.
"""
@njit(cache=True)
def encode_sequence_numba(seq,k,canonical=True):

    n = len(seq)

    if n < k:
        return np.empty(0, dtype=np.uint32)

    result = np.empty(
        n - k + 1,
        dtype=np.uint32,
    )

    mask = np.uint64(
        (1 << (2 * (k - 1))) - 1
    )

    high_shift = 2 * (k - 1)

    forward = np.uint64(0)
    reverse = np.uint64(0)

    valid_length = 0
    result_count = 0

    for i in range(n):

        c = seq[i]

        # -----------------------------------------------------
        # Convert ASCII directly to 2-bit DNA code.
        #
        # A/a -> 0
        # C/c -> 1
        # G/g -> 2
        # T/t -> 3
        # -----------------------------------------------------

        if c == 65 or c == 97:
            base = np.uint64(0)

        elif c == 67 or c == 99:
            base = np.uint64(1)

        elif c == 71 or c == 103:
            base = np.uint64(2)

        elif c == 84 or c == 116:
            base = np.uint64(3)

        else:
            # Invalid base: restart rolling encoders.
            valid_length = 0
            forward = np.uint64(0)
            reverse = np.uint64(0)

            continue

        # -----------------------------------------------------
        # First k-mer of the current valid segment.
        # -----------------------------------------------------

        if valid_length < k:

            forward = (
                (forward << 2)
                | base
            )

            if canonical:

                complement = (
                    np.uint64(3) - base
                )

                reverse |= (
                    complement
                    << (2 * valid_length)
                )

            valid_length += 1

            if valid_length == k:

                if canonical:

                    if forward < reverse:
                        result[result_count] = (
                            np.uint32(forward)
                        )
                    else:
                        result[result_count] = (
                            np.uint32(reverse)
                        )

                else:

                    result[result_count] = (
                        np.uint32(forward)
                    )

                result_count += 1

            continue

        # -----------------------------------------------------
        # Subsequent k-mers.
        # -----------------------------------------------------

        forward = (
            ((forward & mask) << 2)
            | base
        )

        if canonical:

            complement = (
                np.uint64(3) - base
            )

            reverse = (
                (reverse >> 2)
                | (
                    complement
                    << high_shift
                )
            )

            if forward < reverse:
                result[result_count] = (
                    np.uint32(forward)
                )
            else:
                result[result_count] = (
                    np.uint32(reverse)
                )

        else:

            result[result_count] = (
                np.uint32(forward)
            )

        result_count += 1

    return result[:result_count]

"""
Fast Python wrapper around the Numba encoder.

canonical=True:
    Encode canonical k-mers using the minimum of the forward
    and reverse-complement representations.

canonical=False:
    Encode k-mers only in their original sequence orientation.
"""
def encode_kmer_numba(
    seq_bytes: bytes,
    k: int,
    canonical: bool = True,
):

    if len(seq_bytes) < k:
        return np.empty(
            0,
            dtype=np.uint32,
        )

    seq_array = np.frombuffer(
        seq_bytes,
        dtype=np.uint8,
    )

    return encode_sequence_numba(
        seq_array,
        k,
        canonical,
    )