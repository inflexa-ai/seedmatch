import numpy as np
from numpy.typing import NDArray

def hamming_counts(
    codes: NDArray[np.uint8],
    offsets: NDArray[np.int64],
    pattern: NDArray[np.uint8],
    compat: NDArray[np.bool_],
    max_mm: int,
    min_mm: int,
    /,
) -> NDArray[np.int64]: ...
def hamming_ranges(
    codes: NDArray[np.uint8],
    offsets: NDArray[np.int64],
    pattern: NDArray[np.uint8],
    compat: NDArray[np.bool_],
    max_mm: int,
    min_mm: int,
    /,
) -> tuple[NDArray[np.int64], NDArray[np.int64]]: ...
def indel_counts(
    codes: NDArray[np.uint8],
    offsets: NDArray[np.int64],
    pattern: NDArray[np.uint8],
    compat: NDArray[np.bool_],
    max_mm: int,
    /,
) -> NDArray[np.int64]: ...
def indel_ranges(
    codes: NDArray[np.uint8],
    offsets: NDArray[np.int64],
    pattern: NDArray[np.uint8],
    compat: NDArray[np.bool_],
    max_mm: int,
    /,
) -> tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.int64]]: ...
