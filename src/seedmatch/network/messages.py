"""The messages that the genes of a network send to the genes that they point to."""

from enum import StrEnum
from typing import assert_never

import numpy as np
from numpy.typing import NDArray


class Aggregator(StrEnum):
    """How the messages of the senders of a gene join, with the values of `agg` of R."""

    MEAN = "mean"
    DEGREE_IMPORTANCE = "degree_importance"


def sender_weights(aggregator: Aggregator, out_degree: NDArray[np.int64]) -> NDArray[np.float64]:
    """Return the weight of the message of each gene.

    The weight is 1 for `MEAN`, and 1 / `out_degree` for `DEGREE_IMPORTANCE`. A gene with the
    out-degree 0 sends no message, thus its weight is 0.
    """
    match aggregator:
        case Aggregator.MEAN:
            return np.ones(len(out_degree), dtype=np.float64)
        case Aggregator.DEGREE_IMPORTANCE:
            weights = np.zeros(len(out_degree), dtype=np.float64)
            np.divide(1.0, out_degree, out=weights, where=out_degree > 0)
            return weights
        case _:
            assert_never(aggregator)


def mean_messages(
    receivers: NDArray[np.intp],
    senders: NDArray[np.intp],
    values: NDArray[np.float64],
    weights: NDArray[np.float64],
) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
    """Return, for each gene g, the number of its senders and its score.

    Edge i goes from the gene `senders[i]` to the gene `receivers[i]`. Let S(g) be the senders of
    g, f the `values`, and w the `weights`. The score is the sum over S(g) of w(s) * f(s) / |S(g)|,
    plus f(g), all divided by 2. A gene without a sender has the score NaN. The time and the memory
    are O(E + G) for E edges and G genes.
    """
    genes = len(values)
    counts = np.bincount(receivers, minlength=genes).astype(np.int64, copy=False)
    totals = np.bincount(receivers, weights=weights[senders] * values[senders], minlength=genes)
    means = np.full(genes, np.nan)
    np.divide(totals, counts, out=means, where=counts > 0)
    return counts, (means + values) / 2
