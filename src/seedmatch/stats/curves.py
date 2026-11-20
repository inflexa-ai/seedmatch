"""The PR and ROC curves of scores and 0/1 labels."""

from dataclasses import dataclass
from typing import Final, NamedTuple

import numpy as np
import numpy.typing as npt
import pandas as pd
from numpy.typing import ArrayLike

from seedmatch._checks import shown
from seedmatch.errors import InputError
from seedmatch.stats.types import as_sample, as_vector

_LABEL_VALUES: Final = "0, 1, False, or True"
_RECALL_STEPS: Final = 100
# The fuzz of `seq` in R: a quotient that rounds to just below an integer keeps its last row.
_STEP_FUZZ: Final = 1e-10


@dataclass(frozen=True, eq=False)
class PrCurve:
    """The result of `pr_curve`.

    `auc_integral` is the area of Keilwagen, Grosse, and Grau (2014). It loses fewer digits than
    PRROC if the negatives far outnumber the positives (`D-050`), but some error remains: it is
    3.2e-10 for 1 positive in 1,000,000 scores. `auc_davis_goadrich` is the area of Davis and
    Goadrich (2006). `curve` has the columns `recall`, `precision`, and `threshold`. A frame has no
    truth value, so the objects compare by identity.
    """

    auc_integral: float
    auc_davis_goadrich: float
    curve: pd.DataFrame
    n_positive: int
    n_negative: int


@dataclass(frozen=True, eq=False)
class RocCurve:
    """The result of `roc_curve`. `curve` has the columns `fpr`, `tpr`, and `threshold`."""

    auc: float
    curve: pd.DataFrame
    n_positive: int
    n_negative: int


class _Steps(NamedTuple):
    """The counts at each distinct score, from the highest score to the lowest score.

    `tp[k]` and `fp[k]` count the positives and the negatives with a score of at least
    `thresholds[k]`. The `_prev` arrays hold the counts of the step before, and 0 for the first
    step.
    """

    thresholds: npt.NDArray[np.float64]
    tp: npt.NDArray[np.float64]
    fp: npt.NDArray[np.float64]
    tp_prev: npt.NDArray[np.float64]
    fp_prev: npt.NDArray[np.float64]
    n_positive: int
    n_negative: int


def pr_curve(scores: ArrayLike, labels: ArrayLike) -> PrCurve:
    """Return the PR curve of the scores, as `pr.curve` of PRROC 1.4 does for hard labels.

    A label of 1 or True is a positive, and a higher score ranks a gene higher. The areas and the
    rows obey the rules of "PR and ROC curves" in `docs/design.md`. `curve` goes from the recall 1
    to the recall 0. It has a row for each distinct score. It has more rows between two points
    with a large step, at the recall step min(1, P / 100) / P for P positives. The `threshold` of a
    row is the score of its point. The first row and the last row copy the precision and the
    threshold of the row next to them.

    NaN, a label that is not 0, 1, False, or True, labels of another length than the scores, and
    labels without a positive raise `InputError`. An infinite score is valid. The time is
    O(n log n) for one sort of n scores. The memory is O(n), because the curve has at most
    2K + 2 + max(100, P) rows for K distinct scores.
    """
    scores_array, positive = _as_scored_labels(scores, labels, need_negative=False)
    steps = _steps(scores_array, positive)
    return PrCurve(
        auc_integral=_auc_integral(steps),
        auc_davis_goadrich=_auc_davis_goadrich(steps),
        curve=_pr_frame(steps),
        n_positive=steps.n_positive,
        n_negative=steps.n_negative,
    )


def roc_curve(scores: ArrayLike, labels: ArrayLike) -> RocCurve:
    """Return the ROC curve of the scores, as `roc.curve` of PRROC 1.4 does for hard labels.

    The inputs obey the rules of `pr_curve`, and the labels must also have a negative. `auc` is the
    area under the curve with the trapezoid rule. `curve` has a row for each distinct score, after
    the first row (0, 0). The `threshold` of a row is the score of the points that the next row
    adds. Thus the first row has the highest score, and the last row copies the lowest score. The
    time is O(n log n) for one sort of n scores, and the memory is O(n).
    """
    scores_array, positive = _as_scored_labels(scores, labels, need_negative=True)
    steps = _steps(scores_array, positive)
    fpr = np.concatenate(([0.0], steps.fp / steps.n_negative))
    tpr = np.concatenate(([0.0], steps.tp / steps.n_positive))
    area = (fpr[1:] - fpr[:-1]) * (tpr[1:] + tpr[:-1]) / 2
    curve = pd.DataFrame(
        {
            "fpr": fpr,
            "tpr": tpr,
            "threshold": np.concatenate((steps.thresholds, steps.thresholds[-1:])),
        }
    )
    return RocCurve(
        auc=float(area.sum()),
        curve=curve,
        n_positive=steps.n_positive,
        n_negative=steps.n_negative,
    )


def _as_scored_labels(
    scores: ArrayLike, labels: ArrayLike, *, need_negative: bool
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.bool_]]:
    sample = as_sample("scores", scores, finite=False)
    positive = _as_labels(labels)
    if positive.size != sample.size:
        message = (
            f"scores has {sample.size} values and labels has {positive.size}; "
            "expected the same number"
        )
        raise InputError(message)
    if not positive.any():
        message = "labels have no positive; expected at least one label 1 or True"
        raise InputError(message)
    if need_negative and positive.all():
        message = "labels have no negative; expected at least one label 0 or False"
        raise InputError(message)
    return sample, positive


def _as_labels(labels: ArrayLike) -> npt.NDArray[np.bool_]:
    array = as_vector("labels", labels, kinds="biuf", expected=_LABEL_VALUES)
    if array.dtype == np.bool_:
        return array.astype(np.bool_, copy=False)
    if array.dtype.kind in "iuf":
        is_positive: npt.NDArray[np.bool_] = array == 1
        valid = is_positive | (array == 0)
        if valid.all():
            return is_positive
        index = int(valid.argmin())
        raise InputError(_bad_label(index, array[index]))
    for index, value in enumerate(array.tolist()):
        if not _is_label(value):
            raise InputError(_bad_label(index, value))
    return np.array(array.tolist(), dtype=np.bool_)


def _is_label(value: object) -> bool:
    if isinstance(value, bool | np.bool_):
        return True
    return isinstance(value, int | float | np.integer | np.floating) and value in (0, 1)


def _bad_label(index: int, value: object) -> str:
    plain = value.item() if isinstance(value, np.generic) else value
    return f"labels[{index}] is {shown(plain)}; expected {_LABEL_VALUES}"


def _steps(scores: npt.NDArray[np.float64], positive: npt.NDArray[np.bool_]) -> _Steps:
    # The stable sort puts the last zero of the input last in its run, so the threshold 0 has the
    # sign that PRROC shows. A run of -0.0 and 0.0 is one score.
    order = np.argsort(-scores, kind="stable")
    descending = scores[order]
    positives_up_to = np.cumsum(positive[order], dtype=np.int64)
    # The counts of a score are the cumulative counts at the end of its run of equal scores.
    last = np.flatnonzero(np.append(descending[1:] != descending[:-1], True))
    tp = positives_up_to[last]
    fp = last + 1 - tp
    n_positive = int(positives_up_to[-1])
    return _Steps(
        thresholds=descending[last],
        tp=tp.astype(np.float64),
        fp=fp.astype(np.float64),
        tp_prev=np.concatenate(([0.0], tp[:-1])),
        fp_prev=np.concatenate(([0.0], fp[:-1])),
        n_positive=n_positive,
        n_negative=len(scores) - n_positive,
    )


def _slope(steps: _Steps) -> npt.NDArray[np.float64]:
    """Return h of each step, the negatives per positive. A step without a positive gets 1."""
    gained = steps.tp - steps.tp_prev
    return np.divide(steps.fp - steps.fp_prev, gained, out=np.ones_like(gained), where=gained > 0)


def _precision(steps: _Steps) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Return the precision of each step, and the precision of the step before it.

    The step before the first step has no precision, so its entry is NaN.
    """
    precision = steps.tp / (steps.tp + steps.fp)
    return precision, np.concatenate(([np.nan], precision[:-1]))


def _auc_integral(steps: _Steps) -> float:
    positives = float(steps.n_positive)
    h = _slope(steps)
    grew = steps.tp > steps.tp_prev
    a = np.where(grew, 1 + h, 1.0)
    b = np.where(grew, (steps.fp_prev - h * steps.tp_prev) / positives, 0.0)
    width = steps.tp / positives - steps.tp_prev / positives
    area = width / a
    curved = b != 0
    a_curved, b_curved, width_curved = a[curved], b[curved], width[curved]
    before = (steps.tp_prev[curved] + steps.fp_prev[curved]) / positives
    # The difference of two logarithms loses digits when b / a is large (D-050). The value of
    # a * r + b at the step before is (TP + FP) / P, which needs no subtraction.
    area[curved] = (
        width_curved - b_curved / a_curved * np.log1p(a_curved * width_curved / before)
    ) / a_curved
    return float(area.sum())


def _auc_davis_goadrich(steps: _Steps) -> float:
    positives = float(steps.n_positive)
    gained = steps.tp - steps.tp_prev
    precision, precision_prev = _precision(steps)
    interpolated = (gained > 1) & (precision != precision_prev)
    # The first step has no step before it, and PRROC uses its precision for the whole width.
    interpolated[0] = False
    trapezoid = gained / positives * (precision_prev + precision) / 2
    trapezoid[0] = gained[0] / positives * precision[0]
    chosen = np.flatnonzero(interpolated)
    return float(trapezoid[~interpolated].sum() + _integer_points_area(steps, chosen).sum())


def _integer_points_area(steps: _Steps, chosen: npt.NDArray[np.intp]) -> npt.NDArray[np.float64]:
    """Return the areas of the unit trapezoids at each integer TP of the chosen steps."""
    positives = float(steps.n_positive)
    sizes = (steps.tp[chosen] - steps.tp_prev[chosen]).astype(np.int64)
    step_of = np.repeat(chosen, sizes)
    x = (np.arange(sizes.sum()) - np.repeat(np.cumsum(sizes) - sizes, sizes)).astype(np.float64)
    h = _slope(steps)[step_of]
    tp_prev, fp_prev = steps.tp_prev[step_of], steps.fp_prev[step_of]

    def precision_at(offset: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return (tp_prev + offset) / (tp_prev + offset + fp_prev + h * offset)

    return 1 / positives * (precision_at(x + 1) + precision_at(x)) / 2


def _pr_frame(steps: _Steps) -> pd.DataFrame:
    positives = float(steps.n_positive)
    recall = steps.tp / positives
    precision, precision_prev = _precision(steps)
    step_size = min(1.0, positives / _RECALL_STEPS) / positives
    refined = ((steps.tp - steps.tp_prev) / positives > step_size) & (precision != precision_prev)
    refined[0] = False
    chosen = np.flatnonzero(refined)
    width = recall[chosen] - recall[chosen - 1]
    sizes = (width / step_size + _STEP_FUZZ).astype(np.int64) + 1
    step_of = np.repeat(chosen, sizes)
    j = np.arange(sizes.sum()) - np.repeat(np.cumsum(sizes) - sizes, sizes)
    x = np.minimum(j * step_size, np.repeat(width, sizes))
    inner_recall = recall[step_of - 1] + x
    inner_precision = inner_recall / (
        inner_recall + steps.fp_prev[step_of] / positives + _slope(steps)[step_of] * x
    )
    columns = (
        np.concatenate((recall, inner_recall)),
        np.concatenate((precision, inner_precision)),
        np.concatenate((steps.thresholds, steps.thresholds[step_of])),
    )
    # lexsort is stable, and its last key is the main key: recall high to low, then threshold low.
    order = np.lexsort((columns[2], -columns[0]))
    recall, precision, threshold = (column[order] for column in columns)
    return pd.DataFrame(
        {
            "recall": np.concatenate(([1.0], recall, [0.0])),
            "precision": np.concatenate(([precision[0]], precision, [precision[-1]])),
            "threshold": np.concatenate(([threshold[0]], threshold, [threshold[-1]])),
        }
    )
