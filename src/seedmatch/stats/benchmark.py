"""The benchmark of the seed counts against a gold standard of genes that go down."""

from dataclasses import dataclass

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from seedmatch._checks import (
    as_cutoff,
    as_data_frame,
    as_text,
    check_columns,
    check_gene_ids_are_text,
    check_numeric_columns,
    column_values,
)
from seedmatch.errors import InputError
from seedmatch.stats.curves import PrCurve, RocCurve, pr_curve, roc_curve


@dataclass(frozen=True)
class Benchmark:
    """The result of `benchmark_predictions`. The PR curve and the ROC curve use the same labels."""

    pr: PrCurve
    roc: RocCurve


def _is_score(values: NDArray[np.float64]) -> NDArray[np.bool_]:
    return ~np.isnan(values)


def benchmark_predictions(
    de: pd.DataFrame,
    *,
    column: str = "mer7m8",
    padj_cutoff: float = 0.15,
    l2fc_cutoff: float = 0.0,
) -> Benchmark:
    """Return the PR and ROC curves of the seed counts against a gold standard (`benchmark_preds`).

    `de` has the columns `gene_id` (text), `padj`, `log2FoldChange`, and `column`. The score of a
    row is its number in `column`, and a higher score ranks a row higher. The gold standard is the
    set of the gene IDs of the rows with `padj <= padj_cutoff` and `log2FoldChange <
    l2fc_cutoff` (`D-006`). A NaN fails each comparison, thus a row without a `padj` or without a
    `log2FoldChange` does not pass. A row is a positive if its gene is in the set, as `%in%` of R
    does. Thus each row of a gene in the set is a positive, and each other row is a negative. The
    missing gene IDs (None, NaN, and `pd.NA`) are one gene, as `NA` is in R. A row without a
    `gene_id` is a positive if a row without a `gene_id` passes. The PR curve and the ROC curve
    use these labels (`D-006`).

    A cutoff that is a bool, text, or NaN raises `InputError`. A NaN score raises it too, and so
    does a gold standard without a positive or without a negative row (`D-048`). An infinite score
    is valid. A count column that `annotate_seed_counts` makes has no NaN. With
    `shared_genes=False`, a gene without a sequence gets the count 0.

    The time is O(n log n) for n rows, because each curve sorts the scores one time. The memory is
    O(n). Measured on an Apple M2 Max: 1,000,000 rows with distinct scores take 0.33 s, and the
    peak memory above the start is 220 MB. 1,000,000 rows with 8 distinct scores take 0.15 s and
    43 MB.
    """
    table = as_data_frame(de)
    name = as_text("column", column)
    check_columns(table, ("gene_id", "padj", "log2FoldChange", name))
    check_gene_ids_are_text(table)
    check_numeric_columns(table, ("padj", "log2FoldChange", name))
    padj_limit = as_cutoff("padj_cutoff", padj_cutoff)
    l2fc_limit = as_cutoff("l2fc_cutoff", l2fc_cutoff)
    scores = column_values(table, name, _is_score, "a score for each row")
    padj = table["padj"].to_numpy(dtype=np.float64, na_value=np.nan)
    l2fc = table["log2FoldChange"].to_numpy(dtype=np.float64, na_value=np.nan)
    genes = table["gene_id"]
    passing = (padj <= padj_limit) & (l2fc < l2fc_limit)
    # R has one NA, but isin treats None, NaN, and pd.NA as different values.
    missing = genes.isna().to_numpy()
    positive = genes.isin(genes[passing]).to_numpy() | (missing & bool((missing & passing).any()))
    rule = f"padj <= {padj_limit} and log2FoldChange < {l2fc_limit}"
    if not positive.any():
        message = f"no gene has {rule}; expected at least one positive gene"
        raise InputError(message)
    if positive.all():
        message = f"each row is a positive with {rule}; expected at least one negative row"
        raise InputError(message)
    return Benchmark(pr=pr_curve(scores, positive), roc=roc_curve(scores, positive))
