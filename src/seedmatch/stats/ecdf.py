"""The ECDF test of two gene sets in a DE table, and the group of each row (`R/deseq_fc_ecdf.R`)."""

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from enum import StrEnum
from functools import partial
from typing import assert_never

import numpy as np
import numpy.typing as npt
import pandas as pd

from seedmatch._checks import (
    as_items,
    as_member,
    as_text,
    check_columns,
    check_gene_ids_are_text,
    check_numeric_columns,
    shown,
)
from seedmatch.errors import InputError, warn
from seedmatch.stats.dts import check_dts_range, dts_test
from seedmatch.stats.ks import ks_test
from seedmatch.stats.types import Alternative, PValueMethod, TwoSampleResult, as_alternative
from seedmatch.stats.wasserstein import (
    WassersteinAreas,
    check_wasserstein_range,
    check_wasserstein_sample,
    wasserstein_areas,
)
from seedmatch.stats.wilcoxon import wilcoxon_test


class EcdfTest(StrEnum):
    """The test of an ECDF comparison. The values are the `stats_test` names of R.

    R spells the name of the Wilcoxon test `Wilcoxen`.
    """

    KS = "KS"
    WILCOXON = "Wilcoxon"
    WASSERSTEIN = "Wass"
    DTS = "DTS"


@dataclass(frozen=True)
class EcdfTestResult:
    """The result of `ecdf_test`.

    `p_value` and `method` are `None` for `WASSERSTEIN`, which has no p-value. `alternative` is
    `None` for `WASSERSTEIN` and `DTS`. `areas` is set only for `WASSERSTEIN`.
    """

    test: EcdfTest
    statistic: float
    p_value: float | None
    method: PValueMethod | None
    alternative: Alternative | None
    n_target: int
    n_background: int
    areas: WassersteinAreas | None


# Each keyword-only argument is an independent option of `ecdf_stat_test` or of one of its tests.
def ecdf_test(  # noqa: PLR0913
    de: pd.DataFrame,
    *,
    target: Iterable[str],
    background: Iterable[str],
    test: EcdfTest = EcdfTest.KS,
    alternative: Alternative | None = None,
    rng: np.random.Generator | int | None = None,
    n_boot: int | None = None,
    n_interp: int | None = None,
) -> EcdfTestResult:
    """Compare the `log2FoldChange` of two gene sets, as `ecdf_stat_test` of R does.

    `de` needs the columns `gene_id` and `log2FoldChange`, and each of them must occur once
    (`D-051`). `target` and `background` are iterables of gene IDs. A `str` is one ID, not a set,
    thus it raises `InputError`. A sample holds the `log2FoldChange` of each row whose `gene_id` is
    in the set, in the order of `de`. A gene with more than one row gives each row. Each test drops
    the rows with a NaN `log2FoldChange` (`D-036`), and the result counts the rows that remain. A
    set without such a row raises `InputError`.

    `KS` and `WILCOXON` compare the target with the background. `alternative=None` gives `GREATER`,
    the default of R. But `GREATER` has another meaning for each test. For `KS`, it means that the
    target ECDF is above the background ECDF, thus the target values are lower. For `WILCOXON`, it
    means that the target values are higher. To find a down-shift of the target with `WILCOXON`,
    use `Alternative.LESS`. On the test data of the package, `KS` with `GREATER` gives p = 1.04e-4,
    but `WILCOXON` with `GREATER` gives p = 0.9996.

    `DTS` and `WASSERSTEIN` have no direction, and an `alternative` raises `InputError` (`D-035`).
    `DTS` needs `rng`, and `n_boot` is its number of permutations.
    `n_interp` is the number of grid points of `WASSERSTEIN`. For `n_boot` and `n_interp`, `None`
    gives the default of the test: 2000 for `DTS` and 10000 for `WASSERSTEIN`. A value that is not
    `None` for an option of another test raises `InputError`.

    `DTS` and `WASSERSTEIN` raise `InputError` for an infinite `log2FoldChange`. The message
    names the set and the gene ID of the first such row (`D-034` for `DTS`). They also raise
    `InputError` for a pooled range that a float cannot hold (`D-034` and `D-038`), and the message
    names `target` and `background`. `WASSERSTEIN` also needs at least 10 values and at least 2
    distinct values in each set, after the drop of the NaN rows. It takes the target as the seed
    sample and the background as the non-seed sample. Its `statistic` is `areas.negative_censored`,
    as in R.
    """
    chosen = as_member("test", test, EcdfTest)
    _check_options(chosen, alternative, rng, n_boot, n_interp)
    score = _scorer(chosen, alternative, rng, n_boot, n_interp)
    check_columns(de, ("gene_id", "log2FoldChange"))
    check_gene_ids_are_text(de)
    check_numeric_columns(de, ("log2FoldChange",))
    finite = chosen in (EcdfTest.DTS, EcdfTest.WASSERSTEIN)
    x = _fold_changes("target", de, target, chosen, finite=finite)
    y = _fold_changes("background", de, background, chosen, finite=finite)
    scored = score(x, y)
    if isinstance(scored, WassersteinAreas):
        return EcdfTestResult(
            test=chosen,
            statistic=scored.negative_censored,
            p_value=None,
            method=None,
            alternative=None,
            n_target=x.size,
            n_background=y.size,
            areas=scored,
        )
    return EcdfTestResult(
        test=chosen,
        statistic=scored.statistic,
        p_value=scored.p_value,
        method=scored.method,
        alternative=scored.alternative,
        n_target=scored.n_x,
        n_background=scored.n_y,
        areas=None,
    )


def assign_groups(de: pd.DataFrame, groups: Mapping[str, Iterable[str]]) -> pd.Series:
    """Return the group of each row of `de`, as `.assign_groups` of R and its `factor` do.

    The result is a categorical `Series` named `group`, with the index of `de`. Its categories are
    the keys of `groups`, in their order. A row whose `gene_id` is in no group has no value. If a
    gene is in more than one group, the last group wins, as in R. Then the function warns with
    `SeedMatchWarning`, and the message names the first pair of groups that share a gene. A gene
    that occurs two times in one group gives no warning (`D-039`). A repeated `gene_id` column
    raises `InputError` (`D-051`).
    """
    if not isinstance(groups, Mapping):
        message = (
            f"groups {shown(groups)} has type {type(groups).__name__}; "
            "expected a mapping from a group name to gene IDs"
        )
        raise InputError(message)
    check_columns(de, ("gene_id",))
    check_gene_ids_are_text(de)
    members = {
        as_text("group name", name): _as_genes(f"groups[{name!r}]", genes)
        for name, genes in groups.items()
    }
    _warn_for_shared_gene(members)
    gene_ids = de["gene_id"]
    codes = np.full(len(de), -1, dtype=np.int64)
    for code, genes in enumerate(members.values()):
        codes[gene_ids.isin(genes).to_numpy()] = code
    labels = pd.Categorical.from_codes(codes.tolist(), categories=pd.Index(list(members)))
    return pd.Series(labels, index=de.index, name="group")


def _check_options(
    test: EcdfTest, alternative: object, rng: object, n_boot: object, n_interp: object
) -> None:
    """Raise `InputError` for an `alternative` or an option that the test does not use."""
    if test in (EcdfTest.DTS, EcdfTest.WASSERSTEIN) and alternative is not None:
        given = (
            f"Alternative.{alternative.name}"
            if isinstance(alternative, Alternative)
            else shown(alternative)
        )
        message = (
            f"alternative {given} does not apply to {_label(test)}; "
            "expected None, because the test has no direction"
        )
        raise InputError(message)
    options = (
        ("rng", rng, EcdfTest.DTS),
        ("n_boot", n_boot, EcdfTest.DTS),
        ("n_interp", n_interp, EcdfTest.WASSERSTEIN),
    )
    for name, value, owner in options:
        if value is not None and test is not owner:
            message = (
                f"{name} {shown(value)} applies only to {_label(owner)}, not to {_label(test)}; "
                "expected None"
            )
            raise InputError(message)


def _label(test: EcdfTest) -> str:
    return f"EcdfTest.{test.name}"


def _scorer(
    test: EcdfTest,
    alternative: Alternative | None,
    rng: np.random.Generator | int | None,
    n_boot: int | None,
    n_interp: int | None,
) -> Callable[
    [npt.NDArray[np.float64], npt.NDArray[np.float64]], TwoSampleResult | WassersteinAreas
]:
    """Return the function that scores the two samples. `DTS` without `rng` raises `InputError`.

    The function is made before the rows are selected, so that a missing `rng` and an `alternative`
    that is not valid are the first errors.
    """
    directed = Alternative.GREATER if alternative is None else as_alternative(alternative)
    match test:
        case EcdfTest.KS:
            return partial(ks_test, alternative=directed)
        case EcdfTest.WILCOXON:
            return partial(wilcoxon_test, alternative=directed)
        case EcdfTest.DTS:
            if rng is None:
                message = (
                    "rng is None; expected a numpy.random.Generator or an int, "
                    "because the DTS p-value comes from permutations"
                )
                raise InputError(message)
            return partial(_dts, rng=rng, n_boot=n_boot)
        case EcdfTest.WASSERSTEIN:
            return partial(_wasserstein, n_interp=n_interp)
        case _:
            assert_never(test)


def _dts(
    x: npt.NDArray[np.float64],
    y: npt.NDArray[np.float64],
    *,
    rng: np.random.Generator | int,
    n_boot: int | None,
) -> TwoSampleResult:
    check_dts_range("target", "background", x, y)
    if n_boot is None:
        return dts_test(x, y, rng=rng)
    return dts_test(x, y, rng=rng, n_boot=n_boot)


def _wasserstein(
    x: npt.NDArray[np.float64], y: npt.NDArray[np.float64], *, n_interp: int | None
) -> WassersteinAreas:
    check_wasserstein_sample("target", x)
    check_wasserstein_sample("background", y)
    check_wasserstein_range("target", "background", x, y)
    if n_interp is None:
        return wasserstein_areas(x, y)
    return wasserstein_areas(x, y, n_interp=n_interp)


def _as_genes(name: str, genes: Iterable[str]) -> list[str]:
    """Return the distinct gene IDs in the order of their first occurrence."""
    items = as_items(name, genes, as_text, "gene IDs", container="an iterable")
    return list(dict.fromkeys(items))


def _fold_changes(
    name: str, de: pd.DataFrame, genes: Iterable[str], test: EcdfTest, *, finite: bool
) -> npt.NDArray[np.float64]:
    """Return the `log2FoldChange` of the rows in the set, without the NaN values.

    With `finite`, an infinite value raises `InputError` that names the set, the gene ID of the
    first such row, and the test.
    """
    selected = de.loc[de["gene_id"].isin(_as_genes(name, genes)), ["gene_id", "log2FoldChange"]]
    selected = selected.dropna(subset=["log2FoldChange"])
    if selected.empty:
        message = (
            f"{name} has 0 rows with a log2FoldChange in the DE table; expected at least 1 row"
        )
        raise InputError(message)
    values = selected["log2FoldChange"].to_numpy(dtype=np.float64)
    infinite = np.isinf(values)
    if finite and infinite.any():
        position = int(infinite.argmax())
        message = (
            f"{name} gene {shown(selected['gene_id'].iloc[position])} has log2FoldChange "
            f"{float(values[position])}; expected a finite value for {_label(test)}"
        )
        raise InputError(message)
    return values


def _warn_for_shared_gene(members: Mapping[str, list[str]]) -> None:
    owner: dict[str, str] = {}
    for name, genes in members.items():
        for gene in genes:
            first = owner.setdefault(gene, name)
            if first != name:
                message = (
                    f"groups {first!r} and {name!r} share the gene {gene!r}; the last group wins"
                )
                warn(message)
                return
