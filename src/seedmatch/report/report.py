"""`SeedMatchReport` of R: the seed matches of a guide and their effect on a DE table.

R also returns a `gt` table of the report. `report_html` gives its HTML.
"""

from dataclasses import dataclass
from typing import Final

import pandas as pd

from seedmatch._checks import (
    as_cutoff,
    check_columns,
    check_gene_ids_are_text,
    check_numeric_columns,
)
from seedmatch.report.grid import FullSearchOptions, best_match, full_search
from seedmatch.seed import Guide
from seedmatch.sequence import SequenceSet

# The order of `expected.names` (R/SeedMatchReport.R:1-5): the 15mer classes go from D4 to D0.
REPORT_CATEGORIES: Final[tuple[str, ...]] = (
    *(f"Full-D{distance}" for distance in range(5)),
    *(f"18mer-D{distance}" for distance in range(5)),
    *(f"15mer-D{distance}" for distance in range(4, -1, -1)),
    "8mer",
    "7mer-m8",
    "7mer-A1",
    "6mer",
)

_DE_COLUMNS: Final = ("gene_id", "log2FoldChange", "padj", "baseMean")
_DE_VALUE_COLUMNS: Final = _DE_COLUMNS[1:]
_SEARCH_COLUMNS: Final = ("seq_id", "gene_id", "group", "best_match")
_COUNT_ROWS: Final = ("In silico predictions", "Expressed predictions", "Off-target predictions")


@dataclass(frozen=True)
class ReportOptions:
    """The options of `seed_match_report`, with the names of R in snake case.

    The defaults copy `R/SeedMatchReport.R:46-55`. `indels`, `wobble`, `seed_edit_distances`,
    `fixed_pattern`, and `fixed_subject` are the fields of `FullSearchOptions` with the same names.
    A gene is an off-target if its `log2FoldChange` is at most `log2fc_cut` and its `padj` is at
    most `padj_cut`. Each field is checked when the object is made. `wobble=True` needs
    `fixed_pattern=False`, because the wobble target has R at each U of the guide, and a fixed
    pattern matches R only to the letter R (`D-026`). The default of `fixed_pattern` is `True`, as
    in R, thus set it to `False` with the wobble.
    """

    indels: bool = True
    wobble: bool = False
    log2fc_cut: float = -0.7655347
    padj_cut: float = 0.05
    seed_edit_distances: tuple[int, ...] = (0,)
    fixed_pattern: bool = True
    fixed_subject: bool = True

    def __post_init__(self) -> None:
        """Make sure that the options are valid, and store plain floats, tuples, and bools."""
        search = _search_options(self)
        # The dataclass is frozen, so the constructor is the one writer of the normalized values.
        set_field = object.__setattr__
        set_field(self, "log2fc_cut", as_cutoff("log2fc_cut", self.log2fc_cut))
        set_field(self, "padj_cut", as_cutoff("padj_cut", self.padj_cut))
        for name in ("indels", "wobble", "seed_edit_distances", "fixed_pattern", "fixed_subject"):
            set_field(self, name, getattr(search, name))


def _search_options(options: ReportOptions) -> FullSearchOptions:
    """Return the `full_search` options of the report. The non-seed distances are the R default."""
    return FullSearchOptions(
        seed_edit_distances=options.seed_edit_distances,
        indels=options.indels,
        wobble=options.wobble,
        fixed_pattern=options.fixed_pattern,
        fixed_subject=options.fixed_subject,
    )


@dataclass(frozen=True, eq=False)
class SeedMatchReport:
    """The result of `seed_match_report`.

    `report` has one column for each name in `REPORT_CATEGORIES`, then `Total`. `matches` has one
    row for each sequence. A frame has no truth value, so the objects compare by identity.
    """

    report: pd.DataFrame
    matches: pd.DataFrame


def _matches(searched: pd.DataFrame, de: pd.DataFrame) -> pd.DataFrame:
    """Add `detected` and the DE values to the `full_search` table, with one row per sequence."""
    # R takes the first row of a gene (`match`), so a later row never copies a sequence row.
    first_rows = de.drop_duplicates("gene_id")[list(_DE_COLUMNS)]
    values = first_rows.astype(dict.fromkeys(_DE_VALUE_COLUMNS, "float64"))
    detected = searched["gene_id"].isin(values["gene_id"])
    return searched.assign(detected=detected).merge(values, on="gene_id", how="left")


def _gene_table(matches: pd.DataFrame, searched: pd.DataFrame) -> pd.DataFrame:
    """Return one row for each gene: its DE values and its best match over its sequences."""
    counts = searched.drop(columns=list(_SEARCH_COLUMNS))
    # The best match is the first column of the ranking with a count of 1 or more in any sequence
    # of the gene, which is the best match of the maximum of each column.
    by_gene = best_match(counts.groupby(searched["gene_id"]).max())
    # The columns that the merge added are equal in each row of a gene.
    genes = matches[["gene_id", "detected", "log2FoldChange", "padj"]].drop_duplicates("gene_id")
    genes = genes.set_index("gene_id")
    genes["best_match"] = by_gene
    return genes


# Each argument is an independent input of SeedMatchReport, and the options are in ReportOptions.
def seed_match_report(  # noqa: PLR0913
    de: pd.DataFrame,
    sequences: SequenceSet,
    guide: Guide,
    options: ReportOptions = ReportOptions(),
    group: str = "SeedMatchReport",
    *,
    workers: int | None = None,
) -> SeedMatchReport:
    """Count the genes with a seed match, and the genes among them that go down (`SeedMatchReport`).

    `de` needs the columns `gene_id`, `log2FoldChange`, `padj`, and `baseMean`. The `report` has
    the rows `In silico predictions` (the genes with a match), `Expressed predictions` (the genes
    with a match that are in `de`), `Off-target predictions` (the genes with a match that pass
    both cutoffs), and `% off-target` (`100 * off-target / expressed`, and 0 for no expressed
    gene). A gene counts once, in the class of its best match over its sequences (`D-003`).
    `matches` is the `full_search` table with the columns `detected`, `log2FoldChange`, `padj`, and
    `baseMean`. A gene without a row in `de` is not detected and has missing values, and a missing
    value is never an off-target. For a gene with more than one row in `de`, the first row counts.
    `workers` is the width of the thread pool of `full_search`.
    """
    check_columns(de, _DE_COLUMNS)
    check_gene_ids_are_text(de)
    check_numeric_columns(de, _DE_VALUE_COLUMNS)
    searched = full_search(guide, sequences, _search_options(options), group, workers=workers)
    matches = _matches(searched, de)
    genes = _gene_table(matches, searched)

    matched = genes["best_match"].isin(REPORT_CATEGORIES)
    expressed = matched & genes["detected"]
    off_target = (
        matched
        & (genes["log2FoldChange"] <= options.log2fc_cut)
        & (genes["padj"] <= options.padj_cut)
    )
    table = pd.DataFrame(
        {
            name: genes.loc[selected, "best_match"]
            .value_counts()
            .reindex(REPORT_CATEGORIES, fill_value=0)
            for name, selected in zip(_COUNT_ROWS, (matched, expressed, off_target), strict=True)
        }
    )
    table.loc["Total"] = table.sum()
    ratio = table["Off-target predictions"] / table["Expressed predictions"] * 100
    table["% off-target"] = ratio.where(table["Expressed predictions"] > 0, 0.0)

    report = table.T.astype("float64").rename_axis(index="category", columns=None).reset_index()
    report["group"] = group
    return SeedMatchReport(report=report, matches=matches)
