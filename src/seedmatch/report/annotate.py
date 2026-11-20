"""Seed counts as a column of a DE table: `SeedMatchR` of R in `DESEQ2` mode."""

import re
from typing import Final

import numpy as np
import pandas as pd

from seedmatch._checks import (
    as_data_frame,
    as_flag,
    check_gene_ids_are_text,
    check_unique_columns,
    shown,
)
from seedmatch.errors import InputError
from seedmatch.search import SearchOptions, count_matches
from seedmatch.seed import Target
from seedmatch.sequence import SequenceSet
from seedmatch.sequence.alphabet import MASKS

_AMBIGUITY_LETTER: Final = re.compile(
    f"[{''.join(letter for letter, mask in MASKS.items() if mask.bit_count() > 1)}]",
    flags=re.IGNORECASE,
)


# Each keyword-only parameter is an independent option of SeedMatchR, so no parameter object fits.
def annotate_seed_counts(  # noqa: PLR0913
    de: pd.DataFrame,
    sequences: SequenceSet,
    target: Target,
    *,
    column: str | None = None,
    # The defaults of SeedMatchR in R are max.mismatch = 0, with.indels = TRUE, and fixed = TRUE
    # (R/SeedMatchR.R:77-79).
    options: SearchOptions = SearchOptions(indels=True),
    shared_genes: bool = True,
) -> pd.DataFrame:
    """Return a copy of the DE table with the seed count of each gene in a column.

    The count of a gene is the maximum count over its sequences, and 0 for a gene without a
    sequence. The column is named `column`, or the name of the target. A column with that name is
    replaced in place. With `shared_genes`, the result keeps only the rows whose gene has a
    sequence, in the order of the input, and has a new index. With `options.fixed_pattern`, a target
    with an ambiguity letter, for example a wobble target, raises `InputError` (`D-026`). A repeated
    `gene_id` column, or a repeated column with the name of the count column, raises `InputError`
    (`D-051`).
    """
    shared_genes = as_flag("shared_genes", shared_genes)
    if column is not None and not isinstance(column, str):
        message = f"column {shown(column)} has type {type(column).__name__}; expected a str"
        raise InputError(message)
    ambiguous = _AMBIGUITY_LETTER.search(target.sequence) if options.fixed_pattern else None
    if ambiguous is not None:
        letter = ambiguous.group()
        message = (
            f"target {target.name!r} has the ambiguity letter {letter!r}; expected "
            "fixed_pattern=False in the options, because a fixed pattern matches "
            f"{letter!r} only to the letter {letter!r}"
        )
        raise InputError(message)
    if "gene_id" not in as_data_frame(de).columns:
        message = f"DE table has no 'gene_id' column; found {list(de.columns)}"
        raise InputError(message)
    name = target.name if column is None else column
    if name == "gene_id":
        message = "column 'gene_id' would overwrite the gene IDs; expected another column name"
        raise InputError(message)
    check_unique_columns(de, ("gene_id", name))
    check_gene_ids_are_text(de)
    counts = pd.Series(count_matches(target, sequences, options))
    by_gene = counts.groupby(np.asarray(sequences.gene_ids, dtype=object)).max()
    if shared_genes:
        annotated = de.loc[de["gene_id"].isin(by_gene.index)].reset_index(drop=True)
    else:
        annotated = de.copy()
    annotated[name] = by_gene.reindex(annotated["gene_id"], fill_value=0).to_numpy(dtype=np.int64)
    return annotated
