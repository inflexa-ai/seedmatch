"""Seed counts as a column of a DE table: `SeedMatchR` of R in `DESEQ2` mode."""

import numpy as np
import pandas as pd

from seedmatch._checks import check_gene_ids_are_text
from seedmatch.errors import InputError
from seedmatch.search import SearchOptions, count_matches
from seedmatch.seed import Target
from seedmatch.sequence import SequenceSet


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
    sequence, in the order of the input, and has a new index.
    """
    if "gene_id" not in de.columns:
        message = f"DE table has no 'gene_id' column; found {list(de.columns)}"
        raise InputError(message)
    name = target.name if column is None else column
    if name == "gene_id":
        message = "column 'gene_id' would overwrite the gene IDs; expected another column name"
        raise InputError(message)
    check_gene_ids_are_text(de)
    counts = pd.Series(count_matches(target, sequences, options))
    by_gene = counts.groupby(np.asarray(sequences.gene_ids, dtype=object)).max()
    if shared_genes:
        annotated = de.loc[de["gene_id"].isin(by_gene.index)].reset_index(drop=True)
    else:
        annotated = de.copy()
    annotated[name] = by_gene.reindex(annotated["gene_id"], fill_value=0).to_numpy(dtype=np.int64)
    return annotated
