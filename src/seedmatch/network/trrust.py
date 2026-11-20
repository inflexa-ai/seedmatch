"""The reader of the TRRUST table of regulators and targets, and the scores of its genes."""

import os
import warnings
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Final

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from seedmatch._checks import (
    as_items,
    as_member,
    as_text,
    check_columns,
    check_text_columns,
    concatenate_and_clear,
    shown,
)
from seedmatch._io import local_path, open_text
from seedmatch._lines import (
    field_lines,
    frame_of,
    not_utf8,
    parse_data_lines,
    read_chunks,
    vocabulary_codes,
)
from seedmatch.errors import InputError, SeedMatchWarning
from seedmatch.network.messages import Aggregator, mean_messages, sender_weights

TRRUST_MODES: Final = ("Activation", "Repression", "Unknown")
_COLUMNS: Final = ("source", "target", "edge_type", "evidence")
_REQUIRED: Final = ("source", "target", "edge_type")
_LABEL: Final = "TRRUST table"
_CHUNK_LINES: Final = 2**16
_MODE_CODES: Final = {mode: code for code, mode in enumerate(TRRUST_MODES)}
_EXPECTED_MODES: Final = f"{', '.join(TRRUST_MODES[:-1])}, or {TRRUST_MODES[-1]}"


def read_trrust(path: str | os.PathLike[str]) -> pd.DataFrame:
    """Read the TRRUST table of regulators and targets.

    The path is a local path: a URL raises `InputError`, because the core modules do not use the
    network. The file is UTF-8 text, or gzip, bzip2, or xz of UTF-8 text, and the first bytes, not
    the suffix, show the compression. It is the file that `seedmatch-fetch trrust` downloads. The
    file has no header. Each line has 4 fields with a tab between them: the regulator, the target,
    the mode, and the PubMed IDs. A blank line is skipped.

    The columns are `source`, `target`, `edge_type`, and `evidence`, as `R/get_network_matches.R`
    of SeedMatchR names them. The dtype is `object`, equal texts are one object, and the rows keep
    the order of the file. Each line is a row, also line 1. R reads line 1 as the column names
    (`D-057`). A quote character, `#`, and the text `NA` are text (`D-057`). The mode is
    `Activation`, `Repression`, or `Unknown`.

    These inputs raise `InputError`:

    - a line with another number of fields
    - a mode that is not in `TRRUST_MODES`
    - an empty regulator or an empty target
    - a line of more than 1024 characters with its line end
    - a line with a NUL character
    - text that is not valid UTF-8
    - a file without a data line

    R skips a byte order mark at the start of the file, and so does the reader. The line number in
    a message counts the blank lines. The reader reads 2^16 lines at a time, and the time is
    linear in the number of lines. The result holds four pointers (32 bytes) for each row, and
    the distinct strings. The TRRUST human file of `seedmatch-fetch trrust` has 9,396 lines. The
    longest line has 273 characters with its line end. The read takes 0.008 seconds, and the peak
    memory above the start is about 8 MB.
    """
    source = os.fspath(path)
    file = local_path(source)
    rows = _Rows()
    try:
        with open_text(file) as handle:
            for first, chunk in read_chunks(handle, _CHUNK_LINES):
                is_data = field_lines(
                    chunk,
                    first,
                    source,
                    separator="\t",
                    fields=len(_COLUMNS),
                    described="the regulator, the target, the mode, and the PubMed IDs",
                )
                rows.add(chunk, is_data, first, source)
    except UnicodeDecodeError as exc:
        raise not_utf8(source, exc, "TRRUST") from exc
    return rows.table(source)


@dataclass
class _Rows:
    """The rows of the chunks read so far.

    The four columns share one vocabulary, and a column holds the code of each text. Thus the
    chunks keep one string for each distinct text, not one for each row.
    """

    texts: dict[str, int] = field(default_factory=dict)
    codes: dict[str, list[NDArray[np.int32]]] = field(
        default_factory=lambda: {name: [] for name in _COLUMNS}
    )
    data_lines: int = 0

    def add(self, chunk: list[str], is_data: NDArray[np.bool_], first: int, source: str) -> None:
        """Parse the data lines of the chunk. `first` is the number of the first line."""
        numbers = first + np.flatnonzero(is_data)
        if not numbers.size:
            return
        frame = parse_data_lines(chunk, is_data, first, source, separator="\t", names=_COLUMNS)
        _check_rows(frame, numbers, source)
        self.data_lines += len(frame)
        for name in _COLUMNS:
            self.codes[name].append(vocabulary_codes(frame[name], self.texts))

    def table(self, source: str) -> pd.DataFrame:
        """Return the table of the rows. A file without a row raises `InputError`."""
        if not self.data_lines:
            message = f"{source!r}: the file has no line; expected TRRUST lines"
            raise InputError(message)
        texts = np.array(list(self.texts), dtype=object)
        return frame_of({name: texts[concatenate_and_clear(self.codes[name])] for name in _COLUMNS})


def _check_rows(frame: pd.DataFrame, numbers: NDArray[np.int64], source: str) -> None:
    """Raise `InputError` for the first row with an empty gene or a mode that is not valid."""
    regulators = frame["source"].to_numpy() == ""
    targets = frame["target"].to_numpy() == ""
    modes = frame["edge_type"]
    bad = regulators | targets | ~modes.isin(TRRUST_MODES).to_numpy()
    if not bad.any():
        return
    row = int(np.argmax(bad))
    if regulators[row]:
        problem = "the regulator is empty; expected a gene symbol"
    elif targets[row]:
        problem = "the target is empty; expected a gene symbol"
    else:
        problem = f"mode {modes.iloc[row]!r} is not valid; expected {_EXPECTED_MODES}"
    message = f"{source!r}: line {numbers[row]}: {problem}"
    raise InputError(message)


def trrust_scores(
    network: pd.DataFrame,
    seed_genes: Iterable[str],
    *,
    aggregator: Aggregator = Aggregator.MEAN,
) -> pd.DataFrame:
    """Return the repression score and the activation score of each gene of the network.

    `network` has the columns `source`, `target`, and `edge_type` of `read_trrust`, and each of
    them holds only `str` values. The column `evidence` is optional. `seed_genes` names the genes
    with a seed match, by the symbol that TRRUST uses. A gene of `seed_genes` that is not in the
    network has no effect. If no gene of `seed_genes` is in the network, the function warns with
    `SeedMatchWarning`, and each seed flag is 0 (`D-058`).

    An edge is a distinct (`source`, `target`, `edge_type`) triple, thus a triple on more than one
    line counts one time. An edge goes from a regulator to its target. A `Repression` edge is in
    the repression channel, and an `Activation` edge is in the activation channel. An `Unknown`
    edge is in no channel. A pair with two modes is in each channel of its modes.

    The score of a gene in a channel comes from its regulators in the channel. It is the mean of
    the seed flags of the regulators, plus the seed flag of the gene, all divided by 2. A gene
    without a regulator in the channel has the score NaN, and its count is 0. With
    `Aggregator.DEGREE_IMPORTANCE`, the flag of a regulator is divided by its out-degree, the
    number of its distinct targets over each mode, also `Unknown`. R takes the messages from the
    targets and not from the regulators, and it ignores the seed flag of the gene with the mean
    (`D-007`).

    The result has one row for each gene, in the order of the first occurrence of the gene in
    `source` and then in `target`, as R orders the vertices. The columns are `gene`, `seed`
    (`bool`), `repressors` (`int64`), `repression_score` (`float64`), `activators` (`int64`), and
    `activation_score` (`float64`). `repressors` and `activators` are the numbers of regulators in
    the channels.

    These inputs raise `InputError`:

    - a `network` that is not a `DataFrame`, or that lacks one of the three columns
    - a value of these columns that is not a `str`
    - a mode that is not in `TRRUST_MODES`
    - a `network` without a row
    - `seed_genes` that is a `str`, or that holds a value that is not a `str`
    - an `aggregator` that is not an `Aggregator`

    The time is O(R log R) for R rows, because `np.unique` sorts the edge keys. The memory is
    O(R + G) for G genes. On the TRRUST human file of `seedmatch-fetch trrust`, the call takes
    0.003 seconds, and the peak memory above the start is about 1 MB.
    """
    check_columns(network, _REQUIRED, name="network", label=_LABEL)
    check_text_columns(network, _REQUIRED, label=_LABEL)
    if not len(network):
        message = "the TRRUST table has no row; expected at least one edge"
        raise InputError(message)
    modes = _mode_codes(network["edge_type"])
    seeds = frozenset(as_items("seed_genes", seed_genes, as_text, "str"))
    method = as_member("aggregator", aggregator, Aggregator)

    rows = len(network)
    codes, genes = pd.factorize(
        np.concatenate([network[name].to_numpy(dtype=object) for name in ("source", "target")])
    )
    gene_count = len(genes)
    is_seed = np.fromiter((gene in seeds for gene in genes), dtype=np.bool_, count=gene_count)
    if not is_seed.any():
        message = (
            "no gene of seed_genes is a gene of the TRRUST table; each seed flag is 0. "
            f"TRRUST names a gene by its symbol, for example {shown(genes[0])}"
        )
        warnings.warn(message, SeedMatchWarning, stacklevel=2)

    # A key is below 3 * G^2 for G genes, and G is at most 2 * rows. Thus a key fits an int64
    # for a table of fewer than 8e8 rows.
    keys = np.unique((codes[:rows] * gene_count + codes[rows:]) * len(TRRUST_MODES) + modes)
    pairs, edge_modes = np.divmod(keys, len(TRRUST_MODES))
    senders, receivers = np.divmod(pairs, gene_count)
    is_new_pair = np.ones(len(pairs), dtype=np.bool_)
    is_new_pair[1:] = pairs[1:] != pairs[:-1]
    out_degree = np.bincount(senders[is_new_pair], minlength=gene_count)
    weights = sender_weights(method, out_degree.astype(np.int64, copy=False))

    flags = is_seed.astype(np.float64)
    columns: dict[str, NDArray[np.generic]] = {"gene": genes, "seed": is_seed}
    for mode, count_name, score_name in (
        ("Repression", "repressors", "repression_score"),
        ("Activation", "activators", "activation_score"),
    ):
        in_channel = edge_modes == _MODE_CODES[mode]
        columns[count_name], columns[score_name] = mean_messages(
            receivers[in_channel], senders[in_channel], flags, weights
        )
    return frame_of(columns)


def _mode_codes(column: pd.Series) -> NDArray[np.int64]:
    """Return the code of each mode. A mode that is not in `TRRUST_MODES` raises `InputError`."""
    values = column.to_numpy(dtype=object)
    codes = pd.Categorical(values, categories=list(TRRUST_MODES)).codes.astype(np.int64)
    invalid = np.flatnonzero(codes < 0)
    if invalid.size:
        row = int(invalid[0])
        message = (
            f"{_LABEL} column 'edge_type' holds {shown(values[row])} in row {row}; "
            f"expected {_EXPECTED_MODES}"
        )
        raise InputError(message)
    return codes
