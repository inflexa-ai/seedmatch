"""The reader of the STRING protein links, their gene pairs, and the scores of their genes."""

import os
from dataclasses import dataclass, field
from typing import Final

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from seedmatch._checks import (
    as_cutoff,
    as_limit,
    as_member,
    as_text,
    check_columns,
    check_gene_ids_are_text,
    check_integer_columns,
    check_numeric_columns,
    check_text_columns,
    check_type,
    column_values,
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
from seedmatch._rmath import r_quantile
from seedmatch.annotation.reference import Reference
from seedmatch.errors import InputError, warn
from seedmatch.network.messages import Aggregator, mean_messages, pair_keys, sender_weights

STRING_HEADER: Final = "protein1 protein2 combined_score"
STRING_SCORE_MAX: Final = 1000
_PROTEIN_COLUMNS: Final = ("protein1", "protein2")
_SCORE_COLUMNS: Final = ("combined_score",)
_COLUMNS: Final = (*_PROTEIN_COLUMNS, *_SCORE_COLUMNS)
_LINKS_LABEL: Final = "links table"
_GENE_COLUMNS: Final = ("gene1", "gene2")
_NETWORK_LABEL: Final = "network table"
_CHUNK_LINES: Final = 2**16


@dataclass(frozen=True, eq=False)
class StringLinks:
    """The links of the STRING network with a score above the cutoff, and the cutoff.

    `links` has the columns `protein1` and `protein2` (`object`, equal texts are one object) and
    `combined_score` (`int64`), and a `RangeIndex`. The rows keep the order of the file. `cutoff`
    is the percentile of the scores of the links that the file holds with the score threshold or
    more.
    """

    links: pd.DataFrame
    cutoff: float


def read_string_links(
    path: str | os.PathLike[str], *, score_threshold: int = 400, percentile: float = 0.95
) -> StringLinks:
    """Read the STRING protein links, and keep the links with a score above a percentile.

    This is the selection of `get_stringPPI_to_stringGGI` of SeedMatchR, which reads the links
    through STRINGdb 2.22.0 (`R/rstring.R:1064-1067`). The path is a local path: a URL raises
    `InputError`, because the core modules do not use the network. The file is UTF-8 text, or
    gzip, bzip2, or xz of UTF-8 text, and the first bytes, not the suffix, show the compression.
    It is the file that `seedmatch-fetch string` downloads. Line 1 is the header
    `protein1 protein2 combined_score`. Each other line has 3 fields with one space between them:
    two protein IDs and the combined score. A blank line is skipped.

    The reader keeps the links with `combined_score >= score_threshold`, and it takes the
    `percentile` of their scores with `r_quantile`, the quantile of R. Each kept row counts, thus
    a symmetric file counts each pair two times, as R does. The result keeps the links with a
    score above this cutoff, in file order. `score_threshold` is an int from 0 to 1000. The
    default 400 is the default of STRINGdb 2.22.0. `percentile` is a number from 0 to 1. The
    default 0.95 is the default `p` of `get_stringPPI_to_stringGGI`.

    The score is an integer from 0 to 1000. A score with an integer value, for example `400.0`
    or `4e2`, is valid (`D-059`). The protein IDs stay as they are: a quote character and `#` have
    no special meaning (`D-059`). These inputs raise `InputError`:

    - a first line that is not the header
    - a line with another number of fields (`D-059`)
    - a score that is not an integer from 0 to 1000 (`D-059`)
    - an empty protein ID
    - a line of more than 1024 characters with its line end (`D-059`)
    - a line with a NUL character (`D-059`)
    - text that is not valid UTF-8
    - a file without a link, with no link of `score_threshold` or more, or with no link above the
      cutoff. R stops with `non-character argument` for the last input.

    R skips a byte order mark at the start of the file, and so does the reader. The line number
    in a message counts the blank lines.

    The reader reads 2^16 lines at a time, and the time is linear in the number of lines. The work
    memory is one chunk, 10 bytes for each kept row, and the distinct protein IDs. The 10 bytes
    are two 4-byte codes and a 2-byte score. The cutoff needs 8 more bytes for each kept row for
    a moment. The result holds 24 bytes for each link and the distinct strings.

    The STRING rat file of `seedmatch-fetch string 10116` has 14,100,659 lines and 761,435,565
    bytes of text. The longest line has 54 characters with its line end. With the defaults, the
    read takes 8.6 seconds, and the peak memory above the start is about 95 MB. The result has
    98,920 links, the cutoff is 940.0, and the links have 8,809 proteins. With `score_threshold=0`
    and `percentile=0`, the result has 13,894,314 links, the read takes 10.1 seconds, and the peak
    is about 660 MB. Machine: Apple M2 Max, Python 3.12.
    """
    source = os.fspath(path)
    file = local_path(source)
    threshold = _as_threshold(score_threshold)
    probability = _as_percentile(percentile)
    rows = _Rows(threshold)
    try:
        with open_text(file) as handle:
            for first, chunk in read_chunks(handle, _CHUNK_LINES):
                rows.add(chunk, first, source)
    except UnicodeDecodeError as exc:
        raise not_utf8(source, exc, "STRING links") from exc
    return rows.result(source, probability)


def _as_threshold(value: object) -> int:
    threshold = as_limit("score_threshold", value)
    if not 0 <= threshold <= STRING_SCORE_MAX:
        message = (
            f"score_threshold {threshold} is out of range; "
            f"expected an integer from 0 to {STRING_SCORE_MAX}"
        )
        raise InputError(message)
    return threshold


def _as_percentile(value: object) -> float:
    percentile = as_cutoff("percentile", value)
    if not 0.0 <= percentile <= 1.0:
        message = f"percentile {percentile} is out of range; expected a number from 0 to 1"
        raise InputError(message)
    return percentile


@dataclass
class _Rows:
    """The kept links of the chunks read so far.

    A link holds the codes of its two proteins in one vocabulary, and its score. Thus the chunks
    keep one string for each distinct protein ID, not one for each link.
    """

    threshold: int
    proteins: dict[str, int] = field(default_factory=dict)
    first_codes: list[NDArray[np.int32]] = field(default_factory=list)
    second_codes: list[NDArray[np.int32]] = field(default_factory=list)
    scores: list[NDArray[np.int16]] = field(default_factory=list)
    has_header: bool = False
    data_lines: int = 0

    def add(self, chunk: list[str], first: int, source: str) -> None:
        """Parse the data lines of the chunk, and keep the links with a score of the threshold.

        `first` is the number of the first line. A chunk that starts at line 1 starts with the
        header.
        """
        if first == 1:
            _check_header(chunk[0], source)
            self.has_header = True
        is_data = field_lines(
            chunk,
            first,
            source,
            separator=" ",
            fields=len(_COLUMNS),
            described="protein1, protein2, and combined_score",
        )
        if first == 1:
            is_data[0] = False
        numbers = first + np.flatnonzero(is_data)
        if not numbers.size:
            return
        frame, scores = _read_chunk(chunk, is_data, first, source, numbers)
        self.data_lines += len(frame)
        keep = scores >= self.threshold
        if not keep.any():
            return
        frame = frame[keep]
        self.first_codes.append(vocabulary_codes(frame["protein1"], self.proteins))
        self.second_codes.append(vocabulary_codes(frame["protein2"], self.proteins))
        self.scores.append(scores[keep].astype(np.int16))

    def result(self, source: str, percentile: float) -> StringLinks:
        """Return the links above the cutoff. A file without such a link raises `InputError`."""
        if not self.has_header:
            message = f"{source!r}: the file has no line; expected the header {STRING_HEADER!r}"
            raise InputError(message)
        if not self.data_lines:
            message = f"{source!r}: the file has no line after the header; expected STRING links"
            raise InputError(message)
        if not self.scores:
            message = (
                f"{source!r}: no link has a score of {self.threshold} or more; "
                "expected at least one link"
            )
            raise InputError(message)
        scores = concatenate_and_clear(self.scores)
        cutoff = r_quantile(scores, percentile)
        is_above = scores > cutoff
        if not is_above.any():
            message = (
                f"{source!r}: no link has a score above the cutoff {cutoff}, "
                f"the percentile {percentile} of the {len(scores)} scores of "
                f"{self.threshold} or more; expected at least one link"
            )
            raise InputError(message)
        texts = np.array(list(self.proteins), dtype=object)
        links = frame_of(
            {
                "protein1": texts[concatenate_and_clear(self.first_codes)[is_above]],
                "protein2": texts[concatenate_and_clear(self.second_codes)[is_above]],
                "combined_score": scores[is_above].astype(np.int64),
            }
        )
        return StringLinks(links=links, cutoff=cutoff)


def _check_header(line: str, source: str) -> None:
    """Raise `InputError` if line 1 is not the header."""
    text = line.rstrip("\r\n")
    if text != STRING_HEADER:
        message = f"{source!r}: line 1 is {text!r}; expected the header {STRING_HEADER!r}"
        raise InputError(message)


def _read_chunk(
    chunk: list[str],
    is_data: NDArray[np.bool_],
    first: int,
    source: str,
    numbers: NDArray[np.int64],
) -> tuple[pd.DataFrame, NDArray[np.int64]]:
    """Return the table of the data lines and their scores. The first bad line raises `InputError`.

    The read with an integer score is the fast path. If pandas rejects a score, or a score is
    outside the range, the function reads the chunk again with a text score to find the first bad
    line. A score such as `4e2` has an integer value, and it is valid.
    """
    try:
        # pandas casts an infinite score to int64, which warns, before it raises ValueError.
        with np.errstate(invalid="ignore"):
            frame = parse_data_lines(
                chunk,
                is_data,
                first,
                source,
                separator=" ",
                names=_COLUMNS,
                dtype={"protein1": str, "protein2": str, "combined_score": np.dtype(np.int64)},
            )
    except InputError:
        raise
    except (ValueError, OverflowError):
        pass
    else:
        scores = frame["combined_score"].to_numpy()
        if scores.min() >= 0 and scores.max() <= STRING_SCORE_MAX:
            _check_rows(frame, np.zeros(len(frame), dtype=np.bool_), numbers, source)
            return frame, scores
    frame = parse_data_lines(chunk, is_data, first, source, separator=" ", names=_COLUMNS)
    values = pd.to_numeric(frame["combined_score"], errors="coerce").to_numpy(dtype=np.float64)
    is_score = (values == np.floor(values)) & (values >= 0) & (values <= STRING_SCORE_MAX)
    _check_rows(frame, ~is_score, numbers, source)
    return frame, values.astype(np.int64)


def _check_rows(
    frame: pd.DataFrame, is_bad_score: NDArray[np.bool_], numbers: NDArray[np.int64], source: str
) -> None:
    """Raise `InputError` for the first row with an empty protein ID or a bad score."""
    first_empty = frame["protein1"].eq("").to_numpy()
    second_empty = frame["protein2"].eq("").to_numpy()
    bad = first_empty | second_empty | is_bad_score
    if not bad.any():
        return
    row = int(np.argmax(bad))
    if first_empty[row]:
        problem = "protein1 is empty; expected a STRING protein ID"
    elif second_empty[row]:
        problem = "protein2 is empty; expected a STRING protein ID"
    else:
        problem = (
            f"score {frame['combined_score'].iloc[row]!r} is not an integer from 0 to "
            f"{STRING_SCORE_MAX}; expected a STRING combined score"
        )
    message = f"{source!r}: line {numbers[row]}: {problem}"
    raise InputError(message)


@dataclass(frozen=True, eq=False)
class GeneNetwork:
    """The gene pairs of STRING links, and the proteins that have no gene.

    `edges` has the columns `gene1` and `gene2` (`object`) and `combined_score` (`int64`), and a
    `RangeIndex`. It has a row for each order of each gene pair. `unmapped_proteins` holds the
    STRING ID of each protein that the annotation maps to no gene, in the order of their first
    occurrence in `protein1` and then in `protein2`.
    """

    edges: pd.DataFrame
    unmapped_proteins: tuple[str, ...]


def string_gene_network(links: pd.DataFrame, reference: Reference) -> GeneNetwork:
    """Map the proteins of STRING links to genes with the EnsDb of the reference.

    This is the map of `get_stringPPI_to_stringGGI` of SeedMatchR, with the `protein` table of
    the EnsDb in place of BioMart (`D-013`). `links` has the columns `protein1` and `protein2`
    (`str`) and `combined_score` (an integer dtype), as `StringLinks.links` has. `reference` must
    have an EnsDb. An annotation of `load_gtf` has no `protein` table.

    A STRING protein ID has the form `<taxon>.<protein ID>`. The protein ID of the EnsDb is the
    text after the first `.`, also if the ID has a second `.` (`D-061`). A link of two proteins
    with a gene each gives the pair of their genes. A link of two proteins of one gene gives no
    pair, and a link of a protein with itself gives no pair (`D-060`). A link of a protein with no
    gene gives no pair (`D-061`). Thus the network has gene IDs only, and `unmapped_proteins`
    lists the other proteins. A protein of a transcript without a gene is one of them.

    R keeps the first protein of each gene and drops the links of the others. Here the links of
    each protein of a gene join, and a pair gets the highest score of its links (`D-060`).

    The edges have a row for each order of each gene pair. The genes come in the order of the
    first occurrence of one of their proteins, in `protein1` and then in `protein2`. The rows
    come in the order of `gene1`, and then of `gene2`. For a symmetric table in the order of
    `protein1`, as STRING gives it, the genes come in the vertex order of R.

    These inputs raise `InputError`:

    - a `links` that is not a `DataFrame`, or that lacks one of the three columns
    - a protein ID that is not a `str`
    - a score that is not an integer, or that is not from 0 to 1000
    - a `links` without a row
    - a `reference` that is not a `Reference`, or that has no `protein` table
    - a protein ID without a taxon prefix, an empty taxon, or an empty protein ID
    - a protein with more than one gene in the EnsDb
    - links that join no two genes. The message names the annotation.

    The time is O(E log E) for E links, because the function sorts the keys `gene1 * G + gene2`
    for G genes one time. The memory is O(E). The map asks the EnsDb for the proteins of the links
    with one query. The default links of the STRING rat file (`read_string_links`) have 98,920
    links and 8,809 proteins. With the EnsDb 109 of `rnor7`, 628 proteins have no gene, and 96
    genes have more than one protein (at most 4). The result has 84,654 rows and 7,979 genes. This
    is 1,112 rows more than the gene pairs of R. The first call takes 0.2 seconds, a later call
    takes 0.04 seconds, and the peak memory above the start is about 25 MB. With
    `score_threshold=0` and `percentile=0`, the links have 13,894,314 rows. The call takes 3.6
    seconds, and the result has 11,097,016 rows. The peak memory above the start is about 1.8 GB,
    or 130 bytes for each link. Machine: Apple M2 Max, Python 3.12.
    """
    check_columns(links, _COLUMNS, name="links", label=_LINKS_LABEL)
    check_text_columns(links, _PROTEIN_COLUMNS, label=_LINKS_LABEL)
    check_integer_columns(links, _SCORE_COLUMNS, label=_LINKS_LABEL)
    _check_scores(links["combined_score"])
    if not len(links):
        message = "the links table has no row; expected at least one link"
        raise InputError(message)
    check_type("reference", reference, Reference)

    rows = len(links)
    codes, proteins = pd.factorize(
        np.concatenate([links[name].to_numpy(dtype=object) for name in _PROTEIN_COLUMNS])
    )
    ids = [_protein_id(protein) for protein in proteins]
    gene_of = reference.annotation.protein_genes(ids)
    found = [gene_of.get(protein_id) for protein_id in ids]
    is_mapped = np.fromiter((gene is not None for gene in found), dtype=np.bool_, count=len(ids))
    gene_codes, genes = pd.factorize(
        np.array([gene for gene in found if gene is not None], dtype=object)
    )
    gene_of_protein = np.full(len(ids), -1, dtype=np.int64)
    gene_of_protein[is_mapped] = gene_codes

    first = gene_of_protein[codes[:rows]]
    second = gene_of_protein[codes[rows:]]
    is_pair = (first >= 0) & (second >= 0) & (first != second)
    if not is_pair.any():
        message = (
            f"no link joins two genes of the annotation {reference.annotation.source!r}; "
            "expected the STRING links of the species of the annotation"
        )
        raise InputError(message)
    first, second = first[is_pair], second[is_pair]
    scores = links["combined_score"].to_numpy(dtype=np.int64)[is_pair]

    gene_count = len(genes)
    keys = np.concatenate(
        [pair_keys(first, second, gene_count), pair_keys(second, first, gene_count)]
    )
    order = np.argsort(keys)
    keys = keys[order]
    is_start = np.ones(len(keys), dtype=np.bool_)
    is_start[1:] = keys[1:] != keys[:-1]
    starts = np.flatnonzero(is_start)
    best = np.maximum.reduceat(np.concatenate([scores, scores])[order], starts)
    gene1, gene2 = np.divmod(keys[starts], gene_count)
    edges = frame_of({"gene1": genes[gene1], "gene2": genes[gene2], "combined_score": best})
    return GeneNetwork(edges=edges, unmapped_proteins=tuple(proteins[~is_mapped].tolist()))


def _check_scores(scores: pd.Series) -> None:
    """Raise `InputError` for the first score outside 0 to `STRING_SCORE_MAX`.

    The comparison runs on the column as it is, thus a `uint64` score above 2^63 - 1 is found
    before the cast to `int64`.
    """
    is_bad = ((scores < 0) | (scores > STRING_SCORE_MAX)).to_numpy(dtype=np.bool_)
    if is_bad.any():
        row = int(np.argmax(is_bad))
        message = (
            f"{_LINKS_LABEL} column 'combined_score' holds {int(scores.iloc[row])} in row {row}; "
            f"expected an integer from 0 to {STRING_SCORE_MAX}"
        )
        raise InputError(message)


def _protein_id(protein: str) -> str:
    """Return the text after the first `.` of a STRING ID. An ID without a taxon raises."""
    taxon, dot, identifier = protein.partition(".")
    if not (taxon and dot and identifier):
        message = (
            f"{_LINKS_LABEL}: protein ID {protein!r} has no taxon prefix; "
            "expected the form '<taxon>.<protein ID>'"
        )
        raise InputError(message)
    return identifier


def string_scores(
    network: pd.DataFrame,
    de: pd.DataFrame,
    *,
    column: str,
    aggregator: Aggregator = Aggregator.MEAN,
) -> pd.DataFrame:
    """Return the seed count, the neighbors, and the score of each gene of a STRING network.

    This is `create_graph` of SeedMatchR, which sets the feature of each node to its seed count
    (`D-062`). R has no score of a STRING network, thus `neighbors` and `score` are new. `network`
    has the columns `gene1` and `gene2` (`str`), as `GeneNetwork.edges` has. The column
    `combined_score` is optional, and the function does not read it. A row is an edge from `gene1`
    to `gene2`. The network must have a row for each order of each pair, as `string_gene_network`
    gives it. The neighbors of a gene are the `gene2` values of its rows. `de` has the columns
    `gene_id` (`str`) and `column`, a numeric column of seed counts.

    The seed count of a gene is the highest count of its rows in `de`, and it is 0 for a gene
    without a row (`D-062`). A row of `de` with a gene outside the network has no effect. If no
    gene of `de` is a gene of the network, the function warns with `SeedMatchWarning`, and each
    seed count is 0 (`D-062`).

    The score of a gene is the mean of the weighted seed counts of its neighbors, plus its own seed
    count, all divided by 2. With `Aggregator.MEAN`, each weight is 1. With
    `Aggregator.DEGREE_IMPORTANCE`, the weight of a neighbor is 1 divided by its number of
    neighbors. A gene of the network has a neighbor, thus no score is NaN.

    The result has one row for each gene, in the order of the first occurrence of the gene in
    `gene1` and then in `gene2`, and a `RangeIndex`. The columns are `gene_id` (`object`),
    `seed_count` (`float64`), `neighbors` (`int64`), and `score` (`float64`).

    These inputs raise `InputError`:

    - a `network` that is not a `DataFrame`, or that lacks `gene1` or `gene2`
    - a gene that is not a `str`
    - a `network` without a row
    - a row of a gene with itself
    - a pair on more than one row
    - a pair without its reverse row
    - a `column` that is not a `str`, or that `de` does not have
    - a `de` that is not a `DataFrame`, that lacks `gene_id`, or whose `gene_id` is not text
    - a `column` that is not numeric
    - a count that is NaN, infinite, or negative, also in a row of a gene outside the network
      (`D-062`)
    - an `aggregator` that is not an `Aggregator`

    The time is O(E log E + D) for E rows of the network and D rows of `de`. The check of the
    reverse rows sorts the keys one time. The memory is O(E + G + D) for G genes. On the default
    network of the STRING rat file, with 84,654 rows and 7,979 genes, and a `de` of 20,000 rows,
    the call takes 0.015 seconds, and the peak memory above the start is about 5 MB. With
    `score_threshold=0` and `percentile=0`, the network has 11,097,016 rows and 18,223 genes. The
    call takes 1.8 seconds, and the peak memory above the start is about 630 MB, or 57 bytes for
    each row. Machine: Apple M2 Max, Python 3.12.
    """
    check_columns(network, _GENE_COLUMNS, name="network", label=_NETWORK_LABEL)
    check_text_columns(network, _GENE_COLUMNS, label=_NETWORK_LABEL)
    if not len(network):
        message = "the network table has no row; expected at least one gene pair"
        raise InputError(message)
    count_column = as_text("column", column)
    check_columns(de, ("gene_id", count_column))
    check_gene_ids_are_text(de)
    check_numeric_columns(de, [count_column])
    method = as_member("aggregator", aggregator, Aggregator)
    counts = column_values(de, count_column, _is_count, "a finite count of 0 or more")

    rows = len(network)
    codes, genes = pd.factorize(
        np.concatenate([network[name].to_numpy(dtype=object) for name in _GENE_COLUMNS])
    )
    receivers, senders = codes[:rows], codes[rows:]
    _check_pairs(receivers, senders, genes)

    gene_count = len(genes)
    de_genes = pd.Index(de["gene_id"].to_numpy(dtype=object), dtype=object)
    row_genes = pd.Index(genes, dtype=object).get_indexer(de_genes)
    in_network = row_genes >= 0
    if not in_network.any():
        message = (
            "no gene of de is a gene of the network; each seed count is 0. "
            f"Expected the gene IDs of the network, for example {shown(genes[0])}"
        )
        warn(message)
    seed_counts = np.zeros(gene_count, dtype=np.float64)
    np.maximum.at(seed_counts, row_genes[in_network], counts[in_network])

    out_degree = np.bincount(receivers, minlength=gene_count).astype(np.int64, copy=False)
    weights = sender_weights(method, out_degree)
    neighbors, scores = mean_messages(receivers, senders, seed_counts, weights)
    return frame_of(
        {"gene_id": genes, "seed_count": seed_counts, "neighbors": neighbors, "score": scores}
    )


def _is_count(values: NDArray[np.float64]) -> NDArray[np.bool_]:
    return np.isfinite(values) & (values >= 0)


def _check_pairs(
    first: NDArray[np.intp], second: NDArray[np.intp], genes: NDArray[np.object_]
) -> None:
    """Raise `InputError` for a gene with itself, a repeated pair, or a pair without a reverse row.

    The codes in `first` and `second` index `genes`. One sort of the keys finds the repeated
    pairs and the reverse rows. The message names the first row in table order that has the fault.
    """
    is_loop = first == second
    if is_loop.any():
        row = int(np.argmax(is_loop))
        message = (
            f"{_NETWORK_LABEL}: the pair {_pair(genes, first[row], second[row])} joins a gene to "
            "itself; expected two different genes"
        )
        raise InputError(message)

    gene_count = len(genes)
    keys = pair_keys(first, second, gene_count)
    order = np.argsort(keys)
    ordered = keys[order]
    is_repeat = ordered[1:] == ordered[:-1]
    if is_repeat.any():
        in_group = np.zeros(len(keys), dtype=np.bool_)
        in_group[1:] |= is_repeat
        in_group[:-1] |= is_repeat
        row = int(order[in_group].min())
        times = np.searchsorted(ordered, keys[row], side="right") - np.searchsorted(
            ordered, keys[row], side="left"
        )
        message = (
            f"{_NETWORK_LABEL}: the pair {_pair(genes, first[row], second[row])} occurs "
            f"{times} times; expected one row for each pair"
        )
        raise InputError(message)

    reverse_keys = pair_keys(second, first, gene_count)
    position = np.minimum(np.searchsorted(ordered, reverse_keys), len(ordered) - 1)
    has_reverse = ordered[position] == reverse_keys
    if not has_reverse.all():
        row = int(np.argmin(has_reverse))
        message = (
            f"{_NETWORK_LABEL}: the pair {_pair(genes, first[row], second[row])} has no reverse "
            f"row {_pair(genes, second[row], first[row])}; expected both orders of each pair"
        )
        raise InputError(message)


def _pair(genes: NDArray[np.object_], first: np.intp, second: np.intp) -> str:
    return repr((genes[first], genes[second]))
