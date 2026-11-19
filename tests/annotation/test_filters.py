import dataclasses
import re

import numpy as np
import pandas as pd
import pytest

from seedmatch.annotation.filters import (
    STANDARD_CHROMOSOMES,
    FeatureType,
    PerGene,
    TranscriptFilter,
)
from seedmatch.errors import InputError

TEXT_FIELDS = ["gene_ids", "tx_ids", "symbols", "tx_biotypes", "gene_biotypes"]
SET_FIELDS = [*TEXT_FIELDS, "entrez_ids"]
FLAGS = ["standard_chromosomes", "require_tx_id", "protein_coding", "canonical"]


def test_feature_types_have_the_names_of_load_annotations():
    assert {feature.name: feature.value for feature in FeatureType} == {
        "UTR3": "3UTR",
        "UTR5": "5UTR",
        "EXONS": "exons",
        "CDS": "cds",
    }
    assert FeatureType("3UTR") is FeatureType.UTR3


def test_per_gene_modes():
    assert [mode.name for mode in PerGene] == ["EACH_TRANSCRIPT", "LONGEST", "UNION"]


def test_standard_chromosomes_are_1_to_22_x_y_with_and_without_chr():
    names = [*map(str, range(1, 23)), "X", "Y"]
    expected = {*names, *(f"chr{name}" for name in names)}
    assert sorted(STANDARD_CHROMOSOMES) == sorted(expected)
    assert not STANDARD_CHROMOSOMES & {"MT", "chrM", "0", "23", "chr0", "SCAF1.1", "x"}


def test_defaults_are_the_defaults_of_load_annotations():
    # R/load_annotations.R:40-70: build_annotation_filter has canonical = FALSE, but
    # load_annotations passes its own default, canonical = TRUE.
    default = TranscriptFilter()
    assert [getattr(default, name) for name in FLAGS] == [True, True, True, True]
    assert default.support_level is None
    assert [getattr(default, name) for name in SET_FIELDS] == [None] * len(SET_FIELDS)


@pytest.mark.parametrize("name", TEXT_FIELDS)
def test_set_field_rejects_a_str(name):
    message = (
        f"{name} 'ENSSYNG00000000001' is a str; expected a sequence of str, "
        "for example ('ENSSYNG00000000001',)"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        TranscriptFilter(**{name: "ENSSYNG00000000001"})


@pytest.mark.parametrize("name", SET_FIELDS)
def test_set_field_rejects_an_empty_sequence(name):
    with pytest.raises(InputError, match=f"{name} is empty; expected at least one value"):
        TranscriptFilter(**{name: ()})


def test_set_field_rejects_a_value_that_is_not_a_sequence():
    with pytest.raises(InputError, match=r"gene_ids 5 has type int; expected a sequence"):
        TranscriptFilter(gene_ids=5)


@pytest.mark.parametrize("name", TEXT_FIELDS)
@pytest.mark.parametrize("item", [1, None, b"G1", float("nan")])
def test_set_field_rejects_an_item_that_is_not_a_str(name, item):
    message = f"{name}[1] {item!r} has type {type(item).__name__}; expected a str"
    with pytest.raises(InputError, match=re.escape(message)):
        TranscriptFilter(**{name: ("ENSSYNG00000000001", item)})


@pytest.mark.parametrize("level", [0, -1])
def test_support_level_below_1_raises_D031(level):
    message = f"support_level {level} is out of range; expected an integer from 1 to 5"
    with pytest.raises(InputError, match=re.escape(message)):
        TranscriptFilter(support_level=level)


def test_support_level_above_5_raises():
    message = "support_level 6 is out of range; expected an integer from 1 to 5"
    with pytest.raises(InputError, match=re.escape(message)):
        TranscriptFilter(support_level=6)


@pytest.mark.parametrize("level", [True, 1.5, "1", np.float64(1.0), [1]])
def test_support_level_must_be_an_integer(level):
    message = f"support_level {level!r} has type {type(level).__name__}; expected an integer"
    with pytest.raises(InputError, match=re.escape(message)):
        TranscriptFilter(support_level=level)


@pytest.mark.parametrize("level", range(1, 6))
def test_support_level_accepts_1_to_5(level):
    assert TranscriptFilter(support_level=level).support_level == level


def test_support_level_from_numpy_is_a_plain_int():
    level = TranscriptFilter(support_level=np.int64(2)).support_level
    assert level == 2
    assert type(level) is int


@pytest.mark.parametrize("entrez_id", [0, -5])
def test_entrez_id_below_1_raises_D031(entrez_id):
    message = f"entrez_ids[1] {entrez_id} is below 1; expected an Entrez ID of 1 or more"
    with pytest.raises(InputError, match=re.escape(message)):
        TranscriptFilter(entrez_ids=(1006, entrez_id))


def test_entrez_ids_reject_a_str_with_an_example_of_integers():
    message = "entrez_ids '1006' is a str; expected a sequence of int, for example (1006,)"
    with pytest.raises(InputError, match=re.escape(message)):
        TranscriptFilter(entrez_ids="1006")


@pytest.mark.parametrize("item", ["1006", True, 1006.0, None])
def test_entrez_ids_reject_an_item_that_is_not_an_integer(item):
    message = f"entrez_ids[1] {item!r} has type {type(item).__name__}; expected an integer"
    with pytest.raises(InputError, match=re.escape(message)):
        TranscriptFilter(entrez_ids=(1006, item))


@pytest.mark.parametrize("name", FLAGS)
@pytest.mark.parametrize("value", ["yes", 1, None])
def test_flag_rejects_a_value_that_is_not_bool(name, value):
    message = f"{name} {value!r} has type {type(value).__name__}; expected a bool"
    with pytest.raises(InputError, match=re.escape(message)):
        TranscriptFilter(**{name: value})


def test_flag_from_numpy_is_a_plain_bool():
    value = TranscriptFilter(canonical=np.False_).canonical
    assert value is False


@pytest.mark.parametrize("name", TEXT_FIELDS)
def test_text_set_is_stored_as_a_tuple(name):
    assert getattr(TranscriptFilter(**{name: ["a", "b", "a"]}), name) == ("a", "b", "a")


def test_ids_of_a_dataframe_column_make_a_tuple_of_plain_values():
    de = pd.DataFrame({"gene_id": ["G1", "G2"], "entrez": np.array([5, 6], dtype=np.int64)})
    transcript_filter = TranscriptFilter(gene_ids=de["gene_id"], entrez_ids=de["entrez"])
    assert transcript_filter.gene_ids == ("G1", "G2")
    assert transcript_filter.entrez_ids == (5, 6)
    assert {type(value) for value in transcript_filter.entrez_ids} == {int}


def test_filters_with_the_same_values_are_equal_and_hashable():
    first = TranscriptFilter(gene_ids=["a", "b"], entrez_ids=[1])
    second = TranscriptFilter(gene_ids=("a", "b"), entrez_ids=(1,))
    assert first == second
    assert hash(first) == hash(second)


def test_filter_is_frozen():
    with pytest.raises(dataclasses.FrozenInstanceError):
        TranscriptFilter().canonical = False
