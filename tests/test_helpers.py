import gzip

import numpy as np
import pandas as pd
import pytest
from helpers import R_COLUMNS, assert_frame_matches_golden, golden


def write(tmp_path, text):
    path = tmp_path / "case.tsv"
    path.write_text("# generated-by: test\n" + text)
    return str(path)


def test_r_columns_rename_the_r_names_to_the_python_names():
    assert R_COLUMNS["seqs"] == R_COLUMNS["seqnames"] == "seq_id"
    assert R_COLUMNS["best.match"] == "best_match"
    assert R_COLUMNS["IN.DESEQ2"] == "detected"


def test_matching_frame_passes_after_the_rename_and_the_removal_of_names(tmp_path):
    path = write(
        tmp_path,
        "seqs\tnames\tGroup\tbest.match\nt1\tt1\ts\tFull-D0\nt2\tt2\ts\tNo.matches\n",
    )
    got = pd.DataFrame(
        {
            "seq_id": ["t1", "t2"],
            "gene_id": ["g1", "g2"],
            "group": ["s", "s"],
            "best_match": ["Full-D0", "No.matches"],
        }
    )
    assert_frame_matches_golden(got, path)


def test_row_count_difference_fails_with_both_counts(tmp_path):
    path = write(tmp_path, "a\n1\n2\n")
    with pytest.raises(AssertionError, match=r"3 rows.*expected 2"):
        assert_frame_matches_golden(pd.DataFrame({"a": [1, 2, 3]}), path)


def test_golden_column_absent_from_the_frame_fails_with_the_column_name(tmp_path):
    path = write(tmp_path, "a\tTargetSeq\n1\tACG\n")
    with pytest.raises(AssertionError, match=r"no column 'target_seq'"):
        assert_frame_matches_golden(pd.DataFrame({"a": [1]}), path)


def test_first_difference_names_its_row_and_column(tmp_path):
    path = write(tmp_path, "a\tb\n1\t5\n2\t6\n3\t7\n")
    got = pd.DataFrame({"a": [1, 2, 4], "b": [5, 9, 7]})
    with pytest.raises(AssertionError, match=r"row 1 .*column 'b'.*got 9.*expected 6"):
        assert_frame_matches_golden(got, path)


@pytest.mark.parametrize(
    ("value", "passes"),
    [(0.1, True), (0.1 * (1 + 1e-10), True), (0.1 * (1 + 1e-8), False), (0.0, False)],
)
def test_floats_agree_within_a_relative_tolerance(tmp_path, value, passes):
    path = write(tmp_path, "x\n0.10000000000000001\n")
    if passes:
        assert_frame_matches_golden(pd.DataFrame({"x": [value]}), path)
    else:
        with pytest.raises(AssertionError, match=r"column 'x'"):
            assert_frame_matches_golden(pd.DataFrame({"x": [value]}), path)


def test_integers_must_be_equal(tmp_path):
    path = write(tmp_path, "n\n7\n")
    with pytest.raises(AssertionError, match=r"column 'n'"):
        assert_frame_matches_golden(pd.DataFrame({"n": np.array([8], dtype=np.int64)}), path)


def test_r_booleans_and_detection_labels_read_as_booleans(tmp_path):
    path = write(tmp_path, "AllowIndels\tIN.DESEQ2\nTRUE\tDetected\nFALSE\tUndetected\n")
    good = pd.DataFrame({"indels": [True, False], "detected": [True, False]})
    assert_frame_matches_golden(good, path)
    bad = pd.DataFrame({"indels": [True, True], "detected": [True, False]})
    with pytest.raises(AssertionError, match=r"row 1 .*column 'indels'"):
        assert_frame_matches_golden(bad, path)


def test_na_in_the_golden_is_a_missing_value(tmp_path):
    path = write(tmp_path, "x\ty\n1.5\tNA\nNA\tt\n")
    got = pd.DataFrame({"x": [1.5, np.nan], "y": [None, "t"]})
    assert_frame_matches_golden(got, path)
    with pytest.raises(AssertionError, match=r"row 1 .*column 'x'"):
        assert_frame_matches_golden(pd.DataFrame({"x": [1.5, 2.0], "y": [None, "t"]}), path)
    with pytest.raises(AssertionError, match=r"row 0 .*column 'x'"):
        assert_frame_matches_golden(pd.DataFrame({"x": [np.nan, np.nan], "y": [None, "t"]}), path)


def test_a_string_must_be_equal(tmp_path):
    path = write(tmp_path, "s\nFull-D0\n")
    with pytest.raises(AssertionError, match=r"column 's'.*got 'Full-D1'.*expected 'Full-D0'"):
        assert_frame_matches_golden(pd.DataFrame({"s": ["Full-D1"]}), path)


def test_a_hash_inside_a_data_line_is_data(tmp_path):
    path = write(tmp_path, "case\tcolor\nx\t#FF0000\n")
    assert golden(path)["color"].tolist() == ["#FF0000"]


def test_the_colors_of_the_benchmark_plot_golden_start_with_a_hash():
    table = golden("plots/benchmark_plot.tsv")
    first_roc = table[(table["case"] == "benchmark_1") & (table["curve"] == "ROC")].iloc[0]
    assert first_roc["color"].startswith("#")


def test_a_gzipped_golden_is_decompressed_and_keeps_a_hash_in_a_data_line(tmp_path):
    path = tmp_path / "case.tsv.gz"
    path.write_bytes(gzip.compress(b"# generated-by: test\ncase\tcolor\nx\t#FF0000\n"))
    assert golden(str(path))["color"].tolist() == ["#FF0000"]
