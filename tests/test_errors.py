import pytest
from helpers import sirna_de

from seedmatch import SeedMatchWarning, ecdf_test
from seedmatch.plot import plot_ecdf


def test_a_warning_names_the_first_line_outside_the_package_for_each_call_depth():
    de = sirna_de()
    seed = de.loc[de["mer7m8"] >= 1, "gene_id"]
    background = de.loc[de["mer7m8"] == 0, "gene_id"]
    result = ecdf_test(de, target=seed, background=background)
    # assign_groups warns for the shared gene, and plot_ecdf calls it.
    shared = {"Background": [*background, seed.iloc[0]], "mer7m8": seed}
    with pytest.warns(SeedMatchWarning, match="share the gene") as record:
        plot_ecdf(de, shared, result)
    assert record[0].filename == __file__
