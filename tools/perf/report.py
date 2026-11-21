"""Measure the wall time and the peak memory of `seed_match_report`.

The input is the FASTA file of `tools/perf/utr3_fasta.py` and the DE table of the R goldens
(`tests/data/DESEQ2_results.tsv`). Set `NUMBA_CACHE_DIR` to an empty directory for a cold run.

Run it from the repository root:
NUMBA_CACHE_DIR=tmp/perf/numba uv run python tools/perf/report.py tmp/perf/rnor7_utr3.fa
"""

import argparse
import resource
import sys
import time
from importlib.metadata import version
from pathlib import Path

started = time.perf_counter()
import pandas as pd  # noqa: E402

from seedmatch import Guide, read_fasta, report_html, seed_match_report  # noqa: E402

IMPORT_SECONDS = time.perf_counter() - started
GUIDE = "UUAUAGAGCAAGAACACUGUUUU"
DE_TABLE = Path("tests/data/DESEQ2_results.tsv")
VERSIONS = ("seedmatch", "numba", "llvmlite", "numpy", "pandas")


def peak_mib() -> float:
    """Return the peak resident memory of the process. macOS gives bytes, and Linux gives KiB."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / 2**20 if sys.platform == "darwin" else peak / 2**10


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", type=Path)
    parser.add_argument("--calls", type=int, default=3)
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()
    if args.calls < 1:
        parser.error("--calls must be 1 or more")

    print(f"import seedmatch: {IMPORT_SECONDS:.2f} s")
    sequences = read_fasta(args.fasta)
    # The file has a ninth column, mer7m8, thus read_deseq2_table rejects it (D-021).
    de = pd.read_csv(DE_TABLE, sep="\t")
    guide = Guide(GUIDE)
    print(f"{len(sequences)} sequences; peak before the first call: {peak_mib():.0f} MiB")
    for call in range(1, args.calls + 1):
        started = time.perf_counter()
        result = seed_match_report(de, sequences, guide, workers=args.workers)
        seconds = time.perf_counter() - started
        print(f"call {call}: {seconds:.2f} s, peak {peak_mib():.0f} MiB")
    started = time.perf_counter()
    report_html(result, guide)
    print(f"report_html: {time.perf_counter() - started:.3f} s")
    total = result.report.set_index("category").loc["In silico predictions", "Total"]
    print(f"In silico predictions, Total: {total:.0f}")
    versions = ", ".join(f"{name} {version(name)}" for name in VERSIONS)
    print(f"{versions}, Python {sys.version.split()[0]}")


if __name__ == "__main__":
    main()
