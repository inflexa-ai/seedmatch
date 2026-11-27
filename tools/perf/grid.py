"""Time each search of the `full_search` grid, and hash its ranges.

Each line is `name<TAB>indels=on|off<TAB>seconds<TAB>sha256`. The time is the best of `--repeat`
calls of `count_matches` after one call to warm up. The hash is the SHA-256 of the bytes of the
arrays of `match_ranges`. A `total` line ends each grid.

Run it from the repository root:
uv run python tools/perf/grid.py tmp/perf/rnor7_utr3.fa
"""

import argparse
import hashlib
import time
from pathlib import Path

from seedmatch import Guide, read_fasta
from seedmatch.report.grid import FullSearchOptions, _searches
from seedmatch.search import count_matches, match_ranges

GUIDE = "UUAUAGAGCAAGAACACUGUUUU"
GRIDS = (FullSearchOptions(), FullSearchOptions(indels=True))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", type=Path)
    parser.add_argument("--repeat", type=int, default=3)
    args = parser.parse_args()
    if args.repeat < 1:
        parser.error("--repeat must be 1 or more")

    sequences = read_fasta(args.fasta)
    guide = Guide(GUIDE)
    for options in GRIDS:
        indels = "on" if options.indels else "off"
        total = 0.0
        for search in _searches(guide, options):
            count_matches(search.target, sequences, search.options)
            best = float("inf")
            for _ in range(args.repeat):
                started = time.perf_counter()
                count_matches(search.target, sequences, search.options)
                best = min(best, time.perf_counter() - started)
            total += best
            index, starts, ends = match_ranges(search.target, sequences, search.options)
            digest = hashlib.sha256(index.tobytes() + starts.tobytes() + ends.tobytes())
            print(f"{search.name}\tindels={indels}\t{best:.4f}\t{digest.hexdigest()}", flush=True)
        print(f"total\tindels={indels}\t{total:.3f}", flush=True)


if __name__ == "__main__":
    main()
