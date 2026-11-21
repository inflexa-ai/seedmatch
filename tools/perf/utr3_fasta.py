"""Write the canonical 3' UTRs of rnor7 to a FASTA file, with the gene ID as the record ID.

R names each sequence of `load_annotations` with its gene ID, and `SeedMatchReport` joins the
sequences with the DE table by these names. Thus the R script and the Python script read this file.

Run it from the repository root:
uv run python tools/perf/utr3_fasta.py "$SEEDMATCH_TEST_REFS/rnor7" tmp/perf/rnor7_utr3.fa
"""

import argparse
import dataclasses
from pathlib import Path

from seedmatch import FeatureType, Reference, write_fasta


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path, help="the directory of the EnsDb and the 2bit file")
    parser.add_argument("out", type=Path, help="the FASTA file to write")
    args = parser.parse_args()
    with Reference.from_directory(args.reference) as ref:
        utr3 = ref.sequences(FeatureType.UTR3)
    # The canonical transcripts give one sequence for each gene, thus each record ID is unique.
    if len(set(utr3.gene_ids)) != len(utr3):
        message = f"{args.reference}: a gene has more than one canonical 3' UTR"
        raise SystemExit(message)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(dataclasses.replace(utr3, ids=utr3.gene_ids), args.out)
    print(f"{args.out}: {len(utr3)} sequences, {int(utr3.lengths().sum())} bases")


if __name__ == "__main__":
    main()
