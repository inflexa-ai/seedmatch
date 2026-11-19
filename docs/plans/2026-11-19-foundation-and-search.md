# Foundation and search core: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port the SeedMatchR search, full-search, and report paths to `seedmatch`, with output identical to R except for the recorded divergences.

**Architecture:** A `src/` uv project. `SequenceSet` keeps all bases as IUPAC 4-bit masks in one `uint8` buffer. Numba kernels do the mismatch search and the indel search on that buffer. A pure-Python oracle and an R corpus make sure that the kernels obey the Biostrings rules. Higher functions return pandas DataFrames.

**Tech Stack:** Python 3.12, uv, uv_build, ruff, mypy, numpy, pandas, numba, pytest, hypothesis, pytest-socket, pytest-cov. R 4.5 with Biostrings 2.78 and SeedMatchR 2.0.0 for the goldens only.

**Spec:** `docs/design.md`

**Later plans (not in this plan):**

- statistics and the ECDF test
- annotation (2bit, EnsDb, GTF, `Reference`)
- resources and fetch
- the other analyses
- plots, the report HTML, and the documentation set

## Global Constraints

- `requires-python = ">=3.12"`. `.python-version` is `3.12`.
- Runtime dependencies: `numpy>=2.0`, `pandas>=2.2`, `numba>=0.60`. No other runtime dependency in this plan. The sandbox locks numpy 2.5.2, pandas 2.3.3, scipy 1.16.3, numba 0.67.0.
- License: MIT, the license of SeedMatchR. `LICENSE` keeps the SeedMatchR copyright: "Copyright (c) 2025 Eli Lilly and Company".
- No core module imports a network module or `seedmatch.fetch` (ruff `TID251`).
- No function writes a file unless it has an explicit output path parameter.
- Coordinates are 1-based and inclusive.
- Each change of behavior from R gets an entry `## D-NNN: <title>` in `docs/divergences.md`, and at least one test whose name or docstring contains the ID.
- The IDs `D-001` to `D-014` follow the row order of the table in `docs/design.md`, "Divergences (defect fixes)". A new divergence gets the next free ID, from `D-015`.
- A rename of a column or a change of a value type is not a divergence. `docs/r-map.md` records it.
- Errors derive from `seedmatch.errors.SeedMatchError`. Invalid input raises `InputError`. Each message names the bad value and the expected value.
- Documents use ASD-STE100. Code comments obey the admission test of the project rules.
- The word "seam" does not occur in names, comments, or documents.
- Git: `git init`, each `git add`, and each commit need the consent of the user. A commit uses `git commit -s`. Each task ends with a checkpoint, not with a commit.
- The R source of truth is `../SeedMatchR` at commit `2fc0c4f`. Golden files record the R version, the Biostrings version, and the SeedMatchR version and commit.

## Review Focus

1. Lowercase or soft-masked subject bases: they must match as uppercase. A `U` in a subject must raise `InputError`, not match `T`. Owner: Task 3.
2. Runs of `N` in a subject with `fixed_subject=False` (the `full_search` default): each `N` matches each base, as in R. Owner: Task 7.
3. An empty `SequenceSet`, an empty sequence, and a sequence shorter than the pattern: no crash. A short sequence can still have a mismatch match that extends past its ends. Owner: Task 6.
4. A DE table with genes that have no sequence, and a gene with more than one sequence: counts are 0 for an absent gene, and the maximum for a gene with more than one sequence. Owner: Task 8.
5. A guide that is too short for a default `full_search` seed (for example, 18 nt with `18mer`): `InputError` that names the seed and the guide length, before any search. Owner: Task 9.

---

### Task 1: Project scaffold and tooling

**Files:**
- Make: `pyproject.toml`, `.python-version`, `.gitignore`, `LICENSE`, `README.md`
- Make: `src/seedmatch/__init__.py`, `src/seedmatch/errors.py`, `src/seedmatch/py.typed`
- Make: `docs/divergences.md`, `docs/r-map.md`
- Test: `tests/test_package.py`, `tests/helpers.py`

**Interfaces:**
- Produces: `seedmatch.errors.SeedMatchError(Exception)`, `seedmatch.errors.InputError(SeedMatchError, ValueError)`.
- Produces: `tests/helpers.py`, which the tests import as `helpers` (pytest `pythonpath = ["tests"]`). It is not `conftest.py`, because pytest collects each imported function with a `test_` name. It has:
  - `DATA = Path(__file__).parent / "data"`, `GOLDEN = Path(__file__).parent / "golden"`, `GUIDE = "UUAUAGAGCAAGAACACUGUUUU"`
  - `golden(path: str) -> pandas.DataFrame`, which reads `GOLDEN / path` with `sep="\t"`, `comment="#"`, `dtype=str`, `keep_default_na=False`
  - `golden_rows(path: str) -> list[Any]`, the rows of `golden(path)` as named tuples

- [ ] **Step 1: Ask the user for consent to run `git init` in `../seedmatch`.** Do not continue the git part without it.

- [ ] **Step 2: Write `pyproject.toml`**

Values that the spec pins:
- `[project]`: `name = "seedmatch"`, `version = "0.1.0"`, `requires-python = ">=3.12"`, `license = "MIT"`, `license-files = ["LICENSE"]`.
- `dependencies = ["numpy>=2.0", "pandas>=2.2", "numba>=0.60"]`.
- `[dependency-groups] dev = ["pytest>=8.3", "pytest-cov>=6.0", "hypothesis>=6.100", "pytest-socket>=0.7", "mypy>=1.13", "pandas-stubs>=2.2", "ruff>=0.15"]`.
- `[build-system] requires = ["uv_build>=0.9,<0.10"]`, `build-backend = "uv_build"`.
- `[tool.ruff]`: `line-length = 100`, `target-version = "py312"`. `[tool.ruff.lint] select = ["ALL"]`. Ignore `D203`, `D213` (they conflict with `D211`, `D212`), and `COM812`, `ISC001` (they conflict with the formatter). Write the reason next to each ignore.
- `[tool.ruff.lint.flake8-tidy-imports.banned-api]`: ban `urllib.request`, `http.client`, `socket`, `ssl`, `urllib3`, `requests`, `httpx`, and `seedmatch.fetch`, each with `msg = "Network access belongs only in seedmatch.fetch."`.
- `[tool.ruff.lint.per-file-ignores]`: `"src/seedmatch/fetch/**" = ["TID251"]`. `"tests/**" = ["S101", "PLR2004", "D", "ANN", "INP001"]`. `"tools/**" = ["T201", "INP001", "D"]`.
- `[tool.mypy]`: `strict = true`, `files = ["src"]`. An override for `numba.*` with `ignore_missing_imports = true`. An override for `seedmatch.search.kernels` with `disallow_untyped_decorators = false`.
- `[tool.pytest.ini_options]`: `testpaths = ["tests"]`, `pythonpath = ["tests"]`, `addopts = "--disable-socket --strict-markers -ra"`, `markers = ["refs: needs SEEDMATCH_TEST_REFS"]`.
- `[tool.coverage.run]`: `branch = true`, `source = ["seedmatch"]`, `omit = ["*/seedmatch/plot/*"]`. `[tool.coverage.report]`: `fail_under = 95`.

- [ ] **Step 3: Write the other files**

- `.python-version`: `3.12`.
- `.gitignore`: `.venv/`, `__pycache__/`, `.coverage`, `htmlcov/`, `.pytest_cache/`, `.mypy_cache/`, `.ruff_cache/`, `.hypothesis/`, `tools/golden/.rlib/`, `tmp/`.
- `LICENSE`: the MIT License text from `../SeedMatchR/LICENSE.md`, with the copyright lines `Copyright (c) 2026 Inflexa, Inc.` and `Copyright (c) 2025 Eli Lilly and Company`.
- `README.md`: one paragraph on the purpose, then the commands: `uv sync`, `uv run ruff check`, `uv run ruff format --check`, `uv run mypy`, `uv run pytest`, and `NUMBA_DISABLE_JIT=1 uv run pytest --cov` for coverage.
- `docs/divergences.md`: a title, the field list of the spec ("Divergences (defect fixes)"), and no entries yet.
- `docs/r-map.md`: a title and an empty table with the columns "R function", "Python", "Notes".

- [ ] **Step 4: Write the failing tests in `tests/test_package.py`**

```python
def test_import_loads_no_fetch_and_no_plot_modules():
    out = subprocess.run(
        [sys.executable, "-c",
         "import sys, seedmatch; print(sorted(m for m in sys.modules "
         "if m.startswith(('seedmatch.fetch', 'matplotlib', 'pysam'))))"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    assert out == "[]"

def test_each_divergence_id_is_unique_and_has_a_test():
    ids = re.findall(r"^## (D-\d{3}):", DIVERGENCES.read_text(), flags=re.MULTILINE)
    assert len(ids) == len(set(ids))
    test_text = "\n".join(p.read_text() for p in TESTS.rglob("test_*.py"))
    assert [i for i in ids if i not in test_text] == []
```

`DIVERGENCES` is `Path(__file__).parents[1] / "docs/divergences.md"`. `TESTS` is `Path(__file__).parent`.
- [ ] **Step 5: Run `uv sync`, then `uv run pytest tests/test_package.py -v`.** Expected: PASS for the two tests (the package and the record exist, and the record has no entries yet). Then make a temporary entry `## D-999: probe` in `docs/divergences.md`, run the test again, and make sure that it fails with `D-999` in the output. Remove the entry.

- [ ] **Step 6: Run `uv run ruff check`, `uv run ruff format --check`, `uv run mypy`.** Expected: no errors.

- [ ] **Step 7: Checkpoint.** Tell the user that the scaffold is ready to commit.

### Task 2: R golden environment

**Files:**
- Make: `tools/golden/install.R`, `tools/golden/run.sh`, `tools/golden/R/common.R`, `tools/golden/R/check_env.R`, `tools/golden/README.md`

**Interfaces:**
- Produces: `tools/golden/run.sh <script.R>` runs one R script from the repository root with `R_LIBS_USER=tools/golden/.rlib`.
- Produces: `common.R::write_golden(df, path, script)` writes `tests/golden/<path>` as TSV with no row names and no quotes. The first lines are comments: `# generated-by: tools/golden/R/<script>`, `# R: <version>`, `# Biostrings: <version>`, `# SeedMatchR: <version> (2fc0c4f)`.
- Produces: `common.R::load_test_data()` returns `list(res = <DESEQ2_results.tsv as data.frame>, seqs = <Test_seqs.fa as DNAStringSet>, guide = "UUAUAGAGCAAGAACACUGUUUU")` from `tests/data/`.

- [ ] **Step 1: Write `install.R`.** It stops if `git -C ../SeedMatchR rev-parse HEAD` is not `2fc0c4f` (prefix match) or if `git -C ../SeedMatchR status --porcelain` lists a tracked change. It installs `BiocManager`, then the `Imports` of `../SeedMatchR/DESCRIPTION` with `BiocManager::install(..., lib = "tools/golden/.rlib", update = FALSE, ask = FALSE)`, then `install.packages("../SeedMatchR", repos = NULL, type = "source", lib = "tools/golden/.rlib")`.

- [ ] **Step 2: Write `run.sh`, `common.R`, and `check_env.R`.** `check_env.R` runs `SeedMatchReport(res = res, seqs = seqs, guide.seq = guide)$report` on the test data. It stops unless the row `In silico predictions` has `Total` 125 and `15mer-D4` 73. These values come from `../SeedMatchR/tests/testthat/_snaps/SeedMatchReport.md`.

- [ ] **Step 3: Copy the test data.** Copy `../SeedMatchR/tests/testthat/testdata/DESEQ2_results.tsv` and `Test_seqs.fa` to `tests/data/`. Write `tests/data/README.md` with the source path, the commit `2fc0c4f`, and the MIT notice.

- [ ] **Step 4: Run `Rscript tools/golden/install.R`, then `tools/golden/run.sh tools/golden/R/check_env.R`.** Expected: the last line is `golden environment OK: SeedMatchR 2.0.0, Biostrings 2.78.x`.

- [ ] **Step 5: Write `tools/golden/README.md`** with the two commands and the rule that tests never run R.

- [ ] **Step 6: Checkpoint.**

### Task 3: Alphabet, `SequenceSet`, and FASTA

**Files:**
- Make: `src/seedmatch/sequence/__init__.py`, `alphabet.py`, `sequence_set.py`, `fasta.py`
- Modify: `tests/helpers.py`
- Test: `tests/sequence/test_alphabet.py`, `tests/sequence/test_sequence_set.py`

**Interfaces:**
- Produces: `alphabet.MASKS: Final[Mapping[str, int]]` with `A=1, C=2, G=4, T=8`, and each IUPAC code as the union: `M=A|C, R=A|G, W=A|T, S=C|G, Y=C|T, K=G|T, V=A|C|G, H=A|C|T, D=A|G|T, B=C|G|T, N=15`.
- Produces: `alphabet.encode_dna(text: str) -> NDArray[np.uint8]`, `alphabet.decode_dna(codes: NDArray[np.uint8]) -> str`, `alphabet.reverse_complement(text: str) -> str`.
- Produces: `SequenceSet` with `SequenceSet.from_records(records: Iterable[tuple[str, str]], gene_ids: Sequence[str] | None = None) -> SequenceSet`, the attributes `ids: tuple[str, ...]`, `gene_ids: tuple[str, ...]`, `codes: NDArray[np.uint8]`, `offsets: NDArray[np.int64]` (length `n + 1`), and the methods `__len__() -> int`, `sequence(i: int) -> str`, `lengths() -> NDArray[np.int64]`. If `gene_ids` is `None`, `gene_ids` equals `ids`.
- Produces: `fasta.read_fasta(path: str | os.PathLike[str]) -> SequenceSet`. It reads the file line by line. The ID is the header text up to the first white space.
- Produces: `fasta.write_fasta(sequences: SequenceSet, path: str | os.PathLike[str], *, width: int = 80) -> None`. The width copies the default of `Biostrings::writeXStringSet`.
- Produces: `helpers.sirna_seqs() -> SequenceSet`, which is `read_fasta(DATA / "Test_seqs.fa")`.

- [ ] **Step 1: Write the failing tests**

```python
def test_reverse_complement_of_wobble_seed_matches_biostrings():
    assert reverse_complement("YYAYAGAG") == "CTCTRTRR"

@given(st.text(alphabet="ACGTMRWSYKVHDBN"))
def test_reverse_complement_is_an_involution(s):
    assert reverse_complement(reverse_complement(s)) == s

def test_lowercase_bases_encode_as_uppercase():
    assert np.array_equal(encode_dna("acgtn"), encode_dna("ACGTN"))

def test_uracil_in_dna_is_rejected_with_its_position():
    with pytest.raises(InputError, match=r"'U' at position 3"):
        encode_dna("ACUG")

def test_read_fasta_keeps_each_record_of_the_test_data():
    seqs = read_fasta(DATA / "Test_seqs.fa")
    expected = naive_fasta_records(DATA / "Test_seqs.fa")  # ">" lines and joined bodies
    assert list(seqs.ids) == [name for name, _ in expected]
    assert [seqs.sequence(i) for i in range(len(seqs))] == [body for _, body in expected]
    assert seqs.gene_ids == seqs.ids

def test_empty_set_has_no_sequences():
    s = SequenceSet.from_records([])
    assert len(s) == 0 and s.offsets.tolist() == [0]

def test_write_then_read_gives_the_same_set(tmp_path):
    seqs = SequenceSet.from_records([("a", "ACGT" * 50), ("b", ""), ("c", "NNRY")])
    write_fasta(seqs, tmp_path / "x.fa", width=7)
    back = read_fasta(tmp_path / "x.fa")
    assert back.ids == seqs.ids and np.array_equal(back.codes, seqs.codes)
    assert max(len(line) for line in (tmp_path / "x.fa").read_text().splitlines()) == 7
```

`naive_fasta_records` is a ten-line helper in the test file.

- [ ] **Step 2: Run `uv run pytest tests/sequence -v`.** Expected: FAIL with `ModuleNotFoundError: seedmatch.sequence`.

- [ ] **Step 3: Write the modules.** `encode_dna` uses a 256-entry lookup table with `0` for an invalid byte. The error message gives the first invalid letter and its 1-based position. The complement swaps `A/T`, `C/G`, `M/K`, `R/Y`, `V/B`, `H/D`, and keeps `W`, `S`, `N`.

- [ ] **Step 4: Run `uv run pytest tests/sequence -v`.** Expected: PASS.

- [ ] **Step 5: Add the rows for `Biostrings::readDNAStringSet` and `reverseComplement` to `docs/r-map.md`. Checkpoint.**

### Task 4: Guide, seeds, and targets

**Files:**
- Make: `src/seedmatch/seed.py`, `tools/golden/R/seed_targets.R`
- Make: `tests/golden/seed/targets.tsv` (made by the R script)
- Test: `tests/test_seed.py`

**Interfaces:**
- Produces: `Seed(StrEnum)` with `MER7M8 = "mer7m8"`, `MER8 = "mer8"`, `MER6 = "mer6"`, `MER7A1 = "mer7A1"`, `FULL = "Full"`, `MER18 = "18mer"`, `MER15 = "15mer"`.
- Produces: `CustomSeed(name: str, start: int, stop: int)` (frozen dataclass), and `SeedSpec = Seed | CustomSeed`.
- Produces: `Target(name: str, sequence: str, seed_rna: str)` (frozen dataclass). `name` is `seed.value` for a `Seed` and `seed.name` for a `CustomSeed`. `sequence` is the DNA search pattern. `seed_rna` is the guide window as RNA, equal to the reverse complement of `sequence` with `T` as `U`.
- Produces: `Guide(sequence: str)` (frozen dataclass) and `Guide.target(seed: SeedSpec, *, wobble: bool = False) -> Target`.

Target rules (from `R/get_seed.R:92-177`). `g[i..j]` is the 1-based guide window. `rc` is the reverse complement as DNA. With `wobble=True`, each `U` of the guide becomes `Y` before the window is cut.

| Seed | Target |
|---|---|
| `FULL` | `rc(g[2..len])` |
| `MER18` | `rc(g[2..19])` |
| `MER15` | `rc(g[2..16])` |
| `MER7M8` | `rc(g[2..8])` |
| `MER6` | `rc(g[2..7])` |
| `MER8` | `rc(g[2..8]) + "A"` |
| `MER7A1` | `rc(g[2..7]) + "A"` |
| `CustomSeed(s, e)` | `rc(g[s..e])` |

- [ ] **Step 1: Write `seed_targets.R`.** It writes one row for each combination of these guides, each canonical seed, and each `wobble` value. Guides: `UUAUAGAGCAAGAACACUGUUUU`, `UGAGGUAGUAGGUUGUAUAGUU`, `AUUUUUUUUUUUUUUUUUUUA`. It also writes the custom windows `(2, 19)` and `(1, 23)` of the first guide. Columns: `guide`, `seed`, `start`, `stop`, `wobble`, `target`, `seed_rna`. A canonical row has empty `start` and `stop`. Run it with `tools/golden/run.sh`.

- [ ] **Step 2: Write the failing tests**

```python
@pytest.mark.parametrize("row", golden_rows("seed/targets.tsv"), ids=row_id)
def test_target_matches_r(row):
    seed = CustomSeed("custom", int(row.start), int(row.stop)) if row.start else Seed(row.seed)
    t = Guide(row.guide).target(seed, wobble=row.wobble == "TRUE")
    assert (t.sequence, t.seed_rna) == (row.target, row.seed_rna)

def test_guide_rejects_thymine_and_names_uracil():
    with pytest.raises(InputError, match="'T' at position 1.*use 'U'"):
        Guide("TTAUAGAGCAAGAACACUGUUUU")

def test_canonical_seed_needs_eight_nucleotides():
    with pytest.raises(InputError, match="at least 8"):
        Guide("UUAUAGA").target(Seed.MER7M8)

def test_custom_window_must_lie_inside_the_guide():
    with pytest.raises(InputError, match=r"stop 24 .* length 23"):
        Guide("UUAUAGAGCAAGAACACUGUUUU").target(CustomSeed("x", 2, 24))
```

- [ ] **Step 3: Run `uv run pytest tests/test_seed.py -v`.** Expected: FAIL with `ModuleNotFoundError`.

- [ ] **Step 4: Write `seed.py`.** `Guide` accepts `A`, `C`, `G`, `U` in either case and stores uppercase. Each canonical seed needs a guide of at least 8 nt. `MER15` needs at least 16 nt, and `MER18` needs at least 19 nt. The message is `f"{name} needs a guide of at least {n} nt, got {length}"`. `seed_rna` copies `Seed.Seq.RNA` of `R/get_seed.R:69`.

- [ ] **Step 5: Run the tests.** Expected: PASS. Add the `get_seed` row to `docs/r-map.md`. Checkpoint.

### Task 5: Compatibility table, oracle, and R corpus

**Files:**
- Make: `src/seedmatch/search/__init__.py`, `src/seedmatch/search/compat.py`
- Make: `tools/golden/R/iupac_table.R`, `tools/golden/make_corpus.py`, `tools/golden/R/search_corpus.R`
- Make (made by the scripts): `tests/golden/search/iupac_compat.tsv`, `corpus_cases.tsv`, `corpus_matches.tsv`, `corpus_counts.tsv`
- Make: `tests/search/oracle.py`
- Test: `tests/search/test_compat_and_oracle.py`

**Interfaces:**
- Produces: `compat.compat_table(fixed_pattern: bool, fixed_subject: bool) -> NDArray[np.bool_]` with shape `(16, 16)`, indexed `[pattern_mask, subject_mask]`.
- Produces (test helper): `oracle.hamming(pattern: str, subject: str, max_mm: int, min_mm: int, fixed_pattern: bool, fixed_subject: bool) -> list[tuple[int, int]]` and `oracle.indels(pattern: str, subject: str, max_mm: int, fixed_pattern: bool, fixed_subject: bool) -> list[tuple[int, int]]`. Each returns `(start, end)` pairs, 1-based, in R order.

The compatibility rule, from an R probe on Biostrings 2.78.0:

| `fixed_pattern` | `fixed_subject` | Match if |
|---|---|---|
| `True` | `True` | `p == s` |
| `True` | `False` | `p & ~s == 0` (the pattern letter is a subset of the subject letter) |
| `False` | `True` | `s & ~p == 0` |
| `False` | `False` | `p & s != 0` |

The oracle rules are rules 2 and 4 of `docs/design.md`, "Search rules". The R corpus showed that Biostrings 2.78.0 does not obey the "best local match" definition of its documentation. Thus, the oracle obeys the Biostrings C source, and the definition stays only as a test helper.

- [ ] **Step 1: Write `iupac_table.R`.** For the 15 IUPAC letters and the four `fixed` pairs, it writes `pattern`, `subject`, `fixed_pattern`, `fixed_subject`, `match` from `countPattern(p, DNAString(s), fixed = c(pattern = fp, subject = fs)) > 0`.

- [ ] **Step 2: Write `make_corpus.py`.** It uses `random.Random(20261119)` and writes 400 cases (`case_id`, `pattern`, `subject`):
  - pattern length 1 to 12, subject length 0 to 60
  - bases mostly `ACGT`, with runs of `N` and a few IUPAC codes in each of the pattern and the subject
  - in most cases, a copy of the pattern with 0 to 3 random edits at a random place, also at the two ends
  - some subjects made of a repeated short unit

- [ ] **Step 3: Write `search_corpus.R`.** For each case and each option set, it records `matchPattern(pattern, DNAString(subject), ...)` ranges and `countPattern(...)`:
  - `max_mm` 0 to 4, `with.indels` `FALSE` and `TRUE`, the four `fixed` pairs
  - `min_mm` 0 and 1, only with `with.indels = FALSE`
  - `corpus_matches.tsv` columns: `case_id`, `max_mm`, `min_mm`, `indels`, `fixed_pattern`, `fixed_subject`, `start`, `end`
  - `corpus_counts.tsv` has the same option columns and `count`
  - skip an option set for which Biostrings stops with an error, and record the skip in a `# skipped:` comment line

- [ ] **Step 4: Run the three scripts.** Expected: the four TSV files exist and are not empty.

- [ ] **Step 5: Write the failing tests**

```python
def test_compat_table_matches_biostrings():
    for row in golden_rows("search/iupac_compat.tsv"):
        table = compat_table(row.fixed_pattern == "TRUE", row.fixed_subject == "TRUE")
        assert table[MASKS[row.pattern], MASKS[row.subject]] == (row.match == "TRUE"), row

def test_oracle_matches_biostrings_on_corpus():
    for case, opts, expected in corpus_groups():  # groups corpus_matches by case and options
        found = oracle_for(opts)(case.pattern, case.subject, **opts.kwargs())
        assert found == expected, (case.case_id, opts)

def test_biostrings_counts_equal_its_range_counts():
    assert corpus_count_mismatches() == []
```

- [ ] **Step 6: Run `uv run pytest tests/search/test_compat_and_oracle.py -v`.** Expected first: FAIL with `ModuleNotFoundError`. After `compat.py` and `oracle.py` exist: PASS. If the indel oracle disagrees with R on some cases, list the cases and their pattern. Change the oracle only to settle the ties that the definition leaves open (for example, the order of the ranges). Record the settled rule in `docs/design.md`, "Search rules", rule 4. If `test_biostrings_counts_equal_its_range_counts` fails, stop and report it to the user, because the API depends on it.

- [ ] **Step 7: Checkpoint.**

### Task 6: Numba search kernels

**Files:**
- Make: `src/seedmatch/search/kernels.py`
- Test: `tests/search/test_kernels.py`

**Interfaces:**
- Consumes: `SequenceSet.codes`, `SequenceSet.offsets` (Task 3), `compat_table` (Task 5).
- Produces, each `@njit(cache=True, nogil=True)`:
  - `hamming_counts(codes, offsets, pattern, compat, max_mm, min_mm) -> NDArray[np.int64]` (one count for each sequence)
  - `hamming_ranges(codes, offsets, pattern, compat, max_mm, min_mm) -> tuple[NDArray[np.int64], NDArray[np.int64]]` (sequence index, 1-based start)
  - `indel_counts(codes, offsets, pattern, compat, max_mm) -> NDArray[np.int64]`
  - `indel_ranges(codes, offsets, pattern, compat, max_mm) -> tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.int64]]` (sequence index, start, end)

The indel kernel does the steps of rule 4 of `docs/design.md`, "Search rules", as `tests/search/oracle.py::indels` does. The hamming kernel obeys rule 2, with the starts from 2 − m to n. Make the kernels equal to the oracle and to the two R corpora (`tests/search/corpus.py`). Stop the DP of a candidate when each cell of a column is greater than `max_mm`.

The two-pass shape (a count pass, then a fill pass) keeps the output arrays exact without growth.

- [ ] **Step 1: Write the failing tests**

```python
CASES = st.tuples(st.text("ACGTN" * 4 + "MRY", min_size=1, max_size=12),
                  st.text("ACGTN" * 6 + "WSK", max_size=60),
                  st.integers(0, 3), st.booleans(), st.booleans())

@given(CASES, st.integers(0, 1))
def test_hamming_kernel_equals_oracle(case, min_mm):
    pattern, subject, k, fp, fs = case
    assume(min_mm <= k)
    assert kernel_hamming(pattern, subject, k, min_mm, fp, fs) == oracle.hamming(pattern, subject, k, min_mm, fp, fs)

@given(CASES)
def test_indel_kernel_equals_oracle(case):
    pattern, subject, k, fp, fs = case
    assert kernel_indels(pattern, subject, k, fp, fs) == oracle.indels(pattern, subject, k, fp, fs)

def test_kernels_equal_biostrings_on_corpus():
    assert corpus_kernel_mismatches() == []

def test_counts_equal_number_of_ranges_for_each_sequence():
    ...  # many sequences in one SequenceSet, counts == bincount of range sequence indexes
```

`kernel_hamming` and `kernel_indels` are test helpers that build a one-sequence `SequenceSet` and call the kernel.

- [ ] **Step 2: Run `uv run pytest tests/search/test_kernels.py -v`.** Expected: FAIL with `ModuleNotFoundError`.

- [ ] **Step 3: Write `kernels.py`.** No Python object enters a kernel. Each kernel loops over the sequences with the offsets.

- [ ] **Step 4: Run the tests, then `NUMBA_DISABLE_JIT=1 uv run pytest tests/search/test_kernels.py`.** Expected: PASS in the two modes.

- [ ] **Step 5: Checkpoint.**

### Task 7: Search API

**Files:**
- Make: `src/seedmatch/search/api.py`
- Modify: `src/seedmatch/search/__init__.py`, `src/seedmatch/__init__.py`, `docs/divergences.md`, `docs/r-map.md`
- Test: `tests/search/test_api.py`

**Interfaces:**
- Consumes: Tasks 3, 4, 5, 6.
- Produces: `SearchOptions(max_mismatch: int = 0, min_mismatch: int = 0, indels: bool = False, fixed_pattern: bool = True, fixed_subject: bool = True)` (frozen dataclass). `__post_init__` raises `InputError` for a negative value, for `min_mismatch > max_mismatch`, and for `indels` with `min_mismatch > 0`.
- Produces: `count_matches(pattern: str | Target, sequences: SequenceSet, options: SearchOptions = SearchOptions()) -> NDArray[np.int64]`.
- Produces: `find_matches(pattern: str | Target, sequences: SequenceSet, options: SearchOptions = SearchOptions()) -> pandas.DataFrame` with the columns `seq_id`, `gene_id`, `start`, `end`, `width`. The rows are in sequence order, then in the R order of Task 5.
- If `indels` is `True` and `max_mismatch == 0`, the API uses the exact search (search rule 5).

- [ ] **Step 1: Write the failing tests.** The expected values come from R probes on Biostrings 2.78.0.

```python
def s(**named):
    return SequenceSet.from_records(named.items())

def test_mismatch_match_can_extend_past_both_ends():
    got = find_matches("ACGT", s(a="CGTAAAAAAAAACGT", b="TTTTTTACG"), SearchOptions(max_mismatch=1))
    assert got[["seq_id", "start", "end"]].values.tolist() == [["a", 0, 3], ["a", 12, 15], ["b", 7, 10]]

def test_indel_match_stays_inside_the_sequence():
    got = find_matches("ACGT", s(b="TTTTTTACG"), SearchOptions(max_mismatch=1, indels=True))
    assert got[["start", "end"]].values.tolist() == [[7, 9]]

def test_ranges_obey_the_indels_flag_D001():
    seqs = s(a="TTTTACGACGTTTT")
    assert find_matches("ACGTACG", seqs, SearchOptions(max_mismatch=1)).empty
    with_indels = find_matches("ACGTACG", seqs, SearchOptions(max_mismatch=1, indels=True))
    assert with_indels[["start", "end"]].values.tolist() == [[5, 10]]

def test_ranges_obey_min_mismatch_D001():
    assert find_matches("ACGAACG", s(a="TTTTACGAACGTTTT"),
                        SearchOptions(max_mismatch=1, min_mismatch=1)).empty

def test_n_in_subject_matches_each_base_when_subject_is_not_fixed():
    seqs = s(a="ACNTGG")
    assert count_matches("ACGT", seqs).tolist() == [0]
    assert count_matches("ACGT", seqs, SearchOptions(fixed_pattern=False, fixed_subject=False)).tolist() == [1]

def test_overlapping_matches_count_separately():
    assert count_matches("AA", s(a="AAAA")).tolist() == [3]

@given(...)  # random SequenceSet and SearchOptions
def test_count_equals_number_of_ranges(seqs, options):
    counts = count_matches("ACGTAC", seqs, options)
    ranges = find_matches("ACGTAC", seqs, options)
    assert counts.tolist() == ranges.groupby("seq_id", sort=False).size().reindex(seqs.ids, fill_value=0).tolist()

@given(...)  # random SequenceSet and SearchOptions
def test_each_range_is_within_the_edit_limit(seqs, options):
    for row in find_matches("ACGTAC", seqs, options).itertuples():
        text = seqs.sequence(seqs.ids.index(row.seq_id))
        assert distance("ACGTAC", text, row.start, row.end, options) <= options.max_mismatch
```

`distance` is a test helper, independent of the oracle. Without indels, it counts the positions outside the sequence and the incompatible positions. With indels, it is a plain Levenshtein distance with the compatibility table.

```python
def test_indels_with_min_mismatch_is_rejected():
    with pytest.raises(InputError, match="min_mismatch must be 0 when indels"):
        SearchOptions(max_mismatch=1, min_mismatch=1, indels=True)

def test_empty_pattern_is_rejected():
    with pytest.raises(InputError, match="empty pattern"):
        count_matches("", s(a="ACGT"))
```

- [ ] **Step 2: Run `uv run pytest tests/search/test_api.py -v`.** Expected: FAIL with `ImportError`.

- [ ] **Step 3: Write `api.py`.** Encode the pattern with `encode_dna`. Pick the kernel from the options.

- [ ] **Step 4: Run the tests.** Expected: PASS.

- [ ] **Step 5: Write `## D-001: Ranges obey with.indels and min.mismatch` in `docs/divergences.md`.** Use the fields of the spec. Give this evidence:
  - `vmatchPattern` stops with "vmatchPattern() does not support indels yet" for a `DNAStringSet`. Thus `vmatchPattern2` calls the C matcher with "indels" for each call (`R/vmatchPattern2.R:43`).
  - `vmatchPattern2("ACGTACG", DNAStringSet(c(a = "TTTTACGACGTTTT")), max.mismatch = 1, with.indels = FALSE)` returns `5-10`. `vmatchPattern` with the same arguments returns no range.
  - `vmatchPattern2` ignores `min.mismatch`.
  - The `iranges` mode of `SeedMatchR` calls `vmatchPattern` with the default `with.indels = TRUE` (`R/SeedMatchR.R:128-133`). Thus it stops with the same error if `max.mismatch > 0`. Add the rows for `vcountPattern`, `vmatchPattern`, and `vmatchPattern2` to `docs/r-map.md`. Checkpoint.

### Task 8: DE table and seed-count annotation

**Files:**
- Make: `src/seedmatch/de.py`, `src/seedmatch/report/__init__.py`, `src/seedmatch/report/annotate.py`
- Make: `tools/golden/R/de_and_counts.R`, `tests/data/de_filter_cases.tsv`
- Make (made by the R script): `tests/golden/de/filter_res.tsv`, `tests/golden/seedmatchr/deseq2_mode.tsv`
- Modify: `docs/divergences.md`, `docs/r-map.md`, `src/seedmatch/__init__.py`, `tests/helpers.py`
- Test: `tests/test_de.py`, `tests/report/test_annotate.py`

**Interfaces:**
- Produces: `helpers.sirna_de() -> pandas.DataFrame`, which reads `DATA / "DESEQ2_results.tsv"` with `pandas.read_csv(sep="\t")`. It keeps the `mer7m8` column of the file, thus it does not use `read_deseq2_table`.
- Produces: `de.DE_REQUIRED_COLUMNS = ("gene_id", "baseMean", "log2FoldChange", "pvalue", "padj")`.
- Produces: `de.read_deseq2_table(path: str | os.PathLike[str]) -> pandas.DataFrame`. It reads a tab-separated table with a header line and exactly 8 columns, also gzip. It names the columns `gene_id`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`, `symbol` (`R/get_data.R:53-55`).
- Produces: `de.filter_de(de: pandas.DataFrame, *, fdr_cutoff: float = 1.0, fc_cutoff: float = 0.0, basemean_cutoff: float = 30.0) -> pandas.DataFrame`. It removes the rows with a missing `padj`, `pvalue`, or `log2FoldChange`. Then it keeps `padj <= fdr_cutoff`, `abs(log2FoldChange) >= fc_cutoff`, and `baseMean >= basemean_cutoff`. It keeps the row order and resets the index.
- Produces: `annotate.annotate_seed_counts(de: pandas.DataFrame, sequences: SequenceSet, target: Target, *, column: str | None = None, options: SearchOptions = SearchOptions(indels=True), shared_genes: bool = True) -> pandas.DataFrame`. The column name defaults to `target.name`. The count of a gene is the maximum count over the sequences of that gene, and 0 for a gene without a sequence. If `shared_genes` is `True`, the result keeps only the genes with a sequence. A column with the same name is replaced.

- [ ] **Step 1: Write `de_filter_cases.tsv`.** It has 12 rows with these cases:
  - one row each with a missing `padj`, `pvalue`, `log2FoldChange`, or `baseMean`
  - `padj` equal to the cutoff
  - `|log2FoldChange|` equal to the cutoff
  - `baseMean` equal to 30, and `baseMean` equal to 29.999
  - a negative `log2FoldChange`

- [ ] **Step 2: Write `de_and_counts.R`.** It writes the `gene_id` values that `filter_res` keeps for the cutoff sets `(1, 0, 30)` and `(0.05, 0.5, 10)`. It also writes `SeedMatchR(seqs = seqs, sequence = guide, res = res, seed.name = s)` for `s` in `mer7m8`, `mer8`, `mer6`, `mer7A1`, `Full`, as `gene_id`, `seed`, `count`. Use named arguments, because the R order is `seqs`, `sequence`, `res`. The count is the last column of the R result. The test data has a `mer7m8` column already, thus R adds a second column with the same name (D-015). Run it.

- [ ] **Step 3: Write the failing tests**

```python
@pytest.mark.parametrize("cutoffs", [(1, 0, 30), (0.05, 0.5, 10)])
def test_filter_de_matches_r(cutoffs):
    got = filter_de(read_cases(), fdr_cutoff=cutoffs[0], fc_cutoff=cutoffs[1], basemean_cutoff=cutoffs[2])
    assert got["gene_id"].tolist() == golden_filter_ids(cutoffs)

def test_filter_de_names_the_missing_columns():
    with pytest.raises(InputError, match=r"missing columns: \['padj', 'pvalue'\]"):
        filter_de(pd.DataFrame({"gene_id": [], "baseMean": [], "log2FoldChange": []}))

def test_read_deseq2_table_names_the_eight_columns(tmp_path):
    ...  # write a gzip table with other header names, then assert the eight DESeq2 names

@pytest.mark.parametrize("seed", ["mer7m8", "mer8", "mer6", "mer7A1", "Full"])
def test_seed_counts_match_r(seed):
    target = Guide(GUIDE).target(Seed(seed))
    got = annotate_seed_counts(sirna_de(), sirna_seqs(), target)
    assert got[["gene_id", seed]].values.tolist() == golden_counts(seed)

def test_gene_count_is_the_maximum_over_its_sequences_D003():
    seqs = SequenceSet.from_records([("t1", "CTCTATA"), ("t2", "CTCTATACTCTATACTCTATA")], gene_ids=["g", "g"])
    de = pd.DataFrame({"gene_id": ["g"]})
    assert annotate_seed_counts(de, seqs, Guide(GUIDE).target(Seed.MER7M8))["mer7m8"].tolist() == [3]

def test_existing_column_is_replaced_not_duplicated_D015():
    got = annotate_seed_counts(sirna_de(), sirna_seqs(), Guide(GUIDE).target(Seed.MER7M8))
    assert list(got.columns).count("mer7m8") == 1

def test_gene_without_sequence_counts_zero_when_not_shared():
    de = pd.DataFrame({"gene_id": ["absent"]})
    got = annotate_seed_counts(de, sirna_seqs(), Guide(GUIDE).target(Seed.MER7M8), shared_genes=False)
    assert got["mer7m8"].tolist() == [0]
```

- [ ] **Step 4: Run `uv run pytest tests/test_de.py tests/report/test_annotate.py -v`.** Expected: FAIL with `ModuleNotFoundError`.

- [ ] **Step 5: Write `de.py` and `annotate.py`.** Compute the gene counts with one `groupby("gene_id").max()` on the sequence counts, then one `map` on the DE table.

- [ ] **Step 6: Run the tests.** Expected: PASS.

- [ ] **Step 7: Write `## D-003: Gene-level results use the maximum over the sequences of a gene` and `## D-015: A seed-count column replaces a column with the same name`.** R evidence for D-003: `R/SeedMatchR.R:207` (`match` takes the first sequence). R evidence for D-015: `R/SeedMatchR.R:214-217`, and the test data, which already has a `mer7m8` column. Add the rows for `filter_res`, `download_parse_file`, and `SeedMatchR` (`DESEQ2` mode) to `docs/r-map.md`. Checkpoint.

### Task 9: Full search, best match, and ranges

**Files:**
- Make: `src/seedmatch/report/grid.py`
- Make: `tools/golden/R/full_search.R`, `tools/golden/patches/vmatchPattern2.R`
- Make (made by the R script): `tests/golden/seedmatchr/full_search_counts_<variant>.tsv`, `full_search_ranges.tsv`, `dataframe_mode.tsv`
- Modify: `docs/divergences.md`, `docs/r-map.md`, `src/seedmatch/__init__.py`, `tests/helpers.py`
- Test: `tests/report/test_full_search.py`

**Interfaces:**
- Consumes: Tasks 4, 7.
- Produces: `DEFAULT_SEEDS: tuple[Seed, ...] = (Seed.FULL, Seed.MER18, Seed.MER15, Seed.MER8, Seed.MER7M8, Seed.MER7A1, Seed.MER6)`. `NONSEED_SEEDS = frozenset({Seed.FULL, Seed.MER18, Seed.MER15})`.
- Produces: `FullSearchOptions(seeds: tuple[Seed, ...] = DEFAULT_SEEDS, nonseed_edit_distances: tuple[int, ...] = (0, 1, 2, 3, 4), seed_edit_distances: tuple[int, ...] = (0, 1), indels: bool = False, wobble: bool = False, fixed_pattern: bool = False, fixed_subject: bool = False)`. The defaults copy `R/full_search.R:26-35`.
- Produces: `full_search(guide: Guide, sequences: SequenceSet, options: FullSearchOptions = FullSearchOptions(), group: str = "siRNA", *, workers: int | None = None) -> pandas.DataFrame` with the columns `seq_id`, `gene_id`, one column `f"{seed.value}.mm{d}"` for each seed and distance in option order, `group`, `best_match`.
- Produces: `best_match(counts: pandas.DataFrame) -> pandas.Series`. The priority order and the labels copy `R/full_search.R:96-125`: `Full-D0` to `Full-D4`, `18mer-D0` to `18mer-D4`, `15mer-D0` to `15mer-D4`, then `mer8.mm0` as `8mer`, `mer7m8.mm0` as `7mer-m8`, `mer7A1.mm0` as `7mer-A1`, `mer6.mm0` as `6mer`, else `No.matches`. An absent column is skipped.
- Produces: `FullSearchRangesOptions(indel_options: tuple[bool, ...] = (True, False), edit_distances: tuple[int, ...] = (0, 1, 2, 3, 4), seeds: tuple[Seed, ...] = DEFAULT_SEEDS, fixed_pattern: bool = False, fixed_subject: bool = False, max_seed_edit_distance: int = 1)`.
- Produces: `full_search_ranges(guide: Guide, sequences: SequenceSet, options: FullSearchRangesOptions = FullSearchRangesOptions(), group: str = "siRNA", *, workers: int | None = None) -> pandas.DataFrame` with the columns `seq_id`, `gene_id`, `start`, `end`, `width`, `query`, `target_seq`, `edit_distance`, `indels`, `group`. `query` is `f"{seed.value}.mm{d}.indel{'TRUE' if indels else 'FALSE'}"`. The rows follow the R `expand.grid` order: `indels` varies fastest, then the seed, then the distance (`R/full_search.R:156-183`). The grid has no row for a seed outside `NONSEED_SEEDS` with a distance above `max_seed_edit_distance` (`R/full_search.R:162`).
- `workers` is the width of the thread pool. `None` means `os.cpu_count()`.

- [ ] **Step 1: Write `patches/vmatchPattern2.R`.** It defines `vmatchPattern2` with the same arguments and the same `data.frame` output (`start`, `end`, `width`, `names`). If `with.indels` is `FALSE`, it calls `Biostrings::vmatchPattern(pattern, subject, max.mismatch, min.mismatch, with.indels = FALSE, fixed)`. If `with.indels` is `TRUE`, it calls `Biostrings::matchPattern` on each sequence, because `vmatchPattern` does not support indels for a `DNAStringSet`. `full_search.R` assigns it into the SeedMatchR namespace with `assignInNamespace` before it makes the range goldens.

- [ ] **Step 2: Write `full_search.R`.** It writes:
  - `full_search(guide, seqs)` counts with the defaults, with `allow_wobbles = TRUE`, and with `indel.bool = TRUE`, as three files (`<variant>` is `default`, `wobble`, `indels`)
  - `full_search_ranges(guide, seqs)` with the patched function
  - the `data.frame` mode of `SeedMatchR` with the patched function, for the two calls in `../SeedMatchR/tests/testthat/helper.R:47-51` (the subset at line 47, the calls at lines 49 and 51)

  Run it.

- [ ] **Step 3: Write the failing tests**

```python
@pytest.mark.parametrize("variant, options", [
    ("default", FullSearchOptions()),
    ("wobble", FullSearchOptions(wobble=True)),
    ("indels", FullSearchOptions(indels=True)),
])
def test_full_search_counts_match_r(variant, options):
    got = full_search(Guide(GUIDE), sirna_seqs(), options)
    assert_frame_matches_golden(got, f"seedmatchr/full_search_counts_{variant}.tsv")

def test_best_match_skips_absent_columns_D012():
    got = full_search(Guide(GUIDE), sirna_seqs(), FullSearchOptions(nonseed_edit_distances=(0, 1)))
    assert set(got["best_match"]) <= {"Full-D0", "Full-D1", "18mer-D0", "18mer-D1", "15mer-D0",
                                      "15mer-D1", "8mer", "7mer-m8", "7mer-A1", "6mer", "No.matches"}

def test_full_search_ranges_match_patched_r():
    got = full_search_ranges(Guide(GUIDE), sirna_seqs())
    assert_frame_matches_golden(got, "seedmatchr/full_search_ranges.tsv")

def test_short_guide_fails_before_any_search():
    with pytest.raises(InputError, match=r"18mer needs a guide of at least 19 nt, got 18"):
        full_search(Guide("UUAUAGAGCAAGAACACU"), sirna_seqs())
```

`assert_frame_matches_golden(got: pandas.DataFrame, path: str) -> None` is a helper in `tests/helpers.py`. It does these steps:
  - It removes the R column `names`, because `seqnames` has the same values.
  - It renames the R columns with one dictionary, `R_COLUMNS`:
    - `seqs` and `seqnames` to `seq_id`
    - `Category` to `category`, `Group` to `group`, `best.match` to `best_match`
    - `TargetSeq` to `target_seq`, `EditDistance` to `edit_distance`, `AllowIndels` to `indels`
    - `IN.DESEQ2` to `detected`
  - It changes `TRUE` and `FALSE` to booleans, `Detected` and `Undetected` to `True` and `False`, and `NA` to a missing value.
  - It compares `got[<golden columns>]` with the golden, row by row. Integers, booleans, and strings must be equal. Floats must agree within a relative tolerance of `1e-9`.

The same dictionary gives the column rows of `docs/r-map.md`.

- [ ] **Step 4: Run `uv run pytest tests/report/test_full_search.py -v`.** Expected: FAIL with `ModuleNotFoundError`.

- [ ] **Step 5: Write `grid.py`.** Make all the targets of a call first, so that a guide that is too short fails before any search. Run the searches on a `ThreadPoolExecutor` with `max_workers=workers`, and keep the column order of the options.

- [ ] **Step 6: Run the tests.** Expected: PASS.

- [ ] **Step 7: Write `## D-012: The best match uses only the columns that the search made`.** R evidence: `R/full_search.R:98-121` names each column. Add the rows for `full_search`, `full_search_ranges`, `.get_best_match`, and `SeedMatchR` (`data.frame` mode) to `docs/r-map.md`. Checkpoint.

### Task 10: Seed-match report

**Files:**
- Make: `src/seedmatch/report/report.py`, `tools/golden/R/report.R`
- Make (made by the R script): `tests/golden/seedmatchr/report_<variant>.tsv`, `report_matches_<variant>.tsv`
- Modify: `docs/r-map.md`, `src/seedmatch/__init__.py`
- Test: `tests/report/test_report.py`, `tests/test_offline.py`

**Interfaces:**
- Consumes: Tasks 8, 9.
- Produces: `REPORT_CATEGORIES: tuple[str, ...]` in the R column order of `R/SeedMatchReport.R:1-5`, which puts the `15mer` classes from `D4` to `D0`.
- Produces: `ReportOptions(indels: bool = True, wobble: bool = False, log2fc_cut: float = -0.7655347, padj_cut: float = 0.05, seed_edit_distances: tuple[int, ...] = (0,), fixed_pattern: bool = True, fixed_subject: bool = True)`.
- Produces: `SeedMatchReport(report: pandas.DataFrame, matches: pandas.DataFrame)` (frozen dataclass).
- Produces: `seed_match_report(de: pandas.DataFrame, sequences: SequenceSet, guide: Guide, options: ReportOptions = ReportOptions(), group: str = "SeedMatchReport", *, workers: int | None = None) -> SeedMatchReport`.
  - `report` has the column `category` with the rows `In silico predictions`, `Expressed predictions`, `Off-target predictions`, `% off-target`, then one column for each name in `REPORT_CATEGORIES`, then `Total`, then `group`.
  - The counts are gene counts. The best match of a gene is the best over its sequences.
  - `% off-target` is `100 * off-target / expressed`, and 0 if `expressed` is 0. `Total` of `% off-target` uses the totals.
  - `matches` is the `full_search` table with the columns `detected` (bool), `log2FoldChange`, `padj`, `baseMean`.

- [ ] **Step 1: Write `report.R`.** It writes `SeedMatchReport(res = res, seqs = seqs, guide.seq = guide)` (`report_default`, `report_matches_default`) and the same with `allow_wobbles = TRUE` (`wobble`). Run it.

- [ ] **Step 2: Write the failing tests**

```python
@pytest.mark.parametrize("variant, options", [("default", ReportOptions()), ("wobble", ReportOptions(wobble=True))])
def test_report_matches_r(variant, options):
    got = seed_match_report(sirna_de(), sirna_seqs(), Guide(GUIDE), options)
    assert_frame_matches_golden(got.report, f"seedmatchr/report_{variant}.tsv")
    assert_frame_matches_golden(got.matches, f"seedmatchr/report_matches_{variant}.tsv")

def test_report_known_snapshot_values():
    row = seed_match_report(sirna_de(), sirna_seqs(), Guide(GUIDE)).report.set_index("category").loc["In silico predictions"]
    assert (row["Total"], row["15mer-D4"], row["8mer"], row["6mer"]) == (125, 73, 1, 4)

def test_report_counts_genes_not_sequences_D003():
    seqs = SequenceSet.from_records([("t1", "CTCTATAA"), ("t2", "CTCTATAA")], gene_ids=["g", "g"])
    de = pd.DataFrame({"gene_id": ["g"], "log2FoldChange": [-2.0], "padj": [0.01], "baseMean": [100.0]})
    row = seed_match_report(de, seqs, Guide(GUIDE)).report.set_index("category").loc["In silico predictions"]
    assert row["8mer"] == 1

def test_percent_is_zero_when_no_gene_is_expressed():
    de = pd.DataFrame({"gene_id": ["other"], "log2FoldChange": [0.0], "padj": [1.0], "baseMean": [100.0]})
    pct = seed_match_report(de, sirna_seqs(), Guide(GUIDE)).report.set_index("category").loc["% off-target"]
    assert (pct[list(REPORT_CATEGORIES)] == 0).all()
```

`tests/test_offline.py`:

```python
NETWORK_MODULES = ("seedmatch.fetch", "urllib.request", "http.client", "ssl", "urllib3", "requests", "httpx")

def run_core_pipeline():
    seqs = read_fasta(DATA / "Test_seqs.fa")
    de = filter_de(sirna_de(), basemean_cutoff=0)
    annotate_seed_counts(de, seqs, Guide(GUIDE).target(Seed.MER7M8))
    full_search_ranges(Guide(GUIDE), seqs)
    seed_match_report(de, seqs, Guide(GUIDE))

def test_core_pipeline_runs_with_sockets_disabled():
    run_core_pipeline()  # pytest-socket makes each connection attempt raise

def test_core_pipeline_loads_no_network_module():
    code = f"import sys, test_offline; test_offline.run_core_pipeline(); print([m for m in {NETWORK_MODULES!r} if m in sys.modules])"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True,
                         cwd=Path(__file__).parent, env=os.environ | {"PYTHONPATH": str(Path(__file__).parent)}).stdout
    assert out.strip() == "[]"
```

The list leaves out `socket`, because numpy and pandas import it.

- [ ] **Step 3: Run `uv run pytest tests/report/test_report.py tests/test_offline.py -v`.** Expected: FAIL with `ModuleNotFoundError`.

- [ ] **Step 4: Write `report.py`.** Use `full_search` with the report options, `best_match` for each gene in the priority order of Task 9, and one `merge` with the DE table.

- [ ] **Step 5: Run the tests.** Expected: PASS.

- [ ] **Step 6: Add the `SeedMatchReport` row to `docs/r-map.md`. Checkpoint.**

### Task 11: Verification of the plan

- [ ] **Step 1: Run `uv run ruff check`, `uv run ruff format --check`, `uv run mypy`.** Expected: no errors.
- [ ] **Step 2: Run the full suite one time: `uv run pytest`.** Expected: all tests pass, and the `refs` marker selects no test yet.
- [ ] **Step 3: Run `NUMBA_DISABLE_JIT=1 uv run pytest --cov --cov-report=term-missing`.** Expected: branch coverage of at least 95%. Give a reason for each line without coverage, or add a test.
- [ ] **Step 4: Make sure that `docs/divergences.md` has D-001, D-003, D-012, and D-015, and that each has the fields of the spec.** D-002 (the `reduce.features` names) belongs to the annotation plan. D-015 is new, thus add a row for it to the table in `docs/design.md`.
- [ ] **Step 5: Report to the user:** the test results, the coverage, the list of divergences, and the files that are ready to commit.
