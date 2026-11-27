# The Rust search kernels: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move the search kernels of `seedmatch.search` to a Rust extension with bit-parallel algorithms, keep the output of the numba kernels, and remove numba.

**Architecture:** One crate in `rust/` holds the algorithms as plain Rust. The cargo feature `python` adds the PyO3 module `seedmatch.search._kernels`. maturin turns on the feature, and `cargo test` runs without it, thus the Rust tests do not link to Python. The Python API and its checks stay in `src/seedmatch/search/api.py`.

**Tech Stack:** Rust 1.100.0, PyO3 0.30, rust-numpy 0.30, maturin 1.16, Python 3.12, numpy 2.5.2, pytest, and Hypothesis.

**Spec:** `docs/plans/2026-11-27-rust-search-design.md`. The search rules are in `docs/design.md`, section "Search rules".

## Global Constraints

- The output is the output of the numba kernels: the same matches, in the same order, with the same counts.
- No Python runtime dependency is new. `numba` and `llvmlite` go away.
- The crate has `#![forbid(unsafe_code)]`.
- The crate is in `rust/`. The extension module is `seedmatch.search._kernels`.
- The kernels keep the names and the argument order of `src/seedmatch/search/kernels.py`.
- PyO3 has the feature `abi3-py312`. `rust-toolchain.toml` pins Rust 1.100.0 with `rustfmt` and `clippy`.
- The Linux wheels have the tag `manylinux_2_28`. The runners are `ubuntu-latest`, `ubuntu-24.04-arm`, and `macos-latest`.
- Each action of a workflow has a commit SHA and a version comment.
- The crate has no thread pool. Each kernel releases the GIL for its full search.
- The Python coverage floor stays at 95%. The Rust code gets no coverage gate.
- Put scratch files in `tmp/claude/rust-search/`.
- On the development machine, run `. "$HOME/.cargo/env"` before each cargo command.
- Run each changed test file in its task. Run the full Python suite one time, in Task 6.

Git:

- Before Task 1, ask the user which branch to use. Make a branch only if the user asks for one.
- Each commit step needs the consent of the user for that commit. Use `git commit -s` and Conventional Commits.
- Write the commit body in ASD-STE100. Add no agent credit. Do not push.
- Put these git rules in the prompt of each subagent.

## Review Focus

1. A sequence that is longer than one block of the filter, with matches near a block limit. The production block size must give the scalar result. Task 4 has the test.
2. A sequence of low complexity, for example a run of `N` or of `A`, where most positions are candidates. The bit path must give the scalar result. The generator of Tasks 3 and 4 makes such sequences.
3. A `SequenceSet` with an offsets array that is C-contiguous but not aligned to 8 bytes. rust-numpy 0.30 gives no slice for it (`src/array.rs:765` of rust-numpy). The search must give the result of an aligned copy. Task 2 has the test.
4. A pattern of 64 letters with max_mm 62 or 63 on the indel route, where the top bit of the word and the add with overflow matter. Task 4 has the test.
5. A change of a Rust file, then `uv run`. The extension must build again, or a developer gets the old kernels without a message. Task 2 has the check.

## Files

| File | Task | Responsibility |
|---|---|---|
| `src/seedmatch/report/grid.py` | 1, 2 | `_searches` gives the searches of a grid. The module docstring names the search kernels, not numba. |
| `tools/perf/grid.py` | 1 | The time and the SHA-256 of each search of the two grids. |
| `tools/perf/report.py` | 1, 2, 6 | `--hashes`, `VERSIONS` without numba, the docstring without `NUMBA_CACHE_DIR`. |
| `rust-toolchain.toml` | 2 | The Rust version. |
| `rust/Cargo.toml`, `rust/Cargo.lock` | 2 | The crate `seedmatch-kernels`, lib `_kernels`. |
| `rust/src/lib.rs` | 2 | The modules, and `#![forbid(unsafe_code)]`. |
| `rust/src/letters.rs` | 2 | The letter table and its words. |
| `rust/src/hamming.rs` | 2, 3 | Rule 2: the scalar count, then shift-and. |
| `rust/src/indel.rs` | 2, 4 | Rule 4: the scalar DP, then the Myers filter. |
| `rust/src/kernels.rs` | 2 | The four kernels on a set, the offsets check, and the route of rules 4b and 5. |
| `rust/src/python.rs` | 2 | The PyO3 module, behind the feature `python`. |
| `rust/tests/common/mod.rs` | 3 | The random generator and the compat helper of the Rust tests. |
| `rust/tests/hamming_paths.rs` | 3 | Shift-and against the scalar count. |
| `rust/tests/indel_paths.rs` | 4 | The filter and the bit path against the scalar DP. |
| `src/seedmatch/search/_kernels.pyi` | 2 | The types of the extension for mypy. |
| `src/seedmatch/search/kernels.py` | 2 | Deleted. |
| `src/seedmatch/search/api.py` | 2 | The import of `_kernels`, and aligned offsets. |
| `tests/search/test_kernels.py` | 2, 4 | The import of `_kernels`, and the new kernel tests. |
| `tests/search/test_api.py` | 2 | The test of unaligned offsets. |
| `pyproject.toml`, `uv.lock`, `.gitignore` | 2 | The maturin build, without numba. |
| `.github/workflows/*.yml`, `.github/dependabot.yml` | 5 | The CI checks and the wheels. |
| `tools/coverage/serial_pool.py` | 5 | Deleted. |
| `README.md`, `CONTRIBUTING.md`, `docs/design.md`, `docs/performance.md` | 6 | The docs. |

## Interfaces of the crate

Tasks 2, 3, and 4 use these names. A code is a subject mask, and each lookup reads `code & 15`.

```rust
// rust/src/lib.rs
#![forbid(unsafe_code)]
pub mod hamming;
pub mod indel;
pub mod kernels;
pub mod letters;
#[cfg(feature = "python")]
mod python;

// rust/src/letters.rs
pub const MASKS: usize = 16;
pub const WORD: usize = 64;
pub type Compat = [[bool; MASKS]; MASKS]; // [pattern mask][subject mask]
pub struct Letters { /* ok[a][s] for each pattern letter a */ }
impl Letters {
    pub fn new(pattern: &[u8], compat: &Compat) -> Self;
    pub fn len(&self) -> usize;
    pub fn ok(&self, a: usize, code: u8) -> bool;
    pub fn first(&self, code: u8) -> usize; // the first a with ok(a, code), else len()
    pub fn words(&self) -> [u64; MASKS]; // bit a is ok(a, s); only for len() <= WORD
    pub fn reversed_words(&self) -> [u64; MASKS]; // bit a is ok(len() - 1 - a, s); only for len() <= WORD
}

// rust/src/hamming.rs: emit gets the 1-based start of each match, in increasing order.
pub fn scalar(seq: &[u8], letters: &Letters, max_mm: u64, min_mm: u64, emit: &mut impl FnMut(i64));
pub fn bits(seq: &[u8], letters: &Letters, max_mm: u64, min_mm: u64, emit: &mut impl FnMut(i64)); // Task 3
pub fn search(seq: &[u8], letters: &Letters, max_mm: u64, min_mm: u64, emit: &mut impl FnMut(i64));

// rust/src/indel.rs: emit gets the 1-based start and end. The caller makes sure that 1 <= max_mm < len().
pub const BLOCK: usize = 1 << 16;
pub fn scalar(seq: &[u8], letters: &Letters, max_mm: usize, emit: &mut impl FnMut(i64, i64));
pub fn bits(seq: &[u8], letters: &Letters, max_mm: usize, block: usize, emit: &mut impl FnMut(i64, i64)); // Task 4
pub fn search(seq: &[u8], letters: &Letters, max_mm: usize, emit: &mut impl FnMut(i64, i64));
pub fn accepted(seq: &[u8], letters: &Letters, max_mm: usize) -> Vec<usize>; // Task 4: each 0-based j that the exact step accepts
pub fn candidates(seq: &[u8], letters: &Letters, max_mm: usize, block: usize) -> Vec<usize>; // Task 4: each 0-based j with S(j) <= max_mm

// rust/src/kernels.rs
pub struct OffsetsError { pub codes: usize }
// Display: "offsets do not fit codes of length {codes}; expected the offsets of a SequenceSet"
pub fn hamming_counts(codes: &[u8], offsets: &[i64], pattern: &[u8], compat: &Compat, max_mm: u64, min_mm: u64) -> Result<Vec<i64>, OffsetsError>;
pub fn hamming_ranges(codes: &[u8], offsets: &[i64], pattern: &[u8], compat: &Compat, max_mm: u64, min_mm: u64) -> Result<(Vec<i64>, Vec<i64>), OffsetsError>; // index, start
pub fn indel_counts(codes: &[u8], offsets: &[i64], pattern: &[u8], compat: &Compat, max_mm: u64) -> Result<Vec<i64>, OffsetsError>;
pub fn indel_ranges(codes: &[u8], offsets: &[i64], pattern: &[u8], compat: &Compat, max_mm: u64) -> Result<(Vec<i64>, Vec<i64>, Vec<i64>), OffsetsError>; // index, start, end
```

The Python stub `src/seedmatch/search/_kernels.pyi` gives the four functions with positional-only arguments. A limit is an `int`, each array is an `NDArray` of its dtype, and the results are `NDArray[np.int64]` values or tuples of them.

---

### Task 1: The timing tool and the numbers of the numba build

**Files:**
- Modify: `src/seedmatch/report/grid.py` (`full_search`, the loop that makes `searches`)
- Modify: `tools/perf/report.py`
- Make: `tools/perf/grid.py`

**Interfaces:**
- Produces: `_searches(guide: Guide, options: FullSearchOptions) -> list[_Search]` in `grid.py`.
- Produces: the files `tmp/claude/rust-search/before/grid.tsv` and `tmp/claude/rust-search/before/report-*.txt`. Tasks 2, 3, 4, and 6 compare with them.

- [ ] **Step 1: Move the loop of `full_search` that makes `searches` into `_searches(guide, options)`, and call it from `full_search`.**

The order of the searches does not change.

- [ ] **Step 2: Run the tests of the grid**

Run: `uv run pytest tests/report/test_full_search.py -q`
Expected: PASS.

- [ ] **Step 3: Write `tools/perf/grid.py`**

Usage: `uv run python tools/perf/grid.py FASTA [--repeat 3]`. Use the guide `UUAUAGAGCAAGAACACUGUUUU` of `tools/perf/report.py`. The grids are `FullSearchOptions()` and `FullSearchOptions(indels=True)`. For each search of `_searches`:

1. Call `count_matches` one time, then time `repeat` calls in the thread of the caller. Keep the best time.
2. Call `match_ranges` one time. Hash `index.tobytes() + starts.tobytes() + ends.tobytes()` with SHA-256.
3. Print `{name}\tindels={on|off}\t{seconds:.4f}\t{sha256}`.

After each grid, print `total\tindels={on|off}\t{seconds:.3f}`.

- [ ] **Step 4: Add `--hashes` to `tools/perf/report.py`**

With `--hashes`, print the SHA-256 of `result.report` and of `result.matches`. Hash the UTF-8 text of `frame.to_csv(sep="\t", index=False)`.

- [ ] **Step 5: Run the lint**

Run: `uv run ruff check && uv run ruff format --check && uv run mypy`
Expected: no error.

- [ ] **Step 6: Measure the numba build**

Run each command with `NUMBA_CACHE_DIR="$PWD/tmp/perf/numba"`, and save the output in `tmp/claude/rust-search/before/`:

```bash
uv run python tools/perf/grid.py tmp/perf/rnor7_utr3.fa > tmp/claude/rust-search/before/grid.tsv
uv run python tools/perf/report.py tmp/perf/rnor7_utr3.fa --workers 1 --hashes > tmp/claude/rust-search/before/report-workers1-1.txt
uv run python tools/perf/report.py tmp/perf/rnor7_utr3.fa --workers 1 --hashes > tmp/claude/rust-search/before/report-workers1-2.txt
uv run python tools/perf/report.py tmp/perf/rnor7_utr3.fa --hashes > tmp/claude/rust-search/before/report-default-1.txt
uv run python tools/perf/report.py tmp/perf/rnor7_utr3.fa --hashes > tmp/claude/rust-search/before/report-default-2.txt
```

Expected: the grid file has a line for each search and a `total` line for each grid. The total with indels on is near 9.9 seconds. A call with `--workers 1` takes near 8.7 seconds. If `tmp/perf/rnor7_utr3.fa` does not exist, make it with `tools/perf/utr3_fasta.py`, as `docs/performance.md` tells.

- [ ] **Step 7: Commit, after the consent of the user**

```bash
git add src/seedmatch/report/grid.py tools/perf/grid.py tools/perf/report.py
git commit -s -m "chore: add a tool that times each search of the full_search grid"
```

### Task 2: The crate, the build, and the scalar kernels

**Files:**
- Make: `rust-toolchain.toml`, `rust/Cargo.toml`, `rust/Cargo.lock`, `rust/src/lib.rs`, `rust/src/letters.rs`, `rust/src/hamming.rs`, `rust/src/indel.rs`, `rust/src/kernels.rs`, `rust/src/python.rs`, `src/seedmatch/search/_kernels.pyi`
- Delete: `src/seedmatch/search/kernels.py`
- Modify: `src/seedmatch/search/api.py`, `src/seedmatch/report/grid.py:3-4`, `tools/perf/report.py` (`VERSIONS`), `pyproject.toml`, `uv.lock`, `.gitignore`
- Test: `tests/search/test_kernels.py`, `tests/search/test_api.py`

**Interfaces:**
- Consumes: `tools/perf/grid.py` and `tmp/claude/rust-search/before/` from Task 1.
- Produces: each item of "Interfaces of the crate" except the items that name Task 3 or Task 4. `hamming::search` and `indel::search` call `scalar` in this task.

- [ ] **Step 1: Change the tests**

In `tests/search/test_kernels.py`:

- Import the four kernels from `seedmatch.search._kernels`.
- Remove `from numba import config` and `CORPUS_STRIDE`. The corpus test reads each group.
- Remove the comment above `KERNEL_SETTINGS`. Keep `deadline=None` only if the tests fail with the default deadline, and then write that cause in the comment.

Add these tests. `long_case(m, size)` gives a pattern of m letters from `ACGTRYN` and a subject of `size` letters, from `random.Random(m)`. The subject has random letters, one exact copy of the pattern, and one copy with 2 random edits.

```python
LONG_LENGTHS = (63, 64, 65, 80)
ALL_FIXED = list(product((True, False), repeat=2))


@pytest.mark.parametrize("m", LONG_LENGTHS)
@pytest.mark.parametrize("max_mm", [0, 1, 2, 3])
@pytest.mark.parametrize("fixed", ALL_FIXED)
def test_long_patterns_equal_the_oracle(m, max_mm, fixed):
    pattern, subject = long_case(m, 4 * m)
    expected = oracle.hamming(pattern, subject, max_mm, 0, *fixed)
    assert expected
    assert kernel_hamming(pattern, subject, max_mm, 0, fixed) == expected
    assert kernel_indels(pattern, subject, max_mm, fixed) == oracle.indels(pattern, subject, max_mm, *fixed)


KERNEL_CALLS = [
    lambda c, o, p, t: hamming_counts(c, o, p, t, 1, 0),
    lambda c, o, p, t: hamming_ranges(c, o, p, t, 1, 0),
    lambda c, o, p, t: indel_counts(c, o, p, t, 1),
    lambda c, o, p, t: indel_ranges(c, o, p, t, 1),
]
COMPAT = compat_table(fixed_pattern=True, fixed_subject=True)


@pytest.mark.parametrize("call", KERNEL_CALLS)
def test_kernels_reject_offsets_outside_codes(call):
    message = "offsets do not fit codes of length 2; expected the offsets of a SequenceSet"
    with pytest.raises(ValueError, match="^" + re.escape(message) + "$"):
        call(np.array([1, 2], np.uint8), np.array([0, 3], np.int64), encode_dna("AC"), COMPAT)


@pytest.mark.parametrize("call", KERNEL_CALLS)
def test_kernels_read_a_code_above_15_as_its_low_four_bits(call):
    offsets = np.array([0, 3], np.int64)
    high = call(np.array([200, 1, 8], np.uint8), offsets, encode_dna("TA"), COMPAT)
    low = call(np.array([8, 1, 8], np.uint8), offsets, encode_dna("TA"), COMPAT)
    for found, expected in zip(np.atleast_2d(high), np.atleast_2d(low), strict=True):
        np.testing.assert_array_equal(found, expected)


INT64_MAX = 2**63 - 1


def test_limits_near_the_int64_maximum():
    assert hamming_count_list("AC", EMPTY_INSIDE, INT64_MAX, 0) == [3, 1, 4]
    assert hamming_count_list("AC", EMPTY_INSIDE, INT64_MAX, INT64_MAX) == [0, 0, 0]
    assert indel_hits("AC", EMPTY_INSIDE, INT64_MAX) == EMPTY_INSIDE_MISMATCH_HITS
```

In `tests/search/test_api.py`, add:

```python
def test_search_of_a_set_with_unaligned_offsets_equals_the_aligned_search():
    offsets = np.frombuffer(bytearray(1 + 8 * 3), dtype=np.int64, offset=1)
    offsets[:] = [0, 4, 8]
    assert not offsets.flags.aligned
    codes = encode_dna("ACGTACGA")
    unaligned = SequenceSet(("a", "b"), ("a", "b"), codes.copy(), offsets)
    aligned = SequenceSet(("a", "b"), ("a", "b"), codes.copy(), offsets.copy())
    for options in (SearchOptions(max_mismatch=1), SearchOptions(max_mismatch=1, indels=True)):
        np.testing.assert_array_equal(
            count_matches("ACG", unaligned, options), count_matches("ACG", aligned, options)
        )
```

- [ ] **Step 2: Run the tests to see them fail**

Run: `uv run pytest tests/search/test_kernels.py tests/search/test_api.py -q`
Expected: FAIL. `test_kernels.py` stops with `ModuleNotFoundError: No module named 'seedmatch.search._kernels'`. The result of the new test of `test_api.py` on numba does not matter. The test guards Step 4.

- [ ] **Step 3: Write the crate and the build**

- `rust-toolchain.toml` at the root of the repository: the channel `1.100.0`, the profile `minimal`, and the components `rustfmt` and `clippy`.
- `rust/Cargo.toml`:
  - the package `seedmatch-kernels`, version `0.1.0`, edition 2024, `publish = false`, license MIT
  - `[lib] name = "_kernels"` and `crate-type = ["cdylib", "rlib"]`
  - `[features] python = ["dep:pyo3", "dep:numpy", "pyo3/extension-module"]`
  - the optional dependencies `pyo3 = { version = "0.30", features = ["abi3-py312"] }` and `numpy = "0.30"`
- `rust/src/`: the modules of "Interfaces of the crate". `hamming::scalar` is `_hamming_sequence` of `kernels.py`. `indel::scalar` is `_indel_sequence` with `_right_alignment`. Keep each branch of the numba code, and keep the DP property of the design: each value of the row is written before it is read.
- `kernels.rs` reads each offset one time. It returns `OffsetsError` in these conditions:
  - an offset is negative
  - an offset is less than the offset before it
  - an offset is more than `codes.len()`
- `kernels.rs` routes `indel_*`: if max_mm is 0 or max_mm ≥ the pattern length, it gives the rule 2 result with min_mm 0, and the end is start + length − 1. Else it calls `indel::search`.
- A ranges kernel adds each match to its vectors in one pass.
- `python.rs`: one `#[pyfunction]` for each kernel, with `#[pyo3(signature = (..., /))]`. Each one does these steps:
  1. Take `PyReadonlyArray1` for `codes`, `offsets`, and `pattern`, `PyReadonlyArray2<bool>` for `compat`, and `u64` limits.
  2. Make sure that `compat` has the shape (16, 16). Else, raise `ValueError("compat has shape {shape}; expected (16, 16)")`.
  3. Call the kernel inside `py.detach`. Give an `OffsetsError` as `ValueError` with its `Display` text.
  4. Return the vectors with `PyArray1::from_vec`.
- `pyproject.toml`:
  - `[build-system]`: `requires = ["maturin>=1.16,<2"]` and `build-backend = "maturin"`.
  - `[tool.maturin]`: `manifest-path = "rust/Cargo.toml"`, `python-source = "src"`, `module-name = "seedmatch.search._kernels"`, and `features = ["python"]`.
  - `[tool.uv] cache-keys`: `pyproject.toml`, `rust/Cargo.toml`, `rust/Cargo.lock`, `rust/src/**/*.rs`, and `rust-toolchain.toml`.
  - Remove `numba` from `dependencies`, `numba` and `llvmlite` from `constraint-dependencies`, the mypy overrides of `numba.*` and `seedmatch.search.kernels`, and the ruff ignore of `kernels.py`.
- `.gitignore`: add `rust/target/` and `*.so`.
- `api.py`: import `_kernels`. Give the kernels `np.require(sequences.offsets, requirements="A")`. Remove the comment about the numba cache entry (line 45).
- `grid.py:3-4`: "The search kernels release the GIL."
- `tools/perf/report.py`: `VERSIONS = ("seedmatch", "numpy", "pandas")`.
- Delete `src/seedmatch/search/kernels.py`, and write `_kernels.pyi`.
- Run `uv lock`.

- [ ] **Step 4: Build and run the tests**

Run: `uv sync && uv run pytest tests/search tests/report/test_full_search.py -q`
Expected: PASS.

- [ ] **Step 5: Make sure that the build profile is release**

Run: `uv sync --reinstall-package seedmatch -v 2>&1 | grep -i 'profile'`
Expected: `release` profile. If the profile is `dev`, set `profile = "release"` in `[tool.maturin]`, and do Step 4 again.

- [ ] **Step 6: Make sure that a change of a Rust file builds the extension again**

Add a blank line to `rust/src/lib.rs`, then run `uv run python -c "import seedmatch.search._kernels"`.
Expected: uv builds `seedmatch` again. Remove the blank line.

- [ ] **Step 7: Run the Rust checks and the lint**

Run: `cargo fmt --manifest-path rust/Cargo.toml --check && cargo clippy --manifest-path rust/Cargo.toml --all-targets --all-features -- -D warnings && cargo test --manifest-path rust/Cargo.toml && uv run ruff check && uv run ruff format --check && uv run mypy`
Expected: no error.

- [ ] **Step 8: Compare with the numba build**

Run: `uv run python tools/perf/grid.py tmp/perf/rnor7_utr3.fa > tmp/claude/rust-search/scalar/grid.tsv`
Expected: each SHA-256 is equal to its value in `before/grid.tsv`. If a value is different, stop and find the cause.

- [ ] **Step 9: Commit, after the consent of the user**

```bash
git add rust-toolchain.toml rust pyproject.toml uv.lock .gitignore src tests tools/perf/report.py
git commit -s -m "perf: move the search kernels from numba to a Rust extension"
```

The body names seedmatch#27, because the change removes its cause.

### Task 3: Shift-and for rule 2

**Files:**
- Make: `rust/tests/common/mod.rs`, `rust/tests/hamming_paths.rs`
- Modify: `rust/src/hamming.rs`

**Interfaces:**
- Consumes: `Letters`, `Compat`, `hamming::scalar`, and `WORD` from Task 2.
- Produces: `hamming::bits`. `hamming::search` calls `bits` if the pattern has `WORD` letters or fewer. Task 4 uses `rust/tests/common/mod.rs`.

- [ ] **Step 1: Write the test helpers and the failing test**

`rust/tests/common/mod.rs` has:

- a SplitMix64 generator
- `compat(fixed_pattern: bool, fixed_subject: bool) -> Compat`, with the four rules of `src/seedmatch/search/compat.py` and a false row and column for mask 0
- `sequence(rng, len) -> Vec<u8>`. Most sequences have masks from 1 to 15, with more weight on 1, 2, 4, and 8. One in 10 is a run of 15 or a run of 1.
- `plant(rng, seq, pattern, edits)`: put a copy of the pattern with `edits` random edits at a random position. Each letter of the copy is a base of its pattern mask.

`rust/tests/hamming_paths.rs`:

```rust
#[test]
fn bits_equal_scalar_on_random_input() {
    // 20,000 cases from a fixed seed. Each case: m from 1 to 64, a sequence of 0 to 300 letters
    // with 0 to 2 planted copies, each pair of fixed flags, max_mm from 0 to m + 2,
    // and min_mm from 0 to max_mm. One case in 50 uses max_mm = min_mm = u64::MAX.
    for case in cases() {
        let (mut bits, mut scalar) = (Vec::new(), Vec::new());
        hamming::bits(&case.seq, &case.letters, case.max_mm, case.min_mm, &mut |s| bits.push(s));
        hamming::scalar(&case.seq, &case.letters, case.max_mm, case.min_mm, &mut |s| scalar.push(s));
        assert_eq!(bits, scalar, "case {}", case.id);
    }
}
```

- [ ] **Step 2: Run the test to see it fail**

Run: `cargo test --manifest-path rust/Cargo.toml --test hamming_paths`
Expected: FAIL to compile with "cannot find function `bits` in module `hamming`".

- [ ] **Step 3: Write `hamming::bits`, and route `hamming::search`**

Use the update and the virtual letters of the design, section "Rule 2: shift-and with k mismatches". Cap each state index at the pattern length.

- [ ] **Step 4: Run the tests**

Run: `cargo test --manifest-path rust/Cargo.toml && uv sync && uv run pytest tests/search -q`
Expected: PASS.

- [ ] **Step 5: Compare the speed and the output**

Run: `uv run python tools/perf/grid.py tmp/perf/rnor7_utr3.fa > tmp/claude/rust-search/shift-and/grid.tsv`
Expected: each SHA-256 is equal to its value in `before/grid.tsv`. Compare each line with `indels=off` with `scalar/grid.tsv`. If shift-and is slower for a class of search, route this class to `scalar` with a rule on the pattern length and max_mm. Record the numbers for Task 6.

- [ ] **Step 6: Run the Rust checks**

Run: `cargo fmt --manifest-path rust/Cargo.toml --check && cargo clippy --manifest-path rust/Cargo.toml --all-targets --all-features -- -D warnings`
Expected: no error.

- [ ] **Step 7: Commit, after the consent of the user**

```bash
git add rust
git commit -s -m "perf: count the mismatch matches with shift-and"
```

### Task 4: The Myers filter for rule 4

**Files:**
- Make: `rust/tests/indel_paths.rs`
- Modify: `rust/src/indel.rs`
- Test: `tests/search/test_kernels.py`

**Interfaces:**
- Consumes: `Letters::reversed_words`, `Letters::first`, `indel::scalar`, and `rust/tests/common/mod.rs`.
- Produces: `indel::bits`, `indel::accepted`, and `indel::candidates`. `indel::search` calls `bits` with `BLOCK` if the pattern has `WORD` letters or fewer.

- [ ] **Step 1: Write the failing tests**

`rust/tests/indel_paths.rs`:

```rust
const BLOCKS: [usize; 7] = [1, 2, 3, 5, 8, 64, 301];

#[test]
fn bits_equal_scalar_on_random_input() {
    // 20,000 cases from a fixed seed: m from 2 to 64, max_mm from 1 to m - 1, a sequence of 0 to 300
    // letters with 0 to 2 planted copies, each pair of fixed flags, and a block from BLOCKS.
    for case in cases() {
        let (mut bits, mut scalar) = (Vec::new(), Vec::new());
        indel::bits(&case.seq, &case.letters, case.max_mm, case.block, &mut |s, e| bits.push((s, e)));
        indel::scalar(&case.seq, &case.letters, case.max_mm, &mut |s, e| scalar.push((s, e)));
        assert_eq!(bits, scalar, "case {}", case.id);
    }
}

#[test]
fn filter_marks_each_accepted_position() {
    for case in cases() {
        let candidates = indel::candidates(&case.seq, &case.letters, case.max_mm, case.block);
        for j in indel::accepted(&case.seq, &case.letters, case.max_mm) {
            assert!(candidates.binary_search(&j).is_ok(), "case {} position {j}", case.id);
        }
    }
}

#[test]
fn bits_equal_scalar_on_a_sequence_longer_than_one_block() {
    // One sequence of 3 * BLOCK + 17 letters, m = 21, max_mm from 1 to 4. Planted copies with
    // 0 to 4 edits start at BLOCK - 21, BLOCK - 1, BLOCK, 2 * BLOCK - 10, and 3 * BLOCK.
}
```

In `tests/search/test_kernels.py`:

```python
@pytest.mark.parametrize("max_mm", [62, 63])
@pytest.mark.parametrize("fixed", ALL_FIXED)
def test_pattern_of_64_letters_with_max_mm_near_its_length_equals_the_oracle(max_mm, fixed):
    pattern, subject = long_case(64, 100)
    expected = oracle.indels(pattern, subject, max_mm, *fixed)
    assert expected
    assert kernel_indels(pattern, subject, max_mm, fixed) == expected
```

- [ ] **Step 2: Run the tests to see them fail**

Run: `cargo test --manifest-path rust/Cargo.toml --test indel_paths`
Expected: FAIL to compile with "cannot find function `bits` in module `indel`".

Run: `uv run pytest tests/search/test_kernels.py -q -k 64_letters`
Expected: PASS on the scalar path. The test guards Step 3.

- [ ] **Step 3: Write the filter and the bit path**

- Move the work at one position j of `scalar` into one function: the check of `first`, the DP, and the cost. Move the held-candidate logic into one type. `scalar`, `bits`, and `accepted` use the two, thus the exact step is the same code on each path.
- `bits` visits the blocks in order. For each block, it finds the candidates in a bit set of `block` bits, and runs the exact step at each candidate in increasing order. The held candidate carries to the next block.
- `candidates` gives each marked position of each block, in increasing order.

The reversed Myers run of one block, for the block from b0 to b1 − 1 (0-based), with `peq = letters.reversed_words()`:

```text
top = min(n, b1 + m + max_mm - 1)
high = 1 << (m - 1)
pv = u64::MAX; mv = 0; score = m
for i from top - 1 down to b0:
    eq = peq[seq[i] & 15]
    xv = eq | mv
    xh = ((eq & pv).wrapping_add(pv) ^ pv) | eq
    ph = mv | !(xh | pv)
    mh = pv & xh
    if ph & high != 0: score += 1
    if mh & high != 0: score -= 1
    ph <<= 1
    mh <<= 1
    pv = mh | !(xv | ph)
    mv = ph & xv
    if i < b1 and score <= max_mm: mark i
```

The `ph <<= 1` has no `| 1`, thus a text part can start at each position of the reversed run. The design gives the proof of the bound of `top`.

- [ ] **Step 4: Run the tests**

Run: `cargo test --manifest-path rust/Cargo.toml && uv sync && uv run pytest tests/search -q`
Expected: PASS.

- [ ] **Step 5: Compare the speed and the output**

Run: `uv run python tools/perf/grid.py tmp/perf/rnor7_utr3.fa > tmp/claude/rust-search/myers/grid.tsv`
Expected: each SHA-256 is equal to its value in `before/grid.tsv`. The `total` with `indels=on` is less than in `scalar/grid.tsv`. Record the numbers for Task 6.

- [ ] **Step 6: Run the Rust checks**

Run: `cargo fmt --manifest-path rust/Cargo.toml --check && cargo clippy --manifest-path rust/Cargo.toml --all-targets --all-features -- -D warnings`
Expected: no error.

- [ ] **Step 7: Commit, after the consent of the user**

```bash
git add rust tests/search/test_kernels.py
git commit -s -m "perf: filter the indel candidates with the Myers bit-vector algorithm"
```

### Task 5: CI and the wheels

**Files:**
- Modify: `.github/workflows/test.yml`, `.github/workflows/lint.yml`, `.github/workflows/coverage.yml`, `.github/dependabot.yml`
- Make: `.github/workflows/wheels.yml`
- Delete: `tools/coverage/serial_pool.py`

**Interfaces:**
- Consumes: the build of Task 2.
- Produces: the wheels as artifacts of each run of `wheels.yml`.

- [ ] **Step 1: Make sure that rustup reads `rust-toolchain.toml`**

Run: `rustup toolchain install` at the root of the repository.
Expected: rustup installs or keeps `1.100.0` with `rustfmt` and `clippy`. CI uses the same command.

- [ ] **Step 2: Change the workflows**

- `test.yml`: in each job, before "Install the dependencies", add the step "Install the Rust toolchain" with `rustup toolchain install`. In the job `test`, after "Run the tests", add "Run the Rust tests" with `cargo test --manifest-path rust/Cargo.toml --locked`.
- `lint.yml`: add the toolchain step. After "Mypy", add these steps:
  - "Rustfmt": `cargo fmt --manifest-path rust/Cargo.toml --check`
  - "Clippy": `cargo clippy --manifest-path rust/Cargo.toml --locked --all-targets --all-features -- -D warnings`, with `PYO3_PYTHON: ${{ github.workspace }}/.venv/bin/python`
- `coverage.yml`: add the toolchain step. Remove `NUMBA_DISABLE_JIT`, `PYTHONPATH`, `-p serial_pool`, and the comment lines about numba and the plugin. Task 6 measures the time of the run and corrects the comment about 30 minutes.
- Delete `tools/coverage/serial_pool.py`.
- `dependabot.yml`: remove `numba` and `llvmlite`. Add the ecosystem `cargo` for `/rust`, with the schedule, the limit, the commit prefix, and the groups of the `uv` entry.

- [ ] **Step 3: Write `wheels.yml`**

- The triggers: `pull_request` and `push` to `main`, and `push` of the tags `v*`. Use `permissions: contents: read` and the `concurrency` block of `test.yml`.
- The job `linux`: a matrix of `ubuntu-latest` with `x86_64` and `ubuntu-24.04-arm` with `aarch64`. Build with `PyO3/maturin-action`, `manylinux: "2_28"`, and `args: --release --locked --out dist`.
- The job `macos`: `macos-latest` with `aarch64`, and the same `args`.
- Read the README of `maturin-action` at the pinned tag. Make sure that the build uses the toolchain of `rust-toolchain.toml`.
- In each wheel job, after the build:
  1. Make a new environment with `uv venv --python 3.12`.
  2. Install the wheel and the `dev` group.
  3. Run `python -m pytest tests/search` from the root of the repository. `pyproject.toml` puts only `tests` on the path, thus the tests import the installed wheel.
- The job `sdist`: run `maturin sdist --out dist`. Then run `uv build --wheel` on the sdist, to make sure that the sdist holds each file of the build.
- Each job keeps `dist/` with `actions/upload-artifact`.
- Pin each action with a commit SHA and a version comment.

- [ ] **Step 4: Do a check of the workflows**

Run: `uvx --from actionlint-py actionlint`
Expected: no output.

- [ ] **Step 5: Build and do a test of a wheel and of the sdist on this machine**

```bash
uvx maturin build --release --locked --out tmp/claude/rust-search/dist
uv venv --python 3.12 tmp/claude/rust-search/wheel-env
uv pip install --python tmp/claude/rust-search/wheel-env tmp/claude/rust-search/dist/*.whl --group dev
tmp/claude/rust-search/wheel-env/bin/python -m pytest tests/search -q
uvx maturin sdist --out tmp/claude/rust-search/sdist
uv build --wheel tmp/claude/rust-search/sdist/seedmatch-0.1.0.tar.gz --out-dir tmp/claude/rust-search/sdist-wheel
```

Expected: the wheel name has `cp312-abi3` and `macosx` and `arm64`. The tests pass. The sdist gives a wheel. If the sdist has no `rust-toolchain.toml`, add it to `[tool.maturin] include` for the sdist, and do this step again.

- [ ] **Step 6: Commit, after the consent of the user**

```bash
git add .github tools/coverage pyproject.toml
git commit -s -m "ci: build the wheels and run the Rust checks"
```

CI runs only after the user pushes. Give the user the list of the new jobs to watch.

### Task 6: The numbers after the change, and the docs

**Files:**
- Modify: `docs/performance.md`, `docs/design.md`, `README.md`, `CONTRIBUTING.md`, `tools/perf/report.py` (docstring), `.github/workflows/coverage.yml` (comment)

**Interfaces:**
- Consumes: `tmp/claude/rust-search/before/` from Task 1, and the files of Tasks 2 to 4.

- [ ] **Step 1: Measure**

Do the commands of Task 1, Step 6, without `NUMBA_CACHE_DIR`. Write the output to `tmp/claude/rust-search/after/`.

- [ ] **Step 2: Compare with the targets of the design**

Expected:

- Each SHA-256 of `after/` is equal to its value in `before/`.
- Each call with `--workers 1` takes 1.75 seconds or less.
- No search in `after/grid.tsv` takes more time than in `before/grid.tsv`. If one does, measure it again. If it is still slower, record it.
- Call 1 takes at most 0.1 seconds more than call 2.

If a hash is different, stop: the change is wrong. If a time target is not met, give the numbers to the user and stop.

- [ ] **Step 3: Update `docs/performance.md`**

Give the new numbers, the method without `NUMBA_CACHE_DIR`, and the versions: the seedmatch commit, Rust 1.100.0, PyO3, and rust-numpy. Remove the row "First call, empty numba cache". Keep the numbers of R. Give the numba numbers of `before/` as the state before the change. Give the per-search numbers of the grids with indels on.

- [ ] **Step 4: Update `docs/design.md`**

- Line 41: maturin is the backend, and the crate is in `rust/`.
- Line 44: the required dependencies are `numpy`, `pandas`, and `scipy`.
- Line 51: the wheel has a compiled extension. The build images have no Rust toolchain, thus the index must have a wheel for each platform of the sandbox. Name the wheels.
- Line 52: remove it.
- Line 97: "Rust kernels and the search API."
- Line 562: the property tests compare the Rust kernels with the oracle. The Rust tests compare the bit path with the scalar path.
- Line 623: give the complexity of shift-and, of the filter with the DP at the candidates, and of the scalar path for more than 64 letters.
- Line 624: the kernels release the GIL.
- Line 663: mypy with `--strict`, without an override.
- Line 698: the Rust kernels with bit-parallel algorithms, and a link to the design of this change.
- The target of "Performance and bounds": no "after the JIT compilation".
- "Tooling": cargo, rustfmt, clippy, and maturin.

- [ ] **Step 5: Update `README.md`, `CONTRIBUTING.md`, and `tools/perf/report.py`**

- `README.md:18`: the package is Python with a compiled search extension, on numpy, pandas, and scipy.
- `README.md:405-419`: a wheel needs no Rust toolchain. An install from GitHub or from a clone builds from source with the toolchain of `rust-toolchain.toml`.
- `README.md:421`: remove the text about the numba cache.
- `README.md:441-445`: the coverage command is `uv run pytest --cov`. Remove the paragraph about the tracer and the plugin. Add the cargo checks to the list of checks.
- `CONTRIBUTING.md`: in "Set up the project", install rustup. In "Do the checks", add the cargo test, rustfmt, and clippy commands.
- `tools/perf/report.py`: remove `NUMBA_CACHE_DIR` from the docstring and from its example command.

- [ ] **Step 6: Run the STE gate on each changed document**

Run: `"$HOME/.local/share/logos/bin/bun" "$HOME/.local/share/logos/src/ste/cli.ts" --file <path>` for each changed Markdown file.
Expected: no `[hard]` finding.

- [ ] **Step 7: Measure the Rust coverage one time**

Run: `cargo llvm-cov --manifest-path rust/Cargo.toml --summary-only`
Expected: a line coverage number for the pull request. If `cargo llvm-cov` is not installed, ask the user before you install it.

- [ ] **Step 8: Run each check, and the full Python suite one time**

Run: `uv run ruff check && uv run ruff format --check && uv run mypy && cargo fmt --manifest-path rust/Cargo.toml --check && cargo clippy --manifest-path rust/Cargo.toml --all-targets --all-features -- -D warnings && cargo test --manifest-path rust/Cargo.toml`
Expected: no error.

Run: `time uv run pytest --cov --cov-report=term-missing`
Expected: PASS, and a total coverage of 95% or more. The run ends without the plugin `serial_pool`. If it does not end, put the plugin back with the new cause in its docstring. Write the measured time in the comment of `coverage.yml`.

- [ ] **Step 9: Commit, after the consent of the user**

```bash
git add docs README.md CONTRIBUTING.md tools/perf/report.py .github/workflows/coverage.yml
git commit -s -m "docs: record the time of the Rust kernels and remove numba from the docs"
```

After this task, tell the user these facts:

- The pull request closes seedmatch#27.
- `CONTRIBUTING.md` asks for an issue before a large change. Ask the user if you must open one.
- The wheels are not published. `INF-01` of the standard makes the publication necessary for the sandbox.
