# The Rust search kernels: design

Date: 2026-11-27. Status: a draft for the approval of the user.

## Intent

The search kernels of `seedmatch.search` use most of the time of a report. On the canonical 3' UTRs of `rnor7`, `seed_match_report` with `workers=1` takes 8.73 seconds. Of this time, `count_matches` takes 8.69 seconds. The indel searches take the most time, because the indel kernel runs a banded dynamic program (DP) at each position of each sequence.

This change moves the kernels to a Rust extension, and it adds bit-parallel algorithms. A bit-parallel filter finds the positions where an indel match can start. The exact DP then runs only at these positions. The output stays identical, and the Python API does not change. The change also removes numba, llvmlite, the compile time of the first call, and the numba cache.

## Goals

- The output of the numba kernels for each input: the same matches, in the same order, with the same counts.
- A report with `workers=1` that is 5 times faster or more. "Measurement" gives the targets.
- No compile step at run time. The first call in a new process takes the time of a later call.
- No `numba` and no `llvmlite` in the runtime dependencies.
- Prebuilt wheels for Linux x86_64, Linux aarch64, and macOS arm64.
- The same public API: `count_matches`, `match_ranges`, `find_matches`, `SearchOptions`, and each function that calls them.

## Non-goals

- A Rust port of a different part of the package. The profile shows no other slow part. For example, `read_fasta` takes 0.08 seconds of the report.
- A change of a search rule or of a divergence.
- A Rust thread pool. The pool of `seedmatch._pool` stays the only source of parallel work.
- Wheels for free-threaded Python, for Windows, or for macOS x86_64.
- The publication of the wheels to a package index. `INF-01` of the Agent-Ready Package Standard makes it necessary, but it is a separate step with its own consent.

## Anti-goals

- A result that is different from the numba result for one input, also for a rare input. A filter that drops a true match makes the change a failure, also if it is fast.
- A Python fallback for the kernels. Two versions of one rule can give different results.
- An `unsafe` block in the crate.

## Constraints

This change replaces these constraints of `docs/design.md`:

- `docs/design.md:41` names the `uv_build` backend. After this change, the backend is maturin.
- `docs/design.md:51` says that the wheel is pure Python, because the inflexa build images have no Rust toolchain. After this change, the wheel has a compiled extension. The build images still have no Rust toolchain, thus the index must have a prebuilt wheel for each platform of the sandbox.
- `docs/design.md:52` says that numba writes a cache. After this change, the package writes no cache.

These constraints stay:

- The Python dependencies come only from the inflexa lock. The change removes `numba` and `llvmlite`, and it adds no Python dependency. The Rust crates compile into the extension, thus they are not runtime dependencies.
- Python 3.12 or later. The sandbox runs `/usr/bin/python3` of its base image, `rocker/r-ver:4.6.0`.
- The package store has `lock.amd64.json` and `lock.arm64.json`, thus the sandbox runs on Linux x86_64 and on Linux aarch64.
- The coordinates are 1-based and inclusive. The matches come in the order of Biostrings 2.78.0.
- Each action of a workflow has a commit SHA and a version comment.

## The search rules

The kernels obey the search rules of `docs/design.md`. This list gives only the parts that the new algorithms touch. A pattern has m letters, and a sequence has n letters.

- Rule 2: a mismatch match can extend past each end of the sequence. Each pattern letter outside the sequence is a mismatch. The 1-based starts go from 2 − m to n.
- Rule 4: at each position j, i0 is the index of the first pattern letter that matches the letter at j. The kernel aligns the rest of the pattern to the text after j, with the cost d. The cost of j is i0 + d. The held-candidate logic then selects the matches.
- Rule 4b: if max_mm ≥ m, the indel search gives the result of rule 2.
- Rule 5: if max_mm is 0, the search ignores indels.
- An indel match stays inside its sequence.

## Algorithms

### The letter table

The pattern is m masks. `compat_table` gives a table of 16 × 16 values. The kernel makes the letter table from the two: ok[a, s] is true if pattern letter a matches the subject mask s. Each algorithm below reads only this table. Thus each IUPAC letter and each pair of `fixed` flags give the result of today.

For the bit-parallel algorithms, the kernel makes one word of 64 bits for each subject mask s. Bit a of the word is 1 if ok[a, s] is true. Myers names this table Peq, and shift-and names it B.

### Rule 2: shift-and with k mismatches

The kernel uses the shift-and algorithm with substitutions only (Baeza-Yates and Gonnet 1992, Wu and Manber 1992). For k = max_mm, it keeps the state words R0 to Rk. Bit a of Rd is 1 if the pattern letters 0 to a match the last a + 1 text letters with d mismatches or fewer. Each state starts at 0. For each text letter c, the update is:

```text
R0' = ((R0 << 1) | 1) & B[c]
Rd' = (((Rd << 1) | 1) & B[c]) | ((R(d-1) << 1) | 1)     for d = 1 to k
```

If bit m − 1 of Rk is 1 after a letter, a window ends at that letter with k mismatches or fewer. If min_mm is more than 0, a window is a match if bit m − 1 of Rk is 1 and bit m − 1 of R(min_mm − 1) is 0.

Rule 2 counts each pattern letter outside the sequence as a mismatch. Thus the kernel reads m − 1 virtual letters before the sequence and m − 1 virtual letters after it, each with B = 0. The windows of this longer text are the windows of rule 2, from the start 2 − m to the start n. Thus the result is exact, and no window gets a second check.

A window of m letters has m mismatches or fewer. Thus the kernel caps each state index at m, and a large max_mm does not make a large state. The time for a sequence is O((n + 2m) · min(max_mm, m)) word operations.

### Rule 4: a filter, then the exact DP

The current kernel runs the DP at each position j, but most positions cannot start a match. The new kernel first finds the candidates: the positions where a match can start. Then it runs the current DP and the current held-candidate logic only at the candidates, in the increasing order of j.

**The filter.** S(j) is the least edit distance between the pattern and a text part that starts at position j and ends inside the sequence:

```text
S(j) = min over e from j to n − 1 of ED(P, T[j..e])
```

ED uses the letter table: a substitution of pattern letter a with text letter s costs 0 if ok[a, s] is true, and 1 if it is false. A position j is a candidate if S(j) ≤ max_mm.

The kernel computes S(j) with the bit-vector algorithm of Myers (1999). For each end of a text part, Myers gives the least edit distance over each start. Thus the kernel runs Myers on the reversed pattern and on the reversed sequence. The edit distance of two reversed strings is the edit distance of the two strings. Thus each step of the reversed run gives S(j) for one j.

Myers reads only the Peq table, thus each IUPAC mask works. For m ≤ 64, each text letter costs a fixed number of word operations.

**Why the output is identical.** At a position j, the current kernel accepts an alignment of three parts:

1. The pattern letters 0 to i0 − 1 are deletions, with the cost i0.
2. Pattern letter i0 matches text letter j, with the cost 0.
3. The pattern letters after i0 align to the text from j + 1 to j + w, with the cost d.

This is an alignment of the full pattern to T[j..j+w], with the cost i0 + d. Thus ED(P, T[j..j+w]) ≤ i0 + d. The DP keeps j + w inside the sequence, thus S(j) ≤ i0 + d. If the kernel accepts j, then i0 + d ≤ max_mm. Thus S(j) ≤ max_mm, and j is a candidate. As a result, the candidates contain each position that the kernel accepts.

At a position that the current kernel does not accept, the kernel changes nothing: it goes to the next position, and the held candidate stays the same. The new kernel does not visit the positions that are not candidates, and these positions are not accepted. Thus the two kernels accept the same positions, in the same order, with the same costs and ends. The held-candidate logic reads only these values, thus the matches are the same.

The DP at a position reads only the codes and its own row. It writes each value of the row before it reads that value (`kernels.py:169-193`). Thus the DP at a candidate gives the same result when the kernel does not visit the positions before it. The Rust port must keep this property.

The filter does not replace the DP. S(j) can be less than i0 + d, for example if pattern letter 0 is a substitution at j. Thus the DP decides each candidate, and the cost and the end come only from the DP.

**Bounded memory.** The kernel keeps the candidates of one block of positions in a bit set of a fixed size. For the block from b0 to b1 − 1 (0-based), the reversed run reads the text from the position min(n, b1 + m + max_mm − 1) − 1 down to b0. A text part with an edit distance of max_mm or less has m + max_mm letters or fewer. Thus each such part that starts in the block ends at or before the first position of the run. As a result, the decision S(j) ≤ max_mm is exact for each j of the block.

Each block reads m + max_mm − 1 extra letters. The held candidate carries from one block to the next block.

### Patterns of more than 64 letters

A pattern of more than 64 letters does not fit in one word. For such a pattern, the kernel uses the scalar algorithms: the current count with an early stop for rule 2, and the current DP at each position for rule 4. They are Rust ports of the current numba code. `MAX_PATTERN_LENGTH` is 20000, thus these patterns can occur. A siRNA guide has approximately 21 letters, thus a report uses only the bit-parallel path.

### The route of a search

- `indels` off: the kernel of rule 2.
- `indels` on, and max_mm is 0 or max_mm ≥ m: the kernel of rule 2, as now.
- `indels` on, and max_mm is from 1 to m − 1: the kernel of rule 4.
- m ≤ 64: the bit-parallel path. m > 64: the scalar path.

The scalar path exists for each m in the crate, thus the tests can compare the two paths. Shift-and can be slower than the scalar count for a class of search of the default grid. If a measurement shows this, the builder sends this class to the scalar path. The builder then gives the numbers in `docs/performance.md`.

## The extension module

- The crate is in `rust/`: `rust/Cargo.toml`, `rust/Cargo.lock`, and `rust/src/`.
- The extension module is `seedmatch.search._kernels`. The stub file `src/seedmatch/search/_kernels.pyi` gives its types to mypy.
- The change deletes `src/seedmatch/search/kernels.py`, and `api.py` imports `_kernels`.
- The functions keep the names and the arguments of the current kernels:

```text
hamming_counts(codes, offsets, pattern, compat, max_mm, min_mm) -> counts
hamming_ranges(codes, offsets, pattern, compat, max_mm, min_mm) -> (index, starts)
indel_counts(codes, offsets, pattern, compat, max_mm) -> counts
indel_ranges(codes, offsets, pattern, compat, max_mm) -> (index, starts, ends)
```

- The inputs are contiguous numpy arrays: `codes` of uint8, `offsets` of int64, `pattern` of uint8, and `compat` of bool with the shape (16, 16). The limits are int64, as `_checks.as_limit` makes them. The outputs are int64 arrays.
- The extension reads the arrays through rust-numpy, without a copy.
- Each function releases the GIL for the full search, with `Python::detach` of PyO3. Each call makes its own pattern and compat arrays, and `SequenceSet` makes its own arrays read-only.
- A write through a base array of a `SequenceSet` can still change `codes` or `offsets` during a search. The docstring of `SequenceSet` (`sequence_set.py:38-42`) says that such a write can give a wrong result or a crash.
- Thus each function reads each offset one time, and makes sure that the sequence is inside `codes`. If it is not, the function raises `ValueError`. Each table lookup uses `code & 15`. As a result, such a write can give a wrong result, but not a panic.
- A ranges function does one pass and adds each match to its output vectors. The current ranges kernels count first and fill second, thus they do each search two times.
- The crate has `#![forbid(unsafe_code)]`.
- The Python boundary keeps each check of the input, and the kernels trust their caller, as now. The kernels must not panic for an input that the boundary accepts.

## Build and packaging

- The build backend is maturin 1.x.
- `[tool.maturin]` sets `manifest-path = "rust/Cargo.toml"`, `python-source = "src"`, `module-name = "seedmatch.search._kernels"`, and the PyO3 feature `extension-module`.
- PyO3 and rust-numpy have the same minor version, 0.30 at the date of this design. PyO3 has the feature `abi3-py312`. Thus one wheel for each platform works on Python 3.12 and on each later version.
- `rust-toolchain.toml` pins the Rust version. rustup reads it in CI and on a development machine.
- `[tool.uv] cache-keys` lists `pyproject.toml`, `rust/Cargo.toml`, `rust/Cargo.lock`, and `rust/src/**/*.rs`. Without these keys, `uv run` does not build the extension again after a change of the Rust code.
- The wheels:
  - Linux x86_64 and Linux aarch64 with the tag `manylinux_2_28`. Each wheel builds on a runner of its own architecture, `ubuntu-latest` or `ubuntu-24.04-arm`.
  - macOS arm64 on `macos-latest`.
  - An sdist. A build from the sdist needs a Rust toolchain.
- A new workflow, `wheels.yml`, builds the wheels with `PyO3/maturin-action`. It runs on each pull request, on each push to `main`, and on each tag `v*`. Each job installs its wheel in a new environment, and it runs `tests/search` on the installed package, not on the source tree. The workflow keeps the wheels as artifacts of the run, and it does not publish them.

## CI

- `test.yml`: each job gets the toolchain of `rust-toolchain.toml`, and `uv sync --locked` builds the extension. A new step runs `cargo test` in `rust/`. The matrix does not change.
- `lint.yml`: new steps run `cargo fmt --check` and `cargo clippy --all-targets -- -D warnings`.
- `coverage.yml`: remove `NUMBA_DISABLE_JIT`, `PYTHONPATH: tools/coverage`, and `-p serial_pool`, and delete `tools/coverage/serial_pool.py`. The plugin exists only because the kernels ran as Python code on the pool threads, and the tracer then stopped the run. If the run with real threads does not end, keep the plugin, and write the new cause in its docstring.
- `.github/dependabot.yml`: remove `numba` and `llvmlite` from the ignore list. Add the ecosystem `cargo` for `rust/`, with the schedule and the commit prefix of the other entries.

## Tests: the proof that the output is identical

The tests of today stay, and they run on the Rust kernels:

- `tests/search/test_kernels.py`. The Hypothesis property tests compare each kernel with the oracle of `tests/search/oracle.py`. They cover each pair of `fixed` flags, IUPAC patterns, and sets of sequences. Examples are `test_hamming_kernel_equals_oracle`, `test_indel_kernel_equals_oracle`, and `test_ranges_of_a_set_equal_the_oracle_for_each_sequence`.
- `test_kernels_equal_biostrings_on_corpus`. The random corpus and the seed corpus of `tools/golden/make_corpus.py` hold the ranges and the counts of Biostrings 2.78.0. The test compares each group, because the stride of the numba coverage run goes away.
- `tests/search/test_compat_and_oracle.py`, which compares the oracle with R.
- The tests of the report, of the non-canonical searches, and of the e-value, with their golden files from R.

The patterns of these tests have 12 letters or fewer, thus they run the bit-parallel path. New tests cover the parts that these tests do not reach:

- The public kernels against the oracle, with patterns of 63, 64, 65, and 80 letters. The tests use max_mm from 0 to 3, with indels on and off. The patterns of 63 and 64 letters use the bit-parallel path, and the patterns of 65 and 80 letters use the scalar path. The extension has no private function that selects a path.
- A Rust test that compares the bit-parallel path with the scalar path on random input from a fixed seed. The input has each subject mask from 1 to 15 and each pair of `fixed` flags. It has max_mm from 0 to m, min_mm from 0 to max_mm, and sequences of 0 to 300 letters.
- The block size is a parameter of the internal function. Thus the Rust test uses small blocks, and the candidates occur on the two sides of each block limit.
- A Rust test of the filter. For each position that the scalar indel kernel accepts, the test makes sure that the filter marks the position as a candidate.
- A test with max_mm near the largest int64, for rule 2 and for rule 4b. It makes sure that the cap of the state index holds.
- The check of "Measurement": the report of the `rnor7` UTRs is the same before and after the change.

pytest-cov measures the Python code with the floor of 95%, as now. The Rust code gets no coverage gate. The builder measures the line coverage of the Rust tests one time with `cargo llvm-cov`, and gives the number in the pull request.

## Measurement

The numbers come from one machine: an Apple M2 Max with 12 CPUs and Python 3.12, with no other large process. The input is the input of `docs/performance.md`: the canonical 3' UTRs of `rnor7` (18,614 sequences, 17.9 million bases) and the test DE table.

Do these steps on the numba build first, with a filled numba cache. Then do them again after the change:

1. Run a new tool, `tools/perf/grid.py`. For each search of the default grid of `full_search`, it gives the pattern name, max_mm, indels, and the best time of 3 runs on one thread. It does the grid with indels on and with indels off.
2. Run `tools/perf/report.py` with `workers=1` and with the default workers, two times each, as `docs/performance.md` tells.
3. Write the report table to a TSV file, and keep its SHA-256. Also keep the SHA-256 of the `match_ranges` output of each search of the two grids.
4. Measure the time of `import seedmatch` in a new process.

These are the targets:

- Each SHA-256 value after the change is equal to its value before the change. If one value is different, the change is wrong, and no time counts.
- The report with `workers=1` takes 1.75 seconds or less, 5 times less than 8.7 seconds.
- No search of the grids takes more time than before the change.
- The first call in a new process takes at most 0.1 seconds more than a later call.

If a target is not met, the builder gives the numbers to the user. The builder does not change the algorithms in a way that breaks the proof of identical output.

These numbers come from a profile in the session of this design, with scratch scripts on the numba build:

- The grid with indels on takes 9.88 seconds. The largest searches are `18mer` with max_mm 4 (1.50 seconds), `15mer` with max_mm 4 (1.28 seconds), and the full guide with max_mm 4 (0.82 seconds).
- The grid with indels off takes 3.45 seconds.
- The first call compiles the kernels and adds approximately 0.9 seconds. With an empty cache, the report takes 2.86 seconds.

My estimate is that the filter costs 1 to 2 nanoseconds for each letter on this machine. Thus the filter of one search takes approximately 30 milliseconds, and the grid with indels on takes less than 1 second. This is a guess, not a measurement.

`docs/performance.md` gets the new numbers, the method, and the versions. The row "First call, empty numba cache" goes away. The numbers of R stay.

## Files to change

- `pyproject.toml`:
  - the build backend, `[tool.maturin]`, and `[tool.uv] cache-keys`
  - no `numba` in `dependencies`, and no `numba` or `llvmlite` in `constraint-dependencies`
  - no mypy override for `numba.*` or for `seedmatch.search.kernels`, and no ruff ignore for `kernels.py`
- `uv.lock`, without numba and llvmlite.
- `src/seedmatch/search/api.py`: the import of `_kernels`. The comment at line 45 about the numba cache entry goes away.
- `src/seedmatch/report/grid.py:3`: the docstring names the search kernels, not numba.
- `tests/search/test_kernels.py`: the import of `_kernels`, and no `numba.config`.
- `README.md`:
  - "Installation" (lines 405 to 419) says that a wheel needs no Rust toolchain. An install from GitHub or from a clone builds from source with a Rust toolchain.
  - Line 18 says that the package is Python with a compiled search extension.
  - The text of line 421 about the numba cache goes away.
  - The coverage command of line 443 has no `NUMBA_DISABLE_JIT`.
- `docs/design.md`: lines 41, 44, 51, 52, 97, 562, 623, 624, 663, and 698, and the target of "Performance and bounds" without "after the JIT compilation". The tool list adds cargo and maturin. "Calls of the designer" names the Rust kernels and this document.
- `docs/performance.md`: the new numbers and the method.
- `tools/perf/report.py`: no `NUMBA_CACHE_DIR` in the docstring, and `VERSIONS` without numba and llvmlite.
- `CONTRIBUTING.md`: the Rust toolchain and the cargo checks.
- Issue seedmatch#27, "The import error for a numba cache that is not writable does not name NUMBA_CACHE_DIR". The change removes its cause, thus the pull request closes it.

## Delivery order

1. `tools/perf/grid.py`, and the numbers and SHA-256 values of the numba build.
2. The crate, the build, and the scalar ports of the four kernels. The tests of today pass on them. numba goes away in this step.
3. Shift-and for rule 2, with its tests.
4. The Myers filter for rule 4, with the blocks and the tests.
5. The wheels workflow and the CI changes.
6. The numbers after the change, `docs/performance.md`, and the other docs.

Each step ends with ruff, mypy, `cargo fmt`, clippy, and the tests of that step.

## Decisions

### Decided by the user

- A Rust extension for the search kernels, in this repository.
- The bit-parallel algorithms, with the output of the current kernels.
- Prebuilt wheels for Linux x86_64, Linux aarch64, and macOS arm64.
- The same Python API.
- No numba.

### Calls of the designer

- The Myers filter on the reversed text, and the exact DP at each candidate. A forward Myers run can also find the candidates, but it marks the m + max_mm starts before each end. Thus it gives more candidates and more DP.
- Shift-and for rule 2, with virtual letters at the two ends. It is exact, thus it needs no second check.
- One word of 64 bits, and the scalar path for longer patterns. A Myers version with more than one word gives speed only to long patterns, and a report has no long pattern.
- maturin, PyO3, and rust-numpy, with `abi3-py312`. Projects such as `pathway` and `rapier` use rust-numpy with abi3.
- The crate in `rust/`, thus `src/` keeps only the Python package.
- The module name `seedmatch.search._kernels`.
- The tag `manylinux_2_28`, because the manylinux2014 images ended their support with CentOS 7. Each Linux wheel builds on a runner of its own architecture.
- No Rust thread pool, and the GIL released for each call.
- The new CI checks: `cargo test`, `cargo fmt --check`, clippy, and the wheels workflow. Dependabot for cargo.
- No coverage gate for Rust.
- The targets of "Measurement".

### Left to the builder

- The internal names and the module layout of the crate.
- The block size of the filter, within a fixed size.
- The Rust version of `rust-toolchain.toml`, and the versions of maturin, PyO3, and rust-numpy. Use the latest release at the time of the build.
- The seed and the size of the random input of the Rust tests.
- Changes for speed that keep each test green, for example iterators that remove bounds checks.
