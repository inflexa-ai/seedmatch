# Contribute to seedmatch

seedmatch is a Python port of [SeedMatchR](https://github.com/tacazares/SeedMatchR). On the core paths, seedmatch gives the same results as SeedMatchR. Each difference from SeedMatchR is a recorded divergence in [`docs/divergences.md`](docs/divergences.md).

## Before you start

- Search the [issues](../../issues) for a duplicate.
- For a large change, open an issue first. A new capability, a change to the public API, and a new divergence are large changes.

## Set up the project

Install [uv](https://docs.astral.sh/uv/) and [rustup](https://rustup.rs). uv installs the correct Python version and the dependencies. rustup installs the Rust version of `rust-toolchain.toml`, and `uv sync` builds the Rust extension:

```bash
git clone https://github.com/inflexa-ai/seedmatch.git
cd seedmatch
uv sync
```

## Do the checks

Before you push, make sure that these commands pass:

```bash
uv run ruff check
uv run ruff format --check
uv run mypy
uv run pytest
cargo test --manifest-path rust/Cargo.toml
cargo fmt --manifest-path rust/Cargo.toml --check
PYO3_PYTHON="$PWD/.venv/bin/python" cargo clippy --manifest-path rust/Cargo.toml --all-targets --all-features -- -D warnings
```

The feature `python` of the crate builds PyO3, thus clippy needs the Python interpreter of the environment in `PYO3_PYTHON`.

`uv sync` installs the versions of the inflexa sandbox, which `constraint-dependencies` in `pyproject.toml` pins. The CI also runs the tests with the latest versions, for example pandas 3. To do this test on your machine, run these commands:

```bash
uv venv tmp/latest
VIRTUAL_ENV=tmp/latest uv pip install --no-config -e ".[plot]" --group dev pyarrow
tmp/latest/bin/python -m pytest
```

`--no-config` makes uv skip the pins in `pyproject.toml`.

The tests disable the sockets. A test that opens a network connection fails.

### Do the tests on a real reference

The tests with the marker `refs` compare seedmatch with R on the rat reference `rnor7` and on the human reference `hg38`. They also compare the GTF path with the EnsDb path. These tests read the files of the reference from the directory of the variable `SEEDMATCH_TEST_REFS`. [`tools/golden/README.md`](tools/golden/README.md) gives the layout of the directory and the source of each file. To do only these tests, run this command:

```bash
SEEDMATCH_TEST_REFS=~/Development/refs/seedmatch uv run pytest -m refs -v
```

These tests take approximately 4 minutes. They use up to 4 GiB of memory. The tests of `hg38` keep its GTF reference and its EnsDb reference open at the same time. If the variable is set, `uv run pytest` also runs them.

If the variable is not set, each of these tests skips with the reason `SEEDMATCH_TEST_REFS is not set`. If a file is absent, each test that reads the file skips, and the reason gives the path of the file. The CI workflows do not set the variable, thus CI skips these tests.

## Rules for a change

- **Results from R.** The golden files in `tests/golden/` hold the results of SeedMatchR. Do not edit a golden file by hand. To make or to update a golden file, use the scripts in `tools/golden/`. [`tools/golden/README.md`](tools/golden/README.md) gives the commands.
- **Divergences.** If a change makes a result different from SeedMatchR, add an entry to `docs/divergences.md`. Give each of the six fields of the entry. Put the ID of the entry in the name of the test, for example `test_D001_...`.
- **No network in the core.** The core reads each resource from a path. Network access belongs only in `seedmatch.fetch`.
- **Tests.** A test covers each change of behavior. For a bug fix, write the test that fails before the fix.
- **Documentation.** If a user can see the change, update `README.md` or the documents in `docs/`.

## Make a pull request

- Fork the repository and make a branch from `main`.
- Keep each pull request to one change.
- Use [Conventional Commits](https://www.conventionalcommits.org/) for the commit messages, for example `fix(search): count the matches at the end of a sequence`.
- In the description, give what changed and why. Give the link to the issue.

## Make a release

A maintainer makes each release:

1. In a pull request, change `version` in `pyproject.toml`.
2. After the merge, add the tag `v` and the version to the commit on `main`, for example `v0.2.0`. Push the tag.
3. The job `publish` of `.github/workflows/wheels.yml` uploads the wheels and the sdist to PyPI. Before the upload, the job makes sure that each file has the version of the tag.
4. Write the changes of the release in a GitHub release for the tag.

PyPI keeps each version for all time. Thus, a version that has an error gets a new version, not a new upload.

The PyPI project trusts the job `publish` in the GitHub environment `pypi`. Thus the job uses no token.

## Developer Certificate of Origin

seedmatch uses the [Developer Certificate of Origin](https://developercertificate.org/) (DCO), version 1.1. Each commit must have a `Signed-off-by` line. The `-s` option of `git commit` adds the line:

```bash
git commit -s -m "fix(search): count the matches at the end of a sequence"
```

A workflow does a check of each commit for the line. If you forget it, use `git commit --amend -s`, or a rebase for more commits.

Your contribution has the license in [`LICENSE`](LICENSE), MIT.

## Security problems

Do not report a security problem in a public issue. Obey [`SECURITY.md`](SECURITY.md).
