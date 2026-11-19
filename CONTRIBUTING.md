# Contribute to seedmatch

seedmatch is a Python port of [SeedMatchR](https://github.com/tacazares/SeedMatchR). On the core paths, seedmatch gives the same results as SeedMatchR. Each difference from SeedMatchR is a recorded divergence in [`docs/divergences.md`](docs/divergences.md).

## Before you start

- Search the [issues](../../issues) for a duplicate.
- For a large change, open an issue first. A new capability, a change to the public API, and a new divergence are large changes.

## Set up the project

Install [uv](https://docs.astral.sh/uv/). uv installs the correct Python version and the dependencies:

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
```

The tests disable the sockets. A test that opens a network connection fails.

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

## Developer Certificate of Origin

seedmatch uses the [Developer Certificate of Origin](https://developercertificate.org/) (DCO), version 1.1. Each commit must have a `Signed-off-by` line. The `-s` option of `git commit` adds the line:

```bash
git commit -s -m "fix(search): count the matches at the end of a sequence"
```

A workflow does a check of each commit for the line. If you forget it, use `git commit --amend -s`, or a rebase for more commits.

Your contribution has the license in [`LICENSE`](LICENSE), MIT.

## Security problems

Do not report a security problem in a public issue. Obey [`SECURITY.md`](SECURITY.md).
