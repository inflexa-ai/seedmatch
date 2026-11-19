<!-- Keep each pull request to one change. CONTRIBUTING.md gives the rules. -->

## What and why

<!-- What does this change do, and why? -->

Closes #<!-- the issue number, if there is one -->

## Checklist

- [ ] `uv run ruff check`, `uv run ruff format --check`, `uv run mypy`, and `uv run pytest` pass.
- [ ] A test covers each change of behavior.
- [ ] Each new change of behavior from SeedMatchR has an entry in `docs/divergences.md`, and a test name has its ID.
- [ ] If an R golden file changed, a script in `tools/golden/` made it, and the description gives the reason.
- [ ] The documentation shows each change that a user can see.
- [ ] Each commit has a `Signed-off-by` trailer for the DCO (CONTRIBUTING.md).
