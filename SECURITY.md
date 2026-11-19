# Security policy

seedmatch is a Python library from Inflexa, Inc. It reads FASTA files and DESeq2 tables from the paths that the caller gives. The core modules make no network connection, because they run in sandboxes that have no network egress. A ruff rule in `pyproject.toml` stops the import of a network module outside `seedmatch.fetch`, and the tests run with the sockets disabled.

## Supported versions

Before version 1.0, only the latest release gets security fixes.

## Report a vulnerability

Do not report a security vulnerability in a public issue, a pull request, or a discussion.

Use one of these channels:

- **GitHub private vulnerability reporting.** On the repository, open the **Security** tab, then click **Report a vulnerability**. This is the preferred channel.
- **Email** to security@inflexa.ai.

Give this information, if you can:

- the seedmatch version, or the commit SHA
- the Python version, the operating system, and the CPU architecture
- a description of the problem and its security effect
- the steps that reproduce the problem, with a small input file

Remove private data and secrets from each file that you send.

## What to expect

- We acknowledge the report in 3 business days.
- We do a triage of the severity, and we tell you the result.
- We release a fix on a schedule that agrees with the severity. We agree the disclosure date with you.
- If you agree, the advisory gives your name.

## In scope

- An input file that causes code execution, a file write, or a network connection.
- A network connection from a module outside `seedmatch.fetch`.
- A small input that makes the time or the memory grow much more than the size of the input.
- The integrity of a published package, or of a workflow in `.github/workflows/`.

## Out of scope

- A wrong scientific result with no security effect. Use a bug report for it.
- A vulnerability in a dependency, if seedmatch gives no path to exploit it.
- A problem that is possible only on a host that is already compromised.

## Disclosure

When a fix is available, we publish each confirmed vulnerability as a GitHub security advisory.
