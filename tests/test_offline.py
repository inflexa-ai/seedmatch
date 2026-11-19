import os
import subprocess
import sys
from pathlib import Path

from helpers import run_core_pipeline

TESTS = Path(__file__).parent
NETWORK_MODULES = (
    "seedmatch.fetch",
    "urllib.request",
    "http.client",
    "ssl",
    "urllib3",
    "requests",
    "httpx",
)


def test_core_pipeline_runs_with_sockets_disabled():
    run_core_pipeline()  # pytest-socket makes each connection attempt raise


def test_core_pipeline_loads_no_network_module():
    # socket is not in the list, because numpy and pandas import it.
    code = (
        "import sys, helpers; helpers.run_core_pipeline(); "
        f"print([m for m in {NETWORK_MODULES!r} if m in sys.modules])"
    )
    # The arguments are the interpreter path and a constant script, so no input reaches the command.
    done = subprocess.run(  # noqa: S603
        [sys.executable, "-c", code],
        check=False,
        capture_output=True,
        text=True,
        cwd=TESTS,
        env=os.environ | {"PYTHONPATH": str(TESTS)},
    )
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip() == "[]"
