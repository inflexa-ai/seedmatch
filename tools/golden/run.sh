#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "usage: tools/golden/run.sh <script.R>" >&2
  exit 2
fi

if [ ! -f "$1" ]; then
  echo "run.sh: no such file: $1" >&2
  exit 2
fi
# Resolve the path before the cd, thus a path relative to the caller's directory stays valid.
script="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"

cd "$(dirname "$0")/../.."
R_LIBS_USER="$PWD/tools/golden/.rlib" exec Rscript --vanilla "$script"
