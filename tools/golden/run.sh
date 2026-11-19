#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "usage: tools/golden/run.sh <script.R>" >&2
  exit 2
fi

cd "$(dirname "$0")/../.."
R_LIBS_USER="$PWD/tools/golden/.rlib" exec Rscript --vanilla "$1"
