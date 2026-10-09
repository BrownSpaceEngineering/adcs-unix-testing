#!/usr/bin/env bash
# Regenerate include/matlab_reference.h from the team's MATLAB code (headless), then rebuild and
# run the C tests. Usage: tools/matlab/run_matlab_parity.sh [path/to/adcs]
# MATLAB is found via $MATLAB, then PATH, then /Applications/MATLAB_*.app.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
adcs="${1:-$(cd "$repo/.." && pwd)/adcs}"

matlab="${MATLAB:-$(command -v matlab || true)}"
if [ -z "$matlab" ]; then
    matlab="$(ls -d /Applications/MATLAB_*.app/bin/matlab 2>/dev/null | sort | tail -n 1 || true)"
fi
if [ -z "$matlab" ]; then
    echo "MATLAB not found (set MATLAB=/path/to/matlab)" >&2
    exit 2
fi

echo "Using $matlab and $adcs"
"$matlab" -batch "addpath('$here'); generate_matlab_reference('$adcs', '$repo/include/matlab_reference.h')"

cmake -S "$repo" -B "$repo/build" >/dev/null
cmake --build "$repo/build" -j4
"$repo/build/adcs-test"
