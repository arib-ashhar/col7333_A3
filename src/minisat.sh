#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./minisat.sh test
# runs:
#   ./MiniSat_v1.14_linux test.satinput test.satoutput

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <base-filename>"
  echo "Example: $0 test   # uses test.satinput -> test.satoutput"
  exit 1
fi

BASE="$1"
BIN="./MiniSat_v1.14_linux"
IN="${BASE}.satinput"
OUT="${BASE}.satoutput"

# sanity checks
[[ -x "$BIN" ]] || { echo "Error: MiniSAT binary not found or not executable at $BIN"; exit 1; }
[[ -f "$IN"  ]] || { echo "Error: input CNF not found: $IN"; exit 1; }

echo "Running MiniSAT: $BIN $IN $OUT"
"$BIN" "$IN" "$OUT"
echo "Done. Wrote: $OUT"