#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run2.sh test
#     -> reads:  test.satoutput, test.city, test.varmap
#     -> writes: test.metromap
#
#   ./run2.sh path/to/foo.satoutput
#     -> reads:  path/to/foo.satoutput, path/to/foo.city, path/to/foo.varmap
#     -> writes: path/to/foo.metromap

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <base-name|path/to/file.satoutput>"
  exit 1
fi

ARG="$1"

# Determine base and file paths
if [[ "$ARG" == *.satoutput ]]; then
  BASE="${ARG%.satoutput}"
  MODEL="$ARG"
else
  BASE="$ARG"
  MODEL="${BASE}.satoutput"
fi

CITY="${BASE}.city"
VMAP="${BASE}.varmap"
OUT="${BASE}.metromap"
DECODER="./decoder"

# Sanity checks
[[ -f "$MODEL" ]] || { echo "Error: model file not found: $MODEL"; exit 1; }
[[ -f "$CITY"  ]] || { echo "Error: city file not found:  $CITY";  exit 1; }
[[ -f "$VMAP"  ]] || { echo "Error: varmap file not found: $VMAP"; exit 1; }

# Build decoder if missing
if [[ ! -x "$DECODER" ]]; then
  if [[ -f "decoder.cpp" ]]; then
    echo "Building decoder..."
    g++ -std=c++17 -O2 -pipe -Wall -Wextra -o "$DECODER" decoder.cpp
  else
    echo "Error: decoder binary not found and decoder.cpp missing."
    exit 1
  fi
fi

# Run decoder -> write to .metromap
echo "Decoding $MODEL -> $OUT"
"$DECODER" "$MODEL" "$CITY" "$VMAP" > "$OUT"
echo "Done: $OUT"