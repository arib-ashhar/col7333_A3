#!/usr/bin/env bash
set -euo pipefail

ENC_SRC="encoder.cpp"
DEC_SRC="decoder.cpp"
ENC_BIN="encoder"
DEC_BIN="decoder"

if [[ ! -f "$ENC_SRC" ]]; then
  echo "Error: $ENC_SRC not found in $(pwd)" >&2
  exit 1
fi

if [[ ! -f "$DEC_SRC" ]]; then
  echo "Error: $DEC_SRC not found in $(pwd)" >&2
  exit 1
fi

echo "Compiling $ENC_SRC -> $ENC_BIN ..."
g++ -std=c++17 -O2 -pipe -Wall -Wextra -Wpedantic -o "$ENC_BIN" "$ENC_SRC"

echo "Compiling $DEC_SRC -> $DEC_BIN ..."
g++ -std=c++17 -O2 -pipe -Wall -Wextra -Wpedantic -o "$DEC_BIN" "$DEC_SRC"

echo "OK: built ./$ENC_BIN and ./$DEC_BIN"