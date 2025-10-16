#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <base-name>"
  echo "Example: $0 test   # will read test.city"
  exit 1
fi

BASE="$1"
INPUT="${BASE}.city"
BIN="./encoder"

if [[ ! -f "$INPUT" ]]; then
  echo "Error: input file '$INPUT' not found."
  exit 1
fi

# Build if binary missing but source exists
if [[ ! -x "$BIN" ]]; then
  if [[ -f "encoder.cpp" ]]; then
    echo "Info: '$BIN' not found. Building via ./compile.sh ..."
    bash ./compile.sh
  else
    echo "Error: '$BIN' not found and encoder.cpp missing to build it."
    exit 1
  fi
fi

echo "Running $BIN with input '$INPUT' ..."
"$BIN" "$BASE"