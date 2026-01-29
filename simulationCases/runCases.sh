#!/bin/bash

set -euo pipefail

if [[ -z "${1:-}" ]]; then
  echo "Usage: $0 <case-name>" >&2
  exit 1
fi

CASE_NAME="$1"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
CASE_DIR="$SCRIPT_DIR/$CASE_NAME"
CASE_SOURCE="$SCRIPT_DIR/$CASE_NAME.c"

mkdir -p "$CASE_DIR"

cp "$CASE_SOURCE" "$CASE_DIR/"

(
  cd "$REPO_ROOT"
  qcc -I"$REPO_ROOT/src-local" -O2 -Wall -disable-dimensions "$CASE_SOURCE" -o "$CASE_DIR/$CASE_NAME" -lm
)
(
  cd "$CASE_DIR"
  "./$CASE_NAME"
)
