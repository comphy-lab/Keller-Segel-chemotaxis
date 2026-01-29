#!/bin/bash

set -euo pipefail

if [[ -z "${1:-}" ]]; then
  echo "Usage: $0 <case-name>" >&2
  exit 1
fi

CASE_NAME="$1"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CASE_DIR="$SCRIPT_DIR/$CASE_NAME"

rm -rf "$CASE_DIR"
