#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET_DIR="${ROOT_DIR}/third_party/nf-core-ampliseq"

if [[ -d "${TARGET_DIR}/.git" ]]; then
  echo "nf-core/ampliseq already present at ${TARGET_DIR}"
  exit 0
fi

mkdir -p "${ROOT_DIR}/third_party"
git clone --branch dev --single-branch https://github.com/nf-core/ampliseq.git "${TARGET_DIR}"
echo "Cloned nf-core/ampliseq dev branch to ${TARGET_DIR}"
