#!/usr/bin/env bash
# Clone or update the local nf-core/ampliseq checkout used by RUN_AMPLISEQ.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO="https://github.com/nf-core/ampliseq.git"
REF="dev"
UPDATE="false"
TARGET_DIR="${ROOT_DIR}/third_party/nf-core-ampliseq"

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --repo URL       Git remote (default: ${REPO})
  --ref REF        Branch, tag, or commit (default: ${REF})
  --update BOOL    Fetch and checkout REF when true (default: ${UPDATE})
  --target-dir DIR Checkout path (default: third_party/nf-core-ampliseq)
  -h, --help       Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo)
      REPO="$2"
      shift 2
      ;;
    --ref)
      REF="$2"
      shift 2
      ;;
    --update)
      UPDATE="$2"
      shift 2
      ;;
    --target-dir)
      TARGET_DIR="$2"
      shift 2
      ;;
    -h | --help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

case "${UPDATE}" in
  true | false) ;;
  *)
    echo "--update must be true or false (got: ${UPDATE})" >&2
    exit 1
    ;;
esac

fetch_ref() {
  # Update remote-tracking refs (origin/<branch>) and tags.
  git -C "${TARGET_DIR}" fetch --tags --prune origin "+refs/heads/*:refs/remotes/origin/*" "+refs/tags/*:refs/tags/*"
  # Also fetch the explicit ref in case it is a commit SHA or odd ref name.
  git -C "${TARGET_DIR}" fetch --tags origin "${REF}" || true
}

# Resolve REF to a commit after fetch. Prefer origin/<ref> for branches so
# --update actually moves past a stale local branch tip.
resolve_commit() {
  if git -C "${TARGET_DIR}" rev-parse --verify --quiet "refs/remotes/origin/${REF}^{commit}" >/dev/null; then
    git -C "${TARGET_DIR}" rev-parse "refs/remotes/origin/${REF}^{commit}"
    return
  fi
  if git -C "${TARGET_DIR}" rev-parse --verify --quiet "refs/tags/${REF}^{commit}" >/dev/null; then
    git -C "${TARGET_DIR}" rev-parse "refs/tags/${REF}^{commit}"
    return
  fi
  if git -C "${TARGET_DIR}" rev-parse --verify --quiet "${REF}^{commit}" >/dev/null; then
    git -C "${TARGET_DIR}" rev-parse "${REF}^{commit}"
    return
  fi
  echo "Could not resolve ampliseq ref: ${REF}" >&2
  exit 1
}

checkout_ref() {
  local commit
  commit="$(resolve_commit)"
  if git -C "${TARGET_DIR}" show-ref --verify --quiet "refs/remotes/origin/${REF}"; then
    # Keep a local branch named after the remote branch, forced to the fetched tip.
    git -C "${TARGET_DIR}" checkout -q -B "${REF}" "${commit}"
  else
    git -C "${TARGET_DIR}" checkout -q --detach "${commit}"
  fi
}

if [[ ! -d "${TARGET_DIR}/.git" ]]; then
  mkdir -p "$(dirname "${TARGET_DIR}")"
  git clone "${REPO}" "${TARGET_DIR}"
  fetch_ref
  checkout_ref
elif [[ "${UPDATE}" == "true" ]]; then
  current_origin="$(git -C "${TARGET_DIR}" remote get-url origin)"
  if [[ "${current_origin}" != "${REPO}" ]]; then
    git -C "${TARGET_DIR}" remote set-url origin "${REPO}"
  fi
  fetch_ref
  checkout_ref
fi

git -C "${TARGET_DIR}" log -1 --format='ampliseq_commit=%H%nampliseq_ref=%D%nampliseq_subject=%s'
