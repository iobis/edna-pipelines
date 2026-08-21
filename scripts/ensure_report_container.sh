#!/usr/bin/env bash
# Ensure the summary-report Docker image exists for the current Dockerfile.
# Tags as edna-pipelines-report:<dockerfile-sha256-12>, rebuilds only when the
# Dockerfile changes (or the tagged image is missing). Prints the image name.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DOCKERFILE="${1:-${ROOT}/containers/Dockerfile}"
CONTEXT="$(cd "$(dirname "$DOCKERFILE")" && pwd)"

if [[ ! -f "$DOCKERFILE" ]]; then
  echo "Dockerfile not found: ${DOCKERFILE}" >&2
  exit 1
fi

HASH="$(
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$DOCKERFILE" | cut -c1-12
  else
    shasum -a 256 "$DOCKERFILE" | cut -c1-12
  fi
)"
IMAGE="edna-pipelines-report:${HASH}"

if docker image inspect "$IMAGE" >/dev/null 2>&1; then
  echo "Using existing report image ${IMAGE}" >&2
else
  echo "Building report image ${IMAGE} from ${DOCKERFILE}" >&2
  docker build \
    --platform linux/amd64 \
    -t "${IMAGE}" \
    -f "${DOCKERFILE}" \
    "${CONTEXT}"
fi

echo "$IMAGE"
