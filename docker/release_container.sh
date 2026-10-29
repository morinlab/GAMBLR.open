#!/usr/bin/env bash
# release_container.sh
#
# Admin-only: builds and pushes the gamblr-collective image to GHCR, built
# from envs/gamblr_collective.yaml plus the six GAMBLR packages (see
# docker/Dockerfile / docker/install_gamblr.R). Run this AFTER PRs are
# merged to master on all six repos for a normal release build -- refs can
# be overridden via env vars for pre-merge testing (see docker/RELEASING.md).
#
# Requires:
#  - docker, logged in to ghcr.io with a PAT that has write:packages scope:
#      echo $GHCR_PAT | docker login ghcr.io -u <username> --password-stdin
#
# Usage: ./docker/release_container.sh [tag]
#   (run from the GAMBLR.open repo root -- the build context needs both
#   envs/gamblr_collective.yaml and docker/ at their real repo paths)

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."

IMAGE="ghcr.io/morinlab/gamblr-collective"
PKG_VERSION="$(grep '^Version:' DESCRIPTION | awk '{print $2}')"
TAG="${1:-$PKG_VERSION}"

GAMBLR_DATA_REF="${GAMBLR_DATA_REF:-master}"
GAMBLR_HELPERS_REF="${GAMBLR_HELPERS_REF:-master}"
GAMBLR_UTILS_REF="${GAMBLR_UTILS_REF:-master}"
GAMBLR_VIZ_REF="${GAMBLR_VIZ_REF:-master}"
GAMBLR_PREDICT_REF="${GAMBLR_PREDICT_REF:-master}"
GAMBLR_OPEN_REF="${GAMBLR_OPEN_REF:-master}"

# Same reasoning as GAMBLR.data's release_mutations_db.R: the tag defaults
# to DESCRIPTION's Version, but that's only meaningful if what's committed
# here actually matches what's on origin -- an uncommitted/unpushed bump
# would tag an image nobody else's checkout can reproduce.
if ! git diff --quiet -- DESCRIPTION || ! git diff --cached --quiet -- DESCRIPTION; then
  echo "ERROR: DESCRIPTION has uncommitted changes. Commit and push the version bump first." >&2
  exit 1
fi

ALL_MASTER=true
for ref in "$GAMBLR_DATA_REF" "$GAMBLR_HELPERS_REF" "$GAMBLR_UTILS_REF" \
           "$GAMBLR_VIZ_REF" "$GAMBLR_PREDICT_REF" "$GAMBLR_OPEN_REF"; do
  [[ "$ref" == "master" ]] || ALL_MASTER=false
done

echo "Building ${IMAGE}:${TAG}"
echo "  GAMBLR.data@${GAMBLR_DATA_REF}  GAMBLR.helpers@${GAMBLR_HELPERS_REF}  GAMBLR.utils@${GAMBLR_UTILS_REF}"
echo "  GAMBLR.viz@${GAMBLR_VIZ_REF}  GAMBLR.predict@${GAMBLR_PREDICT_REF}  GAMBLR.open@${GAMBLR_OPEN_REF}"

docker build \
  --build-arg GAMBLR_DATA_REF="${GAMBLR_DATA_REF}" \
  --build-arg GAMBLR_HELPERS_REF="${GAMBLR_HELPERS_REF}" \
  --build-arg GAMBLR_UTILS_REF="${GAMBLR_UTILS_REF}" \
  --build-arg GAMBLR_VIZ_REF="${GAMBLR_VIZ_REF}" \
  --build-arg GAMBLR_PREDICT_REF="${GAMBLR_PREDICT_REF}" \
  --build-arg GAMBLR_OPEN_REF="${GAMBLR_OPEN_REF}" \
  -t "${IMAGE}:${TAG}" \
  -f docker/Dockerfile \
  .

docker push "${IMAGE}:${TAG}"

# Only move the floating "latest" tag for genuine master-based release
# builds -- a pre-merge test build (custom refs) should never become what
# "latest" resolves to for everyone else.
if $ALL_MASTER; then
  docker tag "${IMAGE}:${TAG}" "${IMAGE}:latest"
  docker push "${IMAGE}:latest"
  echo "Tagged and pushed both ${TAG} and latest."
else
  echo "Custom refs used -- pushed ${TAG} only, did not touch latest."
fi

echo
echo "Done. Pull via:"
echo "  docker pull ${IMAGE}:${TAG}"
echo "  apptainer pull docker://${IMAGE}:${TAG}"
