#!/usr/bin/env bash
set -e
# Usage:
#   check_version.sh
#
# Reads version from Cargo.toml and checks it against tags
#
# Modified from
# https://github.com/mrc-ide/dust/blob/master/scripts/version_check
# by  @richfitz for the dust package.
#
VERSION=${1:-$(grep '^version' Cargo.toml  | sed 's/.*= *//' | sed 's/"//g')}
TAG="v${VERSION}"

echo "Proposed version number '$VERSION'"

if echo "$VERSION" | grep -Eq "[0-9]+[.][0-9]+[.][0-9]+"; then
    echo "[OK] Version number in correct format"
else
    echo "[ERROR] Invalid format version number '$VERSION' must be in format 'x.y.z'"
    exit 1
fi

EXIT_CODE=0

echo "Updating remote git data"
git fetch --quiet

BRANCH_DEFAULT=$(git remote show origin | awk '/HEAD branch/ {print $NF}')
LAST_TAG=$(git describe --tags --abbrev=0 "origin/${BRANCH_DEFAULT}")

echo "Pushed tag is $LAST_TAG"

if git rev-parse "$TAG" >/dev/null 2>&1; then
    echo "[OK] Version number $VERSION in Cargo.toml matches tag $TAG"
else
    echo "[ERROR] Tag $TAG does not match version $VERSION in Cargo.toml - update version number in Cargo.toml"
    exit 1
fi
