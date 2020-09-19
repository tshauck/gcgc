#!/usr/bin/env sh
#######################################################################
#                     Make the release for GCGC.                      #
#######################################################################

LAST_TAG="$(git tag --sort=-v:refname | head -n 2 | tail -n 1)"
CURRENT_TAG="$(git tag --sort=-v:refname | head -n 1)"

cz ch --start-rev "$LAST_TAG"

echo "Making release with $LAST_TAG and $CURRENT_TAG."
echo "Current CHANGELOG:"
cat ./CHANGELOG.md

if [ -z "$GITHUB_TOKEN" ]; then
    echo "GITHUB_TOKEN needs to be set."
    exit 1
fi
echo "GITHUB_TOKEN is set."
gh release create "$CURRENT_TAG" -F CHANGELOG.md
