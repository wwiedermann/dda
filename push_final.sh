#!/usr/bin/env bash
# Force-pushes the cleaned dda history (dda-rewritten-final.bundle) to GitHub.
# Run from Git Bash:  cd ~/Documents/dda && bash push_final.sh
set -euo pipefail
cd "$(dirname "$0")"

BUNDLE="dda-rewritten-final.bundle"
REMOTE="https://github.com/wwiedermann/dda.git"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo ">> Verifying bundle..."
git bundle verify "$BUNDLE"

echo ">> Unpacking bundle into temp repo..."
git init --bare --quiet "$TMP/final.git"
git -C "$TMP/final.git" fetch --quiet "$PWD/$BUNDLE" \
  "refs/heads/*:refs/heads/*" "refs/tags/*:refs/tags/*"

echo ">> Sanity check: authors in rewritten history"
git -C "$TMP/final.git" log --all --format='%an <%ae>' | sort | uniq -c

echo ">> Force-pushing to $REMOTE ..."
git -C "$TMP/final.git" push "$REMOTE" --force \
  refs/heads/main:refs/heads/main \
  refs/heads/Bagging:refs/heads/Bagging \
  refs/heads/Bagging_backup:refs/heads/Bagging_backup \
  refs/heads/dda_robust:refs/heads/dda_robust \
  refs/tags/dda:refs/tags/dda \
  refs/tags/dda0.1.1:refs/tags/dda0.1.1

echo ">> Done. Pushed refs:"
git -C "$TMP/final.git" for-each-ref --format='  %(refname:short) %(objectname:short)'
echo ">> Now re-sync your local working repo:"
echo "   cd ~/Documents/dda && git fetch origin && git reset --hard origin/main"
