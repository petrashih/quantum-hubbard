#!/usr/bin/env bash
# Run from the root of your local quantum-hubbard clone.
# Usage: bash apply_cleanup.sh /path/to/unzipped/cleanup_files
set -euo pipefail
SRC="${1:?path to unzipped cleanup_files}"
git checkout -b cleanup-roadmap
mkdir -p docs .github/workflows
git mv TODO.md docs/phase1-backlog.md
cp "$SRC/docs/phase1-backlog.md" docs/phase1-backlog.md
cp "$SRC/ROADMAP.md" ROADMAP.md
cp "$SRC/README.md" README.md
cp "$SRC/.github/workflows/tests.yml" .github/workflows/tests.yml
git add -A
git status --short
echo "Review with 'git diff --cached', then commit and push."

