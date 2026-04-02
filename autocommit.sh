#!/bin/bash
cd /Users/sairshaikh/Math/ || exit 1

git add .

# Only commit if there are changes
if ! git diff --cached --quiet; then
    git commit -m "Auto-commit: $(date '+%Y-%m-%d %H:%M:%S')"
    git push
fi