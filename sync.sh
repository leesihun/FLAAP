#!/bin/bash
# Commit all changes, rebase on remote, push. Run from Git Bash: bash sync.sh ["message"]
cd "$(dirname "$0")" || exit 1

git add .

if ! git diff --cached --quiet; then
    MSG="${1:-$(date '+%Y-%m-%d %H:%M:%S')}"
    git commit -m "$MSG"

    echo "Pulling remote changes..."
    if ! git pull --rebase; then
        echo "Warning: Failed to pull from remote. You may need to manually sync."
        exit 1
    fi

    echo "Pushing to remote..."
    if timeout 30 git push 2>&1; then
        echo "Successfully pushed to remote!"
    else
        echo "Warning: Failed to push to remote. Changes are committed locally."
        echo "Run 'git push' manually when ready."
    fi
else
    echo "No changes to commit."
fi
