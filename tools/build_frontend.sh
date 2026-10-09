#!/usr/bin/env bash
# Build the React front-end into lignoforge/web/static (served by `lignoforge gui`).
# Needs Node.js >= 18.  Run before `python -m build` for a release, or after
# changing anything under frontend/.
set -euo pipefail
cd "$(dirname "$0")/../frontend"
npm ci
npm run build
echo "Front-end bundle written to lignoforge/web/static"
