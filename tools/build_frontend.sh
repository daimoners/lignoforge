#!/usr/bin/env bash
# Build the React front-end into lignoforge/web/static (served by `lignoforge gui`).
# Needs Node.js >= 20.19 (or >= 22.12).  Run before `python -m build` for a release,
# or after changing anything under frontend/.  Not needed with a release wheel/sdist
# (they already contain the built bundle).
set -euo pipefail

if ! command -v node >/dev/null 2>&1; then
  echo "error: Node.js not found. Install Node >= 20.19, e.g.  conda install -c conda-forge 'nodejs>=22'" >&2
  exit 1
fi
v=$(node -p "process.versions.node")
major=${v%%.*}; rest=${v#*.}; minor=${rest%%.*}
if [ "$major" -lt 20 ] || { [ "$major" -eq 20 ] && [ "$minor" -lt 19 ]; } || { [ "$major" -eq 21 ]; } \
   || { [ "$major" -eq 22 ] && [ "$minor" -lt 12 ]; }; then
  cat >&2 <<MSG
error: Node.js $v is too old to build the front-end (needs >= 20.19, or >= 22.12).
  Install a newer Node, for example:
    conda install -c conda-forge 'nodejs>=22'        (inside your conda environment)
  or with nvm:
    nvm install 22 && nvm use 22
MSG
  exit 1
fi

cd "$(dirname "$0")/../frontend"
npm ci
npm run build
echo "Front-end bundle written to lignoforge/web/static"
