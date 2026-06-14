#!/bin/bash
# Deploy the serverless (shinylive/webR) demo to GitHub Pages on the gh-pages branch.
#
# WARNING: this force-pushes ~66 MB of WebAssembly assets to the gh-pages branch of the
# PUBLIC repository, where they remain in history. The site is regenerable
# (toolkit/app-lite/build_shinylive.R), so it need not live in the main branch.
#
# After running, set the repo's Pages source to the gh-pages branch:
#   GitHub > Settings > Pages > Build and deployment > Source: "Deploy from a branch" > gh-pages / (root)
# The demo will then be live at:
#   https://pablobernabeu.github.io/Experimental-data-simulation/
#
# Run from the repository root:  bash toolkit/app-lite/deploy_ghpages.sh
set -euo pipefail

SITE="toolkit/build/shinylive-demo"
REMOTE="$(git remote get-url origin)"

if [ ! -f "$SITE/index.html" ]; then
  echo "Building the shinylive site first..."
  Rscript toolkit/app-lite/build_shinylive.R
fi
touch "$SITE/.nojekyll"   # serve static files as-is (don't run Jekyll)

# Build an orphan gh-pages commit in a scratch dir so the main repo history is untouched.
WORK="$(mktemp -d)"
cp -r "$SITE"/. "$WORK"/
( cd "$WORK"
  git init -q
  git checkout -q -b gh-pages
  git add -A
  git -c user.name="pilotr deploy" -c user.email="deploy@local" commit -q -m "Deploy pilotr shinylive demo"
  git push -f "$REMOTE" gh-pages )
rm -rf "$WORK"

echo "Pushed gh-pages to $REMOTE"
echo "Now set Settings > Pages > Source = gh-pages branch, then visit:"
echo "  https://pablobernabeu.github.io/Experimental-data-simulation/"
