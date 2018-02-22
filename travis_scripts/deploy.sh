#!/usr/bin/env bash

SCRIPT_PATH="$(cd "$(dirname "$0")" && pwd -P)"
sudo pip install mkdocs
mkdir mkdocs_build
cd mkdocs_build
# Initialize gh-pages checkout
DATE=`date`
git clone https://github.com/ARTbio/GalaxyKickStart.git
cd GalaxyKickStart
git config credential.helper "store --file=.git/credentials"
echo "https://${GH_TOKEN}:@github.com" > .git/credentials
mkdocs gh-deploy --clean -m "gh-deployed by travis $DATE"
cd $SCRIPT_PATH
