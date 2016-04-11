#!/usr/bin/env bash

# We start by building and deploying the documentation
SCRIPT_PATH="$(cd "$(dirname "$0")" && pwd -P)"
sudo pip install mkdocs
# Initialize gh-pages checkout
mkdir -p docs/html
(
    cd docs/html
    git init
    git config user.name "${GH_USER_NAME}"
    git config user.email "${GH_USER_EMAIL}"
    git remote add upstream "https://${GH_TOKEN}@${GH_REF}"
    git fetch upstream
    git reset upstream/gh-pages
)

cd $SCRIPT_PATH
mkdocs build --clean --site-dir docs/html

# Commit and push the documentation to gh-pages
(
    cd docs/html
    touch .
    git add -A .
    git commit -m "Rebuild pages at ${rev}"
    git push -q upstream HEAD:gh-pages
)
cd $SCRIPT_PATH

# Next we upload to docker
docker tag galaxy_kickstart artbio/galaxy-kickstart-base:$TRAVIS_COMMIT
docker tag galaxy_kickstart artbio/galaxy-kickstart-base:latest
LOGIN="docker login -e=$DOCKER_EMAIL -u=$DOCKER_USERNAME -p=$DOCKER_PASSWORD"
$LOGIN || (sleep 5s && $LOGIN || echo "login failed twice, quitting" && exit 1)
docker push artbio/galaxy-kickstart-base:$TRAVIS_COMMIT || (sleep 5s && docker push artbio/galaxy-kickstart-base:$TRAVIS_COMMIT || echo "push failed twice, quitting" && exit 1)
docker push artbio/galaxy-kickstart-base:latest || (sleep 5s && docker push artbio/galaxy-kickstart-base:latest || echo "push failed twice, quitting" && exit 1)
