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

# Next we upload to docker
docker tag galaxy_kickstart artbio/galaxy-kickstart-base:$TRAVIS_COMMIT
docker tag galaxy_kickstart artbio/galaxy-kickstart-base:latest
LOGIN="docker login -e=$DOCKER_EMAIL -u=$DOCKER_USERNAME -p=$DOCKER_PASSWORD"
$LOGIN || (sleep 5s && $LOGIN || echo "login failed twice, quitting" && exit 1)
docker push artbio/galaxy-kickstart-base:$TRAVIS_COMMIT || (sleep 5s && docker push artbio/galaxy-kickstart-base:$TRAVIS_COMMIT || echo "push failed twice, quitting" && exit 1)
docker push artbio/galaxy-kickstart-base:latest || (sleep 5s && docker push artbio/galaxy-kickstart-base:latest || echo "push failed twice, quitting" && exit 1)
