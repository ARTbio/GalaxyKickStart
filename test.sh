#!/bin/bash
# this script in github actions is accessible at the env variable $GITHUB_WORKSPACE
# if executed avec action checkout

set -euo pipefail

export GALAXY_TRAVIS_USER="galaxy"
export GALAXY_UID="1450"
export GALAXY_GID="1450"
export GALAXY_HOME="/home/galaxy"
export GALAXY_USER="admin@galaxy.org"
export GALAXY_USER_EMAIL="admin@galaxy.org"
export GALAXY_USER_PASSWD="artbio2020"
export GALAXY_VERSION=release_20.01
export BIOBLEND_GALAXY_API_KEY="artbio2020"
export BIOBLEND_TEST_JOB_TIMEOUT="240"
export BIOBLEND_GALAXY_URL="http://127.0.0.1:80"

sudo rm -f /etc/boto.cfg # to do: understand the purpose of this step

cd /opt/hostedtoolcache/Python/3.7.10/x64/lib/python3.7/site-packages/bioblend/_tests/

bioblend-galaxy-tests --color=yes TestGalaxy*.py --no-summary || true
