#!/usr/bin/env bash
docker tag galaxy_kickstart artbio/galaxy-kickstart-base:$TRAVIS_COMMIT
docker tag galaxy_kickstart artbio/galaxy-kickstart-base:latest
docker login -e="$DOCKER_EMAIL" -u="$DOCKER_USERNAME" -p="$DOCKER_PASSWORD"
docker push artbio/galaxy-kickstart-base:$TRAVIS_COMMIT
docker push artbio/galaxy-kickstart-base:latest
