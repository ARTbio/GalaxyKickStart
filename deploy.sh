#!/usr/bin/env bash
docker tag galaxy_kickstart artbio/galaxy-kickstart-base:$TRAVIS_COMMIT
docker tag galaxy_kickstart artbio/galaxy-kickstart-base:latest
LOGIN="docker login -e=$DOCKER_EMAIL -u=$DOCKER_USERNAME -p=$DOCKER_PASSWORD"
$LOGIN || (sleep 5s && $LOGIN || echo "login failed twice, quitting" && exit 1)
docker push artbio/galaxy-kickstart-base:$TRAVIS_COMMIT || (sleep 5s && docker push artbio/galaxy-kickstart-base:$TRAVIS_COMMIT || echo "push failed twice, quitting" && exit 1)
docker push artbio/galaxy-kickstart-base:latest || (sleep 5s && docker push artbio/galaxy-kickstart-base:latest || echo "push failed twice, quitting" && exit 1)
