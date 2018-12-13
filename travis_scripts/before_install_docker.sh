#!/usr/bin/env bash
set -e
docker --version
docker info
ansible-galaxy install -r requirements_roles.yml -p roles
sudo groupadd -r $GALAXY_TRAVIS_USER -g $GALAXY_GID
sudo useradd -u $GALAXY_UID -r -g $GALAXY_TRAVIS_USER -d $GALAXY_HOME -p travis_testing\
  -c "Galaxy user" $GALAXY_TRAVIS_USER
sudo mkdir $GALAXY_HOME
sudo chown -R $GALAXY_TRAVIS_USER:$GALAXY_TRAVIS_USER $GALAXY_HOME
docker build -t galaxy_kickstart -f Dockerfile.galaxykickstart-base .
sudo mkdir /export && sudo chown $GALAXY_UID:$GALAXY_GID /export
sudo mkdir /export2 && sudo chown $GALAXY_UID:$GALAXY_GID /export2
export CID1=`docker run -d --privileged=true -p 80:80 -p 21:21\
  -e NAT_MASQUERADE=true \
  -e NGINX_GALAXY_LOCATION=/subdir \
  -v /tmp/:/tmp/ \
  -v /export:/export \
  galaxy_kickstart`

export CID2=`docker run -d -p 8080:80 -p 8021:21 -p 8800:8800 \
  --privileged=true \
  -e GALAXY_CONFIG_ALLOW_USER_DATASET_PURGE=True \
  -e GALAXY_CONFIG_ALLOW_LIBRARY_PATH_PASTE=True \
  -e GALAXY_CONFIG_ENABLE_USER_DELETION=True \
  -e GALAXY_CONFIG_ENABLE_BETA_WORKFLOW_MODULES=True \
  -v /tmp/:/tmp/ \
  -v /export2/:/export \
  galaxy_kickstart`
docker ps

# test without   -v /tmp/:/tmp/ \
