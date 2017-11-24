#!/usr/bin/env bash
set -e
docker --version
docker info
pip install ansible
ansible-galaxy install -r requirements_roles.yml -p roles
sudo groupadd -r $GALAXY_TRAVIS_USER -g $GALAXY_GID
sudo useradd -u $GALAXY_UID -r -g $GALAXY_TRAVIS_USER -d $GALAXY_HOME -p travis_testing\
  -c "Galaxy user" $GALAXY_TRAVIS_USER
sudo mkdir $GALAXY_HOME
sudo chown -R $GALAXY_TRAVIS_USER:$GALAXY_TRAVIS_USER $GALAXY_HOME
docker build -t metavisitor .
sudo mkdir /export && sudo chown $GALAXY_UID:$GALAXY_GID /export
sudo mkdir /export2 && sudo chown $GALAXY_UID:$GALAXY_GID /export2
export CUSTOM=`docker run -d --privileged=true --tmpfs /var/run/ \
  -p 8181:80 \
  -e NAT_MASQUERADE=true \
  -e NGINX_GALAXY_LOCATION=/subdir \
  -v /export2:/export \
  metavisitor`
export STANDARD=`docker run -d --tmpfs /var/run/ \
  -p 80:80 -p 8021:21 -p 8800:8800 \
  --privileged=true \
  -e GALAXY_CONFIG_ALLOW_USER_DATASET_PURGE=True \
  -e GALAXY_CONFIG_ALLOW_LIBRARY_PATH_PASTE=True \
  -e GALAXY_CONFIG_ENABLE_USER_DELETION=True \
  -e GALAXY_CONFIG_ENABLE_BETA_WORKFLOW_MODULES=True \
  -v /tmp/:/tmp/ \
  -v /export/:/export \
  metavisitor`
docker ps

