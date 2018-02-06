#!/usr/bin/env bash
set -e
docker --version
docker info
pip install ansible==2.2.0.0
ansible-galaxy install -r requirements_roles.yml -p roles -f
sudo groupadd -r $GALAXY_TRAVIS_USER -g $GALAXY_GID
sudo useradd -u $GALAXY_UID -r -g $GALAXY_TRAVIS_USER -d $GALAXY_HOME -p travis_testing\
  -c "Galaxy user" $GALAXY_TRAVIS_USER
sudo mkdir $GALAXY_HOME
sudo chown -R $GALAXY_TRAVIS_USER:$GALAXY_TRAVIS_USER $GALAXY_HOME
docker build -t metavisitor -f Dockerfile.test .
sudo mkdir /export && sudo chown $GALAXY_UID:$GALAXY_GID /export
sudo mkdir /export2 && sudo chown $GALAXY_UID:$GALAXY_GID /export2
export CID1=`docker run -d --privileged=true -p 80:80 -p 21:21\
  -e NAT_MASQUERADE=true \
  -e NGINX_GALAXY_LOCATION=/subdir \
  -v /export:/export \
  -v /tmp/:/tmp/ \
  metavisitor`

docker ps

