#!/usr/bin/env bash
set -e
docker --version
docker info
pip --version
pip install ansible
ansible-galaxy install -r requirements_roles.yml -p roles
mv extra-files/metavisitor/metavisitor_tool_list.yml.fortestonly extra-files/metavisitor/metavisitor_tool_list.yml
printf "galaxy_runs_in_docker: yes\ngalaxy_docker_container_id: standard\n" >> group_vars/all
sudo groupadd -r $GALAXY_TRAVIS_USER -g $GALAXY_GID
sudo useradd -u $GALAXY_UID -r -g $GALAXY_TRAVIS_USER -d $GALAXY_HOME -p travis_testing\
  -c "Galaxy user" $GALAXY_TRAVIS_USER
sudo mkdir $GALAXY_HOME
sudo chown -R $GALAXY_TRAVIS_USER:$GALAXY_TRAVIS_USER $GALAXY_HOME
docker build -t metavisitor .
sudo mkdir /export && sudo chown $GALAXY_UID:$GALAXY_GID /export
sudo mkdir /export2 && sudo chown $GALAXY_UID:$GALAXY_GID /export2
export STANDARD=`docker run -d \
  --name standard \
  -p 80:80 -p 8021:21 -p 8800:8800 --tmpfs /var/run/ \
  --privileged=true \
  -e GALAXY_CONFIG_ALLOW_USER_DATASET_PURGE=True \
  -e GALAXY_CONFIG_ALLOW_LIBRARY_PATH_PASTE=True \
  -e GALAXY_CONFIG_ENABLE_USER_DELETION=True \
  -e GALAXY_CONFIG_ENABLE_BETA_WORKFLOW_MODULES=True \
  -v /tmp/:/tmp/ \
  -v /export/:/export \
  metavisitor`
#export CUSTOM=`docker run -d --tmpfs /var/run/ --tmpfs /tmp/ \
#  --privileged=true \
#  -p 8181:80 \
#  -e NAT_MASQUERADE=true \
#  -e NGINX_GALAXY_LOCATION=/subdir \
#  -v /export2:/export \
#  metavisitor`
#docker ps

