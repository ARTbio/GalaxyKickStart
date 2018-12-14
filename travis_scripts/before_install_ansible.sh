#!/usr/bin/env bash
set -e
sudo /etc/init.d/postgresql stop
sudo apt-get -y --purge remove postgresql libpq-dev libpq5 postgresql-client-common postgresql-common
sudo rm -rf /var/lib/postgresql
ansible-galaxy install -r requirements_roles.yml -p roles
sudo groupadd -r $GALAXY_TRAVIS_USER -g $GALAXY_GID
sudo useradd -u $GALAXY_UID -r -g $GALAXY_TRAVIS_USER -d $GALAXY_HOME -p travis_testing -c "Galaxy user" $GALAXY_TRAVIS_USER
sudo mkdir $GALAXY_HOME
sudo chown -R $GALAXY_TRAVIS_USER:$GALAXY_TRAVIS_USER $GALAXY_HOME

