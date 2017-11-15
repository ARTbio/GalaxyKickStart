#!/usr/bin/env bash
set -e
ansible-playbook -i inventory_files/galaxy-kickstart galaxy.yml --connection=local --sudo
sleep 60s
curl http://localhost:80/api/version| grep version_major
date > $HOME/date.txt && curl --fail -T $HOME/date.txt ftp://localhost:21 --user $GALAXY_USER:$GALAXY_USER_PASSWD
sudo -E su $GALAXY_TRAVIS_USER -c "export PATH=$GALAXY_HOME/.local/bin/:$PATH &&
cd $GALAXY_HOME &&
bioblend-galaxy-tests -v $GALAXY_HOME/.local/lib/python2.7/site-packages/bioblend/_tests/TestGalaxy*.py"
curl --fail $BIOBLEND_GALAXY_URL/api/version
