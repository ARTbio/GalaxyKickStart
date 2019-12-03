#!/usr/bin/env bash
set -e
ansible-playbook -i inventory_files/galaxy-kickstart --skip-tags=install_tools galaxy.yml
ansible-playbook -i inventory_files/galaxy-kickstart --tags=install_tools galaxy.yml
sleep 10s
ansible-playbook -i inventory_files/galaxy-kickstart --tags=install_tools galaxy.yml
sleep 10s
curl http://localhost:80/api/version| grep version_major
date > $HOME/date.txt && curl --fail -T $HOME/date.txt ftp://localhost:21 --user $GALAXY_USER:$GALAXY_USER_PASSWD
# sudo -E su $GALAXY_TRAVIS_USER -c "export PATH=$GALAXY_HOME/.local/bin/:$PATH &&
# cd $GALAXY_HOME &&
# bioblend-galaxy-tests -v $GALAXY_HOME/.local/lib/python2.7/site-packages/bioblend/_tests/TestGalaxy*.py"

cp $TRAVIS_BUILD_DIR/travis_scripts/testingalaxy.sh $GALAXY_HOME/
chmod 777 $GALAXY_HOME/testingalaxy.sh
su $GALAXY_TRAVIS_USER $GALAXY_HOME/script_docker.sh

curl --fail $BIOBLEND_GALAXY_URL/api/version
