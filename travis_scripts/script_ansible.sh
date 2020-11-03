#!/usr/bin/env bash
set -e
export GALAXY_USER="admin@galaxy.org"
export GALAXY_USER_EMAIL="admin@galaxy.org"
export GALAXY_USER_PASSWD="artbio2020"
export GALAXY_HOME=/home/galaxy
export GALAXY_TRAVIS_USER=galaxy
export GALAXY_UID=1450
export GALAXY_GID=1450
export BIOBLEND_GALAXY_API_KEY=artbio2020
export BIOBLEND_GALAXY_URL=http://127.0.0.1:80
export BIOBLEND_TEST_JOB_TIMEOUT=240

ansible-galaxy install -r requirements_roles.yml -p roles
ansible-playbook -i inventory_files/galaxy-kickstart galaxy.yml
sleep 30
ansible-playbook -i inventory_files/galaxy-kickstart galaxy_tool_install.yml

# simple pings to galaxy server
sudo supervisorctl status
curl http://localhost:80/api/version| grep version_major
curl --fail $BIOBLEND_GALAXY_URL/api/version
# test proftpd
date > $HOME/date.txt && curl --fail -T $HOME/date.txt ftp://127.0.0.1:21 --user $GALAXY_USER:$GALAXY_USER_PASSWD

# install bioblend testing, GKS way.
pip --version
sudo rm -f /etc/boto.cfg
pip install --ignore-installed https://github.com/galaxyproject/bioblend/archive/master.zip pytest

chmod a+rx /home/travis/
sudo -E su $GALAXY_TRAVIS_USER -c "source /home/travis/virtualenv/python3.7/bin/activate &&
cd $GALAXY_HOME &&
bioblend-galaxy-tests -v -k 'not download_dataset and \
              not download_history and \
              not export_and_download and \
              not test_show_nonexistent_dataset and \
              not test_invocation and \
              not test_update_dataset_tags and \
              not test_upload_file_contents_with_tags and \
              not test_create_local_user and \
              not test_show_workflow_versions' \
              /home/travis/virtualenv/python3.7/lib/python3.7/site-packages/bioblend/_tests/TestGalaxy*.py"
cd $TRAVIS_BUILD_DIR
