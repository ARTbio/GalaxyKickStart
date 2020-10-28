#!/usr/bin/env bash
set -e
ansible-galaxy install -r requirements_roles.yml -p roles
ansible-playbook -i inventory_files/galaxy-kickstart --skip-tags install_tools galaxy.yml
sleep 60
ansible-playbook -i inventory_files/galaxy-kickstart --tags install_tools galaxy.yml
# simple pings to galaxy server
sudo supervisorctl status
curl http://localhost:80/api/version| grep version_major
curl --fail $BIOBLEND_GALAXY_URL/api/version
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
