#!/usr/bin/env bash
# this script in github actions is accessible at the env variable $GITHUB_WORKSPACE
# if executed avec action checkout

set -o pipefail
set -e

export GALAXY_TRAVIS_USER="galaxy"
export GALAXY_UID="1450"
export GALAXY_GID="1450"
export GALAXY_HOME="/home/galaxy"
export GALAXY_USER="admin@galaxy.org"
export GALAXY_USER_EMAIL="admin@galaxy.org"
export GALAXY_USER_PASSWD="artbio2020"
export GALAXY_VERSION=release_20.01
export BIOBLEND_GALAXY_API_KEY="artbio2020"
export BIOBLEND_TEST_JOB_TIMEOUT="240"
export BIOBLEND_GALAXY_URL="http://127.0.0.1:80"

/etc/init.d/postgresql stop || true
apt-get -y --purge remove postgresql libpq-dev libpq5 postgresql-client-common postgresql-common || true
rm -rf /var/lib/postgresql || true


which pip3
python3 -m pip install -U pip setuptools 
python3 -m pip install ansible==2.7.4

# git clone https://github.com/artbio/galaxykickstart -b actions 
# cd galaxykickstart

ansible-galaxy install -r requirements_roles.yml -p roles

# tests

ansible-playbook -i inventory_files/galaxy-kickstart galaxy.yml
sleep 15
# ansible-playbook -i inventory_files/galaxy-kickstart galaxy_tool_install.yml

# simple pings to galaxy server
sudo supervisorctl status
curl http://localhost:80/api/version| grep version_major
curl --fail $BIOBLEND_GALAXY_URL/api/version
# test proftpd
date > $HOME/date.txt && curl --fail -T $HOME/date.txt ftp://127.0.0.1:21 --user $GALAXY_USER:$GALAXY_USER_PASSWD

# install bioblend testing, GKS way.
# pip --version
# sudo rm -f /etc/boto.cfg
# pip install --ignore-installed https://github.com/galaxyproject/bioblend/archive/master.zip pytest
# 
# chmod a+rx /home/runner/
# sudo -E su $GALAXY_TRAVIS_USER -c "source /home/travis/virtualenv/python3.7/bin/activate &&
# cd $GALAXY_HOME &&
# bioblend-galaxy-tests -v -k 'not download_dataset and \
#               not download_history and \
#               not export_and_download and \
#               not test_show_nonexistent_dataset and \
#               not test_invocation and \
#               not test_update_dataset_tags and \
#               not test_upload_file_contents_with_tags and \
#               not test_create_local_user and \
#               not test_show_workflow_versions' \
#               /home/travis/virtualenv/python3.7/lib/python3.7/site-packages/bioblend/_tests/TestGalaxy*.py"
# cd $TRAVIS_BUILD_DIR
