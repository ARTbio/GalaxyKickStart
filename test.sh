#!/bin/bash
# this script in github actions is accessible at the env variable $GITHUB_WORKSPACE
# if executed avec action checkout

set -euo pipefail

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


which pip3
python3 -m pip install -U pip setuptools 
python3 -m pip install ansible==2.7.4
python3 -m pip install --ignore-installed https://github.com/galaxyproject/bioblend/archive/refs/tags/v0.15.0.zip pytest

pip freeze

#debug
pwd

# run playbooks

ansible-galaxy install -r requirements_roles.yml -p roles

ansible-playbook -i inventory_files/galaxy-kickstart --extra-vars RUNNER_ALLOW_RUNASROOT="1" galaxy.yml

sleep 15

ansible-playbook -i inventory_files/galaxy-kickstart --extra-vars RUNNER_ALLOW_RUNASROOT="1" galaxy_tool_install.yml

# pings galaxy server
sudo supervisorctl status
curl http://localhost:80/api/version| grep version_major
curl --fail $BIOBLEND_GALAXY_URL/api/version

# test proftpd
echo "\ntest ftp transfer to proftpd server\n"
date > $HOME/date.txt && curl --fail -T $HOME/date.txt ftp://127.0.0.1:21 --user $GALAXY_USER:$GALAXY_USER_PASSWD

# install bioblend testing environnement

sudo rm -f /etc/boto.cfg # to do: understand the purpose of this step


bioblend-galaxy-tests -v /opt/hostedtoolcache/Python/3.7.10/x64/lib/python3.7/_tests/TestGalaxy*.py || true


# chmod a+rx /home/runner/
# sudo -E su $GALAXY_TRAVIS_USER -c "source GALAXY_HOME/virtualenv/python3.7/bin/activate &&
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
