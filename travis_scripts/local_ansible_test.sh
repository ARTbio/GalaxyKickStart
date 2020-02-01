set -e

export GALAXY_TRAVIS_USER="galaxy"
export GALAXY_UID="1450"
export GALAXY_GID="1450"
export GALAXY_HOME="/home/galaxy"
export GALAXY_USER="admin@galaxy.org"
export GALAXY_USER_EMAIL="admin@galaxy.org"
export GALAXY_USER_PASSWD="admin"
export GALAXY_VERSION=release_18.09
export BIOBLEND_GALAXY_API_KEY="admin"
export BIOBLEND_TEST_JOB_TIMEOUT="200"
export BIOBLEND_GALAXY_URL="http://localhost:80"

# before_install:
# make sure pip install 'ansible<2.8'
 
ansible-galaxy install -r ~/galaxykickstart/requirements_roles.yml -p ~/galaxykickstart/roles

# tests

sh script_ansible.sh
