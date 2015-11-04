#!/bin/bash
#Test Dependencies
echo "Minimum requirements: OpenSSH client, Ansible >=1.8, Vagrant >=1.7.4, Virtual Box (compatible with vagrant - see vagrant site) and git."
vagrant up
vagrant ssh
#INSTALL_USER=vagrant PLAYBOOK="../galaxy/galaxy.yml" VAGRANT_LOG=info vagrant provision
#INSTALL_USER=vagrant PLAYBOOK="../galaxy/tools.yml" VAGRANT_LOG=info vagrant provision
