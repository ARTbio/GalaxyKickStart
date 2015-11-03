#!/bin/bash
#Test Dependencies
echo "Minimum requirements: OpenSSH client, Ansible >=1.8, Vagrant >=1.7.4, Virtual Box (compatible with vagrant - see vagrant site) and git."
#INSTALL_USER=vagrant PLAYBOOK="galaxy.yml" VAGRANT_LOG=info vagrant up
#INSTALL_USER=vagrant PLAYBOOK="tools.yml" VAGRANT_LOG=info vagrant provision

vagrant up
vagrant ssh
