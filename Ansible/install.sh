#!/bin/bash
cd $HOME/
git clone git://github.com/ansible/ansible.git --recursive
cd ./ansible
source ./hacking/env-setup
echo "[ARTiMED_host_servers]" > $HOME/ansible/ansible_hosts
echo "MACHINE WHERE GALAXY VM SHOULD BE INSTALLED" >> $HOME/ansible/ansible_hosts
export ANSIBLE_INVENTORY=$HOME/ansible/ansible_hosts
ansible-playbook ARTiMED.yml --become-user=root
