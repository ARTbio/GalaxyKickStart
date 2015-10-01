#!/bin/bash
cd $HOME/
apt-get install git
apt-get install python2.7
apt-get install software-properties-common
apt-add-repository ppa:ansible/ansible
apt-get update
apt-get install ansible
source ./hacking/env-setup
echo "[ARTiMED_host_servers]" > $HOME/ansible/ansible_hosts
echo "127.0.0.1" >> $HOME/ansible/ansible_hosts
export ANSIBLE_INVENTORY=$HOME/ansible/ansible_hosts
ansible-playbook ARTiMED.yml --become-user=root
