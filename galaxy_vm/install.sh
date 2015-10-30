#!/bin/bash
#Parsing parameters 
if [ "$1" == "-h" ]
then
	echo "Just one optional parameter: the git repository."
	exit 0
fi

artimed_git_repo=$1

if [ "$artimed_git_repo" == "" ]
then
	artimed_git_repo="https://github.com/fabiorjvieira/ansible-artimed.git"
fi

#install Dependencies
echo "Minimum requirements: OpenSSH client, Ansible >=1.8 and git."
sudo apt-add-repository ppa:ansible/ansible
sudo apt-get update
sudo apt-get install software-properties-common git openssh-client ansible	

#Installation
if git clone --recursive $artimed_git_repo
then
	cd ansible-artimed/galaxy_vm/
	sudo INSTALL_HOSTNAME=$HOSTNAME INSTALL_USER=$USER ansible-playbook -i "localhost," galaxy.yml
	echo "Wait until galaxy provide the web service on http://localhost:8080"
	echo -n "Press any key to install Galaxy tools..."
	read
	sudo INSTALL_HOSTNAME=$HOSTNAME INSTALL_USER=$USER ansible-playbook -i "localhost," tools.yml
fi
