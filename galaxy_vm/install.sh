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

#Test Dependencies
echo "Minimum requirements: OpenSSH client, Ansible >=1.8 and git."

gitVersion=`git --version | grep -i "^git"`
if [ "$gitVersion" == "" ]
then
	echo "git not present, please install it."
	exit 1
fi

ansibleVersion=`ansible --version | grep -i "^ansible" | sed "s/ansible //"`
ansibleVersion=`printf "$ansibleVersion\n1.8.0\n" | sort -n | head -n1`
if [ "$ansibleVersion" != "1.8.0" ]
then
	echo "ansible (>= 1.8) not present, please install it."
	exit 1
fi

#Installation
if git clone --recursive $artimed_git_repo
then
	cd ansible-artimed/galaxy_vm/
	sudo INSTALL_USER=$USER ansible-playbook -i "localhost," galaxy.yml
	echo "Wait until galaxy provide the web service on http://localhost:8080"
	echo -n "Press any key to install Galaxy tools..."
	read
	sudo INSTALL_USER=$USER ansible-playbook -i "localhost," tools.yml
fi
