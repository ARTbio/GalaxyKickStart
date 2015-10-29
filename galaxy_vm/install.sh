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
echo "Minimum requirements: OpenSSH client, Ansible >=1.8, Vagrant >=1.7.4, Virtual Box (compatible with vagrant - see vagrant site) and git."
#sshVersion=`ssh -V | grep -i "openssh"`
#echo $sshVersion
#if [ "$sshVersion" == "" ]
#then
#	echo "ssh not present, please install it."
#	exit 1
#fi

gitVersion=`git --version | grep -i "^git"`
if [ "$gitVersion" == "" ]
then
	echo "git not present, please install it."
	exit 1
fi

vagrantVersion=`vagrant --version | grep -i "^vagrant"`
if [ "$vagrantVersion" == "" ]
then
	echo "vagrant not present, please install it and install the compatible virtualbox."
	exit 1
fi

ansibleVersion=`ansible --version | grep -i "^ansible" | sed "s/ansible //"`
ansibleVersion=`printf "$ansibleVersion\n1.8.0\n" | sort -n | head -n1`
if [ "$ansibleVersion" != "1.8.0" ]
then
	echo "ansible (>= 1.8) not present, please install it."
	exit 1
fi

#On Debians
#sudo apt-get install vim python-pip python-virtualenv python2.7 software-properties-common -y
#sudo apt-add-repository ppa:ansible/ansible -y
#sudo apt-get update
#sudo apt-get install ssh git virtualbox vagrant ansible -y

#Installation
if git clone --recursive $artimed_git_repo
then
	cd ansible-artimed/galaxy_vm/
	PLAYBOOK="galaxy.yml" VAGRANT_LOG=info vagrant up
	echo "Wait until galaxy provide the web service on http://localhost:8080"
	echo -n "Press any key to install Galaxy tools..."
	read
	PLAYBOOK="tools.yml" VAGRANT_LOG=info vagrant provision
fi
