#!/bin/bash
#Parsing parameters 
if [ "$1" == "-h" ]
then
	echo "Just one optional parameter: the git repository."
	exit
fi

artimed_git_repo=$1
artimed_vm_relative_dir="ansible-artimed/"

if [ "$artimed_git_repo" == "" ]
then
	artimed_git_repo="https://github.com/fabiorjvieira/ansible-artimed.git"
fi

#Dependencies
sudo apt-add-repository ppa:ansible/ansible -y
sudo apt-get update
sudo apt-get install ansible vim git python2.7 software-properties-common virtualbox vagrant python-pip python-virtualenv -y

#VM directory
export vm_dir=$HOME/$artimed_vm_relative_dir
mkdir -p $vm_dir
rm -rf $vm_dir/*
cd $vm_dir

#Installation
if git clone $artimed_git_repo $vm_dir/
then
	cd galaxy_vm/
	PLAYBOOK="galaxy.yml" VAGRANT_LOG=info vagrant up
	echo "ssh to the VM and start galaxy. Wait until galaxy provide the web service on port 8080."
	echo -n "Press any key to continue..."
	read
	PLAYBOOK="tools.yml" VAGRANT_LOG=info vagrant provision
fi
