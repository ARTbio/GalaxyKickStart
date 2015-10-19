#!/bin/bash
#Parsing parameters 
if [ "$1" == "-h" ]
then
	echo "Just one optional parameter: the configuration file."
	exit
fi

if [ "$1" != "" ]
then
	while read parameter
	do
		tag=´echo $parameter | cut -f 1 -d "=" ´
		if [ "$tag" == "artimed_git_repo" ]; then artimed_git_repo=´echo $parameter | cut -f 2 -d "=" ´; 
		elif [ "$tag" == "artimed_vm_relative_dir" ]; then artimed_vm_relative_dir=´echo $parameter | cut -f 2 -d "=" ´;
		else echo "Invalid parameter in configuration file: $parameter"; fi
	done < $1
fi

if [ "$artimed_git_repo" == "" ]
then
	artimed_git_repo="https://github.com/fabiorjvieira/ansible-artimed.git"
fi
if [ "$artimed_vm_relative_dir" == "" ]
then
	artimed_vm_relative_dir="ansible-artimed/"
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
fi
