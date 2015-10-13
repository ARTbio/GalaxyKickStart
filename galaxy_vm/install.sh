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
sudo sh -c "echo -e  'y/n' | ssh-keygen -q -f /root/.ssh/id_rsa -t rsa -N ''"
sudo sh -c "cat /root/.ssh/id_rsa.pub > /root/.ssh/authorized_keys"
sudo apt-add-repository ppa:ansible/ansible -y
sudo apt-get update
sudo apt-get install ansible vim git python2.7 software-properties-common virtualbox vagrant -y

#Can be moved to galaxy.yml as a pre task
sudo ansible-galaxy install galaxyprojectdotorg.galaxy --force
sudo ansible-galaxy install galaxyprojectdotorg.postgresql --force

#Can be moved to galaxy.yml as a pre task
#Fix the problem with language environment settings 
sudo sh -c 'echo "" >> /etc/ansible/roles/galaxyprojectdotorg.postgresql/tasks/debian.yml'
sudo sh -c 'echo "- shell: pg_createcluster {{ postgresql_version }} main --start" >> /etc/ansible/roles/galaxyprojectdotorg.postgresql/tasks/debian.yml'
sudo sh -c 'echo "  ignore_errors: yes" >> /etc/ansible/roles/galaxyprojectdotorg.postgresql/tasks/debian.yml'

#VM dir
export vm_dir=$HOME/$artimed_vm_relative_dir
mkdir -p $vm_dir
rm -rf $vm_dir/*
cd $vm_dir

#Installation
if git clone $artimed_git_repo $vm_dir/
then
	cd galaxy_vm/
	VAGRANT_LOG=info vagrant up
fi
