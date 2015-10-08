#!/bin/bash
galaxy_relative_dir=$1
git_repo=$2
artimed=$3

if [ "$git_repo" == "" ]
then
	git_repo="https://github.com/ARTbio/ansible-artimed.git"
fi

if [ "$artimed" == "" ]
then
	artimed="ansible-artimed/galaxy_vm"
fi

if [ "$galaxy_relative_dir" == "" ]
then
	galaxy_relative_dir="galaxy_install"
fi

export vm_dir=$HOME/$galaxy_relative_dir
mkdir -p $vm_dir
cd $vm_dir

sudo ssh-keygen -f /root/.ssh/id_rsa -t rsa -N ''
sudo cat /root/.ssh/id_rsa.pub > /root/.ssh/authorized_keys
sudo apt-get install vim git python2.7 software-properties-common -y
sudo apt-get install ansible -y

if git clone $git_repo
then
	echo "" >> $vm_dir/ansible_hosts
	echo "[ARTiMED_host_servers]" >> $vm_dir/ansible_hosts
	echo "localhost" >> $vm_dir/ansible_hosts

	ssh-keygen -f "$HOME/.ssh/known_hosts" -R [localhost]:2222
	sudo ansible-playbook -vvvv -i $vm_dir/ansible_hosts -e "VM_DIR=$vm_dir/$artimed/ LOCAL_USER=$USER" $vm_dir/$artimed/ARTiMED.yml
fi
