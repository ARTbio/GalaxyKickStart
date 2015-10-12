#!/bin/bash
#wget https://raw.githubusercontent.com/fabiorjvieira/ansible-artimed/master/galaxy_vm/install.sh?token=AIqEJXyH_Z5PNAo_4Ff6Sr3UszgqqTAbks5WH5fewA%3D%3D
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

sudo sh -c "echo -e  'y/n' | ssh-keygen -q -f /root/.ssh/id_rsa -t rsa -N ''"
sudo sh -c "cat /root/.ssh/id_rsa.pub > /root/.ssh/authorized_keys"
sudo apt-add-repository ppa:ansible/ansible -y
sudo apt-get update
sudo apt-get install ansible vim git python2.7 software-properties-common -y

if git clone $git_repo
then
	echo "" >> $vm_dir/ansible_hosts
	echo "[ARTiMED_host_servers]" >> $vm_dir/ansible_hosts
	echo "localhost" >> $vm_dir/ansible_hosts

	ssh-keygen -f "$HOME/.ssh/known_hosts" -R [localhost]:2222
	sudo ansible-playbook -vvvv -i $vm_dir/ansible_hosts -e "VM_DIR=$vm_dir/$artimed/ LOCAL_USER=$USER" $vm_dir/$artimed/ARTiMED.yml
fi
