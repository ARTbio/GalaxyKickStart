#!/bin/bash
#fabio
galaxy_relative_dir=$1
git_repo=$2

cd
export vm_dir=$HOME/$galaxy_relative_dir
mkdir -p $vm_dir
cd $vm_dir
sudo ssh-keygen -f /root/.ssh/id_rsa -t rsa -N ''
sudo cp -v /root/.ssh/id_rsa.pub /root/.ssh/authorized_keys

sudo apt-get install vim git python2.7 software-properties-common -y
sudo apt-get install ansible -y

echo "" >> $vm_dir/ansible_hosts
echo "[ARTiMED_host_servers]" >> $vm_dir/ansible_hosts
echo "127.0.0.1" >> $vm_dir/ansible_hosts
if git clone $git_repo 
then
	ssh-keygen -f "$HOME/.ssh/known_hosts" -R [localhost]:2222
	mv -f $vm_dir/general/Ansible/* $vm_dir/; rm -rf $vm_dir/general/ #to be removed, however git must be corrected
	sudo ansible-playbook -vvvv -i $vm_dir/ansible_hosts -e "VM_DIR=$vm_dir LOCAL_USER=$USER" $vm_dir/ARTiMED.yml
fi
