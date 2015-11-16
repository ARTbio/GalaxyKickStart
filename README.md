# Requirements
  * your Operational System must be a Ubuntu Trusty 64 bits (can be other one of Debian flavours, but it was tested in Ubuntu Trusty 64 bits).
  * your Operational System user must be in sudo group to do this (do not execute as root or with sudo).
  * your Operational System must have at least 4GB of RAM.
  * your Operational System must have Ansible >= 1.8 (www.ansible.com).
  * your Operational System must have ssh client and server running on port 22, and git client.

# Ansible ARTiMED Galaxy instance
To deploy just execute:
```
git clone --recursive -b dev https://github.com/ARTbio/ansible-artimed.git
cd ansible-artimed/galaxy/
hostIP=`hostname -I | cut -d " " -f 1`
INSTALL_HOSTNAME=$hostIP ansible-playbook -i "localhost," galaxy.yml -vvvv
```

This script will ask you 3 things in 3 different times: your login/password on githut, if you agree to do a "ssh localhost" and if you want to install Galaxy tools.
Galaxy will be avaible in http port 80 of the network ip where it was installed.

# Installing Galaxy NGS tools
If you want to install only galaxy tools, execute the first two lines of the previous script and execute: 
```
ansible-playbook -i "localhost," tools.yml -vvvv
```

# Alterative install - Vagrant
Before continue you must install Vagrant (www.vagrantup.com) and a vagrant compatible Virtual Box (www.virtualbox.org).
After copy the file https://github.com/ARTbio/ansible-artimed/blob/master/galaxy_vm/Vagrantfile to one directory and execute:
```
vagrant up
vagrant ssh
```

Beware that vagrant redirect some ports from the guest machine to the host machine. Therefore if this ports are already in use, you must change the ports specified in the Vagrantfile to other ports.
After "ssh" to the virtual machine, execute the same procedure described in the beginning of this text. 
Galaxy will be avaible in http port 10080 of the host network ip where the guest was installed if you did not changed it in the Vagrantfile.
