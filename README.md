# Requirements
  * The target Operational System must be a Ubuntu Trusty 64 bits (can be other one of Debian flavours, but it was tested in Ubuntu Trusty 64 bits).
  * The target Operational System user must be in sudo group to do this (do not execute as root or with sudo).
  * The target Operational System must have at least 4GB of RAM.
  * The target Operational System must have Ansible >= 1.8 (www.ansible.com).
  * The target Operational System must have ssh client and server running on port 22, and git client.
  * The target Operational System must allow your user to ssh to it locally without user interactively input (https://stackoverflow.com/questions/4388385/how-to-ssh-login-without-password).
  

# Ansible Galaxy instance (web and ftp)
To deploy just execute:
```
git clone --recursive -b min https://github.com/ARTbio/ansible-artimed.git
cd ansible-artimed/galaxy/
hostIP=`hostname -I | cut -d " " -f 1`
INSTALL_HOSTNAME=$hostIP ansible-playbook -i "localhost," galaxy.yml -vvvv
```
Galaxy will be avaible in http port 8080 on the network ip where it was installed.

# Installing Galaxy NGS tools
If you want to install galaxy tools, execute the first two lines of the previous script and execute: 
```
GALAXY_USER="galaxy" GALAXY_PORT="8080" ansible-playbook -i "localhost," tools.yml -vvvv
```
Be sure that Galaxy is running and available in http port 8080 with the Operational System user galaxy, otherwise change the previous command accordingly. 

# Alterative install - Vagrant
Before continue you must install Vagrant (www.vagrantup.com) and a vagrant compatible Virtual Box (www.virtualbox.org).
Next, copy the file https://github.com/ARTbio/ansible-artimed/blob/master/galaxy_vm/Vagrantfile to one directory and execute:
```
vagrant up
vagrant ssh
```

Beware that vagrant redirect some ports from the guest machine to the host machine. 
Therefore, if this ports are already in use, you must change the ports specified in the Vagrantfile to other ports.
After "ssh" to the virtual machine, execute the same procedure described in the beginning of this text. 
Galaxy will be available in http port 8080 on the host network ip where the guest was installed if you did not changed it in the Vagrantfile.
