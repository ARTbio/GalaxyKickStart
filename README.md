# Requirements
  * The target Operational System must be a Ubuntu Trusty 64 bits (can be other one of Debian flavours, but it was tested in Ubuntu Trusty 64 bits).
  * The target Operational System must have at least 4GB of RAM.
  * The target Operational System user must be in sudo group to do this (do not execute as root or with sudo).
  * The machine where the Ansible playbook will be execute (not the target) must have Ansible >= 1.8 (www.ansible.com) and a git client.
  
# Ansible Galaxy instance
To deploy just change the [targethost] for the IP of the target machine, [targetuser] for the remote user and execute:
```
git clone --recursive -b dev https://github.com/ARTbio/ansible-artimed.git
cd ansible-artimed
ansible-playbook -u targetuser -i "targethost," galaxy.yml -vvvv
```
Galaxy will be avaible in http port 80 (proxy NGINX) on the network ip where it was installed.

# Installing Galaxy NGS tools
If you want to install galaxy tools, change the [targethost] and [targetuser] for the IP and user of the target machine respectively, and execute: 
```
cd ansible-artimed/roles/artimed_extras/
GALAXY_USER="galaxy" GALAXY_PORT="80" ansible-playbook -u targetuser -i "targethost," tools.yml -vvvv
```
Be sure that Galaxy is running and available in http port 80 with the Operational System user galaxy, otherwise change the previous command accordingly. 

# Alterative install - Vagrant
Before continue you must install Vagrant (www.vagrantup.com) and a vagrant compatible Virtual Box (www.virtualbox.org).
Next, copy the file https://github.com/ARTbio/ansible-artimed/blob/dev/galaxy/Vagrantfile to one directory and execute:
```
vagrant up
```
Beware that vagrant redirect some ports from the guest machine to the host machine. 
Therefore, if this ports are already in use, you must change the ports specified in the Vagrantfile to other ports.
After "ssh" to the virtual machine, execute the same procedure described in the beginning of this text. 
Galaxy will be available in http port 8080 on the host network ip where the guest was installed if you did not changed it in the Vagrantfile. FTP server will be on 2121.
