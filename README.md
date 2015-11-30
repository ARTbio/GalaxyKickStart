# Requirements
  * The target Operating System (OS) must be a Ubuntu Trusty 64 bits ( it might work on other Debian systems, untested).
  * The target instance must have at least 4GB of RAM.
  * The target user must have sudo rights.
  * You need Ansible >= 1.8 (www.ansible.com) on the machine on which you run the playbook.
  
# Ansible Galaxy instance
To deploy you will need ssh access to an account that can do passwordless sudo.
You may need to include the path to the ssh [--private-key path_to_private_key] if it is not the system's default key (optional).
In the below command change the "targethost" for the IP of the target machine, "targetuser" for the remote user and execute:
```
git clone --recursive https://github.com/ARTbio/ansible-artimed.git
cd ansible-artimed
ansible-playbook -u targetuser --private-key path_to_private_key -i "targethost," galaxy.yml -vvvv
```

Galaxy will be avaible on http port 80 (proxy NGINX) on the "targethost" ip.

# Installing Galaxy NGS tools
This procedure assumes Galaxy has already been installed and configured (for instance with the procedures described above).
If you want to install galaxy tools, change the "targethost" and "targetuser" for the IP and user of the target machine respectively, and execute: 
```
cd ansible-artimed/roles/artimed_extras/
GALAXY_USER="galaxy" GALAXY_PORT="80" ansible-playbook -u targetuser -i "targethost," tools.yml -vvvv
```
# Alterative install - Vagrant
Before continue you must install Vagrant (www.vagrantup.com) and a vagrant compatible Virtual Box (www.virtualbox.org).
Execute the first script of this readme and execute:
```
vagrant up
```
Beware that vagrant redirect some ports from the guest machine to the host machine. 
Therefore, if this ports are already in use, you must change the ports specified in the Vagrantfile to other ports.
After "ssh" to the virtual machine, execute the same procedure described in the beginning of this text. 
Galaxy will be available in http port 8080 on the host network ip where the guest was installed if you did not changed it in the Vagrantfile. FTP server will be on 2121.

# Troubleshooting
The installation of postgresql might fail due to non-standard locale settings that can be propagated by ssh (found on ubuntu systems).
If you are using Ubuntu on your ansible machine, make sure that you deactivate `SendEnv LANG LC_*` in /etc/ssh_config.

Alternatively, execute the file https://gist.github.com/fabiorjvieira/8672f445baf887eb5318 on the target machine and re-execute the Galaxy installation procedure.
It will configure the language environment variables and reinstall postgresql.
