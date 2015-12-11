# Requirements
  * The target Operating System must be a Ubuntu Trusty 64 bits ( it might work on other Debian systems, untested).
  * The target instance must have at least 4GB of RAM.
  * You need ssh access to an account that can do passwordless sudo shell commands.
  * You need git and Ansible >= 1.8 (www.ansible.com) on the machine on which you run the playbook.
  
# Ansible Galaxy instance and NGS tools
To deploy Galaxy first, clone this repository by executing:
```
git clone --recursive https://github.com/ARTbio/ansible-artimed.git
```
Second, change "targethost" (located in the second line of the file ansible-artimed/hosts) by the IP of the target host.
Finally, customize any other variables, if necessary, and execute:
```
cd ansible-artimed
ansible-playbook -i host galaxy.yml
```
If you may need to include the path to your ssh private key change the value of the variable "ansible_ssh_private_key_file" in the ansible-artimed/hosts.
This procedure will install Galaxy on the target host (refered by the IP on the second line of the file hosts).
Galaxy will be available on http port 80 (proxy NGINX) on the target host IP.

# Installing only Galaxy NGS tools
This procedure assumes Galaxy has already been installed and configured (for instance with the procedures described above).
To install only NGS tools on a Galaxy instance, change the value of the variable "install_galaxy" to "False" (located in the file ansible-artimed/hosts) and execute:
```
ansible-playbook -i host galaxy.yml
```
Note that the file https://github.com/ARTbio/ansible-artimed/blob/master/roles/artimed_extras/files/artimed_tool_list.yaml contains the default list of NGS tools to be installed.
Therefore, if you want change it, please see the file https://github.com/galaxyproject/ansible-galaxy-tools/blob/master/files/tool_list.yaml.sample for instructions.
If you want to provide your own list of tools, change the value of the variable "galaxy_tools" in ansible-artimed/hosts.

# Alternative install - Vagrant
Before continue you must install Vagrant (www.vagrantup.com) and a vagrant compatible Virtual Box (www.virtualbox.org).
Execute the first script of this readme and execute:
```
vagrant up
```
Beware that vagrant redirect some ports from the guest machine to the host machine. 
Therefore, if this ports are already in use, you must change the ports specified in the Vagrantfile to other ports.
After "ssh" to the virtual machine, execute the same procedure described in the beginning of this text. 
Galaxy will be available in http port 8080 on the host network IP where the guest was installed if you did not changed it in the Vagrantfile. FTP server will be on 2121.

# Troubleshooting
The installation of postgresql might fail due to non-standard locale settings that can be propagated by ssh (found on ubuntu systems).
If you are using Ubuntu on your ansible machine, make sure that you deactivate `SendEnv LANG LC_*` in /etc/ssh_config.

#Important variables
The files ansible-artimed/hosts contains various important variables and the respectives default values, so you can change it as necessary.
Here we list the function of each one of them by the order of apperence:
- ansible_ssh_user - The login name used to access the target host.
- ansible_ssh_private_key_file - The ssh private key used to access the target host.
- install_galaxy - True for install a Galaxy instace.
- install_tools - True for install the NGS tools.
- run_dm - True for run the data maanager procedure.
- galaxy_user_name - The Operating System user name for galaxy process.
- galaxy_server_dir - The home of Operating System user for galaxy process.
- galaxy_admin - The admin galaxy user.
- galaxy_api - The api key for tool installation.
- galaxy_tools - The files that constants the list of tools to be installed.
- galaxy_data_managers - The reference genomes and indexes to be load and build.
- galaxy_data - The persistent directory where the galaxy config and database directories will be installed or will be recovered.
- galaxy_database - The persistent directory where postgresql will be installed or will be recovered.
- galaxy_db - Connection string for galaxy-postgresql.
