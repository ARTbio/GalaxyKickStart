# Requirements
  * The target Operating System must be a Ubuntu Trusty 64 bits (it might work on other Debian systems, untested yet).
  * The target instance must have at least 4GB of RAM.
  * The target user must have sudo rights.
  * You need git and Ansible >= 1.8 (www.ansible.com) on the control machine on which you run the playbook.
  
# Ansible Galaxy instance, tools and reference genomes
To deploy Galaxy, first clone this repository by executing:
```
git clone --recursive https://github.com/ARTbio/ansible-artimed.git
```
Then, if you are targeting a remote machine, adapt the inventory (a sample can be found in the https://github.com/ARTbio/ansible-artimed/blob/master/hosts file) and execute:
```
cd ansible-artimed
#update/remove "-i hosts" if you want to use or own inventory file
ansible-playbook -i hosts galaxy.yml 
```
If the ssh key for the target is not the default ssh key, add its path to the "ansible_ssh_private_key_file" in the inventory.
This procedure will install Galaxy on the target.
Galaxy will be available on http port 80 (proxy NGINX) on the target IP with the admin user as artimed@gmail.com (password: artimed).
Note that the file https://github.com/ARTbio/ansible-artimed/blob/master/roles/artimed_extras/files/artimed_tool_list.yaml contains the default list of NGS tools to be installed.
Also note that file https://github.com/ARTbio/ansible-artimed/blob/master/roles/artimed_extras/files/data_managers.yaml contains the default list of reference genomes and data managers to run.

# Installing only tools
This procedure assumes Galaxy has already been installed and configured (for instance with the procedures described above).
To install only NGS tools on a Galaxy instance, set the value of the variable "install_galaxy" and "install_dms" to "False" in the inventory, and execute the last command of the previous procedure.
If you want to provide your own list of tools or change the default list, please see the file https://github.com/galaxyproject/ansible-galaxy-tools/blob/master/files/tool_list.yaml.sample for instructions and change the value of the variable "galaxy_tool_list" in the inventory.

# Load only reference genomes
This procedure assumes Galaxy has already been installed, configured and with data managers installed (for instance with the procedures described above).
To load reference genomes on a Galaxy instance, change the value of the variables "install_galaxy" and "install_tools" to "False" in the inventory, and execute the last command of the previous procedure.

# Alternative install - Vagrant
Before continue you must install Vagrant (www.vagrantup.com) and a vagrant compatible Virtual Box (www.virtualbox.org).
Bring up a vagrant machine by using:
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
The file ansible-artimed/hosts contains various important variables that may be adapated, so you can change it as necessary.
These variables are:
- ansible_ssh_user - The login name used to access the target.
- ansible_ssh_private_key_file - The ssh private key used to access the target.
- install_galaxy - True for install a Galaxy instance.
- install_tools - True for install the NGS tools.
- run_dms - True for run the data manager procedure.
- galaxy_user_name - The Operating System user name for galaxy process.
- galaxy_server_dir - The home of Operating System user for galaxy process.
- galaxy_admin - The admin galaxy user.
- galaxy_admin_pw - The admin galaxy password.
- galaxy_admin_api_key - The api key for tool installation and download reference genomes throught galaxy data managers.
- galaxy_tool_list - The files that constants the list of tools to be installed.
- galaxy_data_managers - The reference genomes and indexes to be load and build.
- galaxy_data - The persistent directory where the galaxy config and database directories will be installed or will be recovered.
- galaxy_database - The persistent directory where postgresql will be installed or will be recovered.
- galaxy_db - Connection string for galaxy-postgresql.
- galaxy_changeset_id - The release of Galaxy to be installed (master, dev or release_xx_xx).
