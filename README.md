# Requirements
  * The target Operating System (OS) must be a Ubuntu Trusty 64 bits ( it might work on other Debian systems, untested).
  * The target instance must have at least 4GB of RAM.
  * The target user must have sudo rights.
  * You need Ansible >= 1.8 (www.ansible.com) on the machine on which you run the playbook.
  
# Ansible Galaxy instance
To deploy you will need ssh access to an account that can do passwordless sudo.
In the below command change the "targethost" for the IP of the target machine, "targetuser" for the remote user and execute:
```
git clone --recursive https://github.com/ARTbio/ansible-artimed.git
cd ansible-artimed
ansible-playbook -u targetuser -i "targethost," galaxy.yml -vvvv
```
If you may need to include the path to the ssh [--private-key path_to_private_key] if it is not the system's default key (optional).
```
ansible-playbook -u targetuser --private-key path_to_private_key -i "targethost," galaxy.yml -vvvv
```
Galaxy will be available on http port 80 (proxy NGINX) on the "targethost" IP.

# Installing Galaxy NGS tools
This procedure assumes Galaxy has already been installed and configured (for instance with the procedures described above).
If you want to install galaxy tools, change the "targethost" and "targetuser" for the IP and user of the target machine, respectively, and execute: 
```
cd ansible-artimed
ansible-playbook -u targetuser -i "targethost," -e "INSTALL_GALAXY=False INSTALL_TOOLS=True" galaxy.yml -vvvv
```

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

Alternatively, execute the file https://gist.github.com/fabiorjvieira/8672f445baf887eb5318 on the target machine and re-execute the Galaxy installation procedure.
It will configure the language environment variables and reinstall postgresql.

#Important variables
The two playbooks have some parameters those default values can be modified by using ansible-playbook parameter "-e" (see https://docs.ansible.com/ansible/playbooks_variables.html#passing-variables-on-the-command-line).

For galaxy.yml, the parameters are:
- INSTALL_HOSTNAME - The public network address (IP or full domain name) where the galaxy server will be installed.
- GALAXY_USER - The Operating System user name for galaxy process.
- GALAXY_ADMIN - The admin galaxy user.
- FTP_PORT - The ftp port for the proftpd service.
- GALAXY_KEY - The api key for tool installation.
- GALAXY_DATA - The persistent directory where the galaxy config and database directories will be installed or will be recovered. 
- GALAXY_DATABASE - The persistent directory where postgresql will be installed or will be recovered.
- GALAXY_DB_CONN - Connection string for galaxy-postgresql.
- GALAXY_TOOLS - The file that constants the list of tools to be installed (see file https://github.com/ARTbio/ansible-artimed/blob/master/roles/artimed_extras/files/artimed_tool_list.yaml).
- INSTALL_GALAXY - Installs galaxy if True.
- INSTALL_TOOLS - Installs galaxy tools if True.

#Format of the galaxy tools
The file https://github.com/ARTbio/ansible-artimed/blob/master/roles/artimed_extras/files/artimed_tool_list.yaml contains a default list of NGS tools.
If you want change it or do your on list, just follow the format of this file, that is:
- Do not change the first 3 lines of the file if you want to personalize it or copy the 3 lines of the file to your own file.
- Each tool is represented by a group of 4 lines: the name, owner, the panel id and name are the id and name of section where the tool will be listed, respectively.
- The name should should be preceded by the tag "- name" placed in the first column of the line.
- The owner should be preceded by the tag "owner" placed in the third column of the line.
- The panel id should be preceded by the tag "tool_panel_section_id" placed in the third column of the line.
- The panel name should be preceded by the tag "tool_panel_section_name" placed in the third column of the line.
