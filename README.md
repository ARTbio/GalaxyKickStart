# Ansible ARTiMED virtual machine
Ansible playbooks for ARTiMED Virtual Machine
Deploys ARTiMED Vagrant box (the virtual machine) on Debian distributions. It includes Galaxy with postgresql database. To deploy just download the ansible-artimed/galaxy_vm/install.sh file and run:
```
bash install.sh;
```

The install.sh file will:
 - install the requirements SO packages (git, pip, virtualenv, ansible, vagrant, virtualbox, ...) on Debian distributions.
 - create the instalation directory for the Vagrant box ($HOME/ansible-artimed/galaxy_vm/).
 - git clone this galaxy server playbook for vagrant provision procedure.
 - do a vagrant up to deploy the box.
 
If you want to skip these steps please download the Vagrantfile, galaxy.yml and the roles provided in this repo (ansible-artimed/galaxy_vm directory), install the necessary packages listed above and execute a "PLAYBOOK='galaxy.yml' vagrant up". Note that Vagrantfile, galaxy.yml and roles directory must be on the same directory to work.

# Running Galaxy
To run galaxy you have to "ssh" to the box, execute:
```
sudo su galaxy; 
cd $HOME/galaxy/;
sh run.sh;
```

# Installing NGS tools
To install Galaxy tools execute "PLAYBOOK='tools.yml' vagrant provision" on host machine. Note that Galaxy must be running on guest machine. 

If you want to redo the installation process "cd" to the box directory (where the Vagrantfile is) and do the following steps:
```
vagrant package --output backup.box #backup your box;
vagrant destroy #it will destroy completelly your box, so do not miss the previous command;
vagrant up;
```

FYI, this repository has submodules, so to clone it with the roles "git clone --recursive" command must be used.
Galaxy is configured to run on the default port (8080), to monitor all network cards and the admin user is artimed@gmail.com. Galaxy database is labeled as galaxy with owner galaxy.
