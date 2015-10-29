# Ansible ARTiMED virtual machine
Ansible playbooks for ARTiMED Virtual Machine
Deploys ARTiMED Vagrant box (the virtual machine) on Debian distributions. It includes Galaxy with postgresql database. To deploy just download the ansible-artimed/galaxy_vm/install.sh file and run:
```
#On Debian distros do:
git clone --recursive https://github.com/ARTbio/ansible-artimed.git
cd ansible-artimed/galaxy_vm
bash install.sh;
```

The install.sh file will:
 - install the requirements SO packages (git, pip, virtualenv, ansible, vagrant, virtualbox, ...) on Debian distributions.
 - create the instalation directory for the Vagrant box ($HOME/ansible-artimed/galaxy_vm/).
 - git clone this galaxy server playbook for vagrant provision procedure.
 - do a vagrant up to deploy the box.
 
If you want to skip these steps please clone this repository recussivelly (execute a "git clone --recursive"), install the necessary S.O. packages (ssh, git, ansible, vagrant and virtualbox principally) and execute a "PLAYBOOK='galaxy.yml' vagrant up", that is:
```
#Make sure that you have at least ssh, git, ansible, vagrant and virtualbox in your host system 
#Open one shell on host machine and execute:
git clone --recursive https://github.com/ARTbio/ansible-artimed.git
cd ansible-artimed/galaxy_vm
PLAYBOOK='galaxy.yml' vagrant up
echo "Run galaxy (see "Running Galaxy" next in this tutorial) and "
echo "wait until galaxy provide the web service on port 8080."
echo "Press any key to continue..."
read
PLAYBOOK='tools.yml' vagrant provision
```

# Running Galaxy
To run galaxy open another shell on host machine in the same directory (ansible-artimed/galaxy_vm) and execute:
```
vagrant ssh -c "sudo -i -u galaxy sh /home/galaxy/galaxy/run.sh"
```

# Re-doing
If you want to redo the installation process "cd" to the box directory (where the Vagrantfile is) and do the following steps:
```
vagrant package --output backup.box #backup your box;
vagrant destroy #it will destroy completelly your box, so do not miss the previous command;
vagrant up;
```

FYI, this repository has submodules, so to clone it with the roles "git clone --recursive" command must be used.
Galaxy is configured to run on the default port (8080), to monitor all network cards and the admin user is artimed@gmail.com. Galaxy database is labeled as galaxy with owner galaxy.
