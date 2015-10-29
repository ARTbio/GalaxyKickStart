# Ansible ARTiMED virtual machine
Ansible playbooks for ARTiMED Virtual Machine
Deploys ARTiMED Vagrant box (the virtual machine). It includes Galaxy with postgresql database. To deploy just download the ansible-artimed/galaxy_vm/install.sh file and run:
```
<<<<<<< HEAD
#Download ansible-artimed/galaxy_vm/install.sh and execute:
=======
#On Debian distros, download ansible-artimed/galaxy_vm/install.sh and execute:
>>>>>>> 135668bac7d2c5ebb5c8365d48c23414ce53424d
bash install.sh;
cd ansible-artimed/galaxy_vm
```

The install.sh file will:
 - install verify the requirements (git, pip, virtualenv, ansible, vagrant, virtualbox, ...).
 - create the instalation directory for the Vagrant box ($HOME/ansible-artimed/galaxy_vm/).
 - git clone this galaxy server playbook for vagrant provision procedure.
 - do a vagrant up to deploy the box.
 - start galaxy
 - install the tools listed in into Galaxy
 
<<<<<<< HEAD
# Running from the host machine Galaxy
Galaxy is installed as a SO service, however if you want to run galaxy from the host machine as a shell application, open another shell on host machine on the same install directory (ansible-artimed/galaxy_vm) and execute:
```
vagrant ssh -c "sudo service galaxy stop"
vagrant ssh -c "sudo -i -u galaxy sh /home/galaxy/galaxy/run.sh"
```

=======
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

>>>>>>> 135668bac7d2c5ebb5c8365d48c23414ce53424d
# Re-doing
If you want to redo the installation process "cd" to the box directory (where the Vagrantfile is) and do the following steps:
```
vagrant package --output backup.box #backup your box;
vagrant destroy #it will destroy completelly your box, so do not miss the previous command;
vagrant up;
```

FYI, this repository has submodules, so to clone it with the roles "git clone --recursive" command must be always used.
Galaxy is configured to run on the default port (8080), to monitor all network cards and the admin user is artimed@gmail.com. 
Galaxy database is labeled as galaxy with owner galaxy.
Minimum requirements: Ansible >=1.8, Vagrant >=1.7.4, Virtual Box (compatible with vagrant - see vagrant site) and git. 
