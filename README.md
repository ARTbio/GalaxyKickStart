# ansible-artimed
Ansible playbooks for ARTiMED Virtual Machine
Deploys ARTiMED Vagrant box (the virtual machine). It includes Galaxy with postgresql database. To deploy just download the install.sh file and run:
```
bash install.sh
```

The install.sh file will:
 - install git, ansible (at least version 1.8), vagrant and virtualbox SO packages on Debian distributions.
 - install ansible roles for galaxy servers (galaxyprojectdotorg.galaxy and galaxyprojectdotorg.postgresql).
 - create the instalation directory for the Vagrant box ($HOME/ansible-artimed/galaxy_vm/).
 - git clone this galaxy server playbook for vagrant provision procedure.
 - do a vagrant up to deploy the box.
 
If you want to skip these steps please download the Vagrantfile and galaxy.yml provided in this repo and install the necessary packages and roles. Vagrantfile and galaxy.yml must be on the same directory to work.

To run galaxy you have to "ssh" to the box, cd to the galaxy directory and do "sh run.sh". Galaxy should be running on the default port (8080), listening all network cards and the admin user is artimed@gmail.com. Galaxy database is labeled as galaxy with owner galaxy.

If you want to redo the process "cd" to the box directory (where the Vagrantfile is) and do the following steps:
```
vagrant package --output backup.box #backup your box
vagrant destroy #it will destroy completelly your box, so do not miss the previous command
vagrant up
```

Note that, if your language environment variables are not correctly configured you must include one procedure on the postgresql role file /etc/ansible/roles/galaxyprojectdotorg.postgresql/tasks/debian.yml by executing:
```
sudo sh -c 'echo "" >> /etc/ansible/roles/galaxyprojectdotorg.postgresql/tasks/debian.yml'
sudo sh -c 'echo "- shell: pg_createcluster {{ postgresql_version }} main --start" >> /etc/ansible/roles/galaxyprojectdotorg.postgresql/tasks/debian.yml'
sudo sh -c 'echo "  ignore_errors: yes" >> /etc/ansible/roles/galaxyprojectdotorg.postgresql/tasks/debian.yml'
```
