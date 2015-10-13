# ansible-artimed
Ansible playbooks for ARTiMED Virtual Machine
Deploys ARTiMED Vagrant box (the virtual machine). It includes Galaxy with postgresql database. To deploy just download the install.sh file and run:
```
bash install.sh
```

The install.sh file will:
 - install git, ansible, vagrant and virtualbox SO packages on Debian distributions.
 - install ansible roles for galaxy servers (galaxyprojectdotorg.galaxy and galaxyprojectdotorg.postgresql).
 - create the instalation directory for the Vagrant box.
 - git clone this galaxy server playbook for vagrant provision procedure.
 - do a vagrant up to deploy the box.
 
If you want to skip these steps please download the Vagrantfile and galaxy.yml provided in this repo and install the necessary packages and roles. Vagrantfile and galaxy.yml must be on the same directory to work.
To run galaxy you have to "ssh" to the box, cd to the galaxy directory and do "sh run.sh".

Note that, if your language environment variables are not correctly configured you must include one procedure on the postgresql role file /etc/ansible/roles/galaxyprojectdotorg.postgresql/tasks/debian.yml by executing:
```
sudo echo "" >> /etc/ansible/roles/galaxyprojectdotorg.postgresql/tasks/debian.yml
sudo echo "- shell: pg_createcluster {{ postgresql_version }} main --start" >> /etc/ansible/roles/galaxyprojectdotorg.postgresql/tasks/debian.yml 
sudo echo "  ignore_errors: yes" >> /etc/ansible/roles/galaxyprojectdotorg.postgresql/tasks/debian.yml
```
