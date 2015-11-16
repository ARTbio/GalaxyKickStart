# Requirements
  * your Operational System must a Debian like system (Debian and Ubuntu basically).
  * your Operational System user must be in sudo group to do this (do not execute as root or with sudo).
  * your Operational System must have at least 4GB of RAM.

# Ansible ARTiMED Galaxy instance
Deploys a Galaxy instance on the host machine in Debian flavors. 
It includes Galaxy with postgresql database and some extras (proftpd and nginx).
To deploy just download the https://github.com/ARTbio/ansible-artimed/blob/master/galaxy/install.sh file in your HOME directory and run:
```
cd $HOME
bash install.sh;
```
This script will ask you 3 things in 3 different times: your login/password on githut, if you agree to "ssh localhost" and if you want to install Galaxy tools.

# Using Galaxy instance
Galaxy is provided in http://localhost and the nginx proxy user name and password are artimed and artimed.
The administrator email of Galaxy is artimed@gmail.com (you have to register it in Galaxy).

If you restart the machine where Galaxy was installed, please start Galaxy service with:
```
sudo service galaxy start
```

# Installing Galaxy NGS tools
If you want to install only galaxy tools, donwload https://github.com/ARTbio/ansible-artimed/blob/master/galaxy/tools.yml and https://github.com/ARTbio/ansible-artimed/blob/master/galaxy/artimed_tool_list.yaml and execute:
```
cd $HOME/ansible-artimed/galaxy/
GALAXY_USER=$galaxy_user GALAXY_PORT=$galaxy_port INSTALL_USER=$USER ansible-playbook -i "localhost," tools.yml -vvvv
```

