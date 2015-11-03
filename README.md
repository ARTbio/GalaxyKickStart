# Requirements
  * your Operational System must a Debian like system (Debian and Ubuntu basically).
  * your Operational System user must be sudo to do this.
  * your Operational System must have apt packages candidates for OpenSSH client and git.
  * your Operational System must have at least 4GB of RAM.

# Ansible ARTiMED Galaxy instance
Deploys a Galaxy instance on the host machine in Debian flavors. 
It includes Galaxy with postgresql database and some extras (proftpd and nginx).
To deploy just download the https://github.com/ARTbio/ansible-artimed/blob/master/galaxy/install.sh file and run:
```
#Download https://github.com/ARTbio/ansible-artimed/blob/master/galaxy/install.sh and execute:
bash install.sh;
```
