# GalaxyKickStart

GalaxyKickStart is an [Ansible](http://www.ansible.com/) playbook designed for installing,
testing, deploying and  maintaining production-grade [Galaxy
server](https://galaxyproject.org/admin/get-galaxy/) instances. GalaxyKickStart
playbook code is available in [GitHub](https://github.com/ARTbio/GalaxyKickStart). 

In the basic configuration, the deployed Galaxy servers include:

- a `postgresql` server as database backend 
- a `uwsgi` Web Server Gateway Interface
- a `Slurm` job manager
- a `nginx` web proxy
- a `proftpd` FTP server
- the service manager `supervisor`
