# GalaxyKickStart

GalaxyKickStart is an Ansible playbook designed to help you get one or more production-ready
 [Galaxy servers](https://usegalaxy.org/) based on Ubuntu within minutes, and to maintain these servers.

# Required ansible version >= 2.1.2.0

Optionally, instances can be pre-loaded with tools and workflows.

The playbook has been tested on 

- Cloud Machines
- Vagrant Boxes
- Physical Servers 
- Docker.

GalaxyKickStart has been developed at the [ARTbio platform](http://artbio.fr) and contains roles developed
by the [Galaxy team](https://github.com/galaxyproject/).

List of roles included in this playbook
------

[ansible-postgresql-objects role](https://github.com/ARTbio/ansible-postgresql-objects/tree/ansible_2.2)
[ensure_postgresql_up role](https://github.com/mvdbeek/ensure_postgresql_up.git)  
[galaxy-extras role](https://github.com/galaxyproject/ansible-galaxy-extras)  
[galaxy-tools role](https://github.com/galaxyproject/ansible-galaxy-tools)  
[galaxy-os role](https://github.com/galaxyproject/ansible-galaxy-os)  
[galaxy role](https://github.com/galaxyproject/ansible-galaxy)  
[ansible-trackster role](https://github.com/galaxyproject/ansible-trackster)
[miniconda role](https://github.com/ARTbio/ansible-miniconda-role.git)
