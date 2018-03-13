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

- [ensure_postgresql_up](https://github.com/ARTbio/ensure_postgresql_up.git)
- [galaxy-extras role](https://github.com/galaxyproject/ansible-galaxy-extras)
- [galaxy-tools role](https://github.com/galaxyproject/ansible-galaxy-tools)
- [galaxy-os role](https://github.com/galaxyproject/ansible-galaxy-os)
- [galaxy role](https://github.com/galaxyproject/ansible-galaxy)
- [galaxy-trackster role](https://github.com/galaxyproject/ansible-trackster)
- [natefoo-postgresql_objects](https://github.com/natefoo/ansible-postgresql-objects)
- [miniconda-role](https://github.com/uchida/ansible-miniconda-role.git)
