# GalaxyKickStart

GalaxyKickStart is an Ansible [playbook](https://github.com/ARTbio/GalaxyKickStart)
designed to deploy one or more production-ready [Galaxy servers](https://usegalaxy.org/)
based on Ubuntu within ~30 minutes, and to maintain these servers.

GalaxyKickStart can also install tools and workflows in the deployed Galaxy servers.


# Requirements
The Galaxykickstart playbook is tested in

- Ubuntu **16.04, 18.04 and 20.04**,
- with **ansible >= 9.2.6**,
- and a target machine [pre-installed](https://phoenixnap.com/kb/how-to-install-python-3-ubuntu)
with **Python 3.7**.


# Target environments
The GalaxyKickStart playbook is primarily tested (ci) in virtual machines using GitHub
Actions workflows.

ARTbio is also using the playbook to install and maintain its bare metal Galaxy servers or
its virtual servers in [Google Cloud Platform](https://cloud.google.com/).

Finally, GalaxyKickStart can be used to build your production-ready Docker image and
a Galaxy Docker image is freely available using `docker pull artbio/galaxykickstart:18.04`.

The GalaxyKickStart Ansible playbook is maintained by [ARTbio platform](http://artbio.fr)
and makes use of roles which have been developed by the [Galaxy team](https://github.com/galaxyproject/).

To ensure the GalaxyKickStart stability, these roles (listed below) are forked and maintained
separately in ARTbio GitHub repositories (in the `galaxykickstart` branches).

List of dependency roles included in this playbook:
------

- [galaxy-os role](https://github.com/ARTbio/ansible-galaxy-os)
- [natefoo-postgresql_objects](https://github.com/ARTbio/ansible-postgresql-objects)
- [ensure_postgresql_up](https://github.com/ARTbio/ensure_postgresql_up.git)
- [galaxy role](https://github.com/ARTbio/ansible-galaxy)
- [miniconda-role](https://github.com/ARTbio/ansible-miniconda-role.git)
- [galaxy-extras role](https://github.com/ARTbio/ansible-galaxy-extras)
- [galaxy-trackster role](https://github.com/galaxyproject/ansible-trackster)
- [galaxy-tools role](https://github.com/ARTbio/ansible-galaxy-tools)
