[![Build Status](https://travis-ci.org/ARTbio/GalaxyKickStart.svg?branch=master)](https://travis-ci.org/ARTbio/GalaxyKickStart)

# GalaxyKickStart

GalaxyKickStart is an Ansible playbook designed to help you get one or more
production-ready  [Galaxy servers](https://usegalaxy.org/) based on Ubuntu
within minutes, and to maintain these servers.
Optionally, instances can be pre-loaded with tools and workflows.

Detailed usage instructions are available in the
[Documentation](https://artbio.github.io/GalaxyKickStart/).

### Required Ansible version >= 2.7

The playbook has been tested on

- Cloud Machines
- Physical Servers
- Docker

GalaxyKickStart has been developed at the [ARTbio platform](http://artbio.fr)
and contains roles developed by the [Galaxy
team](https://github.com/galaxyproject/).

List of roles included in this playbook
------
- [ensure_postrgesql_up](https://github.com/ARTbio/ensure_postgresql_up.git)
- [natefoo-postgresql_objects](https://github.com/ARTbio/ansible-postgresql-objects)
- [galaxy-os role](https://github.com/ARTbio/ansible-galaxy-os)
- [galaxy role](https://github.com/ARTbio/ansible-galaxy)
- [miniconda-role](https://github.com/ARTbio/ansible-miniconda-role.git)
- [galaxy-extras role](https://github.com/ARTbio/ansible-galaxy-extras)
- [galaxy-trackster role](https://github.com/galaxyproject/ansible-trackster)
- [galaxy-tools role](https://github.com/ARTbio/ansible-galaxy-tools)


# Troubleshooting
### Installation of postgresql might fails due to non-standard locale ###
If you are using Ubuntu on your Ansible machine, make sure that you deactivate
`SendEnv LANG LC_*` in `/etc/ssh_config`. This will allow locale settings to
be propagated by ssh.
