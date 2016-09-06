[![Build Status](https://travis-ci.org/ARTbio/GalaxyKickStart.svg?branch=master)](https://travis-ci.org/ARTbio/GalaxyKickStart)

# GalaxyKickStart

GalaxyKickStart is an Ansible playbook designed to help you get one or more
production-ready  [Galaxy servers](https://usegalaxy.org/) based on Ubuntu
within minutes, and to maintain these servers.
Optionally, instances can be pre-loaded with tools and workflows.

Detailed usage instructions are available in the
[Documentation](https://artbio.github.io/GalaxyKickStart/).

### Required Ansible version >= 2.1.2.0

The playbook has been tested on

- Cloud Machines
- Vagrant Boxes
- Physical Servers
- Docker

GalaxyKickStart has been developed at the [ARTbio platform](http://artbio.fr)
and contains roles developed by the [Galaxy
team](https://github.com/galaxyproject/).

List of roles included in this playbook
------
- [ensure_postrgesql_up](https://github.com/mvdbeek/ensure_postgresql_up)
- [galaxy-extras role](https://github.com/galaxyproject/ansible-galaxy-extras)
- [galaxy-tools role](https://github.com/galaxyproject/ansible-galaxy-tools)
- [galaxy-os role](https://github.com/galaxyproject/ansible-galaxy-os)
- [galaxy role](https://github.com/galaxyproject/ansible-galaxy)
- [galaxy-trackster role](https://github.com/galaxyproject/ansible-trackster)
- [natefoo-postgresql_objects](https://github.com/natefoo/ansible-postgresql-objects)

# Troubleshooting
### Installation of postgresql might fails due to non-standard locale ###
If you are using Ubuntu on your Ansible machine, make sure that you deactivate
`SendEnv LANG LC_*` in `/etc/ssh_config`. This will allow locale settings to
be propagated by ssh.
