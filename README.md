[![Ansible Testing](https://github.com/ARTbio/GalaxyKickStart/actions/workflows/ci.yaml/badge.svg)](https://github.com/ARTbio/GalaxyKickStart/actions/workflows/ci.yaml.yaml)
[![Readthedocs Published](https://github.com/ARTbio/GalaxyKickStart/actions/workflows/readthedocs.yaml/badge.svg)](https://github.com/ARTbio/GalaxyKickStart/actions/workflows/readthedocs.yaml)
# GalaxyKickStart

GalaxyKickStart is an Ansible playbook designed to help you get one or more
production-ready [Galaxy servers](https://usegalaxy.org/) based on Ubuntu
within minutes, and to maintain these servers.
Optionally, GalaxyKickStart can install tools and workflows in the deployed
Galaxy server.

Detailed usage instructions are available in the
[Documentation](https://artbio.github.io/GalaxyKickStart/).

### Required Ansible version >= 2.9.2

### The playbook is tested on

- Cloud Machines
- Physical Servers
- Docker Images

GalaxyKickStart is developed by the [ARTbio platform](http://artbio.fr)
and uses roles developed by the [Galaxy team](https://github.com/galaxyproject/),
including:

------
- [ensure_postrgesql_up](https://github.com/ARTbio/ensure_postgresql_up.git)
- [natefoo-postgresql_objects](https://github.com/ARTbio/ansible-postgresql-objects)
- [galaxy-os role](https://github.com/ARTbio/ansible-galaxy-os)
- [galaxy role](https://github.com/ARTbio/ansible-galaxy)
- [miniconda-role](https://github.com/ARTbio/ansible-miniconda-role.git)
- [galaxy-extras role](https://github.com/ARTbio/ansible-galaxy-extras)
- [galaxy-trackster role](https://github.com/galaxyproject/ansible-trackster)
- [galaxy-tools role](https://github.com/ARTbio/ansible-galaxy-tools)


### Troubleshooting
Installation of postgresql might fails due to non-standard locale ###
If you are using Ubuntu on your Ansible machine, make sure that you deactivate
`SendEnv LANG LC_*` in `/etc/ssh_config`. This will allow locale settings to
be propagated by ssh.
