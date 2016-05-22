[![Build Status](https://travis-ci.org/ARTbio/galaxy-kickstart.svg?branch=master)](https://travis-ci.org/ARTbio/galaxy-kickstart)


# GalaxyKickstarter

GalaxyKickstarter is an Ansible playbook designed to help you get one or more production-ready
 [Galaxy servers](https://usegalaxy.org/) based on Ubuntu within minutes, and to maintain these servers.

Optionally, instances can be pre-loaded with tools and workflows.

The playbook has been tested on 

- Cloud Machines
- Vagrant Boxes
- Physical Servers 
- Docker.

Detailed instructions are available in the [Documentation](https://artbio.github.io/galaxy-kickstart/)


GalaxyKickstarter has been developed at the [ARTbio platform](http://artbio.fr) and contains roles developed
by the [Galaxy team](https://github.com/galaxyproject/).

List of roles included in this playbook
------

[galaxy-extras role](https://github.com/galaxyproject/ansible-galaxy-extras)  
[galaxy-tools role](https://github.com/galaxyproject/ansible-galaxy-tools)  
[galaxy-os role](https://github.com/galaxyproject/ansible-galaxy-os)  
[galaxy role](https://github.com/galaxyproject/ansible-galaxy)  
 
# Troubleshooting
The installation of postgresql might fail due to non-standard locale settings that can be propagated by ssh (found on ubuntu systems).
If you are using Ubuntu on your ansible machine, make sure that you deactivate `SendEnv LANG LC_*` in /etc/ssh_config.
