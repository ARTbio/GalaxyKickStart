# Getting Started

Make sure that you have a recent version of Ansible installed;
the playbook has been tested with version 2.1.2.0.

# Getting the playbook

[//]: # (TODO: Once we do releases, we include the submodules and hence users can just download the playbook without git)

GalaxyKickStart is hosted on
[github](https://github.com/ARTbio/GalaxyKickStart.git) and uses a number of
dependent Ansible roles that need to be downloaded as part of the installation
step:

```
git clone https://github.com/ARTbio/GalaxyKickStart.git
ansible-galaxy install -r requirements_roles.yml -p roles
```

The playbook (here `galaxy.yml`) should be in the GalaxyKickStart folder.
```bash
ls GalaxyKickStart/
CONTRIBUTORS.md  docs  extra-files  galaxy.yml  group_vars  hosts
LICENSE.txt  mkdocs.yml  pre-commit.sh  README.md  roles  Vagrantfile
```

# Deploying galaxy-kickstart on remote machines.
----

Inside the repository you will find a hosts file.
This is an example inventory.

```
[artimed]
localhost ansible_ssh_user="root" ansible_ssh_private_key_file="~/.ssh/id_rsa"
...
```

Here `[artimed]` is a group, that contains a machine called localhost.
The variables defined in `group_vars/artimed` will be applied to this host.
Ansible will connect by ssh to this machine, using the ssh key in `~/.ssh/id_rsa`.

If you would like to run this playbook on a remote machine by ssh (currently needs to be a debian-type machine),
create a new inventory, and change `localhost` to the IP address of that machine.
`ansible_ssh_user=<user>` controls under which username to connect to this machine.
This user needs to have sudo rights.

Then, run the plabook by typing:
```
ansible-playbook --inventory-file inventory_files/<your_inventory_file> galaxy.yml
```

You can put multiple machines in your inventory.
If you run the playbook a second time, the process will be much faster, since steps that have already been executed will be skipped.
Whenever you change a variable (see [customizations](customizations.md)), you need to run the playbook again.

## Deploying galaxy-kickstart on specified clouds

Inside the repository you will find a [file](https://github.com/ARTbio/GalaxyKickStart/tree/master/inventory_files/cloud) called
`inventory_files/cloud`. This file serves as an example hosts file
for how to deploy galaxy-kickstart on Google Compute Engine(gce), Amazon Web Services(aws), and Jetstream. *Please note
that the `ansible_ssh_user` variable in the file changes for each remote target*.



Specifications for each remote target:

* Jetstream
    * Image needed to deploy galaxy-kickstart:
        `Ubuntu 14.04.3 Development (jetstream image id: 3c3db94e-377b-4583-83d7-082d1024d74a)`
    *  Inventory: `<remote host IP> anisble_ssh_user="root" ansible_ssh_private_key_file="<path/to/your/private/key>"`

* GCE
    * Image needed to deploy galaxy-kickstart: `Ubuntu 14.04 LTS`
    * Inventory: ` <remote host IP> anisble_ssh_user="ubuntu" ansible_ssh_private_key_file="<path/to/your/private/key>"`

* AWS
    * Image needed to deploy galaxy-kickstart: `Ubuntu Server 14.04 LTS (HVM), SSD Volume Type - ami-2d39803a`
    * Inventory: `<target Amazon Web Services IP address> ansible_ssh_user="ubuntu" ansible_ssh_private_key_file="<path/to/your/aws/private/key>"`
