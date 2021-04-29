# Getting Started

## Requirements on the target machine (where ansible will deploy GalaxyKickStart)

- Ubuntu 16.04, 18.04 or 20.04.
  
    Other systems (eg, debian) may work but they are not tested for the GalaxyKickStart
    development.

- Python >= 3.6

    Note that if this requirement is not satisfied, Ansible will try to  install Python 3

## Requirements on the Ansible machine

- **Ansible**

    Ansible uses ssh to send its commands. Thus, Ansible can be installed **remotely** (ie, on a
    machine that *will not* contain the Galaxy server at the end of the deployment), or **locally**
    (ie on the machine that *will* contain the Galaxy server, also called the target machine in
    this tutorial). In the latest case, ssh is used locally on the localhost 127.0.0.1 to chanel
    the commands sent by Ansible.
    
    Whether remote or local, the **Ansible version must be >= 2.9.6**

- git
  
    To clone the GalaxyKickStart GitHub repository


### Install Ansible with pip

A simple way to install the latest Ansible version is using [pip](https://pip.pypa.io/en/stable/quickstart/):

##### Install pip

```
sudo apt-get update -y && sudo apt-get -y install python-pip python-dev
```

- Ensure you have recent pip version installed (sudo -i && pip install upgrade pip maybe necessary)

pip upgrade:
```
sudo pip install -U pip
```
```
$ pip --version
pip 9.0.1 from /usr/local/lib/python2.7/site-packages (python 2.7)
```

- Then

```
sudo pip install ansible==2.7.4
```
### Install Ansible with apt

Alternatively, Ansible may be installed with the Apt package manager (Ubuntu):

```
sudo -i
apt-get install software-properties-common
apt-add-repository ppa:ansible/ansible
apt-get update
apt-get install ansible
```

# Getting the playbook

[//]: # (TODO: Once we do releases, we include the submodules and hence users can just download the playbook without git)

GalaxyKickStart is hosted on
[github](https://github.com/ARTbio/GalaxyKickStart.git) and uses a number of
dependent Ansible roles that need to be downloaded as part of the installation
step:

```
git clone https://github.com/ARTbio/GalaxyKickStart.git
cd GalaxyKickStart
ansible-galaxy install -r requirements_roles.yml -p roles
```

The playbook (here `galaxy.yml`) should be in the GalaxyKickStart folder.
```bash
ls
CONTRIBUTORS.md		Vagrantfile		docs			inventory_files		roles
Dockerfile		ansible.cfg		extra-files		mkdocs.yml		scripts
LICENSE.txt		deploy.sh		galaxy.yml		pre-commit.sh		startup.sh
README.md		dockerfiles		group_vars		requirements_roles.yml	templates
```

# Deploying galaxy-kickstart on remote machines.
----

Inside the `inventory_files` folder, you will find a number of inventory files.
This is an example of inventory taken from the `artimed` inventory file.

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

Inside the repository you will find a [file](https://github.com/ARTbio/GalaxyKickStart/tree/master/inventory_files/cloud)
called `inventory_files/cloud`. This file serves as an example hosts file for
how to deploy galaxy-kickstart on Google Compute Engine (GCE),  Amazon Web
Services(aws), and Jetstream (OpenStack). *Please note that the `ansible_user`
variable in the file changes for each remote target*. If you are wanting to use
this playbook on a cloud other than the ones listed  below, you will need to
update the inventory to add a new section header for the respective target. If
this happens to be a cloud setup, make sure to add the section header under
`[cloud_setup:children]`.

Specifications for each remote target:

* OpenStack
    * Image needed to deploy on [Jetstream](http://jetstream-cloud.org/):
        `Ubuntu 14.04.3 Development (jetstream image id: d7fe3289-943f-4417-90e8-8a6b171677ca)`
    *  Inventory: `<remote host IP> anisble_ssh_user="root" ansible_ssh_private_key_file="<path/to/your/private/key>"`

* GCE
    * Image needed to deploy galaxy-kickstart: `Ubuntu 14.04 LTS`
    * Inventory: ` <remote host IP> anisble_ssh_user="ubuntu" ansible_ssh_private_key_file="<path/to/your/private/key>"`

* AWS
    * Image needed to deploy galaxy-kickstart: `Ubuntu Server 14.04 LTS (HVM), SSD Volume Type - ami-2d39803a`
    * Inventory: `<target Amazon Web Services IP address> ansible_ssh_user="ubuntu" ansible_ssh_private_key_file="<path/to/your/aws/private/key>"`

## Deploying galaxy-kickstart behind a proxy

See [How can I set up GalaxyKickStart behind a proxy?](faq.md)
