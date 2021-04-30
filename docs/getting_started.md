# Getting Started

### 1. Requirements on the target machine (where ansible will deploy GalaxyKickStart)

- **Ubuntu 16.04, 18.04 or 20.04**.
??? tips "Note"    
    Other systems (eg, debian) may work but they are not tested for the GalaxyKickStart
    development.

- **Python >= 3.6**

??? tips "Note"    
    If this requirement is not satisfied, Ansible will try to install Python 3 on the
    target machine

### 2. Requirements on the Ansible machine

- **Ansible**
Whether used remotely or locally, the **Ansible version must be >= 2.9.6**
??? tips "Note"    
    Ansible uses ssh to send its commands. Thus, Ansible can be installed **remotely** (ie, on a
    machine that *will not* contain the Galaxy server at the end of the deployment), or **locally**
    (ie on the machine that *will* contain the Galaxy server, also called the target machine in
    this tutorial). In the latest case, ssh is used locally on the localhost 127.0.0.1 to chanel
    the commands sent by Ansible.
    
    Ansible may be installed using **[pip](https://pip.pypa.io/en/stable/installing/)**
    ```
    pip install ansible==2.9.2
    ```
    or **[apt](https://help.ubuntu.com/community/AptGet/Howto)**
    ```
    sudo apt-get install software-properties-common
    sudo apt-add-repository ppa:ansible/ansible
    sudo apt-get update
    sudo apt-get install ansible
    ```

    
- **git**
    
??? tips "Note"    
    To clone the GalaxyKickStart GitHub repository


### 3. Getting the playbook

[//]: # (TODO: Once we do releases, we include the submodules and hence users can just download the playbook without git)

GalaxyKickStart is hosted on
[github](https://github.com/ARTbio/GalaxyKickStart.git) and uses a number of
dependent Ansible roles that need to be downloaded as part of the installation
step:

```
git clone https://github.com/artbio/galaxyKickstart.git
cd GalaxyKickStart
ansible-galaxy install -r requirements_roles.yml -p roles
```

The playbooks scripts `galaxy.yml` and `galaxy_tool_install.yml` are in the galaxykickstart folder.
```bash
CONTRIBUTORS.md			dockerfiles			group_vars			slurm_slave_node.yml
Dockerfile			docs				inventory_files			startup.sh
Dockerfile.galaxykickstart-base	extra-files			mkdocs.yml			templates
LICENSE.txt			galaxy.yml			requirements_roles.yml
README.md			galaxyToSlurmCluster.yml	roles
ansible.cfg			galaxy_tool_install.yml		scripts
```

### 3. Deploying galaxy-kickstart on remote machines.
----

Inside the `inventory_files` folder, you will find a number of inventory files.
This is an example of inventory taken from the `artimed` inventory file.

```
[artimed]
localhost ansible_connection=local

# change to the line below if remote target host
# <remote host IP> ansible_ssh_user="root" ansible_ssh_private_key_file="<path/to/your/private/key>"
```

Here `[artimed]` is a group, that contains a machine called localhost.

The variables defined in `group_vars/artimed` will be applied to this host,
in addition to, or overwriting, the variables defined in group_vars which apply to any host.

In this example, Ansible is acting locally on the localhost target.

If  instead you want to use Ansible remotely, replace `localhost ansible_connection=local`
by `<remote host IP> ansible_ssh_user="root" ansible_ssh_private_key_file="<path/to/your/private/key>"`

Ansible will connect by ssh to the target `remote host IP`, using the ssh key in `path/to/your/private/key`.

The user specified by `ansible_ssh_user=<user>` may be an other user than `root` but needs
in any case to to have sudo rights.

Then, run the plabook by typing:
```
ansible-playbook --inventory-file inventory_files/<your_inventory_file> galaxy.yml
```
Typically, you can test using:
```
ansible-playbook -i inventory_files/galaxy-kickstart galaxy.yml
```


You can put multiple machines in your inventory file.
If you run the playbook a second time, the process will be much faster, since steps that
have already been executed will be skipped. Whenever you change a variable (see
[customizations](customizations.md)) in group_vars/<your inventory> or in group_vars/all,
you will need to run the playbook again.

### 4. Deploying galaxy-kickstart on specified clouds

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

* GCE
    * Image needed to deploy galaxykickstart: `Ubuntu 18.04 LTS` > `Ubuntu 20.04 LTS` > `Ubuntu 16.04 LTS`
    * Inventory: ` <remote host IP> anisble_ssh_user="ubuntu" ansible_ssh_private_key_file="<path/to/your/private/key>"`

* AWS
    * Image needed to deploy galaxykickstart: `Ubuntu Server 18.04 LTS (HVM), SSD Volume Type - ami-013f17f36f8b1fefb (64 bits x86) / ami-02ed82f3a38303e6f (64 bits Arm)`
    * Inventory: `<target Amazon Web Services IP address> ansible_ssh_user="ubuntu" ansible_ssh_private_key_file="<path/to/your/aws/private/key>"`

* Jetstream (OpenStack)
    * Image needed to deploy galaxykickstart on [Jetstream](http://jetstream-cloud.org/):
        `Ubuntu 18.04 LTS Development + GUI support + Docker (jetstream image id: 15ff25f6-6ac5-4c12-b6ce-c08615ba32be)`
    *  Inventory: `<remote host IP> ansible_ssh_user="root" ansible_ssh_private_key_file="<path/to/your/private/key>"`

### 5. Deploying galaxy-kickstart behind a proxy

See [How can I set up GalaxyKickStart behind a proxy?](faq.md)
