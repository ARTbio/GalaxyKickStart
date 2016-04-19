# Installing Metavisitor with GalaxyKickstarter and Ansible
----
Here, a `Deployment Machine` will install a Metavisitor Galaxy server on `Target Machine`. Note that `Deployment Machine` and `Target Machine` can both be local or remote machines, and that they can be the same machine.

# Requirements
----

- On the `Deployment Machine`, [git](https://git-scm.com/downloads) and [ansible](https://docs.ansible.com/ansible/intro_installation.html) need to be installed.
- The `Target Machine` has to be accessible through ssh connection by the user (you) with `root` privileges. This implies that a correct ssh private key file is available on your `Deployment Machine`, for instance `~/.ssh/id_rsa`. This key will be used for secure transactions managed by ansible between the `Deployment Machine` and the `Target Machine`.

# Getting the ansible playbook
----
This is done on the `Deployment Machine` by cloning the [GalaxyKickstarter (ansible-artimed) repository](https://github.com/ARTbio/ansible-artimed.git) hosted by [the ARTbio organization](https://github.com/ARTbio):

In your terminal, type:
```
git clone --recursive https://github.com/ARTbio/ansible-artimed.git
```
[//]: # (TODO: Once we do releases, we include the submodules and hence users can just download the playbook without git)

Importantly, GalaxyKickstarter makes use of submodules, so care
needs to be taken to also download these submodules. This is why `--recursive` is included in the git command line.

At completion of the git cloning, you will have a new `ansible-artimed` folder, which contains everything need for deployment with ansible, including the playbook file (here `galaxy.yml`). You can verify this by typing in terminal:

`ls -la ansible-artimed`

# Adapting the ansible-artimed folder to your deployment
----

There are only few things to change in the `ansible-artimed` folder before running ansible.

## Adapt the ansible inventory file

In the ansible-artimed folder, there is a `hosts` file called the "inventory file".
For deploying Metavisitor, you need to edit this file so that it just contains

```

[metavisitor]

<ip address> ansible_ssh_user="root" ansible_ssh_private_key_file="<path/to/the/ssh/private/key>"

```

The `<ip address>` is the address of the `Target Machine`. The `<path/to/the/ssh/private/key>` is the path *on the `Deployment Machine`* to your ssh key, to be recognized by the `Target Machine`.

Thus, a practical exemple of the final content on the inventory file `hosts` is:

```

[metavisitor]

192.54.201.126 ansible_ssh_user="root" ansible_ssh_private_key_file="~/.ssh/id_rsa"

```

where `192.54.201.126` is the ip address of the `Target machine` and `~/.ssh/id_rsa` the path to the private ssh key.

## Adapt the group_vars/all file for persisting data, if needed.

In cases where your `Target machine` has volumes where you wish the Galaxy data to be persisted in, you have to edit the `ansible-artimed/group_vars/all` file, to indicate the path to this volume on the `Target machine`.

If you don't understand the previous statement, no worries, just don't do anything and skip this step.

For others, find the lines

```
#persistent data
galaxy_persistent_directory: /export # for IFB it's /root/mydisk, by default, /export
```

in the `ansible-artimed/group_vars/all` file, and change /export to the path of your persistent volume.

Note that if `/export` is not changed, nothing will happen and the deployed galaxy server and all associated data files will be in the `/home/galaxy/galaxy` folder of the `Target Machine`.

# Deploying Metavisitor Galaxy on the `Target Machine`
----

You are almost done.

Navigate with your terminal to your `ansible-artimed` folder and type the following command to run the ansible playbook for deploying metavisitor Galaxy on the `Target Machine`:

```
ansible-playbook --inventory-file=hosts galaxy.yml
```


If everything is ok, you may be asked to authorize the access to the `Target Machine` by typing yes in the terminal, and you will see ansible orchestrating the serveur deployment on the `Target Machine` in this terminal.

When the process is finished, you should be able to access the `Target Machine` by typing its IP address in your web browser.

By default the admin login/password is `admin@galaxy.org` / `admin`. You should change the password for safety.

# Re-deploying Metavisitor Galaxy on the `Target Machine`
----
If you are experimented in using ansible, you may customize your Metavisitor Galaxy instance deployed with GalaxyKickstarter by editing the content of `ansible-artimed`.

In that case, when your changes are done, just run again the command

```
ansible-playbook --inventory-file=hosts galaxy.yml
```
When you run the playbook a second time, the process will be much faster, since steps that have already been executed are skipped.
Whenever you change a variable (see [customizations](customizations.md)), you need to run the playbook again.


You can put multiple machines in your inventory: a simple way to do this is just copying the line the required number of times with the appropriate ip addresses:

```
[metavisitor]
192.54.201.126 ansible_ssh_user="root" ansible_ssh_private_key_file="~/.ssh/id_rsa"
192.54.201.127 ansible_ssh_user="root" ansible_ssh_private_key_file="~/.ssh/id_rsa"
192.54.201.128 ansible_ssh_user="root" ansible_ssh_private_key_file="~/.ssh/id_rsa"

```


