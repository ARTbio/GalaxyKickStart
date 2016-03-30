# Customising the playbook

We strongly encourage users to read the [ansible inventory](https://docs.ansible.com/ansible/intro_inventory.html) documentation first.

Most settings should be editable without modifying the playbook directly,
instead variables can be set in group_vars and host vars.

The playbook comes with an example inventory file `hosts`.
```
[artimed]
localhost ansible_ssh_user="root" ansible_ssh_private_key_file="~/.ssh/id_rsa"
[travis_bioblend]
localhost ansible_connection=local
[aws]
# Put you aws IP and key here to make FTP work in the default VPC.
# If you want further group-specific variables, put the host in these groups as well [e.g artimed].
```
`[artimed]`, `[travis_bioblend]` and `[aws]` are predefined groups. Any host (here we only have localhost) that
is added to one or multiple groups will have the corresponding group variables applied.
Group variables are defined in `group_vars/[name of the group]` and default variables are found in   
`group_vars/all`.
All variables defined in `group_vars/all` are overwritten in `group_vars/[name of the group]`.  

For instance the variable `proftpd_nat_masquerade` is set to `false` in `group_vars/all`, while hosts in the `[aws]` group
apply the `[aws]` group variables which set `proftpd_nat_masquerade` to true, so that hosts in the aws group will have
this aws-specific setting applied. Any combination of groups may be used.

If you want to apply any of the changes you made to the variables you need to run the playbook again, making sure that
the host you are targeting is in the right group. The simplest way to do so is to use an inventory file that only contains
the group and the host you wish to target. If this is for example the group metavisitor, and you target the host localhost,
your inventory file should look like this:

```
[metavisitor]
localhost
```
You can then run the playbook as usual:
```
ansible-playbook --inventory-file=<your_inventory_file> galaxy.yml
```


[//]: # (TODO: Write-up extra-files, tools, workflows, which variables win.)

# Important variables

We aimed for this playbook to be reusable. We therefore made most variables configurable.
The group_vars/all file contains the variables we have chosen as defaults. You may override them either in this file
or you can use ansible group variables to selectively set the variables for certain hosts/groups. See the [ansible documentation
about group variables](http://docs.ansible.com/ansible/intro_inventory.html#splitting-out-host-and-group-specific-data) for details.

These most important variables are:

- ansible_ssh_user - The login name used to access the target.

- ansible_ssh_private_key_file - The ssh private key used to access the target.

- install_galaxy - True for install a Galaxy instance.

- install_tools - True for install the NGS tools.

- run_data_manager - True for run the data manager procedure.

- galaxy_user_name - The Operating System user name for galaxy process.

- galaxy_server_dir - The home of Operating System user for galaxy process.

- galaxy_admin - The admin galaxy user.

- galaxy_admin_pw - The admin galaxy password.

- default_admin_api_key - The api key for tool installation and download reference genomes throught galaxy data managers. To be removed in production.

- galaxy_tool_list - The files that constants the list of tools to be installed.

- galaxy_data_managers - The reference genomes and indexes to be load and build.

- galaxy_data - The persistent directory where the galaxy config and database directories will be installed or will be recovered.

- galaxy_database - The persistent directory where postgresql will be installed or will be recovered.

- galaxy_db - Connection string for galaxy-postgresql.

- galaxy_changeset_id - The release of Galaxy to be installed (master, dev or release_xx_xx).