# Customising the playbook

We strongly encourage users to read the [ansible inventory](https://docs.ansible.com/ansible/intro_inventory.html) documentation first.  
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

[//]: # (TODO: Write-up extra-files, tools, workflows, which variables win.)