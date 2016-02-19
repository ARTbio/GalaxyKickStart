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
ansible-playbook --inventory-file=<your_inventory> galaxy.yml
```

You can put multiple machines in your inventory.  
If you run the playbook a second time, the process will be much faster, since steps that have already been executed will be skipped.