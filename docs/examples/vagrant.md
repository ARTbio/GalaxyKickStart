# Deploying galaxy-kickstart on local virtual machine (VM) using vagrant.
----

GalaxyKickStart is designed to be flexible and powerful, but for demonstration purposes we start a simple vagrant box
that runs this playbook. Following these instructions will not change the host system.
Alternatively, see [examples/docker](docker.md) for running the playbook in docker,
or [getting started](../getting_started.md) for running the playbook on local or remote machines.

# Requirements

To follow the examples [ansible](https://docs.ansible.com/ansible/intro_installation.html), [vagrant](https://www.vagrantup.com/downloads.html)
and [git](https://git-scm.com/downloads) need to be installed.


# Running the playbook on a Virtual Machine

The Vagrantfile describes a Virtual Machine (VM) that is based on Ubuntu 14.04 (codename trusty).
```
VAGRANTFILE_API_VERSION = "2"
   Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|
      config.vm.box = "ubuntu/trusty64"
      config.vm.network "forwarded_port", guest: 80, host: 8080
      config.vm.network "forwarded_port", guest: 21, host: 2121

      config.vm.provider "virtualbox" do |v|
         v.memory = 4096
      end

      config.vm.provision "ansible" do |ansible|
         ansible.extra_vars = {
            ntp_server: "pool.ntp.org",
            ansible_ssh_user: 'vagrant'
         }
         ansible.verbose = 'vvvv'
         ansible.playbook = "galaxy.yml"
      end
   end
```
By default, port 8080 will be forwarded to port 80, and port 2121 will be forwarded to port 21 (for FTP),
and 4096 MB of memory will be attributed to the VM.
Enter the playbook directory `cd ansible-artimed` and type `vagrant up` to download a VM image and run the `galaxy.yml` playbook.

This will take a while. Once finished, you should find a running Galaxy Instance on http://localhost:8080 .
If you would like to see the internals of the VM, you can log into the machine by typing `vagrant ssh`.

`vagrant up` makes use of the ansible provisioner and is equivalent of starting a vagrant machine without the ansible provisioner
and running ansible through an ssh connection to the vagrant machine (which listens by default on port 2222)
The hosts inventory file contains an example for directly pointing ansible to the vagrant machine.
Uncomment the vagrant specific lines and comment or remove the remaining lines:

```
#[artimed]
#localhost ansible_ssh_user="root" ansible_ssh_private_key_file="~/.ssh/id_rsa"
#[travis_bioblend]
#localhost ansible_connection=local
# Uncomment the 2 lines below to point ansible to a local vagrant machine.
[all]
localhost ansible_user="vagrant" ansible_port=2222 ansible_private_key_file=.vagrant/machines/default/virtualbox/private_key
#[aws]
# Put you aws IP and key here to make FTP work in the default VPC.
# If you want further group-specific variables, put the host in these groups as well [e.g artimed].
```

To run the playbook again, type

```
ansible-playbook --inventory-file=<your_inventory> galaxy.yml
```

# Cleaning up

The VM image and various config files have been written to the `.vagrant` folder. Type `vagrant halt` to stop the running instance
and `vagrant destroy` to remove the VM, and then delete the `.vagrant` folder.