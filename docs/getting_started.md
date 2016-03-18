# Getting Started

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

# Deploying galaxy-kickstart on local virtual machine (VM) using vagrant.
----

GalaxyKickstarter is designed to be flexible and powerful, but for demonstration purposes we start a simple vagrant box
that runs this playbook. Following these instructions will not change the host system.
Alternatively, see [examples/docker](examples/docker.md) for running the playbook in docker.
More advanced examples are shown in the examples section.

# Requirements

To follow the examples [ansible](https://docs.ansible.com/ansible/intro_installation.html), [vagrant](https://www.vagrantup.com/downloads.html) 
and [git](https://git-scm.com/downloads) need to be installed.

# Getting the playbook

[//]: # (TODO: Once we do releases, we include the submodules and hence users can just download the playbook without git)


GalaxyKickstarter is hosted on [github](https://github.com/ARTbio/ansible-artimed.git) and makes use of submodules, so care
needs to be taken to also download the submodules. Cloning the repository for the first time can be done like this 
(note the `--recursive`):

```
git clone --recursive https://github.com/ARTbio/ansible-artimed.git
```

The playbook (here `galaxy.yml`) should be in the ansible-artimed folder.
```bash
ls ansible-artimed/
CONTRIBUTORS.md  docs  extra-files  galaxy.yml  group_vars  hosts  
LICENSE.txt  mkdocs.yml  pre-commit.sh  README.md  roles  Vagrantfile
```

# Running the playbook on a Virtual Machine

The Vagrantfile describes a Virtual Machine (VM) that is based on Ubuntu trusty.
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
This will take a while. Once finished, you should find a running Galaxy Instance on http://localhost:8080. You can log to this instance using admin@galaxy.org and admin as a password
If you would like to see the internals of the VM, you can log into the machine by typing `vagrant ssh`.

# Cleaning up

The VM image and various config files have been written to the `.vagrant` folder. Type `vagrant halt` to stop the running instance
and `vagrant destroy` to remove the VM, and then delete the `.vagrant` folder.

