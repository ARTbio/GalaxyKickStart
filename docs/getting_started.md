# Getting Started

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
This will take a while. Once finished, you should find a running Galaxy Instance on http://localhost:8080 .
If you would like to see the internals of the VM, you can log into the machine by typing `vagrant ssh`.

# Cleaning up

The VM image and various config files have been written to the `.vagrant` folder. Type `vagrant stop` to stop the running instance
and `vagrant destroy` to remove the VM, and then delete the `.vagrant` folder.

