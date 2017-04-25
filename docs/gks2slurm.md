# What is GKS2slurm ?
GKS2slurm is a playbook that is played to install a multinode slurm cluster over a GalaxyKickStart single-node installation.
The playbook GKS2slurm `galaxyToSlurmCluster.yml` was tested with multiple virtual machines (VMs) in `Stratuslab`, `Google Cloud Engine (GCE)` and `Amazon Web Services (AWS)` clouds.

# Installation of a Galaxy slurm cluster with GKS2slurm

## Step 1: Install a Galaxy server with GalaxyKickStart


- Report to the [Getting Started](getting_started.md) of this manual for the basics of GalaxyKickStart installation

- install any GalaxyKickStart "flavor" by configuring the inventory file (in inventory_files folder) and the group_vars file (in the group_vars folder) of your choice.
Flavors currently available are `kickstart`, `artimed` and `metavisitor` but other will come soon. Alternatively, you can build you own flavor by customizing a group_vars, extrafiles file and inventory file, which will install your Galaxy tools and workflows.

    **in Step 1, the most important thing to keep track with is to configure your target machine with an extra volume**

    Indeed GKS2slurm has be designed so that the Galaxy slurm cluster can accumulate large amount of data in the long term, which can be more easily shared with the cluster nodes and more importantly backed up.

    Thus in addition of all the adaptations you will do for your own purpose (tools, workflows, etc), edit the `group_vars/all` file and adapt the `galaxy_persistent_directory` variable to your extra volume which should be already formatted and mounted:
    
    Change
    
```
#persistent data
galaxy_persistent_directory: /export # for IFB it's /root/mydisk, by default, /export
```
    
To
    
```
#persistent data
galaxy_persistent_directory: /pathto/mounted/extravolume
```

- Having configured your GalaxyKickStart installation, import the extra roles (if not already done)
```
ansible-galaxy install -r requirements_roles.yml -p roles
```
and run the galaxy.yml playbook:
```
ansible-playbook --inventory-file inventory_files/<your_inventory_file> galaxy.yml
```

## Step 2: Check the single node Galaxy installation

If the playbook was run successfully, connect to your Galaxy instance through http and check that you can login (admin@galaxy.org:admin), and that tools and workflows are correctly installed.

## Step 3: Moving your single node configuration to a multinode slurm configuration

- Start as many compute nodes you want for the slurm cluster and gather information from each node:
    - IP address (all slurm nodes should must be accessible in the same network, ie nodes can be ping-ed from any nodes)
    - hostname
    - number of CPUs
    - memory (in MB)

### Step 3-1
Adapt the inventory file `slurm-kickstart` in the inventory_files folder.
```
[slurm_master]
# adapt the following line to IP address and ssh user of the slurm master node
192.54.201.102 ansible_ssh_user=root ansible_ssh_private_key_file="~/.ssh/mysshkey"

[slurm_slave]
# adapt the following lnes to IP addresses and ssu users of the slum slave nodes
192.54.201.98 ansible_ssh_user=root ansible_ssh_private_key_file="~/.ssh/mysshkey"
192.54.201.99 ansible_ssh_user=root ansible_ssh_private_key_file="~/.ssh/mysshkey"
192.54.201.101 ansible_ssh_user=root ansible_ssh_private_key_file="~/.ssh/mysshkey"
```

### Step 3-2
Adapt the group_vars file `slurm_master` in the `group_vars` folder.
This is done using the information gathered in step 3

```
# nfs sharing
cluster_ip_range: "0.0.0.0/24" # replace by your ip network range

# slave node specifications, adapt to your set of slave nodes
slave_node_dict:
  - {hostname: "slave-1", CPUs: "2", RealMemory: "7985"}
  - {hostname: "slave-2", CPUs: "2", RealMemory: "7985"}
  - {hostname: "slave-3", CPUs: "2", RealMemory: "7985"}
```
### Step 3-3
Adapt the group_vars file `slurm_slave` in the `group_vars` folder

```
# adapt the following variable to the master slurm node IP address
master_slurm_node_ip: "192.54.201.102"
```

### Step 3-4
Run the playbook `galaxyToSlurmCluster.yml` playbook.
from the GalaxyKickStart directory:
```
ansible-playbook -i inventory_files/slurm-kickstart galaxyToSlurmCluster.yml
```

Note that if you configure multiple slave nodes without prior ssh key authentification, you can run the same command with the variable ANSIBLE_HOST_KEY_CHECKING put to False:

```
ANSIBLE_HOST_KEY_CHECKING=False ansible-playbook -i inventory_files/slurm-kickstart galaxyToSlurmCluster.yml
```

# Checking slurm installation

Connect to your master node as root and type `sinfo`
Refer to slurm documentation for more investigation/control


