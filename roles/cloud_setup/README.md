This role is primarily used to build Galaxy Standalone cloud images.

Due to the way the overall playbook and included roles are written, using
this role is a multi-step process. Each of the runs can be invoked with the
following command:

```ansible-playbook --inventory-file inventory_files/my_inventory galaxy.yml```

### Run #1: Build the core image and install the tools
If installing the Main toolset, this process takes about 3 hours per server.
 * Comment out tasks in `cleanup.py` that stop processes and remove ssh keys
 * Launch a new VM and update inventory with the target instance addresses
 * Update the toolset to be installed (check variable `galaxy_tools_tool_list_files` in `cloud_setup`)

### Run #2: Delete Galaxy bootstrap user
This step is required because `galaxy-tools` role flushes all the handlers so
if not split up, the cleanup step (_Run #3_) will run before any of the tools
are installed.
 * In `group_vars/cloud_setup`, set `galaxy_tools_delete_bootstrap_user: yes`
 * To speed up the build, in `galaxy.yml`, comment out all the roles except
   `galaxyprojectdotorg.galaxy-tools`

### Run #3: Cleanup the image for bundling
 * Revert changes from _Run #2_ and uncomment tasks from `cleanup.py`
 * In `galaxy.yml`, comment out all the roles except `cloud_setup`
 * Create an image via the Dashboard or API
