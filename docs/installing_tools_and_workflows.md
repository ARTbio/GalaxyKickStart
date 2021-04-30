# Installing tools and workflows
----

This playbook includes the [ansible-galaxy-tools](https://github.com/galaxyproject/ansible-galaxy-tools)
role which can be used to install tools and workflows into galaxy instances using the
[bioblend](https://bioblend.readthedocs.org/en/latest/) API.

Importantly, in the latest GalaxyKickStart version (`v20.05`), tool installation is performed
using a separate run of `ansible-playbook`, typically:

```
# these steps have should have already been performed

# ansible-galaxy install -r requirements_roles.yml -p roles/
# ansible-playbook -i inventory_files/galaxy-kickstart galaxy.yml

ansible-playbook -i inventory_files/galaxy-kickstart galaxy_tool_install.yml
```
!!! note ""
    Note that this is during the galaxy_tool_install.yml ansible play that
    a Galaxy admin user account is created with the credentials admin@galaxy.org:artbio2020

### Creating a tool_list.yml file
Before running the `galaxy_tool_install.yml` playbook script as shown above, you need to
prepare a `tool_list.yml` file with a list of tools in yaml format and with the following
example content:

```
tools:
- name: blast_to_scaffold
  owner: drosofff
  revisions:
  tool_panel_section_label: Metavisitor
  tool_shed_url: https://toolshed.g2.bx.psu.edu/
- name: blastx_to_scaffold
  owner: drosofff
  revisions:
  tool_panel_section_label: Metavisitor
  tool_shed_url: https://toolshed.g2.bx.psu.edu/
- name: bowtie2
  owner: devteam
  revisions:
  - 019c2a81547a
  tool_panel_section_label: Metavisitor
  tool_shed_url: https://toolshed.g2.bx.psu.edu/
```

If the `revisions` key is empty, the latest available in the Galaxy toolshed revision will
be installed.  

The `tool_panel_section_label` key sets the tool panel section where the tools will show up
in the Galaxy tool bar.

Another example of a a tool list can be found in [here](https://github.com/ARTbio/GalaxyKickStart/blob/master/extra-files/metavisitor/metavisitor_tool_list.yml)

### Obtaining a tool_list.yml file from a running Galaxy server

You can also retrieve the list of tools running in a specific Galaxy instance using the
[ephemeris](https://github.com/galaxyproject/ephemeris) script `workflow-to-tools`.

1. First Install ephememeris using pip
    ```
    pip install ephemeris
    ```
2. This will bring several scripts in your pip environment:
    ```
    run-data-managers --help
    shed-tools install --help
    shed-tools update --help
    workflow-install --help
    setup-data-libraries --help
    get-tool-list --help
    workflow-to-tools --help
    ```
3. use `get-tool-list` to retrieve a yml list of tools using the command:
    ```
    get-tool-list -g https://usegalaxy.org -u <main galaxy username> -p <user password> --get_data_managers -o main_tools_list.yml
    ``` 
### Obtaining a tool_list.yml file from a workflow.ga galaxy file.

You can also retrieve a list of tools used in one or more workflow galaxy files (.ga extension).

These .ga files can be obtain in Galaxy server instances (menu "download workflow file") and
in some other repositories

1. use `workflow-to-tools` to retrieve a yml list of tools using the command:
    ```
    workflow-to-tools -w <Galaxy-Workflow-File.ga> -l <menu_label> -o <tool_list.yml>
    ```

# Adding a tool_list.yml file to a group_variable files

Group variable files are in the group_vars directory.

If you would like to install tools, you need to reference the tool_list.yml in the group variable file.
We typically place additional files in the `extra-files/<hostname>/<hostname>_tool_list.yml` file.

If you would like to add tools to a group that is called metavisitor edit `group_vars/metavisitor` and add these lines:
```
install_tools: true
galaxy_tools_tool_list: "extra-files/metavisitor/metavisitor_tool_list.yml"
```

# Installing workflows

You can also make sure that workflows are available after running the playbook.
As with tools, place the workflows in `extra-files/<hostname>/<hostname><workflow_name>.ga`
Add these lines to the corresponding group_var file:
```
galaxy_tools_install_workflows: true
galaxy_tools_workflows:
  - "extra-files/metavisitor/Galaxy-Workflow-create_model.ga"
  - "extra-files/metavisitor/Galaxy-Workflow-separate_host_and_virus_reads.ga"
  - "extra-files/metavisitor/Galaxy-Workflow-standart_metavisitor_workflow_(input__clipped_dataset).ga"
  - "extra-files/metavisitor/Galaxy-Workflow-Metavisitor_Test_case_1-1_Guided.ga"
  - "extra-files/metavisitor/Galaxy-Workflow-Metavisitor_Test_case_1-2_Guided.ga"
  - "extra-files/metavisitor/Galaxy-Workflow-Metavisitor_Test_case_1-3_Guided.ga"
  - "extra-files/metavisitor/Galaxy-Workflow-Meta-visitor__test_case_Nora_virus,_REMAPPING.ga"
```

# Running the playbook

As per usual, run the playbook with an inventory file that maps your target machine to the metavisitor group.
If the target is localhost, your inventory file should look ike this:

```
[metavisitor]
localhost
```

then run the playbook like so:

```
ansible-playbook --inventory-file=<your_inventory_file> galaxy.yml
```
