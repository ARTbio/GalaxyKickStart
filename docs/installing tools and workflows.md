# Installing tools
----

This playbook includes the [ansible-galaxy-tools](https://github.com/galaxyproject/ansible-galaxy-tools) role which can be used
to install tools and workflows into galaxy instances using the [bioblend](https://bioblend.readthedocs.org/en/latest/) api.  


# Creating a tool_list.yml file
To install tools, you will need to prepare a list of tools in yaml format.
A an example of a a tool list can be found in [here](https://github.com/ARTbio/ansible-artimed/blob/master/extra-files/metavisitor/metavisitor_tool_list.yml)
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
when the revision is empty, the latest available revision will be installed.  
tool_panel_section_label will determine the tool panel section where the tools will be found.

# Obtaining a tool_list.yml file 

We can also obtain a tool list from a runnning galaxy instance.
Note that for server running a galaxy release <16.04, you need a galaxy API keys and bioblend.
A script is included in the extra-files directory.
```
python get_tool_yml_from_gi.py --galaxy <my_galaxy_url> --api-key <my_admin_api_key> --output-file <my_tool_list.yml>
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