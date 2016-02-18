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
Download the script and run it.
```
python get_tool_from_
```