# Extra Files / Cloud Setup

The files described below are the ones the user can modify in order to meet
their specifications.

```
├── data_manager_tasks.yml
├── data_manager_tool_list.yml
├── dependency_resolvers_conf.xml
├── fix_hostname.sh
├── tool_sheds_conf.xml
├── usegalaxy_main_tool_list.yml
└── vimrc
```

* `vimrc`: The VIM configurations used for the remote target.

* `usegalaxy_main_tool_list.yml`: This is largely a replica of the Galaxy Main
tool list, but edited to include the tools that install cleanly. The Galaxy
Main tool list can be obtained using this ephemeris [script](https://g
ithub.com/galaxyproject/ephemeris/blob/master/ephemeris/get_tool_list_from_ga
laxy.py). This file can be changed in `group_vars/cloud_setup` if the user
wishes to use a custom tool list.

If making a new tool list, it is important to keep in mind the following: if a
tool is being installed into a tool section that exists natively in Galaxy
code base, it is necessary to use `tool_panel_section_id` of the existing tool
section. For example:

```
- name: tabular_to_fastq
  owner: devteam
  revisions:
  - b334cd1095ea
  tool_panel_section_id: convert  # Use ID since the group is included w/ Galaxy
  tool_shed_url: https://toolshed.g2.bx.psu.edu/
  install_resolver_dependencies: True
```

On the other hand, if you the tool is being added to a new section, use
`tool_panel_section_label`. It is OK to specify the same label for multiple
tools.

```
- name: fastx_nucleotides_distribution
  owner: devteam
  revisions:
  - 63d6e2daad48
  tool_panel_section_label: 'NGS: QC and manipulation'
  tool_shed_url: https://toolshed.g2.bx.psu.edu/
  install_resolver_dependencies: True
```

 * `data_manager_tool_list.yml`: List of data managers installed.

 * `data_manager_tasks.yml`: Run data manager tasks with dbkeys (organisms).

 * `fix_hostname.sh`: Custom script to fix the hostname on Jetstream instances.

 * `tool_sheds_conf.xml`: Two tool sheds are available right now through this
  file (Galaxy Main and Test ToolSheds); users can reference other ToolSheds
  via this file.
