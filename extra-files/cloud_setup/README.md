# Extra Files / Cloud Setup

The files described below are the ones the user can modify in order to meet their specifications.

```
├── data_manager_tasks.yml
├── data_manager_tool_list.yml
├── dependency_resolvers_conf.xml
├── fix_hostname.sh
├── tool_sheds_conf.xml
├── usegalaxy_main_tool_list.yml
└── vimrc
```

* `vimrc` : The VIM configurations used for the remote target.

* `usegalaxy_main_tool_list.yml` : This is a replica of the Galaxy main tool list for the most part, but edited to fit the size of the instance. The Galaxy main tool list can be obtained using this ephemeris [script](https://github.com/galaxyproject/ephemeris/blob/master/ephemeris/get_tool_list_from_galaxy.py). This file can be changed in group_vars/cloud_setup if the user wishes to use a custom tool list. While making a new tool list, it is important to keep in mind the tags in the `yaml` file.

This tool uses the `tool_panel_section_id` because it needs to be placed in a section which is shipped with the base Galaxy instance i.e default tool sections.
```
- name: tabular_to_fastq
  owner: devteam
  revisions:
  - b334cd1095ea
  tool_panel_section_id: convert  # Use ID since the group is included w/ Galaxy
  tool_shed_url: https://toolshed.g2.bx.psu.edu/
  install_resolver_dependencies: True
```

This tool uses `tool_panel_section_label` as we are making a new section in the tool panel.

```
- name: fastx_nucleotides_distribution
  owner: devteam
  revisions:
  - 63d6e2daad48
  tool_panel_section_label: 'NGS: QC and manipulation'
  tool_shed_url: https://toolshed.g2.bx.psu.edu/
  install_resolver_dependencies: True
```

 * `dependency_resolvers_conf.xml` : Modified `dependency_resolvers_conf` file with conda moved up as the first option.

 * `data_manager_tool_list.yml` : List of data managers installed.

 * `data_manager_tasks.yml` : Run data manager tasks with dbkeys (organisms).

 * `fix_hostname.sh` : Custom script to fix the hostname on Jetstream instances.

 *  `tool_sheds_conf.xml` : Two tool sheds are available right now through this file, if the user has a local toolshed, they can add it to this list.
