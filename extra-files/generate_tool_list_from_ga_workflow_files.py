#!/usr/bin/env python

## Christophe Antoniewski <drosofff@gmail.com>
## https://github.com/ARTbio/ansible-artimed/tree/master/extra-files/generate_tool_list_from_ga_workflow_files.py
## Usage example:
## generate_tool_list_from_ga_workflow_files.py -w workflow1 workflow2 -o mytool_list.yml -l my_panel_label

import yaml
import json
import numpy as np
from argparse import ArgumentParser


def _parse_cli_options():
    """
    Parse command line options, returning `parse_args` from `ArgumentParser`.
    """
    parser = ArgumentParser(usage="usage: python %(prog)s <options>")
    parser.add_argument('-w', '--workflow',
                        dest="workflow_files",
                        required=True,
                        nargs='+',
                        help='A comma-separated of galaxy workflow description files in json format', )
    parser.add_argument('-o', '--output-file',
                        required=True,
                        dest='output_file',
                        help='tool_list.yml output file')
    parser.add_argument('-l', '--panel_label',
                        dest='panel_label',
                        default='Tools from workflows',
                        help='Name of the panel where the tools will show up.'
                             ' If not specified, workflow file base name is taken')
    return parser.parse_args()


def get_workflow_dictionary(json_file):
    with open(json_file, "r") as File:
        mydict = json.load(File)[u'steps']
    return mydict


def translate_workflow_dictionary_to_tool_list(tool_dictionary):
    extract_dict = {}
    extract_dict[u'tools'] = []
    for step in tool_dictionary:
        if tool_dictionary[step][u'tool_id'] and 'toolshed' in tool_dictionary[step][u'tool_id']:
            extract_dict[u'tools'].append(tool_dictionary[step][u'tool_shed_repository'])
    tool_list = []
    for itemlist in extract_dict[u'tools']:
        sub_dic = {'name': itemlist['name'], 'owner': itemlist['owner'], 'revision': itemlist['changeset_revision'],
                      'tool_panel_section_label': options.panel_label, 'tool_shed_url': itemlist['tool_shed']}
        tool_list.append(sub_dic)
    return tool_list

def print_yaml_tool_list(tool_dictionary, output_file):
    with open(output_file, 'w') as F:
        F.write(yaml.safe_dump(tool_dictionary, default_flow_style=False))
    return

if __name__ == "__main__":
    options = _parse_cli_options()
    intermediate_tool_list = []
    for workflow in options.workflow_files:
        workflow_dictionary = get_workflow_dictionary (workflow)
        intermediate_tool_list += translate_workflow_dictionary_to_tool_list (workflow_dictionary)
    intermediate_tool_list = list(np.unique(np.array(intermediate_tool_list)))
    convert_dic = {}
    convert_dic['tools'] = intermediate_tool_list
    print_yaml_tool_list(convert_dic, options.output_file)
