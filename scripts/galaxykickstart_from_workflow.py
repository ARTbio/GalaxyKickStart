import os
import sys
import textwrap
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from shutil import copyfile
from generate_tool_list_from_ga_workflow_files import *



def _parse_cli_options():
    """
    Parse command line options, returning `parse_args` from `ArgumentParser`.
    """
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                            usage="python %(prog)s <options>",
                            epilog="example:\n"
                                   "python %(prog)s -w workflow1 workflow2 -l my_panel_label\n"
                                   "Christophe Antoniewski <drosofff@gmail.com>\n"
                                   "https://github.com/ARTbio/ansible-artimed/tree/master/scritps/galaxykickstart_from_workflow.py")
    parser.add_argument('-w', '--workflow',
                        dest="workflow_files",
                        required=True,
                        nargs='+',
                        help='A comma-separated of galaxy workflow description files in json format', )
    parser.add_argument('-l', '--panel_label',
                        dest='panel_label',
                        default='Tools from workflows',
                        help='The name of the panel where the tools will show up in Galaxy.'
                             'If not specified: "Tools from workflows"')
    return parser.parse_args()

def makedir (path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

def make_inventory (path="../inventory_files/"):
    file_path = path + "gks_workflows"
    try:
        f = os.fdopen(os.open(file_path, os.O_CREAT | os.O_WRONLY | os.O_EXCL), 'w')
    except OSError:
        sys.exit("gks_workflows inventory file already exists in ../inventory_files")
    print >> f, textwrap.dedent('''\
            [gks_workflows]
            localhost ansible_connection=local

            # if remote target, replace the above line by
            # <remote host IP> ansible_ssh_user="root" ansible_ssh_private_key_file="<path/to/your/private/key>"
            ''')
    f.close()

def make_gks_workflows_groupvars (path="../group_vars/gks_workflows",
                                  template="galaxykickstart_from_workflow_templates/group_vars_template.yml"):
    if os.path.exists(path):
        sys.exit("../group_vars/gks_workflows already exists.")
    copyfile(template, path)

def make_gks_workflows_extra_files (workflow_files, panel_label, tool_list_file, extra_files_dir="../extra-files/gks_workflows"):
    makedir(extra_files_dir)
    copyfile("galaxykickstart_from_workflow_templates/tool_sheds_conf.xml.sample", extra_files_dir + "/tool_sheds_conf.xml.sample")
    generate_tool_list_from_workflow(workflow_files, panel_label, tool_list_file)
    copyfile(tool_list_file, extra_files_dir + "/" + tool_list_file)
    for file in workflow_files:
        file = os.path.basename(file)
        copyfile (file, extra_files_dir + "/" + file)





def __main__(workflow_file_list, panel_label, tool_list_file):
    '''
    creates inventory files in inventory_files folder
    creates a group_vars "gks_workflows" in group_vars folder
    creates a tool_list.yml files from galaxy workflow files and
    copy these files in and extra-files/gks_workflows folder
    '''
    make_inventory()
    make_gks_workflows_groupvars()
    make_gks_workflows_extra_files (workflow_file_list, panel_label, tool_list_file)

if __name__ == "__main__":
    options = _parse_cli_options()
    __main__(options.workflow_files, options.panel_label, "tool_list.yml")