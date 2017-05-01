import os
import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import shutil
from generate_tool_list_from_ga_workflow_files import generate_tool_list_from_workflow
from string import Template

INVENTORY_FILE_TEMPLATE = '''
[GKSfromWorkflow]
localhost ansible_connection=local

# if remote target, replace the above line by
# <remote host IP> ansible_ssh_user="root" ansible_ssh_private_key_file="~/.ssh/somekey"
'''

GROUP_VAR_FILE_TEMPLATE = '''
galaxy_tools_tool_list_files:
  - "extra-files/GKSfromWorkflow/GKSfromWorkflow_tool_list.yml"

galaxy_tools_workflows:
$workflow_list

galaxy_web_processes: 2

additional_files_list:
  - { src: "extra-files/galaxy-kickstart/welcome.html", dest: "{{ galaxy_server_dir }}/static/" }
  - { src: "extra-files/galaxy-kickstart/galaxy-kickstart_logo.png", dest: "{{ galaxy_server_dir }}/static/images/" }
  - { src: "extra-files/tool_sheds_conf.xml", dest: "{{ galaxy_config_dir }}" }
'''

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
                        help='A space-separated list of galaxy workflow description files in json format', )
    parser.add_argument('-l', '--panel_label',
                        dest='panel_label',
                        default='Tools from workflows',
                        help='The name of the panel where the tools will show up in Galaxy.'
                             'If not specified: "Tools from workflows"')
    return parser.parse_args()


def makedir (path):
    if os.path.isdir(path):
        shutil.rmtree(path)
    os.makedirs(path)


def make_inventory (file_path="../inventory_files/GKSfromWorkflow"):
    if not os.path.exists (file_path):
        open(file_path, "w").write(INVENTORY_FILE_TEMPLATE)
    else:
        print("GKSfromWorkflow inventory file exists, file unchanged")

def make_groupvars (workflow_file_list, file_path="../group_vars/GKSfromWorkflow"):
    if os.path.exists(file_path):
        print("The GKSfromWorkflow group_vars file already existed and has been overwritten")
    internal_workflow_list = []
    for workflow in workflow_file_list:
        workflow = os.path.basename(workflow)
        workflow = '  - "extra-files/GKSfromWorkflow/' + workflow + '"'
        internal_workflow_list.append(workflow)
    workflow_list = "\n".join(internal_workflow_list)
    template_params = {"workflow_list": workflow_list}
    config_contents = Template(GROUP_VAR_FILE_TEMPLATE).substitute(template_params)
    open(file_path, "w").write(config_contents)
    
def make_extra_files (workflow_files, panel_label, tool_list_file, extra_files_dir="../extra-files/GKSfromWorkflow"):
    if os.path.exists(extra_files_dir):
        print("The extra-files/GKSfromWorkflow directory already existed and has been overwritten")
    makedir(extra_files_dir)
    generate_tool_list_from_workflow(workflow_files, panel_label, tool_list_file)
    os.rename(tool_list_file, extra_files_dir + "/" + tool_list_file)
    for workflow in workflow_files:
        workflow_basename = os.path.basename(workflow)
        shutil.copyfile (workflow, extra_files_dir + "/" + workflow_basename)


def create_gks_flavor (workflow_file_list, panel_label, tool_list_file):
    """
    creates inventory files in inventory_files folder
    creates a group_vars "gks_workflows" in group_vars folder
    creates a tool_list.yml files from galaxy workflow files and
    copy these files in and extra-files/gks_workflows folder
    """
    make_inventory()
    make_groupvars(workflow_file_list)
    make_extra_files (workflow_file_list, panel_label, tool_list_file)

def main():
    options = _parse_cli_options()
    create_gks_flavor (options.workflow_files, options.panel_label, "GKSfromWorkflow_tool_list.yml")

if __name__ == "__main__":
    main()