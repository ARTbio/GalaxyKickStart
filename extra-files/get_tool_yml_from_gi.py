#!/usr/bin/env python

## Example usage: python get_tool_yml_from_gi.py -g https://mississippi.snv.jussieu.fr/ -a <api_key> -o tool_list.yml

import sys
from operator import itemgetter
from argparse import ArgumentParser
from distutils.version import StrictVersion

try:
    import bioblend
    if StrictVersion( bioblend.get_version() ) < StrictVersion( "0.6.0" ):
        print("This script needs bioblend >= 0.6.0.")
        sys.exit(1)
    import bioblend.galaxy
except ImportError:
    print("bioblend is not installed, this script needs bioblend >= 0.6.0.")
    raise(ImportError)
import yaml

class GiToToolYaml:
    def __init__(self, url,
                 api_key,
                 output_file,
                 get_deleted=False,
                 get_packages=False,
                 get_latest_installed=False,
                 include_tool_panel_section_id=False,
                 skip_tool_panel_section_name=True):

        self.url = url
        self.api_key = api_key
        self.output_file = output_file
        self.get_deleted = get_deleted
        self.get_packages = get_packages
        self.get_latest_installed = get_latest_installed
        self.include_tool_panel_section_id = include_tool_panel_section_id
        self.skip_tool_panel_section_name = skip_tool_panel_section_name
        self.gi = self.get_instance()
        self.galaxy_client = bioblend.galaxy.config.ConfigClient( self.gi )
        self.galaxy_config = self.galaxy_client.get_config()
        if not self.galaxy_config[u'is_admin_user']:
            raise ValueError("Current API key is not an admin api-key.")
        self.repository_list = self.get_tools_from_gi()
        self.filtered_repository_list = self.filter_tools()
        self.revisions_to_list()
        self.tool_panel_list = self.get_tool_panel_list_from_gi()
        self.tool_panel_list_filtered = self.filter_tool_panel_list()
        self.add_tool_panel_section()
        self.tool_list = self.repo_to_tool_list()
        self.remove_none()
        self.filtered_tool_list = self.remove_suites()
        self.merge_dicts_in_list()
        self.write_to_yaml()

    def get_instance(self):
        return bioblend.galaxy.GalaxyInstance(url=self.url, key=self.api_key)

    def get_tools_from_gi(self):
        """
        """
        return self.gi.toolShed.get_repositories()

    def get_tool_panel_list_from_gi(self):
        return self.gi.tools.get_tool_panel()

    def parse_panel_section(self, tool, tool_panel_list_filtered):
        if tool[u'model_class'] == u'Tool':
            try:
                tool_panel_section_id = tool["panel_section_id"]
                tool_panel_section_label = tool["panel_section_name"]
                split_repo_url = tool["id"].split("/repos/")
                sub_dir = split_repo_url[1].split("/")
                tool_shed = "https://" + split_repo_url[0] + "/"
                owner = sub_dir[0]
                name = sub_dir[1]
                tool_panel_list_filtered.append(
                    {"tool_panel_section_id": tool_panel_section_id,
                    "tool_panel_section_label": tool_panel_section_label,
                    "tool_shed": tool_shed,
                    "owner": owner,
                    "name": name}
                )
            except IndexError:  # when tools are not coming from the toolshed
                pass

    def filter_tool_panel_list(self):
        tool_panel_list = self.tool_panel_list
        tool_panel_list_filtered = []
        for section in tool_panel_list:
            elems = section.get("elems", None)
            if elems:
                for tool in elems:
                    self.parse_panel_section(tool, tool_panel_list_filtered)
            else:
                try:
                    self.parse_panel_section(section, tool_panel_list_filtered)
                except Exception:
                    print("could not parse section, continuing anyway")
                    continue

        return tool_panel_list_filtered

    def get_tool_panel_section_id(self, required_values):
        required_fields = ("tool_shed",
                           "owner",
                           "name")
        for tool in self.tool_panel_list_filtered:
            tool_values = itemgetter(*required_fields)(tool)
            if tool_values == required_values:
                return tool["tool_panel_section_id"]

    def get_tool_panel_section_label(self, required_values):
        required_fields = ("tool_shed",
                           "owner",
                           "name")
        for tool in self.tool_panel_list_filtered:
            tool_values = itemgetter(*required_fields)(tool)
            if tool_values == required_values:
                return tool["tool_panel_section_label"]

    def filter_tools(self):
        filtered_repository_list = []
        for repo in self.repository_list:
            if (repo['deleted'] is True) and (self.get_deleted is False):
                continue
            if self.get_latest_installed is True and not (self.is_latest_installed(repo)):
                print("skip because is not latest rev.")
                continue
            if self.get_packages is False:
                if repo["name"].startswith("package_"):
                    continue
            repo["tool_shed"] = "https://" + repo["tool_shed"] + "/"
            filtered_repository_list.append(repo)
        return filtered_repository_list

    def is_latest_installed(self, current_repository):
        unique_repo = ("name", "owner", "tool_shed")
        current_repo_comb = itemgetter(*unique_repo)(current_repository)
        current_rev = int(current_repository["ctx_rev"])
        for repo in self.repository_list:
            repo_comb = itemgetter(*unique_repo)(repo)
            repo_rev = int(repo["ctx_rev"])
            if (repo_comb == current_repo_comb) and (current_rev < repo_rev):
                return False
            else:
                continue
        return True

    def add_tool_panel_section(self):
        required_fields = ("tool_shed",
                           "owner",
                           "name")
        for repo in self.filtered_repository_list:
            required_values = itemgetter(*required_fields)(repo)
            tool_panel_section_id = self.get_tool_panel_section_id(required_values)
            tool_panel_section_label = self.get_tool_panel_section_label(required_values)
            if self.include_tool_panel_section_id:
                repo["tool_panel_section_id"] = tool_panel_section_id
            if not self.skip_tool_panel_section_name:
                repo["tool_panel_section_label"] = tool_panel_section_label

    def revisions_to_list(self):
        for repo in self.filtered_repository_list:
            repo["installed_changeset_revision"] = [repo["installed_changeset_revision"]]

    def repo_to_tool_list(self):
        filtered_repository_list = self.filtered_repository_list
        required_fields = ["name",
                           "owner",
                           "tool_panel_section_id",
                           "tool_panel_section_label",
                           "installed_changeset_revision",
                           "tool_shed"]
        yaml_categories = ["name",
                           "owner",
                           "tool_panel_section_id",
                           "tool_panel_section_label",
                           "revisions",
                           "tool_shed_url"]
        if not self.include_tool_panel_section_id:
            required_fields.remove("tool_panel_section_id")
            yaml_categories.remove("tool_panel_section_id")
        if self.skip_tool_panel_section_name:
            required_fields.remove("tool_panel_section_label")
            yaml_categories.remove("tool_panel_section_label")
        tool_list = []
        for repo in filtered_repository_list:
            values = itemgetter(*required_fields)(repo)
            tool_list.append(dict(zip(yaml_categories, values)))
        return tool_list

    def merge_dicts_in_list(self):
        tool_list = self.filtered_tool_list
        for current_tool in tool_list:
            for tool in tool_list:
                if current_tool == tool:
                    continue
                if (tool["name"] == current_tool["name"]) and (tool["owner"] == current_tool["owner"]):
                    current_tool["revisions"].extend(tool["revisions"])
                    tool_list.remove(tool)

    def remove_suites(self):
        """
        Since there appears to be no way to find out
        into which section toolsuites have been installed,
        filter them out by the presence of "suite" in name
        from tool_list.
        """
        filtered_tool_list = []
        for tool in self.tool_list:
            if "suite" in tool["name"]:
                continue
            filtered_tool_list.append(tool)
        return filtered_tool_list

    def remove_none(self):
        for tool in self.tool_list:
            try:
                if tool["tool_panel_section_id"] is None:
                    del tool["tool_panel_section_id"]
                    del tool["tool_panel_section_name"]
            except KeyError:
                continue


    def write_to_yaml(self):
	tool_dict={"tools":self.filtered_tool_list}
        with open(self.output_file, "w") as output:
            output.write(yaml.safe_dump(tool_dict, default_flow_style=False))


def _parse_cli_options():
    """
    Parse command line options, returning `parse_args` from `ArgumentParser`.
    """
    parser = ArgumentParser(usage="usage: python %(prog)s <options>")
    parser.add_argument("-g", "--galaxy",
                        dest="galaxy_url",
                        required=True,
                        help="Target Galaxy instance URL/IP address (required "
                             "if not defined in the tools list file)", )
    parser.add_argument("-a", "--api-key",
                        required=True,
                        dest="api_key",
                        help="Galaxy user API key")
    parser.add_argument("-o", "--output-file",
                        required=True,
                        dest="output",
                        help="tool_list.yml output file")
    parser.add_argument("-d", "--get_deleted",
                        dest="get_deleted",
                        type=bool,
                        default=False,
                        help="Include deleted repositories in tool_list.yml ?")
    parser.add_argument("-p", "--get_packages",
                        dest="get_packages",
                        type=bool,
                        default=False,
                        help="Include packages in tool_list.yml?")
    parser.add_argument("-l", "--get_latest",
                        action="store_true",
                        default=False,
                        help="Include only latest revision of a tool version in tool_list.yml ?")
    parser.add_argument("-include_id", "--include_tool_panel_id",
                        action="store_true",
                        default=False,
                        help="Include tool_panel_id in tool_list.yml ? "
                             "Use this only if the tool panel id already exists. See "
                             "https://github.com/galaxyproject/ansible-galaxy-tools/blob/master/files/tool_list.yaml.sample")
    parser.add_argument("-skip_name", "--skip_tool_panel_name",
                        action="store_true",
                        default=False,
                        help="Do not include tool_panel_name in tool_list.yml ?")
    return parser.parse_args()


if __name__ == "__main__":
    options = _parse_cli_options()
    get_tool_yml_instance = GiToToolYaml(url=options.galaxy_url,
                                         api_key=options.api_key,
                                         output_file=options.output,
                                         get_deleted=options.get_deleted,
                                         get_packages=options.get_packages,
                                         get_latest_installed=options.get_latest,
                                         include_tool_panel_section_id=options.include_tool_panel_id,
                                         skip_tool_panel_section_name=options.skip_tool_panel_name)
