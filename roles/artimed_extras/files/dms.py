import sys
sys.path.insert(0, '../../galaxyprojectdotorg.galaxy-tools/files/')

import install_tool_shed_tools

if __name__ == "__main__":
    global log
    log = _setup_global_logger()
    options = _parse_cli_options()
    if options.dbkeys_list_file:
        run_data_managers(options)
    else:
        log.error("Must provide the tool list file or individual tools info; "
                  "look at usage.")