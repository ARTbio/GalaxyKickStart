#this code is a simplification of the original python script https://github.com/galaxyproject/ansible-galaxy-tools/blob/master/files/install_tool_shed_tools.py
#This is temporary, because it is only a manner of uncomment some lines in the original to make the original work 

import sys
sys.path.insert(0, '../../galaxyprojectdotorg.galaxy-tools/files/')

from install_tool_shed_tools import _setup_global_logger
from install_tool_shed_tools import ArgumentParser
from install_tool_shed_tools import load_input_file
from install_tool_shed_tools import galaxy_instance
from install_tool_shed_tools import dt
from install_tool_shed_tools import ConnectionError
from install_tool_shed_tools import time
#from install_tool_shed_tools import run_data_managers
#from install_tool_shed_tools import _parse_cli_options

#this is because "-d" option is commented in the original
def _parse_cli_options():
    """
    Parse command line options, returning `parse_args` from `ArgumentParser`.
    """
    parser = ArgumentParser(usage="usage: python %(prog)s <options>")
    parser.add_argument("-d", "--dbkeysfile",
                        dest="dbkeys_list_file",
                        help="Reference genome dbkeys to install (see "
                             "dbkeys_list.yaml.sample)",)
    parser.add_argument("-g", "--galaxy",
                        dest="galaxy_url",
                        help="Target Galaxy instance URL/IP address (required "
                             "if not defined in the tools list file)",)
    parser.add_argument("-a", "--apikey",
                        dest="api_key",
                        help="Galaxy admin user API key (required if not "
                             "defined in the tools list file)",)
    return parser.parse_args()

#this is because log is global
def run_data_managers(options):
    """
    Run Galaxy Data Manager to download, index, and install reference genome
    data into Galaxy.

    :type options: OptionParser object
    :param options: command line arguments parsed by OptionParser
    """
    dbkeys_list_file = options.dbkeys_list_file
    kl = load_input_file(dbkeys_list_file)  # Input file contents
    dbkeys = kl['dbkeys']  # The list of dbkeys to install
    dms = kl['data_managers']  # The list of data managers to run
    galaxy_url = options.galaxy_url or kl['galaxy_instance']
    api_key = options.api_key or kl['api_key']
    gi = galaxy_instance(galaxy_url, api_key)

    istart = dt.datetime.now()
    errored_dms = []
    dbkey_counter = 0
    for dbkey in dbkeys:
        dbkey_counter += 1
        dbkey_name = dbkey.get('dbkey')
        dm_counter = 0
        for dm in dms:
            dm_counter += 1
            dm_tool = dm.get('id')
            # Initate tool installation
            start = dt.datetime.now()
            log.debug('[dbkey {0}/{1}; DM: {2}/{3}] Installing dbkey {4} with '
                      'DM {5}'.format(dbkey_counter, len(dbkeys), dm_counter,
                                      len(dms), dbkey_name, dm_tool))
            tool_input = dbkey
            try:
                response = gi.tools.run_tool('', dm_tool, tool_input)
                jobs = response.get('jobs', [])
                # Check if a job is actually running
                if len(jobs) == 0:
                    log.warning("\t(!) No '{0}' job found for '{1}'".format(dm_tool,
                                dbkey_name))
                    errored_dms.append({'dbkey': dbkey_name, 'DM': dm_tool})
                else:
                    # Monitor the job(s)
                    log.debug("\tJob running", extra={'same_line': True})
                    done_count = 0
                    while done_count < len(jobs):
                        done_count = 0
                        for job in jobs:
                            job_id = job.get('id')
                            job_state = gi.jobs.show_job(job_id).get('state', '')
                            if job_state == 'ok':
                                done_count += 1
                            elif job_state == 'error':
                                done_count += 1
                                errored_dms.append({'dbkey': dbkey_name, 'DM': dm_tool})
                        log.debug("", extra={'same_line': True})
                        time.sleep(10)
                    log.debug("\tDbkey '{0}' installed successfully in '{1}'".format(
                              dbkey.get('dbkey'), dt.datetime.now() - start))
            except ConnectionError, e:
                response = None
                end = dt.datetime.now()
                log.error("\t* Error installing dbkey {0} for DM {1} (after {2}): {3}"
                          .format(dbkey_name, dm_tool, end - start, e.body))
                errored_dms.append({'dbkey': dbkey_name, 'DM': dm_tool})
    log.info("All dbkeys & DMs listed in '{0}' have been processed.".format(dbkeys_list_file))
    log.info("Errored DMs: {0}".format(errored_dms))
    log.info("Total run time: {0}".format(dt.datetime.now() - istart))

#this is because in the original "run_data_managers" call is commented
if __name__ == "__main__":
    global log
    log = _setup_global_logger()
    options = _parse_cli_options()
    if options.dbkeys_list_file:
        run_data_managers(options)
    else:
        log.error("Must provide the tool list file or individual tools info; "
                  "look at usage.")