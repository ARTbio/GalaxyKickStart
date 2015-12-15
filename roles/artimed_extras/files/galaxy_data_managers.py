import datetime as dt
import logging
import time
import yaml
from argparse import ArgumentParser

from bioblend.galaxy import GalaxyInstance
from bioblend.galaxy.toolshed import ToolShedClient
from bioblend.toolshed import ToolShedInstance
from bioblend.galaxy.client import ConnectionError

# Omit (most of the) logging by external libraries
logging.getLogger('bioblend').setLevel(logging.ERROR)
logging.getLogger('requests').setLevel(logging.ERROR)
logging.captureWarnings(True)  # Capture HTTPS warngings from urllib3

MTS = 'https://toolshed.g2.bx.psu.edu/'  # Main Tool Shed

class ProgressConsoleHandler(logging.StreamHandler):
    """
    A handler class which allows the cursor to stay on
    one line for selected messages
    """
    on_same_line = False

    def emit(self, record):
        try:
            msg = self.format(record)
            stream = self.stream
            same_line = hasattr(record, 'same_line')
            if self.on_same_line and not same_line:
                stream.write('\r\n')
            stream.write(msg)
            if same_line:
                stream.write('.')
                self.on_same_line = True
            else:
                stream.write('\r\n')
                self.on_same_line = False
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)

def _setup_global_logger():
    formatter = logging.Formatter('%(asctime)s %(levelname)-5s - %(message)s')
    progress = ProgressConsoleHandler()
    file_handler = logging.FileHandler('/tmp/galaxy_tool_install.log')
    console = logging.StreamHandler()
    console.setFormatter(formatter)

    logger = logging.getLogger('test')
    logger.setLevel(logging.DEBUG)
    logger.addHandler(progress)
    logger.addHandler(file_handler)
    return logger

def load_input_file(tool_list_file='tool_list.yaml'):
    """
    Load YAML from the `tool_list_file` and return a dict with the content.
    """
    with open(tool_list_file, 'r') as f:
        tl = yaml.load(f)
    return tl

def galaxy_instance(url=None, api_key=None):
    """
    Get an instance of the `GalaxyInstance` object. If the arguments are not
    provided, load the default values using `load_input_file` method.
    """
    if not (url and api_key):
        tl = load_input_file()
        url = tl['galaxy_instance']
        api_key = tl['api_key']
    return GalaxyInstance(url, api_key)

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
    parser.add_argument("-t", "--toolsfile",
                        dest="tool_list_file",
                        help="Tools file to use (see tool_list.yaml.sample)",)
    parser.add_argument("--name",
                        help="The name of the tool to install (only applicable "
                             "if the tools file is not provided).")
    parser.add_argument("--owner",
                        help="The owner of the tool to install (only applicable "
                             "if the tools file is not provided).")
    parser.add_argument("--section",
                        dest="tool_panel_section_id",
                        help="Galaxy tool panel section ID where the tool will "
                             "be installed (the section must exist in Galaxy; "
                             "only applicable if the tools file is not provided).")
    parser.add_argument("--toolshed",
                        dest="tool_shed_url",
                        help="The Tool Shed URL where to install the tool from. "
                             "This is applicable only if the tool info is "
                             "provided as an option vs. in the tools file.")
    return parser.parse_args()

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

if __name__ == "__main__":
    global log
    log = _setup_global_logger()
    options = _parse_cli_options()
    if options.dbkeys_list_file:
        run_data_managers(options)
    else:
        log.error("Must provide the tool list file or individual tools info; "
                  "look at usage.")
