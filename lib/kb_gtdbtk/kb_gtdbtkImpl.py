# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import subprocess

from pathlib import Path

from kb_gtdbtk.core.api_translation import get_gtdbtk_params
from kb_gtdbtk.core.sequence_downloader import download_sequence
from kb_gtdbtk.core.kb_client_set import KBClients
from kb_gtdbtk.core.gtdbtk_runner import run_gtdbtk
from kb_gtdbtk.core.kb_report_generation import generate_report
#END_HEADER


class kb_gtdbtk:
    '''
    Module Name:
    kb_gtdbtk

    Module Description:
    A KBase module: kb_gtdbtk
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.1.1"
    GIT_URL = "https://github.com/mrcreosote/kb_gtdbtk.git"
    GIT_COMMIT_HASH = "c3241b7286ad0047a114429f7b7c51c9442d92fe"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = Path(config['scratch'])
        self.cpus = 32  # bigmem 32 cpus & 90,000MB RAM
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass

    def run_kb_gtdbtk(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_kb_gtdbtk

        # TODO remove unused files
        # TODO put actual params in spec
        params = get_gtdbtk_params(params)

        logging.info("Get Genome Seqs\n")
        fasta_path = self.shared_folder / 'fastas'
        fasta_path.mkdir(parent=True, exist_ok=True)

        cli = KBClients(self.callback_url, ctx['token'])

        path_to_filename = download_sequence(params.upa, fasta_path, cli)
        for path, fn in path_to_filename.items():
            print(fn, path)

        logging.info("Run gtdbtk classifywf\n")

        output_path = Path(self.shared_folder) / 'output'
        temp_output = Path(self.shared_folder) / 'temp_output'
        output_path.mkdir(parent=True, exist_ok=True)
        temp_output.mkdir(parent=True, exist_ok=True)

        def runner(args):
            env = dict(os.environ)
            env['TEMP_DIR'] = self.shared_folder / 'tmp'
            # should print to stdout/stderr
            subprocess.run(args, check=True, env=env)

        run_gtdbtk(
            runner, path_to_filename, output_path, temp_output, params.min_perc_aa, self.cpus)

        output = generate_report(cli, output_path, params.workspace_id)

        #END run_kb_gtdbtk

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_kb_gtdbtk return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
