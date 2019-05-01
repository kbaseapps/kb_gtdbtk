# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

from installed_clients.KBaseReportClient import KBaseReport
from .utils.misc_utils import load_fastas
from .utils.misc_utils import create_html_report
from .utils.gtdbtk_utils import GTDBTkUtils
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
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        self.config = config
        self.config['callback_url'] = self.callback_url
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

        # TODO: some parameter checking
        try:
            ref = params.get('inputObjectRef')
        except KeyError:
            print("Must provide a ws reference to object with sequences")

        min_perc_aa = params.get('min_perc_aa', 10)

        try:
            workspace_id = params.get('workspace_id')
        except KeyError:
            print("Must provide a workspace id")

        # get the fasta file from the input ref
        # TODO: handle sets
        logging.info("Get Genome Seqs\n")
        fasta_paths = load_fastas(self.config, self.shared_folder, ref)
        print(fasta_paths)

        logging.info("Run gtdbtk classifywf\n")
        gtdbtku = GTDBTkUtils(self.config, self.callback_url, workspace_id, self.cpus)
        results = gtdbtku.gtdbtk_classifywf(fasta_paths, min_perc_aa)
        logging.info(results)
        output = create_html_report(self.callback_url, self.shared_folder, params['workspace_name'])
        
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
