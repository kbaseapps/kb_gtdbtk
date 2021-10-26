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
from kb_gtdbtk.core.genome_obj_update import check_obj_type_genome, update_genome_objs_class
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
    VERSION = "0.1.7"
    GIT_URL = "https://github.com/kbaseapps/kb_gtdbtk"
    GIT_COMMIT_HASH = "2b1fb26322612fd24d388f0a086f793b7297c627"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = Path(config['scratch'])
        self.ws_url = config['workspace-url']
        self.cpus = 32  # bigmem 32 cpus & 251 GB RAM
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def run_kb_gtdbtk(self, ctx, params):
        """
        Run GTDB-tk Classify (deprecated method name)
        :param params: instance of type "GTDBtk_Classify_Params" (Parameters
           for the GTDB-tk Classify (classify_wf) run. Required:
           input_object_ref: A reference to the workspace object to process.
           workspace_id: The integer workspace ID where the results will be
           saved. Optional: min_perc_aa: the minimum sequence alignment as a
           percent, default 10.) -> structure: parameter "input_object_ref"
           of String, parameter "workspace_id" of Long, parameter
           "min_perc_aa" of Double, parameter "overwrite_tax" of type "bool"
        :returns: instance of type "ReportResults" (The results of the
           GTDB-tk run. report_name: The name of the report object in the
           workspace. report_ref: The UPA of the report object, e.g.
           wsid/objid/ver.) -> structure: parameter "report_name" of String,
           parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_kb_gtdbtk
        logging.info("run_kb_gtdbtk() is deprecated.  Use run_kb_gtdbtk_classify_wf() instead.\n")
        output = dict()
        #END run_kb_gtdbtk

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_kb_gtdbtk return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def run_kb_gtdbtk_classify_wf(self, ctx, params):
        """
        Run GTDB-tk Classify
        :param params: instance of type "GTDBtk_Classify_Params" (Parameters
           for the GTDB-tk Classify (classify_wf) run. Required:
           input_object_ref: A reference to the workspace object to process.
           workspace_id: The integer workspace ID where the results will be
           saved. Optional: min_perc_aa: the minimum sequence alignment as a
           percent, default 10.) -> structure: parameter "input_object_ref"
           of String, parameter "workspace_id" of Long, parameter
           "min_perc_aa" of Double, parameter "overwrite_tax" of type "bool"
        :returns: instance of type "ReportResults" (The results of the
           GTDB-tk run. report_name: The name of the report object in the
           workspace. report_ref: The UPA of the report object, e.g.
           wsid/objid/ver.) -> structure: parameter "report_name" of String,
           parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_kb_gtdbtk_classify_wf

        params = get_gtdbtk_params(params)

        logging.info("Get Genome Seqs\n")
        fasta_path = self.shared_folder / 'fastas'
        fasta_path.mkdir(parents=True, exist_ok=True)

        cli = KBClients(self.callback_url, self.ws_url, ctx['token'])

        path_to_filename = download_sequence(params.ref, fasta_path, cli)
        for path, fn in path_to_filename.items():
            print(fn, path)

        logging.info("Run gtdbtk classify_wf\n")

        output_path = self.shared_folder / 'output'
        temp_output = self.shared_folder / 'temp_output'
        output_path.mkdir(parents=True, exist_ok=True)
        temp_output.mkdir(parents=True, exist_ok=True)

        def runner(args):
            env = dict(os.environ)
            env['TEMP_DIR'] = str(self.shared_folder / 'tmp')
            # should print to stdout/stderr
            subprocess.run(args, check=True, env=env)

        classification = run_gtdbtk(
            runner, path_to_filename, output_path, temp_output, params.min_perc_aa, self.cpus)

        objects_created = None
        if check_obj_type_genome (params.ref, cli):
            objects_created = update_genome_objs_class (params.ref, classification, params.overwrite_tax, cli)
        
        output = generate_report(cli, output_path, params.workspace_id, objects_created)

        #END run_kb_gtdbtk_classify_wf

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_kb_gtdbtk_classify_wf return value ' +
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
