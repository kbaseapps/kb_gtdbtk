# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import sys
import subprocess

from datetime import datetime
from pprint import pprint, pformat
from pathlib import Path

from kb_gtdbtk.core.api_translation import get_gtdbtk_params
from kb_gtdbtk.core.sequence_downloader import download_sequence
from kb_gtdbtk.core.kb_client_set import KBClients
from kb_gtdbtk.core.gtdbtk_runner import run_gtdbtk
from kb_gtdbtk.core.krona_runner import run_krona_import_text
from kb_gtdbtk.core.kb_report_generation import generate_report
from kb_gtdbtk.core.genome_obj_update import copy_gtdb_species_reps, get_obj_type, check_obj_type_genome, check_obj_type_assembly, update_genome_assembly_objs_class, process_tree_files, save_gtdb_tree_objs
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
    VERSION = "1.4.0"
    GIT_URL = "https://github.com/kbaseapps/kb_gtdbtk"
    GIT_COMMIT_HASH = "a9b1858968e9d181f989febe0ae4921ec9cdce10"

    #BEGIN_CLASS_HEADER

    ### now_ISOish()
    #
    def now_ISOish(self):
        now_timestamp = datetime.now()
        now_secs_from_epoch = (now_timestamp - datetime(1970,1,1)).total_seconds()
        now_timestamp_in_iso = datetime.fromtimestamp(int(now_secs_from_epoch)).strftime('%Y-%m-%d_%T')
        return now_timestamp_in_iso

    
    ### log()
    #
    def log(self, target, message):
        message = '['+self.now_ISOish()+'] '+message
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = Path(config['scratch'])
        self.ws_url = config['workspace-url']
        self.hs_url = config['handle-service-url']
        self.cpus = config['cpus']  # bigmem 32 cpus & 251 GB RAM.  new gtdb-tk needs less mem.
        self.genome_upas_map_file = config['genome_upas_map_file']
        
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
           percent, default 10.) -> structure: parameter "workspace_id" of
           Long, parameter "input_object_ref" of String, parameter
           "output_tree_basename" of String, parameter "copy_proximals" of
           type "bool", parameter "save_trees" of type "bool", parameter
           "min_perc_aa" of Double, parameter "db_ver" of Long, parameter
           "keep_intermediates" of type "bool", parameter "overwrite_tax" of
           type "bool", parameter "dendrogram_report" of type "bool"
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
           percent, default 10.) -> structure: parameter "workspace_id" of
           Long, parameter "input_object_ref" of String, parameter
           "output_tree_basename" of String, parameter "copy_proximals" of
           type "bool", parameter "save_trees" of type "bool", parameter
           "min_perc_aa" of Double, parameter "db_ver" of Long, parameter
           "keep_intermediates" of type "bool", parameter "overwrite_tax" of
           type "bool", parameter "dendrogram_report" of type "bool"
        :returns: instance of type "ReportResults" (The results of the
           GTDB-tk run. report_name: The name of the report object in the
           workspace. report_ref: The UPA of the report object, e.g.
           wsid/objid/ver.) -> structure: parameter "report_name" of String,
           parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_kb_gtdbtk_classify_wf

        ### Step 00: Init
        params = get_gtdbtk_params(params)
        console = []
        method_name = 'run_kb_gtdbtk_classify_wf'
        self.log(console, 'Running ' + method_name + ' with params='
)
        self.log(console, "\n" + pformat(params))
        
        
        self.log(console, "Get Genome Seqs\n")
        fasta_path = self.shared_folder / 'fastas'
        fasta_path.mkdir(parents=True, exist_ok=True)

        cli = KBClients(self.callback_url, self.ws_url, self.hs_url, ctx['token'])

        path_to_filename = download_sequence(params.ref, fasta_path, cli)
        for path, fn in path_to_filename.items():
            print(fn, path)

        output_path = self.shared_folder / 'output'
        temp_output = self.shared_folder / 'temp_output'
        output_path.mkdir(parents=True, exist_ok=True)
        temp_output.mkdir(parents=True, exist_ok=True)

        
        ### Step 01: run GTDB-Tk Classify WF
        def runner(args):
            self.log(console, "Run gtdbtk classify_wf\n")

            env = dict(os.environ)
            env['TEMP_DIR'] = str(self.shared_folder / 'tmp')
            # should print to stdout/stderr
            subprocess.run(args, check=True, env=env)

        (classification, summary_tables) = run_gtdbtk (runner,
                                                       path_to_filename,
                                                       output_path,
                                                       temp_output,
                                                       params.min_perc_aa,
                                                       params.db_ver,
                                                       params.keep_intermediates,
                                                       self.cpus)


        ### Step 02: Make Krona plot
        self.log(console, "Format Krona plot")
        run_krona_import_text(runner, output_path, temp_output)

        
        ### Step 03: Save Genome and/or Assembly objects with updated lineage
        objects_created = None
        top_query_obj_type = get_obj_type (params.ref, cli)
        if check_obj_type_assembly (top_query_obj_type) or check_obj_type_genome (top_query_obj_type):
            self.log(console, "Update Genome and Assembly objects and lineage files")
            if params.db_ver == 207:
                taxon_assignment_field = 'GTDB_R07-RS207'
            else:
                taxon_assignment_field = 'GTDB_R08-RS214'
            objects_created = update_genome_assembly_objs_class (params.workspace_id,
                                                                 params.ref,
                                                                 classification,
                                                                 params.overwrite_tax,
                                                                 str(params.db_ver),
                                                                 taxon_assignment_field,
                                                                 cli)

        
        ### Step 04: copy over GTDB Species Rep Genomes to calling WS and make GenomeSets
        if params.copy_proximals and check_obj_type_genome (top_query_obj_type):
            self.log(console, "Create Proximal GenomeSets and copy Species Representative Genomes")
            objects_created.extend (copy_gtdb_species_reps (params.workspace_id,
                                                            params.ref,
                                                            self.genome_upas_map_file,
                                                            summary_tables,
                                                            cli))
        

        ### Step 05: process trees
        self.log(console, "Process Trees")
        file_links = process_tree_files (params.ref,
                                         output_path,
                                         summary_tables,
                                         classification,
                                         str(params.db_ver),
                                         params.dendrogram_report,
                                         cli)


        ### Step 06: copy tree genomes and save tree object
        if params.save_trees and check_obj_type_genome (top_query_obj_type):
            self.log(console, "Save Tree object and copy Species Representative Genomes")
            objects_created.extend (save_gtdb_tree_objs (params.workspace_id,
                                                         params.ref,
                                                         output_path,
                                                         params.output_tree_basename,
                                                         self.genome_upas_map_file,
                                                         cli))
        
        
        ### Step 07: make report
        self.log(console, "Generate Report")
        output = generate_report(cli,
                                 output_path,
                                 params.workspace_id,
                                 objects_created,
                                 file_links)

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
