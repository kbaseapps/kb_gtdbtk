import os
import subprocess
import logging
import pandas as pd
import json
import tempfile

from from .misc_utils import mkdir_p


class GTDBTkUtils():
    '''
    Utilities for running GTDB-Tk
    '''

    def __init__(self, config, callback_url, workspace_id, cpus):
        self.shared_folder = config['scratch']
        self.callback_url = callback_url
        self.cpus = cpus
        self.gtdbtk = '/bin/bash -c "source activate py2 && GTDBTK_DATA_PATH=/data gtdbtk'
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)

    def gtdbtk_classifywf(self, output_path, min_perc_aa, upa_to_obj_info):
        '''
        Run the classify workflow on the fasta files
        '''
        with tempfile.NamedTemporaryFile(
                prefix='gtdb_tk_file_input',
                suffix='tmp',
                delete=False,
                dir=self.shared_folder) as tf:
            for val in upa_to_obj_info.values():
                tf.write(val['path'] + '\t' + val['assembly_name'] + '\n')

        gtdbtk_cmd = " ".join(
            [self.gtdbtk,
             "classify_wf",
             "--out_dir", output_path,
             "--batchfile", tf.name,
             "--cpus", str(self.cpus), 
             "--min_perc_aa", str(min_perc_aa), '"'])
        logging.info("Starting Command:\n" + gtdbtk_cmd)

        env = dict(os.environ)
        env['TMPDIR'] = os.paths.join(self.shared_folder, 'tmp')
        mkdir_p(env['TMPDIR'])
        # should figure a way of getting this to run without shell=True, security risk
        # https://docs.python.org/3.7/library/subprocess.html#security-considerations
        output = subprocess.check_output(gtdbtk_cmd, shell=True, env=env).decode('utf-8')
        logging.info(output)

        self._process_output_files(output_path)
        return output
    
    def _process_output_files(self, out_dir):

        for path in (os.path.join(out_dir, 'gtdbtk.ar122.summary.tsv'),
                     os.path.join(out_dir, 'gtdbtk.bac120.summary.tsv'),
                     os.path.join(out_dir, 'gtdbtk.bac120.markers_summary.tsv'),
                     os.path.join(out_dir, 'gtdbtk.ar122.markers_summary.tsv'),
                     os.path.join(out_dir, 'gtdbtk.filtered.tsv')):
            try:
                summary_df = pd.read_csv(path, sep='\t', encoding='utf-8')
                outfile = path + '.json'
                summary_json = '{"data": ' + summary_df.to_json(orient='records') + '}'
                with open(outfile, 'w') as out:
                    out.write(summary_json)
            except Exception as exc:
                # should throw an exception rather than continuing
                logging.info(exc)

        return
