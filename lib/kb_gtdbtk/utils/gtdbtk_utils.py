import os
import subprocess
import json
import logging
import pandas as pd
import json
import tempfile

from .misc_utils import mkdir_p


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

    def gtdbtk_classifywf(self, output_path, min_perc_aa, id_to_obj_info):
        '''
        Run the classify workflow on the fasta files
        '''
        # TODO only use upas for data, update output to include names
        with tempfile.NamedTemporaryFile(
                mode='w',
                prefix='gtdb_tk_file_input_',
                suffix='.tmp',
                delete=False,
                dir=self.shared_folder) as tf:
            for id_, val in id_to_obj_info.items():
                tf.write(val['path'] + '\t' + id_ + '\n')

        gtdbtk_cmd = " ".join(
            [self.gtdbtk,
             "classify_wf",
             "--out_dir", output_path,
             "--batchfile", tf.name,
             "--cpus", str(self.cpus), 
             "--min_perc_aa", str(min_perc_aa), '"'])
        logging.info("Starting Command:\n" + gtdbtk_cmd)

        env = dict(os.environ)
        env['TMPDIR'] = os.path.join(self.shared_folder, 'tmp')
        mkdir_p(env['TMPDIR'])
        # should figure a way of getting this to run without shell=True, security risk
        # https://docs.python.org/3.7/library/subprocess.html#security-considerations
        output = subprocess.check_output(gtdbtk_cmd, shell=True, env=env).decode('utf-8')
        logging.info(output)

        self._process_output_files(output_path, id_to_obj_info)
        return output
    
    def _process_output_files(self, out_dir, id_to_obj_info):

        for path in (os.path.join(out_dir, 'gtdbtk.ar122.summary.tsv'),
                     os.path.join(out_dir, 'gtdbtk.bac120.summary.tsv'),
                     os.path.join(out_dir, 'gtdbtk.bac120.markers_summary.tsv'),
                     os.path.join(out_dir, 'gtdbtk.ar122.markers_summary.tsv'),
                     os.path.join(out_dir, 'gtdbtk.filtered.tsv')):
            try:
                summary_df = pd.read_csv(path, sep='\t', encoding='utf-8')
                outfile = path + '.json'
                summary_json = '{"data": ' + summary_df.to_json(orient='records') + '}'
                sj = json.loads(summary_json)
                for item in sj['data']:
                    item['Name'] = id_to_obj_info[item['Name']]['assembly_name']

                with open(outfile, 'w') as out:
                    out.write(json.dumps(sj))
            except Exception as exc:
                # should throw an exception rather than continuing
                # also some exceptions are expected depending on whether bac or arch
                logging.info(exc)

        return
