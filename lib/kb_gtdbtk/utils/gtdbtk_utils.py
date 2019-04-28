import os
import subprocess
import logging
import pandas as pd
import json


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

    def gtdbtk_classifywf(self, fasta_paths):
        '''
        Run the classify workflow on the fasta files
        '''
        out_dir = os.path.join(self.shared_folder, "output")
        gtdbtk_cmd = " ".join([self.gtdbtk, "classify_wf", "--out_dir", out_dir,
                              "--genome_dir", self.shared_folder, "-x", "fa",
                               "--cpus", str(self.cpus), '"'])
        logging.info("Starting Command:\n" + gtdbtk_cmd)
        output = subprocess.check_output(gtdbtk_cmd, shell=True).decode('utf-8')
        logging.info(output)

        output = output + self._process_output_files(out_dir)
        return output
    
    def _process_output_files(self, out_dir):

        for path in (os.path.join(out_dir, 'gtdbtk.ar122.summary.tsv'),
                     os.path.join(out_dir, 'gtdbtk.bac120.summary.tsv')):
            try:
                summary_file = open(path, 'r')
                file_content = summary_file.read()
                summary_file.close()
                summary_df = pd.read_csv(path, sep='\t', encoding='utf-8')
                outfile = path + '.json'
                summary_json = '{"data": ' + summary_df.to_json(orient='records') + '}'
                with open(outfile, 'w') as out:
                    out.write(summary_json)
            except Exception as exc:
                logging.info(exc)

        return file_content
