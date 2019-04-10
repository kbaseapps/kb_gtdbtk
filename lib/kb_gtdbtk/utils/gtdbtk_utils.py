import os
import subprocess


class GTDBTkUtils():
    '''
    Utilities for running GTDB-Tk
    '''

    def __init__(self, config, callback_url, workspace_id):
        self.shared_folder = config['scratch']
        self.callback_url = callback_url
        self.gtdbtk = '/bin/bash -c "source activate py2 && gtdbtk'

    def gtdbtk_classifywf(self, fasta_paths):
        '''
        Run the classify workflow on the fasta files
        '''
        out_dir = os.path.join(self.shared_folder, "output")
        gtdbtk_cmd = " ".join([self.gtdbtk, "classifywf", "--out_dir", out_dir, "--genome_dir" fasta_paths, '"'])
        print("Starting Command:\n", gtdbtk_cmd)
        output = subprocess.run(gtdbtk_cmd, stdout=subprocess.PIPE, shell=True)
        print(output)