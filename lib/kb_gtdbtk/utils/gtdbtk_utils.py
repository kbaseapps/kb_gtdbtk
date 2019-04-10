import os
import subprocess


class GTDBTkUtils():
    '''
    Utilities for running GTDB-Tk
    '''

    def __init__(self, config, callback_url, workspace_id):
        self.shared_folder = config['scratch']
        self.callback_url = callback_url
        self.gtdbtk = 'gtdbtk'

    def gtdbtk_classifywf(self, fasta_paths):
        '''
        Run the classify workflow on the fasta files
        '''
        gtdbtk_cmd = [self.gtdbtk, "test"]
        pipe = subprocess.Popen(gtdbtk_cmd, stdout=subprocess.PIPE, shell=True)
        output = pipe.communicate()[0]
        print(output)
