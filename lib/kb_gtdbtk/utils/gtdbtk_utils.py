import os
import subprocess


class GTDBTkUtils():
    '''
    Utilities for running GTDB-Tk
    '''

    def __init__(self, config, callback_url, workspace_id, cpus):
        self.shared_folder = config['scratch']
        self.callback_url = callback_url
        self.cpus = cpus
        self.gtdbtk = '/bin/bash -c "source activate py2 && GTDBTK_DATA_PATH=/data gtdbtk'

    def gtdbtk_classifywf(self, fasta_paths):
        '''
        Run the classify workflow on the fasta files
        '''
        out_dir = os.path.join(self.shared_folder, "output")
        gtdbtk_cmd = " ".join([self.gtdbtk, "classify_wf", "--out_dir", out_dir,
                              "--genome_dir", self.shared_folder, "-x", "fa",
                               "--cpus", str(self.cpus), '"'])
        print("Starting Command:\n", gtdbtk_cmd)
        output = subprocess.check_output(gtdbtk_cmd, shell=True).decode('utf-8')
        print(output)

        for path in ('gtdbtk.ar122.summary.tsv', 'gtdbtk.bact120.summary.tsv'):
            try:
                summary_file = open(path, 'r')
                output.append(summary_file.read())
                summary_file.close()
            except:
                pass

        return output
