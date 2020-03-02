'''
Run GTDB-tk against a set of sequence files.
'''

import logging
import json
import pandas as pd
import tempfile

from pathlib import Path
from shutil import copyfile
from typing import Dict, List, Callable


def run_gtdbtk(
        gtdbtk_runner: Callable[[List[str]], None],
        sequences: Dict[str, Dict[str, str]],
        output_dir: Path,
        temp_dir: Path,
        min_perc_aa: float,
        cpus: int) -> None:
    '''
    Run GTDB-tk on a set of sequences in FASTA format. Expects the 'gtdbtk' command to be on the
    system path.

    Any temporary files generated are not deleted unless the 3rd party GTDB-tk code deletes them.

    :param gtdbtk_runner: a callable that takes a list of arguments for GTDB-tk and
        excutes the program with those arguments.
    :param sequences: Information about the fasta files. A mapping from a unique, filename safe
        ID for each file to another mapping containing the keys 'path', which maps to the path
        to the FASTA file, and 'assembly_name', which is the name of the assembly to use in
        the ouput.
    :param output_dir: an extant directory in which to place the output. The output is JSON
        versions of GTDB-tk's TSV output using the 'assembly_name' from the sequences map
        as the sequence name.
    :param temp_dir: an extant temporary directory to use for processing. Any files or
        directories in this directory may be deleted or overwritten.
    :param min_perc_aa: The mimimum sequence alignment in percent.
    :param cpus: the number of CPUs GTDB-tk should use.
    '''
    # TODO input checking
    # TODO test logging, need to install an interceptor. Tested manually for now
    # TODO disallow numbers in string form as IDs - pandas will load it as a number -> KeyError
    # TODO FOR REAL abstract the unique ID part here rather than in the downloader

    # all this complication is due to GTDB-tk choking on many legal, but uncommon, file name
    # characters such as |. Essentially here we provide safe file names and identitifers
    # (which GTDB-tk will use to create temporary files) and then remap to the original,
    # potentially unsafe names.
    with tempfile.NamedTemporaryFile(
            mode='w',
            prefix='gtdb_tk_file_input_',
            suffix='.tmp',
            delete=False,
            dir=temp_dir) as tf:
        for id_ in sorted(sequences):
            tf.write(sequences[id_]['path'] + '\t' + id_ + '\n')

    temp_output = temp_dir / 'output'
    temp_output.mkdir(parents=True, exist_ok=True)

    gtdbtk_cmd = [
        'gtdbtk',
        'classify_wf',
        '--out_dir', str(temp_output),
        '--batchfile', tf.name,
        '--cpus', str(cpus),
        '--min_perc_aa', str(min_perc_aa)]

    logging.info('Starting Command:\n' + ' '.join(gtdbtk_cmd))
    gtdbtk_runner(gtdbtk_cmd)
    _process_output_files(temp_output, output_dir, sequences)


def _process_output_files(temp_output, out_dir, sequences):

    for file_ in ('gtdbtk.ar122.summary.tsv',
                  'gtdbtk.bac120.summary.tsv',
                  'gtdbtk.bac120.markers_summary.tsv',
                  'gtdbtk.ar122.markers_summary.tsv'
                  # skip filtered for now, unused
                  # 'gtdbtk.filtered.tsv'
                  ):
        tmppath = temp_output / file_
        if not tmppath.is_file():
            logging.info('No such file, skipping: ' + str(tmppath))
        else:
            path = out_dir / file_
            copyfile(tmppath, path)
            summary_df = pd.read_csv(path, sep='\t', encoding='utf-8')
            outfile = str(path) + '.json'
            summary_json = '{"data": ' + summary_df.to_json(orient='records') + '}'
            sj = json.loads(summary_json)
            for item in sj['data']:
                key = 'Name' if 'Name' in item else 'user_genome'
                item[key] = sequences[item[key]]['assembly_name']

            with open(outfile, 'w') as out:
                out.write(json.dumps(sj))

    return
