'''
Run GTDB-tk against a set of sequence files.
'''

import logging
import json
import os
import pandas as pd
import tempfile

from datetime import datetime
from pathlib import Path
from shutil import copyfile,copytree,rmtree
from typing import Dict, List, Callable


# timestamp
def now_ISOish():
    now_timestamp = datetime.now()
    now_secs_from_epoch = (now_timestamp - datetime(1970,1,1)).total_seconds()
    #now_timestamp_in_iso = datetime.fromtimestamp(int(now_secs_from_epoch)).strftime('%Y-%m-%d_%T')
    now_timestamp_in_isoish = datetime.fromtimestamp(int(now_secs_from_epoch)).strftime('%Y%m%d_%H%M%S')
    return now_timestamp_in_isoish


# main func
def run_gtdbtk(
        gtdbtk_runner: Callable[[List[str]], None],
        sequences: Dict[Path, str],
        output_dir: Path,
        temp_dir: Path,
        min_perc_aa: float,
        db_ver: int,
        keep_intermediates: int,
        cpus: int) -> None:
    '''
    Run GTDB-tk on a set of sequences in FASTA format. Expects the 'gtdbtk' command to be on the
    system path.

    Any temporary files generated are not deleted unless the 3rd party GTDB-tk code deletes them.

    :param gtdbtk_runner: a callable that takes a list of arguments for GTDB-tk and
        excutes the program with those arguments.
    :param sequences: Information about the fasta files. A mapping from a file path to a display
        name for the file, often simply the file name.
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

    # all this complication is due to GTDB-tk choking on many legal, but uncommon, file name
    # characters such as |. Essentially here we provide safe file names and identitifers
    # (which GTDB-tk will use to create temporary files) and then remap to the original,
    # potentially unsafe names.
    timestamp = now_ISOish()
    temp_links = temp_dir / 'links' / timestamp
    temp_links.mkdir(parents=True, exist_ok=True)
    id_to_name = {}
    with tempfile.NamedTemporaryFile(
            mode='w',
            prefix='gtdb_tk_file_input_',
            suffix='.tmp',
            delete=False,
            dir=temp_dir) as tf:
        for i, path in enumerate(sorted(sequences)):
            id_ = f'id{i}'
            id_to_name[id_] = sequences[path]
            os.symlink(path, temp_links / id_)
            tf.write(str(temp_links / id_) + '\t' + id_ + '\n')

    temp_output = temp_dir / 'output' / timestamp
    temp_output.mkdir(parents=True, exist_ok=True)

    gtdbtk_cmd = [
        'gtdbtk',
        'classify_wf',
        '--out_dir', str(temp_output),
        '--batchfile', tf.name,
        '--cpus', str(cpus),
        '--min_perc_aa', str(min_perc_aa)]

    if db_ver == 207:
        gtdbtk_cmd += ['--skip_ani_screen']
    else:
        gtdbtk_cmd += ['--mash_db', '/data/r'+str(db_ver)+'/mash.db']
        
    if keep_intermediates == 1:
        gtdbtk_cmd += ['--keep_intermediates']

    # set refdata location
    os.environ['GTDBTK_DATA_PATH'] = '/data/'+'r'+str(db_ver)
        
    logging.info('Starting Command:\n' + ' '.join(gtdbtk_cmd))
    gtdbtk_runner(gtdbtk_cmd)

    return _process_output_files(temp_output, output_dir, id_to_name)
    

def _process_output_files(temp_output, out_dir, id_to_name):

    classification = dict()
    summary_tables = dict()
    trimmed_tree_files = dict()
    
    # copy over all created output
    """
    for file_ in os.listdir (temp_output):
        tmppath = temp_output / file_
        if not tmppath.is_file():
            continue
        path = out_dir / file_
        copyfile(tmppath, path)        
    """
    sub_out_dir = Path(out_dir / 'runtime_output')
    if os.path.isdir(sub_out_dir):  # only occurs during unit tests
        rmtree(sub_out_dir)
    copytree(temp_output, sub_out_dir, symlinks=True)

    # save id to name mapping as a file
    id_map_buf = []
    for id_ in sorted(id_to_name.keys()):
        id_map_buf.append("\t".join([id_,id_to_name[id_]]))
    id_map_path = os.path.join(out_dir, 'id_to_name.map')
    with open (id_map_path, 'w') as file_h:
        file_h.write("\n".join(id_map_buf)+"\n")
    
    # make json files for html tables
    tree_files = ['gtdbtk.ar53.classify.tree',
                  'gtdbtk.bac120.classify.tree']
    bb_tree_file = ['gtdbtk.backbone.bac120.classify.tree']

    base_files = ['gtdbtk.ar53.summary.tsv',
                  'gtdbtk.bac120.summary.tsv',
                  'gtdbtk.ar53.markers_summary.tsv',
                  'gtdbtk.bac120.markers_summary.tsv']
    file_folder = {'gtdbtk.ar53.summary.tsv': 'classify',
                   'gtdbtk.bac120.summary.tsv': 'classify',
                   'gtdbtk.ar53.markers_summary.tsv': 'identify',
                   'gtdbtk.bac120.markers_summary.tsv': 'identify',
                   'gtdbtk.ar53.classify.tree': 'classify',
                   'gtdbtk.bac120.classify.tree': 'classify',
                   'gtdbtk.backbone.bac120.classify.tree': 'classify'
                  }
    extra_bac_tree_files = []
    for i in range(10000):
        subtree_file = 'gtdbtk.bac120.classify.tree.'+str(i)+'.tree'
        extra_bac_tree_files.append(subtree_file)
        file_folder[subtree_file] = 'classify'

    for file_ in tree_files + bb_tree_file + base_files + extra_bac_tree_files:
                  # skip filtered for now, unused
                  # 'gtdbtk.filtered.tsv'
        tmppath = temp_output / file_folder[file_] / file_
        if not tmppath.is_file():
            #logging.info('No such file, skipping: ' + str(tmppath))
            continue
        else:
            path = out_dir / file_
            copyfile(tmppath, path)
            if not file_.endswith('.tsv'):
                continue
            summary_df = pd.read_csv(path, sep='\t', encoding='utf-8')
            outfile = str(path) + '.json'
            summary_json = '{"data": ' + summary_df.to_json(orient='records') + '}'
            sj = json.loads(summary_json)
            for item in sj['data']:

                # no blank fields.  messes up datatables in index.html
                for key in item.keys():
                    #print ("RESULTS: k: '"+key+"' val: '"+str(item[key])+"'")  # DEBUG
                    if not item.get(key):
                        item[key] = '-'  # note: this resets data in sj

                # reset id to assembly name
                # Note: field 'Name' was changed to 'name'
                #key = 'Name' if 'Name' in item else 'user_genome'
                if 'name' in item:
                    key = 'name'
                elif 'user_genome' in item:
                    key = 'user_genome'
                else:
                    continue
                this_id = item[key]
                if this_id not in id_to_name:
                    raise ValueError ("missing "+this_id+" in id_to_name dict")
                item[key] = id_to_name[this_id]  # note: this resets data in sj

                # store classification by assembly name
                if 'classification' in item:
                    classification[id_to_name[this_id]] = item['classification']
                
            # rewrite with updated vals
            with open(outfile, 'w') as out:
                out.write(json.dumps(sj))

            # return rest of summary table data
            summary_tables[file_] = sj

    return (classification, summary_tables)
