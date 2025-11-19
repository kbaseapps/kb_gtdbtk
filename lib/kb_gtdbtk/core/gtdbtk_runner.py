'''
Run GTDB-tk against a set of sequence files.
'''

from dataclasses import dataclass
import logging
import json
import os
import shutil
import pandas as pd
import tempfile

from datetime import datetime
from pathlib import Path
from shutil import (
    copyfile,
    copytree,
    rmtree
)
from typing import (
    Callable,
    Dict,
    List,
    Optional,
    Tuple
)
from kb_gtdbtk.core.string_util import now_ISOish

def get_mash_db_path(root_dir: Path, db_ver: int) -> Path:
    """
    Returns the expected path to the MASH database, based on the GTDB-tk database version
    and the root refdata path. If this DB doesn't exist, this raises a RuntimeError.
    """
    # refdata mounted mash db.  Must be generated during docker image registration init as /data is read-only at app runtime
    mash_db_path = root_dir / f"r{db_ver}" / "mash" / "gtdb_ref_sketch.msh"
    if not mash_db_path.exists():
        raise RuntimeError(f"GTDB ref genomes MASH DB not found in expected path {mash_db_path}. This must be generated during refdata initialization.")
    return mash_db_path


# main func
def run_gtdbtk(
        gtdbtk_runner: Callable[[List[str]], None],
        sequences: Dict[Path, str],
        output_dir: Path,
        temp_dir: Path,
        min_perc_aa: float,
        db_ver: int,
        keep_intermediates: int,
        cpus: int,
        data_root_dir: Optional[Path]=None) -> Tuple[dict, dict]:
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
    :param db_ver: The version of the GTDB-tk reference data.
    :param keep_intermediates: If 1, passes the --keep_intermediates flag to gtdbtk
    :param cpus: the number of CPUs GTDB-tk should use.
    :param data_root_dir: The filesystem path where ref data is stored. This will typically be
        /data , but might be somewhere else for testing with mocked files.
    '''
    # TODO input checking
    # TODO test logging, need to install an interceptor. Tested manually for now

    # all this complication is due to GTDB-tk choking on many legal, but uncommon, file name
    # characters such as |. Essentially here we provide safe file names and identifiers
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
            id_to_name[id_] = str(sequences[path]).replace('.gz','')
            os.symlink(path, temp_links / id_)
            tf.write(str(temp_links / id_) + '\t' + id_ + '\n')

    if data_root_dir is None:
        data_root_dir = Path("/data")

    # set refdata location
    os.environ['GTDBTK_DATA_PATH'] = str(data_root_dir / f"r{db_ver}")

    # set output dirs
    # temp_output = main gtdbtk output directory
    # temp_trees_output = only used if a second gtdbtk run is done with --skip_ani_screen
    temp_output = temp_dir / 'output' / timestamp
    temp_trees_output = temp_dir / 'output_trees' / timestamp
    temp_output.mkdir(parents=True, exist_ok=True)

    gtdbtk_cmd = [
        'gtdbtk',
        'classify_wf',
        '--out_dir', str(temp_output),
        '--batchfile', tf.name,
        '--cpus', str(cpus),
        '--min_perc_aa', str(min_perc_aa)
    ]
    if keep_intermediates == 1:
        gtdbtk_cmd += ['--keep_intermediates']

    # refdata mounted mash db.  Must be generated during docker image registration init as /data is read-only at app runtime
    mash_db_path = get_mash_db_path(data_root_dir, db_ver)
    gtdbtk_cmd += ['--mash_db', str(mash_db_path)]

    # run first pass
    logging.info('Starting Command:\n' + ' '.join(gtdbtk_cmd))
    gtdbtk_runner(gtdbtk_cmd)

    # Not all queries may be placed into trees (ANI step may filter)
    if not _all_ids_in_trees (temp_output, id_to_name):
        logging.info('Not all queries placed in trees.  Running second pass with --skip_ani_screen ...')
        temp_trees_output.mkdir(parents=True, exist_ok=True)

        gtdbtk_cmd = [
            'gtdbtk',
            'classify_wf',
            '--out_dir', str(temp_trees_output),
            '--batchfile', tf.name,
            '--cpus', str(cpus),
            '--min_perc_aa', str(min_perc_aa),
            '--skip_ani_screen',
            '--no_mash'
        ]
        if keep_intermediates == 1:
            gtdbtk_cmd += ['--keep_intermediates']
        # run first pass
        logging.info('Starting Command:\n' + ' '.join(gtdbtk_cmd))
        gtdbtk_runner(gtdbtk_cmd)

    return _process_output_files(temp_output, temp_trees_output, output_dir, id_to_name)


def _all_ids_in_trees(temp_output: Path, id_to_name: Dict[str, str]) -> bool:
    """
    Looks at the primary gtdbtk output file - either gtdbtk.ar53.summary.tsv or
    gtdbtk.bac120.summary.tsv (under the temp_output directory). These files have
    20 columns. First is the user_genome id, and 13th is the taxonomy of the tree it was
    placed in. This checks the following:
    1. if all user_genomes exist in the output.
    2. if they're all mapped to trees - that pplacer_taxonomy is not "N/A"
    Returns True if all genomes submitted appear with trees, False otherwise.
    """
    all_ids_found = True

    ids_found = dict()
    for file_ in ['gtdbtk.ar53.summary.tsv', 'gtdbtk.bac120.summary.tsv']:
        summary_path = temp_output / 'classify' / file_
        if not summary_path.is_file():
            continue
        else:
            with open (summary_path, 'r') as summary_h:
                for summary_line in summary_h:
                    if summary_line.startswith('user_genome'):
                        continue
                    [user_genome, classification, fastani_reference, fastani_reference_radius, fastani_taxonomy, fastani_ani, fastani_af, closest_placement_reference, closest_placement_radius, closest_placement_taxonomy, closest_placement_ani, closest_placement_af, pplacer_taxonomy, classification_method, note, other_related_references, msa_percent, translation_table, red_value, warnings] = summary_line.rstrip().split("\t")
                    if pplacer_taxonomy != 'N/A':
                        ids_found[user_genome] = True

    for qid in list(id_to_name.keys()):
        if qid not in ids_found:
            all_ids_found = False
            break

    return all_ids_found


def _process_output_files(
        temp_output: Path,
        temp_trees_output: Path,
        out_dir: Path,
        id_to_name: Dict[str, str]) -> Tuple[dict, dict]:
    """
    Process GTDB-tk output files and consolidate results into a structured format.

    This function handles the post-processing of GTDB-tk classification workflow output.
    It consolidates results from both a primary run and an optional secondary run
    (performed with --skip_ani_screen), merges summary data, converts TSV files to JSON
    format, and remaps internal sequence identifiers to their original assembly names.

    Key operations:
    1. Copies the entire temp_output directory structure to the output directory
    2. Saves the mapping between internal IDs and original assembly names
    3. Copies tree files, preferring those from the secondary run if available
    4. Merges summary TSV files from both runs, with secondary run data taking
       precedence for N/A fields
    5. Converts summary TSV files to JSON format for web UI consumption
    6. Remaps internal IDs back to original assembly names in output data
    7. Filters out blank/null fields by replacing them with '-' for consistency

    :param temp_output: Path to the primary GTDB-tk output directory. Contains
        results from the initial classify_wf run with all markers and ANI screening.
        Expected subdirectories: 'classify' and 'identify'
    :param temp_trees_output: Path to the secondary GTDB-tk output directory. Contains
        results from the secondary classify_wf run (if performed) with --skip_ani_screen
        flag. May be empty if all sequences were placed in trees during the primary run.
        Expected subdirectories: 'classify' and 'identify'
    :param out_dir: Path to the output directory where final results will be written.
        Should be an existing, writable directory. Output includes:
        - 'runtime_output/': Copy of temp_output directory structure
        - 'id_to_name.map': Tab-separated file mapping internal IDs to assembly names
        - Tree files (*.tree): Phylogenetic trees from both runs
        - Summary TSV files: Classification and marker summary data
        - Summary JSON files: JSON-converted versions of summary TSV files
    :param id_to_name: Mapping from internal sequence identifiers (e.g., 'id0', 'id1')
        to original assembly names. Used to remap internal IDs in output data back to
        user-supplied assembly names.

    :return: A tuple of (classification, summary_tables) where:
        - classification (dict): Mapping from assembly name to classification string.
          Keys are original assembly names (from id_to_name values).
          Values are classification strings from the 'classification' field in summary files.
        - summary_tables (dict): Mapping from filename to parsed JSON data.
          Keys are TSV filenames (e.g., 'gtdbtk.bac120.summary.tsv').
          Values are dictionaries with a 'data' key containing a list of record objects
          with assembly names and classifications.
    """

    classification = dict()
    summary_tables = dict()

    # copy over all created output
    sub_out_dir = Path(out_dir / 'runtime_output')
    if os.path.isdir(sub_out_dir):  # should only occur during unit tests
        rmtree(sub_out_dir)
    copytree(temp_output, sub_out_dir, symlinks=True)

    # save id to name mapping as a file
    id_map_buf = [f"{id_}\t{id_to_name[id_]}" for id_ in sorted(id_to_name)]
    with open (out_dir / "id_to_name.map", 'w') as file_h:
        file_h.write("\n".join(id_map_buf)+"\n")

    # make json files for html tables
    tree_files = [
        'gtdbtk.ar53.classify.tree',
        'gtdbtk.bac120.classify.tree'
    ]
    bb_tree_file = ['gtdbtk.backbone.bac120.classify.tree']
    extra_bac_tree_files = ['gtdbtk.backbone.bac120.classify.tree']

    base_files = [
        'gtdbtk.ar53.summary.tsv',
        'gtdbtk.bac120.summary.tsv',
        'gtdbtk.ar53.markers_summary.tsv',
        'gtdbtk.bac120.markers_summary.tsv',
        'gtdbtk.bac120.tree.mapping.tsv'
    ]
    # These are the output files we care about - they're either in the classify/ or identify/
    # subdirectories. We'll update these below based on file structure (and whether or not
    # there was a --skip_ani_screen run)
    file_folder = {
        'gtdbtk.ar53.summary.tsv': 'classify',
        'gtdbtk.bac120.summary.tsv': 'classify',
        'gtdbtk.ar53.markers_summary.tsv': 'identify',
        'gtdbtk.bac120.markers_summary.tsv': 'identify',
        'gtdbtk.bac120.tree.mapping.tsv': 'classify',
        'gtdbtk.ar53.classify.tree': 'classify',
        'gtdbtk.bac120.classify.tree': 'classify',
        'gtdbtk.backbone.bac120.classify.tree': 'classify'
    }

    # TODO: rewrite to sift through available files on the filesystem instead of
    # assuming an upper limit
    # make up a set of possible extra subtree files.
    for i in range(10000):
        subtree_file = 'gtdbtk.bac120.classify.tree.'+str(i)+'.tree'
        extra_bac_tree_files.append(subtree_file)
        file_folder[subtree_file] = 'classify'

    # copy all tree files to the output directory
    # these may be in the temp_trees_output or the temp_output directory
    # prefer the temp_trees_output directory
    for file_ in tree_files + bb_tree_file + extra_bac_tree_files:
        treepath = temp_trees_output / file_folder[file_] / file_
        tmppath = temp_output / file_folder[file_] / file_
        out_path = out_dir / file_
        if treepath.is_file():
            copyfile(treepath, out_path)
        elif tmppath.is_file():
            copyfile(tmppath, out_path)

    # merge summary tsv files
    for file_ in base_files:
        treepath = temp_trees_output / file_folder[file_] / file_
        tmppath = temp_output / file_folder[file_] / file_
        out_path = out_dir / file_
        _merge_summary_tsv_files(tmppath, treepath, out_path)

    # load results
    for file_ in base_files:
                  # skip filtered for now, unused
                  # 'gtdbtk.filtered.tsv'
        path = out_dir / file_
        if not path.is_file():
            #logging.info('No such file, skipping: ' + str(tmppath))
            continue
        else:
            if not file_.endswith('summary.tsv'):
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

def _merge_summary_tsv_files(std_path: Path, tree_path: Path, output_path: Path):
    """
    Merges two TSV files from separate runs of GTDBtk.
    The std_path represents the standard run, using a mash DB.
    The tree_path is from a run using --skip_ani_screen, which generates file in a different path.
    These two TSV files have slightly different outputs, with data in tree_path being generally more
    informative. So fields for each genome in the summary file should have tree_path override any N/A values
    found in std_path. E.g.:

    std_path:
    user_genome    field1    field2    field3
    some_genome    N/A       N/A       X

    tree_path:
    user_genome    field1    field2    field3
    some_genome    Y         N/A       X

    merged:
    user_genome    field1    field2    field3
    some_genome    Y         N/A       X
    """

    std_file_exists = std_path.is_file()
    tree_file_exists = tree_path.is_file()
    if not std_file_exists and not tree_file_exists:
        return
    if std_file_exists and not tree_file_exists:
        shutil.copy(std_path, output_path)
        return
    if tree_file_exists and not std_file_exists:
        shutil.copy(tree_path, output_path)
        return
    std_file_data = _load_summary_tsv_file(std_path)
    tree_file_data = _load_summary_tsv_file(tree_path)
    if std_file_data["header"] != tree_file_data["header"]:
        logging.warning(f"While merging summary files: header of {std_path} (from using a mash DB) does not match {tree_path} (skipping the ANI screen). These cannot be merged, saving only the ANI-skipped file.")
        shutil.copy(tree_path, output_path)
        return
    out_buf = ["\t".join(std_file_data["header"])]
    num_cols = len(std_file_data["header"])
    for qid in std_file_data["id_order"]:
        row = ["N/A"] * num_cols
        if qid in std_file_data["data"]:
            row = std_file_data["data"][qid]
        if qid in tree_file_data["data"]:
            for field_i, value in enumerate(tree_file_data["data"][qid]):
                if row[field_i] == "N/A":
                    row[field_i] = value
        out_buf.append("\t".join(row))

    # write merged summaries
    with open(output_path, "w") as summary_file:
        summary_file.write("\n".join(out_buf)+"\n")


def _load_summary_tsv_file(filepath: Path) -> dict[str, list|dict]:
    """
    Loads a TSV file, includes the header (assumes that each summary file from GTDBtk has
    a header) and all lines.
    Summary files all have a unique identifier (generally the genome or assembly name)
    as the first element of each line. This, then, creates the following structure:
    {
        "header": [],
        "data": {"id": [row]},
        "id_order": []
    }
    """
    header = []
    data = {}
    id_order = []
    with open(filepath, "r") as infile:
        for line_num, line in enumerate(infile):
            row = line.rstrip("\n\r").split("\t")
            if line_num == 0:
                header = row
            else:
                identifier = row[0]
                data[identifier] = row
                id_order.append(identifier)

    return {
        "header": header,
        "data": data,
        "id_order": id_order
    }
