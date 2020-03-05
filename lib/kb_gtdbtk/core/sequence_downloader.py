'''
Downloads sequence data from KBase services.
'''

import os
from typing import Dict
from pathlib import Path
from shutil import copyfile

from kb_gtdbtk.core.kb_client_set import KBClients

# DEV NOTES: This is tested by mocks, and so type changes will not be detected by the tests.
# Type changes should be tested manually.


def download_sequence(
        upa: str,
        destination_dir: Path,
        clients: KBClients,
        ) -> Dict[Path, str]:
    '''
    Download sequence data from KBase.

    Handles types:
    KBaseSets.GenomeSet
    KBaseSearch.GenomeSet
    KBaseGenomes.Genome
    KBaseGenomes.ContigSet
    KBaseGenomeAnnotations.Assembly
    KBaseSets.AssemblySet
    KBaseMetagenomes.BinnedContigs


    :param upa: The KBase UPA (e.g. X/Y/Z format) for the object from which sequence data will be
        downloaded.
    :param destination_dir: Where to store the files; must be an extant directory.
    :param clients: The KBase clients to use for the download operation.
    :returns: A mapping from a file path to a display name for the file (often, but not always,
        the filename).
    '''
    # TODO check input args
    # TODO this code should be folded into AssemblyFileUtils. Similar code exists there but
    # apparently it has permission issues for reference paths.
    # TODO may need to have type version specific code.

    # TODO log server side stack trace for exceptions from clients
    obj_data = clients.dfu().get_objects({'object_refs': [upa]})['data'][0]
    # normalize upa just in case it's a ref vs. an upa
    upa = str(obj_data['info'][6]) + '/' + str(obj_data['info'][0]) + '/' + str(
        obj_data['info'][4])
    obj_type = obj_data['info'][2].split('-')[0]

    dd = destination_dir

    if 'KBaseSets.GenomeSet' == obj_type:
        upas = [gsi['ref'] for gsi in obj_data['data']['items']]
    elif 'KBaseSearch.GenomeSet' == obj_type:
        upas = [gse['ref'] for gse in obj_data['data']['elements'].values()]
    elif 'KBaseGenomes.Genome' == obj_type:
        upas = [upa]
    elif 'KBaseGenomes.ContigSet' == obj_type or 'KBaseGenomeAnnotations.Assembly' == obj_type:
        # in this case we use the assembly file util to get the fasta file
        # didn't need to pull the full object above, could just pull the info to get the type
        # could optimize out
        faf = clients.au().get_assembly_as_fasta(
            {'ref': upa, 'filename': str(_upa_to_path(dd, upa))})
        return {Path(faf['path']): faf['assembly_name']}
    elif 'KBaseSets.AssemblySet' == obj_type:
        id_to_assy_info = {}
        for item_upa in obj_data['data']['items']:
            faf = clients.au().get_assembly_as_fasta(
                {'ref': upa + ';' + item_upa['ref'],
                 'filename': str(_upa_to_path(dd, item_upa['ref']))})
            id_to_assy_info[Path(faf['path'])] = faf['assembly_name']
        return id_to_assy_info
    elif 'KBaseMetagenomes.BinnedContigs' == obj_type:
        return _handle_binned_contigs(upa, clients, dd)
    else:
        raise ValueError(f'{obj_type} type is not supported')

    return _process_genomes(upa, upas, clients, dd)


def _process_genomes(upa, upas, clients, dd):
    id_to_assy_info = {}
    for genome_upa in upas:
        # this could be sped up by batching the get_objects call
        # does assy file util not take bulk calls?
        # maybe doesn't matter since Shock doesn't handle bulk calls
        # this could be optimized in other ways, e.g. putting the assembly ref in the object
        # metadata so the whole object doesn't need to be fetched
        if upa != genome_upa:  # for single genomes, upa and genome_upa will be the same
            genome_upa = upa + ';' + genome_upa
        # why not DFU here?
        # for a single genome this is pulling the same data again. Could optimize here.
        genome_data = clients.ws().get_objects2(
            {'objects': [{'ref': genome_upa}]})['data'][0]['data']
        target_upa = genome_data.get('contigset_ref') or genome_data.get('assembly_ref')
        faf = clients.au().get_assembly_as_fasta({
            'ref': genome_upa + ';' + target_upa,
            'filename': str(_upa_to_path(dd, target_upa))
            })
        id_to_assy_info[Path(faf['path'])] = faf['assembly_name']

    return id_to_assy_info


def _handle_binned_contigs(upa, clients, target_dir):
    # any CoaC issues here are in MetagenomeUtils
    ret = {}
    bin_file_dir = clients.mgu().binned_contigs_to_file(
        {'input_ref': upa, 'save_to_shock': 0})['bin_file_directory']
    for (dirpath, dirnames, filenames) in os.walk(bin_file_dir):
        for fasta_file in filenames:
            # not sure why we're normalizing the extension, but that's how the old code did it
            # *shrug*
            fasta_fixed_ext = os.path.splitext(fasta_file)[0] + '.fa'
            fasta_path = target_dir / fasta_fixed_ext
            copyfile(os.path.join(bin_file_dir, fasta_file), fasta_path)
            # Should I verify that the bins have contigs?
            # is it possible to have empty bins?
            ret[Path(fasta_path)] = fasta_fixed_ext
        break  # not sure why this is necessary
    return ret


def _upa_to_path(target_dir: Path, upa) -> Path:
    return target_dir / _file_safe_upa(upa)


def _file_safe_upa(upa):
    return upa.replace('/', '_')
