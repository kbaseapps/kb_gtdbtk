'''
Downloads sequence data from KBase services.
'''

import os
from typing import Dict
from pathlib import Path
from shutil import copyfile

from kb_gtdbtk.core.kb_client_set import KBClients


def check_obj_type_genome(
        upa: str,
        clients: KBClients,
        ):
    '''
    Check to see if query is a Genome or GenomeSet

    Handles types:
    KBaseSets.GenomeSet
    KBaseSearch.GenomeSet
    KBaseGenomes.Genome

    :param upa: The KBase UPA (e.g. X/Y/Z format) for the object from which sequence data will be
        downloaded.
    :param clients: The KBase clients to use for the download operation.
    :returns: True or False
    '''
    [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I,
     WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple    

    obj_data = clients.dfu().get_objects({'object_refs': [upa]})['data'][0]
    # normalize upa just in case it's a ref vs. an upa
    upa = str(obj_data['info'][WSID_I]) + '/' + str(obj_data['info'][OBJID_I]) + '/' + str(
        obj_data['info'][VERSION_I])
    obj_type = obj_data['info'][TYPE_I].split('-')[0]

    if 'KBaseSets.GenomeSet' == obj_type:
        return True
    elif 'KBaseSearch.GenomeSet' == obj_type:
        return True
    elif 'KBaseGenomes.Genome' == obj_type:
        return True
    else:
        return False


def update_genome_objs_class(
        upa: str,
        classification: Dict[str, str],
        clients: KBClients,
        ):
    '''
    Adjust taxon_assignments field of Genome objects for Genome or GenomeSet

    Handles types:
    KBaseSets.GenomeSet
    KBaseSearch.GenomeSet
    KBaseGenomes.Genome

    :param upa: The KBase UPA (e.g. X/Y/Z format) for the object from which sequence data will be
        downloaded.
    :param classification: a dict with name key to GTDB classification
    :param clients: The KBase clients to use for the download operation.
    :returns: null
    '''
    [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I,
     WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple    

    obj_data = clients.dfu().get_objects({'object_refs': [upa]})['data'][0]
    # normalize upa just in case it's a ref vs. an upa
    upa = str(obj_data['info'][WSID_I]) + '/' + str(obj_data['info'][OBJID_I]) + '/' + str(
        obj_data['info'][VERSION_I])
    obj_type = obj_data['info'][TYPE_I].split('-')[0]

    if 'KBaseSets.GenomeSet' == obj_type:
        upas = [gsi['ref'] for gsi in obj_data['data']['items']]
    elif 'KBaseSearch.GenomeSet' == obj_type:
        upas = [gse['ref'] for gse in obj_data['data']['elements'].values()]
    elif 'KBaseGenomes.Genome' == obj_type:
        upas = [upa]
    else:
        raise ValueError(f'{obj_type} type is not supported')

    return _process_genome_objs(upa, upas, classification, clients)


def _process_genome_objs(upa, upas, classification, clients):
    [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I,
     WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple    

    updated_genome_refs = dict()
    genomeset_query = False
    top_obj = clients.dfu().get_objects({'object_refs': [upa]})['data'][0]
    top_wsid = top_obj['info'][WSID_I]
    
    for genome_upa in upas:
        original_upa = genome_upa
        
        if upa != genome_upa:  # for single genomes, upa and genome_upa will be the same
            #genome_upa = upa + ';' + genome_upa
            genomeset_query = True
        # for a single genome this is pulling the same data again. Could optimize here.
        genome_obj = clients.dfu().get_objects({'object_refs': [genome_upa]})['data'][0]
        genome_name = genome_obj['info'][NAME_I]
        assembly_upa = genome_obj['data'].get('contigset_ref') or genome_obj['data'].get('assembly_ref')
        assembly_name = clients.ws().get_object_info_new({'objects': [{'ref':assembly_upa}]})[0][NAME_I]

        if assembly_name not in classification:
            print ("missing classification for "+assembly_name)
            updated_genome_refs[genome_upa] = genome_upa
            continue
        
        this_classification = classification[assembly_name]

        # set taxon_assignments
        if not genome_obj['data'].get('tason_assignments'):
            genome_obj['data']['taxon_assignments'] = {'GTDB_R06RS202': classification[assembly_name]}
        else:
            genome_obj['data']['taxon_assignments']['GTDB_R06RS202'] = classification[assembly_name]
            
        # set taxonomy (if missing)
        if not genome_obj['data'].get('tasonomy') \
           or genome_obj['data']['taxonomy'].startswith('Unconfirmed') \
           or genome_obj['data']['taxonomy'].startswith('Unknown'):
            genome_obj['data']['taxonomy'] = classification[assembly_name]

        # save updated genome obj
        this_wsid = genome_obj['info'][WSID_I]
        if genomeset_query:
            this_wsid = top_wsid
        updated_obj_info = clients.dfu().save_objects(
            { 'id': this_wsid,
              'objects': [{ 'type': 'KBaseGenomes.Genome',
                            'name': genome_name,
                            'data': genome_obj['data']
                            }]})[0]
        updated_genome_refs[original_upa] = '/'.join([str(updated_obj_info[WSID_I]),
                                                      str(updated_obj_info[OBJID_I]),
                                                      str(updated_obj_info[VERSION_I])])

    # update refs in genomeset
    if genomeset_query:
        new_genomeset_data = dict()
        genomeset_obj = clients.dfu().get_objects({'object_refs': [upa]})['data'][0]
        genomeset_name = genomeset_obj['info'][NAME_I]
        genomeset_type = genomeset_obj['info'][TYPE_I].split('-')[0]

        if genomeset_type == 'KBaseSets.GenomeSet':
            if genomeset_obj['data'].get('description'):
                new_genomeset_data['description'] = genomeset_obj['data']['description']
            else:
                new_genomeset_data['description'] = ''

            items = []
            for gsi in obj_data['data']['items']:
                label = ''
                if gsi.get('label'):
                    label = gsi['label']
                if not updated_genome_refs.get(gsi['ref']):
                    raise ValueError ('unable to find '+gsi['ref']+' in updated_genome_refs')
                new_ref = updated_genome_refs[gsi['ref']]
                items.append({'label': label, 'ref': new_ref})
            new_genomeset_data['itens'] = items
        elif genomeset_type == 'KBaseSearch.GenomeSet':
            if genomeset_obj['data'].get('description'):
                new_genomeset_data['description'] = genomeset_obj['data']['description']
            else:
                new_genomeset_data['description'] = ''
            
            elements = dict()
            old_elements = genomeset_obj['data']['elements']
            for genome_id in old_elements.keys():
                if not updated_genome_refs.get(old_elements[genome_id]['ref']):
                    raise ValueError ('unable to find '+old_elements[genome_id]['ref']+' in updated_genome_refs')
                new_ref = updated_genome_refs[old_elements[genome_id]['ref']]
                elements[genome_id] = {'ref': new_ref}
            new_genomeset_data['elements'] = elements
        else:
            raise ValueError(f'{genomeset_type} type is not supported')

        updated_obj_info = clients.dfu().save_objects(
            { 'id': top_wsid,
              'objects': [{ 'type': 'KBaseSearch.GenomeSet',
                            'name': genomeset_name,
                            'data': new_genomeset_data
                            }]})[0]
        

    return
