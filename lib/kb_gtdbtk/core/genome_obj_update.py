'''
Downloads sequence data from KBase services.
'''

from typing import Dict

from kb_gtdbtk.core.kb_client_set import KBClients


# global indices for KBase obj info list
[OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I,
 WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple    


# get_obj_type()
#
def get_obj_type(
        upa: str,
        clients: KBClients
        ):
    '''
    Get obj type from info

    :param upa: string with KBase obj ref
    :returns: string with obj type
    '''
    obj_info = clients.ws().get_object_info_new({'objects': [{'ref':upa}]})[0]
    obj_type = obj_info[TYPE_I].split('-')[0]
    return obj_type


# check_obj_type_genome()
#
def check_obj_type_genome(
        obj_type: str
        ):
    '''
    Check to see if query is a Genome or GenomeSet

    Handles types:
    KBaseSets.GenomeSet
    KBaseSearch.GenomeSet
    KBaseGenomes.Genome

    :param obj_type: string with KBase obj type
    :returns: True or False
    '''
    genome_obj_types = ['KBaseSets.GenomeSet',
                        'KBaseSearch.GenomeSet',
                        'KBaseGenomes.Genome'
                        ]

    if obj_type in genome_obj_types:
        return True
    else:
        return False


# check_obj_type_assembly()
#
def check_obj_type_assembly(
        obj_type: str
        ):
    '''
    Check to see if query is an Assembly or AssemblySet

    Handles types:
    KBaseSets.AssemblySet
    KBaseGenomeAnnotations.Assembly

    :param obj_type: string with KBase obj type
    :returns: True or False
    '''

    assembly_obj_types = ['KBaseSets.AssemblySet',
                        'KBaseGenomeAnnotations.Assembly'
                        ]

    if obj_type in assembly_obj_types:
        return True
    else:
        return False


# update_genome_assembly_objs_class()
#
def update_genome_assembly_objs_class(
        upa: str,
        classification: Dict[str, str],
        overwrite_tax: int,
        gtdb_ver: str,
        taxon_assignment_field: str,
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

    top_obj = clients.dfu().get_objects({'object_refs': [upa]})['data'][0]
    # normalize upa just in case it's a ref vs. an upa
    upa = upa_from_info(top_obj['info'])
    obj_data = top_obj['data']
    obj_type = top_obj['info'][TYPE_I].split('-')[0]

    if check_obj_type_genome (obj_type):
        if 'KBaseSets.GenomeSet' == obj_type:
            upas = [gsi['ref'] for gsi in obj_data['items']]
        elif 'KBaseSearch.GenomeSet' == obj_type:
            upas = [gse['ref'] for gse in obj_data['elements'].values()]
        elif 'KBaseGenomes.Genome' == obj_type:
            upas = []
        else:
            raise ValueError(f'{obj_type} type is not supported')

        return process_genome_objs(top_obj, upa, upas, classification, overwrite_tax, gtdb_ver, taxon_assignment_field, clients)

    elif check_obj_type_assembly (obj_type):
        if 'KBaseSets.AssemblySet' == obj_type:
            upas = [asi['ref'] for asi in obj_data['items']]
        elif 'KBaseGenomeAnnotations.Assembly' == obj_type:
            upas = []
        else:
            raise ValueError(f'{obj_type} type is not supported')

        return process_assembly_objs(top_obj, upa, upas, classification, overwrite_tax, gtdb_ver, taxon_assignment_field, clients)

    else:
        raise ValueError(f'{obj_type} type is not supported')


# process_genome_objs()
#
def process_genome_objs(top_obj, upa, upas, classification, overwrite_tax, gtdb_ver, taxon_assignment_field, clients):
    objects_created = []
    updated_genome_refs = dict()
    genomeset_query = False
    any_genome_updated = False
    top_wsid = top_obj['info'][WSID_I]
    
    if len(upas) > 0:
        genomeset_query = True
        genomeset_obj = top_obj
    else:
        upas = [upa]

    for genome_upa in upas:
        original_upa = genome_upa

        if not genomeset_query:
            genome_obj = top_obj
        else:
            genome_obj = clients.dfu().get_objects({'object_refs': [genome_upa]})['data'][0]
        genome_name = genome_obj['info'][NAME_I]
        assembly_upa = genome_obj['data'].get('contigset_ref') or genome_obj['data'].get('assembly_ref')
        assembly_obj = clients.dfu().get_objects({'object_refs': [assembly_upa]})['data'][0]
        assembly_name = assembly_obj['info'][NAME_I]

        if assembly_name not in classification:
            print ("missing classification for "+assembly_name)
            updated_genome_refs[genome_upa] = genome_upa
            continue
        
        # set std_lineage GTDB field in genome and assembly objs
        this_classification = classification[assembly_name]
        this_taxon_id = get_taxon_id (this_classification)
        std_lineages = get_std_lineages (this_classification, gtdb_ver, this_taxon_id)

        # update and save assembly and give genome obj new assembly upa
        new_assembly_ref = update_and_save_assembly (assembly_obj, std_lineages, clients)

        # update genome obj with std_lineages and new assembly obj ref
        genome_obj['data']['assembly_ref'] = new_assembly_ref
        genome_obj['data']['std_lineages'] = std_lineages
        any_genome_updated = True
        
        # set taxon_assignments
        if not genome_obj['data'].get('taxon_assignments'):
            genome_obj['data']['taxon_assignments'] = {taxon_assignment_field: classification[assembly_name]}
        else:
            genome_obj['data']['taxon_assignments'][taxon_assignment_field] = classification[assembly_name]
            
        # set taxonomy (if missing or force overwrite)
        this_genome_tax_written = False
        if overwrite_tax == 1 \
           or not genome_obj['data'].get('taxonomy') \
           or genome_obj['data']['taxonomy'].startswith('Unconfirmed') \
           or genome_obj['data']['taxonomy'].startswith('Unknown'):
            this_genome_tax_written = True
            any_genome_updated = True
            genome_obj['data']['taxonomy'] = classification[assembly_name]

        # save updated genome obj
        this_wsid = genome_obj['info'][WSID_I]
        if genomeset_query:
            this_wsid = top_wsid
        new_ref = save_genome_obj (this_wsid, genome_name, genome_obj['data'], clients)
        updated_genome_refs[original_upa] = new_ref
        desc = 'Taxonomy unchanged, taxon_assignment added GTDB'
        if this_genome_tax_written:
            desc = 'Taxonomy and taxon_assignment updated with GTDB'
        objects_created.append({'ref': new_ref, 'description': desc})
        
        
    # update refs in genomeset
    if genomeset_query:
        new_genomeset_ref = update_and_save_genomeset (top_wsid, genomeset_obj, updated_genome_refs, clients)
        desc = 'Taxonomy unchanged, taxon_assignment added GTDB'
        if any_genome_updated:
            desc = 'Taxonomy and taxon_assignment updated with GTDB'
        objects_created.append({'ref': new_genomeset_ref, 'description': desc})
        
    return objects_created


# process_assembly_objs()
#
def process_assembly_objs(top_obj, upa, upas, classification, overwrite_tax, gtdb_ver, taxon_assignment_field, clients):
    objects_created = []
    updated_assembly_refs = dict()
    assemblyset_query = False
    any_assembly_updated = False
    top_wsid = top_obj['info'][WSID_I]

    if len(upas) > 0:
        assemblyset_query = True
        assemblyset_obj = top_obj
    else:
        upas = [upa]
    
    for assembly_upa in upas:
        original_upa = assembly_upa
        
        if not assemblyset_query:
            assembly_obj = clients.dfu().get_objects({'object_refs': [assembly_upa]})['data'][0]
        else:
            assembly_obj = top_obj
        assembly_name = assembly_obj['info'][NAME_I]

        if assembly_name not in classification:
            print ("missing classification for "+assembly_name)
            updated_assembly_refs[assembly_upa] = assembly_upa
            continue
        else:
            any_assembly_updated = True
        

        # set std_lineage GTDB field in assembly obj
        this_classification = classification[assembly_name]
        this_taxon_id = get_taxon_id (this_classification)
        std_lineages = get_std_lineages (this_classification, gtdb_ver, this_taxon_id)

        # update and save assembly
        new_assembly_ref = update_and_save_assembly (assembly_obj, std_lineages, clients)
        updated_assembly_refs[original_upa] = new_assembly_ref
        desc = 'Added GTDB lineage'
        objects_created.append({'ref': new_assembly_ref, 'description': desc})
        
        
    # update refs in assemblyset
    if assemblyset_query and any_assembly_updated:
        new_assemblyset_ref = update_and_save_assemblyset (top_wsid,
                                                           assemblyset_obj,
                                                           updated_assembly_refs,
                                                           clients)
        desc = 'Added GTDB lineage'
        objects_created.append({'ref': new_assemblyset_ref, 'description': desc})
        
    return objects_created


# get_taxon_id ()
#
def get_taxon_id (this_classification):
    #  Note: sometimes lineage does not resolve species, genus, etc.
    #        in this case the string is 's__', 'g__', etc.
    #        for taxon_id we want the first resolved taxon level, hence the skip if len==3

    this_taxon_id = None
    for taxon_id in reversed(this_classification.split(';')):
        if len(taxon_id) == 3:
            continue
        this_taxon_id = taxon_id
        break
    if this_taxon_id == None:
        raise ValueError ("unable to fund taxon+id for {}".format(this_classification))

    return this_taxon_id

    
# get_std_lineages ()
#
def get_std_lineages (this_classification, gtdb_ver, this_taxon_id):
    return { 'gtdb': { 'lineage': this_classification,
                       'source_ver': gtdb_ver,
                       'taxon_id': this_taxon_id
                       #'source_id': None
                     }
    }


# save_genome_obj ()
#
def save_genome_obj (this_wsid, genome_name, genome_obj_data, clients):
    updated_obj_info = clients.dfu().save_objects(
        { 'id': this_wsid,
          'objects': [{ 'type': 'KBaseGenomes.Genome',
                        'name': genome_name,
                        'data': genome_obj_data
          }]})[0]
    new_ref = upa_from_info(updated_obj_info)
    return new_ref


# update_and_save_assembly ()
#
def update_and_save_assembly (assembly_obj, std_lineages, clients):
    assembly_obj['data']['std_lineages'] = std_lineages
    this_wsid = assembly_obj['info'][WSID_I]
    assembly_name = assembly_obj['info'][NAME_I]

    updated_assembly_obj_info = clients.dfu().save_objects({ 'id': this_wsid,
                                                             'objects': [{ 'type': 'KBaseGenomeAnnotations.Assembly',
                                                                           'name': assembly_name,
                                                                           'data': assembly_obj['data']
                                                             }]})[0]
    new_assembly_ref = upa_from_info(updated_assembly_obj_info)
    return new_assembly_ref


# update_and_save_genomeset ()
#
def update_and_save_genomeset (top_wsid, genomeset_obj, updated_genome_refs, clients):
    new_genomeset_data = dict()
    genomeset_name = genomeset_obj['info'][NAME_I]
    genomeset_type = genomeset_obj['info'][TYPE_I].split('-')[0]

    # repoint to new versions of genome objs
    if genomeset_type == 'KBaseSets.GenomeSet':
        if genomeset_obj['data'].get('description'):
            new_genomeset_data['description'] = genomeset_obj['data']['description']
        else:
            new_genomeset_data['description'] = ''

        items = []
        for gsi in genomeset_obj['data']['items']:
            label = ''
            if gsi.get('label'):
                label = gsi['label']
            if not updated_genome_refs.get(gsi['ref']):
                raise ValueError ('unable to find '+gsi['ref']+' in updated_genome_refs')
            new_ref = updated_genome_refs[gsi['ref']]
            items.append({'label': label, 'ref': new_ref})
        new_genomeset_data['items'] = items
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

    # save updated genome set obj
    updated_obj_info = clients.dfu().save_objects(
        { 'id': top_wsid,
          'objects': [{ 'type': 'KBaseSearch.GenomeSet',
                        'name': genomeset_name,
                        'data': new_genomeset_data
          }]})[0]
    new_ref = upa_from_info(updated_obj_info)
    return new_ref


# update_and_save_assemblyset ()
#
def update_and_save_assemblyset (top_wsid, assemblyset_obj, updated_assembly_refs, clients):
    new_assemblyset_data = dict()
    assemblyset_name = assemblyset_obj['info'][NAME_I]
    assemblyset_type = assemblyset_obj['info'][TYPE_I].split('-')[0]

    # repoint to new versions of assembly objs
    if assemblyset_type == 'KBaseSets.AssemblySet':
        if assemblyset_obj['data'].get('description'):
            new_assemblyset_data['description'] = assemblyset_obj['data']['description']
        else:
            new_assemblyset_data['description'] = ''

        items = []
        for asi in assemblyset_obj['data']['items']:
            label = ''
            if asi.get('label'):
                label = asi['label']
            if not updated_assembly_refs.get(asi['ref']):
                raise ValueError ('unable to find '+asi['ref']+' in updated_assembly_refs')
            new_ref = updated_assembly_refs[asi['ref']]
            items.append({'label': label, 'ref': new_ref})
        new_assemblyset_data['items'] = items
    else:
        raise ValueError(f'{assemblyset_type} type is not supported')

    # save updated assembly set obj
    updated_obj_info = clients.dfu().save_objects(
        { 'id': top_wsid,
          'objects': [{ 'type': 'KBaseSets.AssemblySet',
                        'name': assemblyset_name,
                        'data': new_assemblyset_data
          }]})[0]
    new_ref = upa_from_info(updated_obj_info)
    return new_ref


# upa_from_info ()
#
def upa_from_info (obj_info):
    return '/'.join([str(obj_info[WSID_I]),
                     str(obj_info[OBJID_I]),
                     str(obj_info[VERSION_I])])
