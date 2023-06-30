'''
Downloads sequence data from KBase services.
'''

from typing import Dict
import os
import json
import re
import pandas as pd
import subprocess
import ete3

from kb_gtdbtk.core.kb_client_set import KBClients


# global indices for KBase obj info list
[OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I,
 WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple    

# global for tree html
tree_key_taxon_seen = dict()


# get_obj_info ()
def get_obj_info (
        wsid: int,
        objname: str,
        clients: KBClients
        ):
    '''
    Get obj type from info

    :param wsid: int with KBase workspace id
    :param objname: string with KBase obj name
    :returns: string with obj type
    '''

    obj_info = None
    try:
        obj_info = clients.ws().get_object_info3({'objects': [{'wsid':wsid,'name':objname}]})['infos'][0]
    except:
        pass
    return obj_info


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
    obj_info = clients.ws().get_object_info3({'objects': [{'ref':upa}]})['infos'][0]
    obj_type = obj_info[TYPE_I].split('-')[0]
    return obj_type


# get_obj_name()
#
def get_obj_name(
        upa: str,
        clients: KBClients
        ):
    '''
    Get obj type from info

    :param upa: string with KBase obj ref
    :returns: string with obj type
    '''
    obj_info = clients.ws().get_object_info3({'objects': [{'ref':upa}]})['infos'][0]
    obj_name = obj_info[NAME_I]
    return obj_name


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

    
# get_names_list_from_upas_list ()
#
def get_names_list_from_upas_list (upas_list, clients):
    obj_names_list = []
    for upa in upas_list:
        obj_name = get_obj_name(upa,clients)
        obj_names_list.append(obj_name)

    return obj_names_list


# pause()
#
def pause(patience, sleep_interval, target_files):
    found_output = False
    
    time_spent = 0
    while time_spent <= patience:
        found_output = True
        for this_file in target_files:
            if not os.path.exists(this_file) or not os.path.getsize(this_file) > 0:
                found_output = False
                print ("waiting for {} seconds...".format(sleep_interval))
                time.sleep(sleep_interval)
                time_spent += sleep_interval
        if found_output:
            break
    if not found_output:
        raise ValueError ("Never produced output in {} seconds".format(patience))

    return

    
# _format_gtdbtk_tree_to_itol ()
#
def _format_gtdbtk_tree_to_itol (in_tree_path):
    print ("making ITOL format tree "+str(in_tree_path))
    gtdbtk_bin = os.path.join ('/opt', 'conda3', 'bin', 'gtdbtk')
    #itol_tree_path = str(in_tree_path).replace('.tree','-ITOL.tree')
    itol_tree_path = re.sub('.tree$', '-ITOL.tree', str(in_tree_path))

    convert_cmd = [gtdbtk_bin, 'convert_to_itol',
                   '--input_tree', in_tree_path,
                   '--output_tree', itol_tree_path
                   ]
    env = dict(os.environ)
    subprocess.run(convert_cmd, check=True, env=env)
    #pause(600, 10, [itol_tree_path])

    return itol_tree_path

    
# _trim_tree()
#
def _trim_tree (in_tree_path, leaflist_file, leaflist_outfile, lineage_outfile):
    print ("trimming tree "+str(in_tree_path))
    trim_bin = os.path.join ('/kb', 'module', 'bin', 'trim_tree_to_target_leaves.py')
    
    out_tree_paths = []
    
    # just proximal sp rep hits
    #out_tree_path = str(in_tree_path).replace('.tree','-proximals.tree')
    out_tree_path = re.sub('.tree$', '-proximals.tree', str(in_tree_path))
    out_tree_paths.append(out_tree_path)
    trim_cmd = [trim_bin,
                '--intree', str(in_tree_path),
                '--outtree', str(out_tree_path),
                '--leaflist', str(leaflist_file),
                '--targetleafoutfile', str(leaflist_outfile)
                ]
    print ("RUNNING: "+" ".join(trim_cmd))
    env = dict(os.environ)
    subprocess.run(trim_cmd, check=True, env=env)
    #pause(600, 10, [out_tree_path, leaflist_outfile])

    
    # with sister context branches.  Note that leaflist_outfile is updated with context sp reps
    #out_tree_path = str(in_tree_path).replace('.tree','-trimmed.tree')
    out_tree_path = re.sub('.tree$', '-trimmed.tree', str(in_tree_path))
    out_tree_paths.append(out_tree_path)
    trim_cmd = [trim_bin,
                '--intree', str(in_tree_path),
                '--outtree', str(out_tree_path),
                '--leaflist', str(leaflist_file),
                '--targetleafoutfile', str(leaflist_outfile),
                '--gtdblineageoutfile', str(lineage_outfile)
                ]
    trim_cmd.append('--sisters')
    print ("RUNNING: "+" ".join(trim_cmd))
    env = dict(os.environ)
    subprocess.run(trim_cmd, check=True, env=env)
    #pause(600, 10, [out_tree_path, leaflist_outfile, lineage_outfile])
    
    return out_tree_paths


# _write_tree_image_file()
#
def _write_tree_image_file (trimmed_tree_path, query_leaflist_file, leaflist_file, lineage_file):
    print ("making images for trimmed tree "+str(trimmed_tree_path))
    write_image_bin = os.path.join ('/kb', 'module', 'bin', 'make_tree_images.py')

    out_img_base = trimmed_tree_path
    trimmed_tree_image_paths = [ out_img_base+'-rectangle.PNG',
                                 out_img_base+'-rectangle.PDF',
                                 out_img_base+'-circle.PNG',
                                 out_img_base+'-circle.PDF',
                                 out_img_base+'-circle-ultrametric.PNG',
                                 out_img_base+'-circle-ultrametric.PDF'
                               ]

    title = os.path.basename(out_img_base)
    write_image_cmd = [write_image_bin,
                       '--intree', trimmed_tree_path,
                       '--title', title,
                       '--outimgbase', out_img_base,
                       '--queryleaflist', query_leaflist_file,
                       '--leaflist', leaflist_file,
                       '--gtdblineagefile', lineage_file
                       ]
    env = dict(os.environ)
    subprocess.run(write_image_cmd, check=True, env=env)
    #pause(600, 10, [out_img_base+'-circle.PNG', out_img_base+'-circle.PDF'])

    return trimmed_tree_image_paths


# _write_gtdb_tree_html_file ()
#
def _write_gtdb_tree_html_file (out_dir, files_for_html):
    tree_html_path = os.path.join (out_dir, 'gtdb_trees.html')  # if updated, the index.html file should also be updated accordingly

    tax_order = ['p', 'c', 'o', 'f', 'g']  # not handling domain or species
    
    html_buf = []
    table_buf = []
    for file_for_html in files_for_html:
        tree_newick_path = file_for_html['newick_path']
        tree_png_file = file_for_html['png_file']
        taxon_colors_path = file_for_html['taxon_colors_path']
        lineage_path = file_for_html['lineage_path']
        
        # add tree image
        tree_img_width = 500
        tree_img_height = 500
        table_buf += ['<tr>']
        table_buf += ['<td align=left valign=top border=0><img src="{}" border=0 width={} height={}></td>'.format(tree_png_file, tree_img_width, tree_img_height)]
        table_buf += ['<td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td>']

        # get taxon colors
        taxon_color = dict()
        with open (taxon_colors_path, 'r') as taxon_colors_h:
            for line in taxon_colors_h:
                (taxon, color) = line.rstrip().split("\t")
                taxon_color[taxon] = color

        # get full parent mapping from leaves
        all_leaf_lineages = get_all_leaf_lineages (lineage_path)
        parents = get_parents (all_leaf_lineages)
        
        
        # determine lineage structure
        lineages = dict()
        tree = ete3.Tree(tree_newick_path, quoted_node_names=True, format=1)
        tree.ladderize()
        for leaf_node in tree.get_leaves():
            """
            #last_taxon = None
            #for ancestor_node in reversed(leaf_node.get_ancestors()):  # phylum -> genus
            """
            for ancestor_node in leaf_node.get_ancestors():  # genus -> phylum
                if ancestor_node.name:
                    taxon = ancestor_node.name
                    if taxon not in parents:  # handle phylum
                        lineages[taxon] = dict()
                    else:
                        this_taxon = taxon
                        if this_taxon not in lineages:  # handle genus
                            lineages[this_taxon] = dict()
                        level_limit = 5
                        level_i = 0
                        while this_taxon in parents:
                            if level_i > level_limit:  # avoid accidental infinite
                                break
                            level_i += 1
                            parent_taxon = parents[this_taxon]
                            if parent_taxon not in lineages:
                                lineages[parent_taxon] = dict()
                            lineages[parent_taxon][this_taxon] = True
                            this_taxon = parent_taxon

                    """
                    if last_taxon:
                        if last_taxon not in lineages:
                            lineages[last_taxon] = dict()
                        lineages[last_taxon][ancestor_node.name] = True
                    last_taxon = ancestor_node.name
                    """
                    
        # build key
        indent_cnt = 0
        key_box_size = '15px'
        key_font_size = key_box_size
        table_buf += ['<td align=left valign=middle><table border=0>']
        #tax_level_order = ['d', 'p', 'c', 'o', 'f', 'g']
        tax_level_order = ['p', 'c', 'o', 'f', 'g']
        first_row = True
        for tax_level in tax_level_order:  # phylum -> genus
            for taxon in sorted (lineages.keys()):
                if taxon[0] != tax_level:
                    continue
                if taxon not in tree_key_taxon_seen:
                    if first_row:
                        first_row = False
                    else:
                        table_buf += ['<tr><td>&nbsp;</td></tr>']
                    table_buf += add_tree_key_row (taxon, lineages, taxon_color, indent_cnt, key_box_size, key_font_size)
        table_buf += ['</table>']
        table_buf += ['</td></tr>']

        
    # build html
    html_buf += ['<html><head><title>GTDB-Tk Trees</title></head>']
    html_buf += ['<body><table border=0>']
    html_buf += table_buf
    html_buf += ['</table></body></html>']
    with open (tree_html_path, 'w') as tree_h:
        tree_h.write("\n".join(html_buf)+"\n")
        
    return tree_html_path


# get_all_leaf_lineages ()
#
def get_all_leaf_lineages (lineage_file):

    print ("reading target lineages from file {} ...".format(lineage_file))

    all_leaf_lineages = dict()
    
    if lineage_file.lower().endswith('.gz'):
        f = gzip.open(lineage_file, 'rt')
    else:
        f = open(lineage_file, 'r')

    for line in f:
        line = line.rstrip()

        (leaf_name, this_lineage) = line.split("\t")

        base_leaf_id = leaf_name.split(' ')[0].lstrip('"')
        all_leaf_lineages[base_leaf_id] = this_lineage
    f.close()

    return all_leaf_lineages


# get_parents ()
#
def get_parents (all_leaf_lineages):
    parent_taxon_map = dict()

    for leaf_id in all_leaf_lineages.keys():
        lineage = all_leaf_lineages[leaf_id].split(';')

        for tax_i,taxon in enumerate(lineage):
            if tax_i > 0:
                parent_taxon_map[taxon] = lineage[tax_i-1]

    return parent_taxon_map


# add_tree_key_row ()
#
def add_tree_key_row (taxon, lineages, taxon_color, indent_cnt, key_box_size, key_font_size):
    key_buf = []

    if taxon in tree_key_taxon_seen:
        return []
    elif taxon[0] == 's':
        return []
    else:
        tree_key_taxon_seen[taxon] = True
        
    key_buf += ['<tr><td>']
    key_buf += ['<table border=0><tr>']
    for indent_i in range(indent_cnt):
        key_buf += ['<td>&nbsp;</td><td>&nbsp;</td>']
    if taxon in taxon_color:
        key_buf += ['<td style="width:{};height:{};background-color:#{}"></td>'.format(key_box_size, key_box_size, taxon_color[taxon])]
    key_buf += ['<td><p style="font-size:{}">'.format(key_font_size)+taxon+'</p></td>']
    key_buf += ['</tr></table>'] 
    key_buf += ['</td></tr>']

    if taxon not in lineages:
        return key_buf

    indent_cnt += 1
    for child_taxon in sorted(lineages[taxon].keys()):
        key_buf += add_tree_key_row (child_taxon, lineages, taxon_color, indent_cnt, key_box_size, key_font_size)

    return key_buf


# process_tree_files()
#
def process_tree_files (top_upa,
                        out_dir,
                        summary_tables,
                        classification,
                        clients):
    upload_files = []
    file_links = []

    id_map_path = os.path.join(out_dir, 'id_to_name.map')
    tree_files = ['gtdbtk.ar53.classify.tree',
                  'gtdbtk.bac120.classify.tree']
    extra_bac_tree_files = []
    for i in range(10000):
        subtree_file = 'gtdbtk.bac120.classify.tree.'+str(i)+'.tree'
        extra_bac_tree_files.append(subtree_file)

    # add species hits to id_map
    top_obj = clients.dfu().get_objects({'object_refs': [top_upa]})['data'][0]
    query_assembly_to_genome_name = get_query_assembly_to_genome_name (top_obj, clients)
    (all_sp_reps, sp_reps_by_query) = get_sp_rep_hits(summary_tables, query_assembly_to_genome_name)

    # change id map to genome names
    id_map_buf = []
    with open(id_map_path, 'r') as id_map_h:
        for line in id_map_h:
            (qid, assembly_name) = line.rstrip().split("\t") 
            new_name = assembly_name
            if assembly_name in query_assembly_to_genome_name:
                new_name = query_assembly_to_genome_name[assembly_name]
            this_lineage = '-'
            if assembly_name in classification:
                this_lineage = classification[assembly_name]
            id_map_buf.append("\t".join([qid,new_name,this_lineage]))
    new_id_map_path = str(id_map_path).replace('.map','-genomes.map')
    with open(new_id_map_path, 'w') as id_map_h:
        id_map_h.write("\n".join(id_map_buf)+"\n")
            
    # add sp rep hits too 
    for sp_rep_id in sorted(all_sp_reps.keys()):
        id_map_buf.append("\t".join([sp_rep_id,sp_rep_id]))
    id_map_with_sp_rep_hits_path = os.path.join(out_dir, 'id_to_name-with_proximal_sp_reps.map')
    with open(id_map_with_sp_rep_hits_path, 'w') as id_map_h:
        id_map_h.write("\n".join(id_map_buf)+"\n")

    # make itol format files
    for tree_file in tree_files + extra_bac_tree_files:
        in_tree_path = out_dir / tree_file
        if os.path.isfile(in_tree_path):
            itol_tree_path = _format_gtdbtk_tree_to_itol (in_tree_path)
            itol_tree_file = os.path.basename(itol_tree_path)
            upload_files.append({ 'path': in_tree_path,
                                  'name': tree_file,
                                  'description': tree_file+' - whole tree GTDB formatted Newick'})
            upload_files.append({ 'path': itol_tree_path,
                                  'name': itol_tree_file,
                                  'description': itol_tree_file+' - whole tree ITOL formatted Newick'})

            
    # trim tree files and make tree image files
    files_for_html = []
    for tree_file in tree_files + extra_bac_tree_files:
        in_tree_path = out_dir / tree_file
        if os.path.isfile(in_tree_path):

            #new_id_map_with_sp_rep_hits_path = str(in_tree_path).replace('.tree', '.id_to_name-with_all_sp_reps-newleafnames.map')
            new_id_map_with_sp_rep_hits_path = re.sub('.tree$', '.id_to_name-with_proximal_sp_reps-newleafnames.map', str(in_tree_path))
            lineage_path = re.sub('.tree$', '-lineages.map', str(in_tree_path))
            
            trimmed_tree_paths = _trim_tree (in_tree_path, id_map_with_sp_rep_hits_path, new_id_map_with_sp_rep_hits_path, lineage_path)

            for trimmed_tree_path in trimmed_tree_paths:

                trimmed_tree_file = os.path.basename (trimmed_tree_path)
                upload_files.append({ 'path': trimmed_tree_path,
                                      'name': trimmed_tree_file,
                                      'description': trimmed_tree_file+' - Newick'
                                    })

                # proximals failing.  fix later
                if 'proximals' in trimmed_tree_path:
                    continue

                lineage_file = os.path.basename (lineage_path)
                upload_files.append({ 'path': lineage_path,
                                      'name': lineage_file,
                                      'description': lineage_file+' - GTDB lineage'
                                    })

                trimmed_tree_image_paths = _write_tree_image_file (trimmed_tree_path,
                                                                   new_id_map_path,
                                                                   new_id_map_with_sp_rep_hits_path,
                                                                   lineage_path)

                for trimmed_tree_image_path in trimmed_tree_image_paths:
                    trimmed_tree_image_file = os.path.basename (trimmed_tree_image_path)
                    upload_files.append({ 'path': trimmed_tree_image_path,
                                          'name': trimmed_tree_image_file,
                                          'description': trimmed_tree_file+' - Image'
                                        })

                    if '-trimmed.tree-circle-ultrametric.PNG' in trimmed_tree_image_file:
                        taxon_colors_path = str(trimmed_tree_image_path).replace('-circle-ultrametric.PNG','-taxon_colors.map')
                        files_for_html.append({'newick_path': trimmed_tree_path,
                                               'png_file': trimmed_tree_image_file,
                                               'taxon_colors_path': taxon_colors_path,
                                               'lineage_path': lineage_path})

                        
    # upload files and make file links for report
    #
    file_links = []
    for f in upload_files:
        upload_ret = clients.dfu().file_to_shock({'file_path': f['path'],
                                                  'make_handle': 0})
        file_links.append({'shock_id': upload_ret['shock_id'],
                           'name': f['name'],
                           'description': f['description']})


    # Make GTDB Tree html to go in html report
    #
    tree_html_file = _write_gtdb_tree_html_file (out_dir, files_for_html)
        
        
    return file_links


# get_genome_to_upa_map()
#
def get_genome_to_upa_map(genome_upas_map_file):
    genome_id_to_upa_map = dict()
    this_file = genome_upas_map_file
    if this_file.lower().endswith('.gz'):
        f = gzip.open(this_file, 'rt')
    else:
        f = open(this_file, 'r')
    for line in f:
        line = line.rstrip()
        (genome_id, upa) = line.split("\t")
        genome_id_to_upa_map[genome_id] = upa
    f.close()

    return genome_id_to_upa_map


# _save_tree_obj_and_copy_genomes()
#
def  _save_tree_obj_and_copy_genomes(this_tree_path,
                                     tree_name,
                                     obj_name,
                                     tree_short_desc,
                                     top_upa,
                                     workspace_id,
                                     genome_id_to_upa_map,
                                     clients):

    print ("SAVING TREE OBJ {}".format(obj_name))

    new_objects_created = []
    tree_full_desc = tree_name+' '+tree_short_desc
    
    # load tree and get leaf naems and genome ids
    tree = ete3.Tree (this_tree_path, quoted_node_names=True, format=1)
    tree.ladderize()
    newick_buf = tree.write(features=[])
    if not newick_buf.endswith(';'):
        newick_buf += ';'
    leaf_list = []
    genome_ids = []
    for leaf_name in tree.get_leaf_names():
        leaf_list.append(leaf_name)
        genome_ids.append(leaf_name.split(' ')[0])

    # get query upas
    query_upas = dict()
    top_obj = clients.dfu().get_objects({'object_refs': [top_upa]})['data'][0]
    query_upas_list = get_upas_from_set(top_obj)
    query_names_list = get_names_list_from_upas_list (query_upas_list, clients)
    for query_i,query_name in enumerate(query_names_list):
        query_upas[query_name] = query_upas_list[query_i]
    
    # make local copies of genomes
    copied_genome_refs = copy_gtdb_genome_objs (genome_ids, genome_id_to_upa_map, 'GTDB_SP_REP-', workspace_id, clients)
    
    ws_refs = dict()
    default_node_labels = dict()
    for genome_i,genome_id in enumerate(genome_ids):
        leaf_name = leaf_list[genome_i]
        default_node_labels[leaf_name] = leaf_name
        if genome_id in copied_genome_refs:
            genome_ref = copied_genome_refs[genome_id]
        else:
            genome_ref = query_upas[genome_id]
        ws_refs[leaf_name] = dict()
        ws_refs[leaf_name]['g'] = [genome_ref]
    tree_data = { 'name': tree_name,
                  'description': tree_full_desc,
                  'type': 'SpeciesTree',
                  'tree': newick_buf,
                  'leaf_list': leaf_list,
                  'default_node_labels': default_node_labels,
                  'ws_refs': ws_refs
    }
    # save tree
    genome_refs_list = []
    for genome_id in sorted(query_upas.keys()):
        genome_refs_list.append(query_upas[genome_id])
    for genome_id in sorted(copied_genome_refs.keys()):
        genome_refs_list.append(copied_genome_refs[genome_id])
    extra_provenance_input_refs = [top_upa]
    extra_provenance_input_refs.extend(genome_refs_list)
    try:
        tree_out_obj_info = clients.dfu().save_objects({
            'id': workspace_id,
            'objects':[{
                'type': 'KBaseTrees.Tree',
                'data': tree_data,
                'name': obj_name,
                'meta': {},
                'extra_provenance_input_refs': extra_provenance_input_refs
            }]})[0]
    except Exception as e:
        raise ValueError('Unable to save tree '+obj_name+' object to workspace '+str(workspace_id)+': ' + str(e))

    output_tree_ref = upa_from_info (tree_out_obj_info)
    new_objects_created.append({'ref': output_tree_ref, 'description': tree_short_desc})

    # DEBUG
    print ("SAVED TREE OBJ {} ref: {}".format(obj_name,output_tree_ref))
    
    return new_objects_created


# save_gtdb_tree_objs()
#
def save_gtdb_tree_objs (workspace_id,
                         top_upa,
                         out_dir,
                         output_tree_basename,
                         genome_upas_map_file,
                         clients):
    new_objects_created = []

    print ("SAVING GTDB TREE OBJECTS")
    
    # get upas by genome id
    genome_id_to_upa_map = get_genome_to_upa_map(genome_upas_map_file)

    # read trees and collect contained genomes
    tree_files = ['gtdbtk.ar53.classify.tree',
                  'gtdbtk.bac120.classify.tree']
    extra_bac_tree_files = []
    for i in range(10000):
        subtree_file = 'gtdbtk.bac120.classify.tree.'+str(i)+'.tree'
        extra_bac_tree_files.append(subtree_file)

    for tree_file in tree_files + extra_bac_tree_files:
        in_tree_path = out_dir / tree_file
        if os.path.isfile(in_tree_path):

            # proximal tree
            #proximal_tree_path = str(in_tree_path).replace('.tree', '-proximal.tree')
            proximal_tree_path = re.sub('.tree$', '-proximals.tree', str(in_tree_path))
            tree_name = tree_file+'-proximals.tree'
            obj_name = output_tree_basename+'.'+tree_name
            tree_short_desc = 'with proximal GTDB species reps'

            new_objects_created.extend(_save_tree_obj_and_copy_genomes(proximal_tree_path,
                                                                       tree_name,
                                                                       obj_name,
                                                                       tree_short_desc,
                                                                       top_upa,
                                                                       workspace_id,
                                                                       genome_id_to_upa_map,
                                                                       clients))
            # trimed tree
            #trimmed_tree_path = str(in_tree_path).replace('.tree', '-trimmed.tree')
            trimmed_tree_path = re.sub('.tree$', '-trimmed.tree', str(in_tree_path))
            tree_name = tree_file+'-trimmed.tree'
            obj_name = output_tree_basename+'.'+tree_name
            tree_short_desc = 'trimmed with sister context'

            new_objects_created.extend(_save_tree_obj_and_copy_genomes(trimmed_tree_path,
                                                                       tree_name,
                                                                       obj_name,
                                                                       tree_short_desc,
                                                                       top_upa,
                                                                       workspace_id,
                                                                       genome_id_to_upa_map,
                                                                       clients))

    return new_objects_created


# update_genome_assembly_objs_class()
#
def update_genome_assembly_objs_class(
        primary_wsid: int,
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
    KBaseGenomeAnnotations.Assembly
    KBaseSets.AssemblySet

    :param primary_wsid: The ID of the calling workspace.  Important for saving objects from other workspaces
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

        return process_genome_objs(primary_wsid, top_obj, upa, upas, classification, overwrite_tax, gtdb_ver, taxon_assignment_field, clients)

    elif check_obj_type_assembly (obj_type):
        if 'KBaseSets.AssemblySet' == obj_type:
            upas = [asi['ref'] for asi in obj_data['items']]
        elif 'KBaseGenomeAnnotations.Assembly' == obj_type:
            upas = []
        else:
            raise ValueError(f'{obj_type} type is not supported')

        return process_assembly_objs(primary_wsid, top_obj, upa, upas, classification, overwrite_tax, gtdb_ver, taxon_assignment_field, clients)

    else:
        raise ValueError(f'{obj_type} type is not supported')


# process_genome_objs()
#
def process_genome_objs(primary_wsid, top_obj, upa, upas, classification, overwrite_tax, gtdb_ver, taxon_assignment_field, clients):
    objects_created = []
    updated_genome_refs = dict()
    genomeset_query = False
    any_genome_updated = False
    #top_wsid = top_obj['info'][WSID_I]
    
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
        new_assembly_ref = update_and_save_assembly (primary_wsid, assembly_obj, std_lineages, clients)

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
        #this_wsid = genome_obj['info'][WSID_I]
        #if genomeset_query:
        #    this_wsid = top_wsid
        new_ref = save_genome_obj (primary_wsid, genome_name, genome_obj['data'], clients)
        updated_genome_refs[original_upa] = new_ref
        desc = 'Taxonomy unchanged, taxon_assignment added GTDB'
        if this_genome_tax_written:
            desc = 'Taxonomy and taxon_assignment updated with GTDB'
        objects_created.append({'ref': new_ref, 'description': desc})
        
        
    # update refs in genomeset
    if genomeset_query:
        new_genomeset_ref = update_and_save_genomeset (primary_wsid, genomeset_obj, updated_genome_refs, clients)
        desc = 'Taxonomy unchanged, taxon_assignment added GTDB'
        if any_genome_updated:
            desc = 'Taxonomy and taxon_assignment updated with GTDB'
        objects_created.append({'ref': new_genomeset_ref, 'description': desc})
        
    return objects_created


# process_assembly_objs()
#
def process_assembly_objs(primary_wsid, top_obj, upa, upas, classification, overwrite_tax, gtdb_ver, taxon_assignment_field, clients):
    objects_created = []
    updated_assembly_refs = dict()
    assemblyset_query = False
    any_assembly_updated = False
    #top_wsid = top_obj['info'][WSID_I]

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
        new_assembly_ref = update_and_save_assembly (primary_wsid, assembly_obj, std_lineages, clients)
        updated_assembly_refs[original_upa] = new_assembly_ref
        desc = 'Added GTDB lineage'
        objects_created.append({'ref': new_assembly_ref, 'description': desc})
        
        
    # update refs in assemblyset
    if assemblyset_query and any_assembly_updated:
        new_assemblyset_ref = update_and_save_assemblyset (primary_wsid,
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


# copy_gtdb_genome_objs ()
#
def copy_gtdb_genome_objs (genome_ids, genome_id_to_upa_map, new_obj_name_prefix, dst_ws_id, clients):
    genome_refs = dict()

    for genome_id in genome_ids:
        if genome_id not in genome_id_to_upa_map:
            continue
        src_upa = genome_id_to_upa_map[genome_id]
        src_obj_name = get_obj_name(src_upa, clients)

        dst_obj_name = src_obj_name
        if new_obj_name_prefix:
            dst_obj_name = new_obj_name_prefix+dst_obj_name
            
        existing_obj_info = get_obj_info (dst_ws_id, dst_obj_name, clients)
        if existing_obj_info:
            genome_obj_info = existing_obj_info
        else:
            genome_obj_info = _copy_genome_obj(src_upa,
                                               dst_ws_id,
                                               dst_obj_name,
                                               clients)
        genome_refs[genome_id] = upa_from_info(genome_obj_info)

    return genome_refs

            
# save_genome_obj ()
#
def save_genome_obj (primary_wsid, genome_name, genome_obj_data, clients):

    genome_obj_data = fix_unowned_shock_handles ('genome', genome_obj_data, clients)

    updated_obj_info = clients.dfu().save_objects(
        { 'id': primary_wsid,
          'objects': [{ 'type': 'KBaseGenomes.Genome',
                        'name': genome_name,
                        'data': genome_obj_data
          }]})[0]
    new_ref = upa_from_info(updated_obj_info)
    return new_ref


# update_and_save_assembly ()
#
def update_and_save_assembly (primary_wsid, assembly_obj, std_lineages, clients):
    assembly_obj['data']['std_lineages'] = std_lineages
    #this_wsid = assembly_obj['info'][WSID_I]
    assembly_name = assembly_obj['info'][NAME_I]

    assembly_obj['data'] = fix_unowned_shock_handles ('assembly', assembly_obj['data'], clients)
    
    updated_assembly_obj_info = clients.dfu().save_objects({ 'id': primary_wsid,
                                                             'objects': [{ 'type': 'KBaseGenomeAnnotations.Assembly',
                                                                           'name': assembly_name,
                                                                           'data': assembly_obj['data']
                                                             }]})[0]
    new_assembly_ref = upa_from_info(updated_assembly_obj_info)
    return new_assembly_ref


# update_and_save_genomeset ()
#
def update_and_save_genomeset (primary_wsid, genomeset_obj, updated_genome_refs, clients):
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
        { 'id': primary_wsid,
          'objects': [{ 'type': 'KBaseSearch.GenomeSet',
                        'name': genomeset_name,
                        'data': new_genomeset_data
          }]})[0]
    new_ref = upa_from_info(updated_obj_info)
    return new_ref


# update_and_save_assemblyset ()
#
def update_and_save_assemblyset (primary_wsid, assemblyset_obj, updated_assembly_refs, clients):
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
        { 'id': primary_wsid,
          'objects': [{ 'type': 'KBaseSets.AssemblySet',
                        'name': assemblyset_name,
                        'data': new_assemblyset_data
          }]})[0]
    new_ref = upa_from_info(updated_obj_info)
    return new_ref


# get_sp_rep_hits()
#
def get_sp_rep_hits (summary_tables, query_assembly_to_genome_name):
    sp_reps_by_query = dict()
    all_sp_reps = dict()

    single_sp_rep_fields = ['fastani_reference', 'closest_placement_reference']
    multi_sp_rep_field = 'other_related_references(genome_id,species_name,radius,ANI,AF)'

    for summary_file in ['gtdbtk.bac120.summary.tsv', 'gtdbtk.ar53.summary.tsv']:
        if summary_file not in summary_tables:
            continue
        sj = summary_tables[summary_file]
        for item in sj['data']:
            # reset id to assembly name
            # Note: field 'Name' was changed to 'name'
            #key = 'Name' if 'Name' in item else 'user_genome'
            if 'name' in item:
                key = 'name'
            elif 'user_genome' in item:
                key = 'user_genome'
            else:
                continue
            #this_assembly_id = item[key].replace('_assembly','')
            this_assembly_id = item[key]
            #print ("DEBUG: adding hits for query assembly ID {}".format(this_assembly_id))  # DEBUG
            this_genome_id = query_assembly_to_genome_name[this_assembly_id]
            sp_reps_by_query[this_genome_id] = []

            #print ("DEBUG: adding hits for query genome ID {}".format(this_genome_id))  # DEBUG
            
            # single value id
            for sp_rep_f in single_sp_rep_fields:
                if item.get(sp_rep_f) and item.get(sp_rep_f) != '-':
                    sp_rep_id = item[sp_rep_f]
                    all_sp_reps[sp_rep_id] = True
                    if sp_rep_id not in sp_reps_by_query[this_genome_id]:
                        sp_reps_by_query[this_genome_id].append(sp_rep_id)

            # multiple hits
            sp_rep_f = multi_sp_rep_field
            if item.get(sp_rep_f) and item.get(sp_rep_f) != '-':
                for sp_rep_hit in item[sp_rep_f].split(';'):
                    sp_rep_id = sp_rep_hit.split(',')[0]
                    sp_rep_id = sp_rep_id.strip()
                    all_sp_reps[sp_rep_id] = True
                    if sp_rep_id not in sp_reps_by_query[this_genome_id]:
                        sp_reps_by_query[this_genome_id].append(sp_rep_id)

    return (all_sp_reps, sp_reps_by_query)


# get_upas_from_set()
#
def get_upas_from_set(top_obj):

    # normalize upa just in case it's a ref vs. an upa
    top_upa = upa_from_info(top_obj['info'])
    top_name = top_obj['info'][NAME_I]
    obj_data = top_obj['data']
    obj_type = top_obj['info'][TYPE_I].split('-')[0]

    if check_obj_type_genome (obj_type):
        if 'KBaseSets.GenomeSet' == obj_type:
            upas = [gsi['ref'] for gsi in obj_data['items']]
        elif 'KBaseSearch.GenomeSet' == obj_type:
            upas = [gse['ref'] for gse in obj_data['elements'].values()]
        elif 'KBaseGenomes.Genome' == obj_type:
            upas = [top_upa]
        else:
            raise ValueError(f'{obj_type} type is not supported')

    return upas


# get_query_assembly_to_genome_name ()
#
def get_query_assembly_to_genome_name (top_obj, clients):
    query_assembly_to_genome_name = dict()

    obj_type = top_obj['info'][TYPE_I].split('-')[0]
    
    if check_obj_type_genome (obj_type):
        upas = get_upas_from_set (top_obj)

        for query_upa in upas:
            (this_wsid_str, this_objid_str, this_ver) = query_upa.split('/')
            query_genome_obj = clients.ws().get_objects2({'objects': [{'wsid':int(this_wsid_str),'objid':int(this_objid_str)}]})['data'][0]
            genome_name = query_genome_obj['info'][NAME_I]
            genome_upa = upa_from_info (query_genome_obj['info'])
            assembly_upa = query_genome_obj['data'].get('contigset_ref') or query_genome_obj['data'].get('assembly_ref')
            assembly_info = clients.ws().get_object_info3({'objects': [{'ref':assembly_upa}]})['infos'][0]
            assembly_name = assembly_info[NAME_I]
            query_assembly_to_genome_name[assembly_name] = genome_name

    return query_assembly_to_genome_name


# _copy_genome_obj()
#
def _copy_genome_obj(src_upa,
                     dst_ws_id,
                     dst_obj_name,
                     clients):

    (src_wsid_str, src_objid_str, src_ver) = src_upa.split('/')
    
    dst_genome_obj_info = clients.ws().copy_object({'from':{'wsid':int(src_wsid_str), 'objid':int(src_objid_str)}, 'to':{'wsid':dst_ws_id,'name':dst_obj_name}})

    return dst_genome_obj_info


# copy_gtdb_species_reps()
#
def copy_gtdb_species_reps (primary_wsid, top_upa, genome_upas_map_file, summary_tables, clients):
    new_objects_created = []

    # get upas by genome id
    genome_id_to_upa_map = get_genome_to_upa_map(genome_upas_map_file)

    # get query assembly name to genome name mapping
    #  Note: redundant calls: could make global map earlier
    #
    top_obj = clients.dfu().get_objects({'object_refs': [top_upa]})['data'][0]
    top_name = top_obj['info'][NAME_I]
    upas = get_upas_from_set (top_obj)
    query_assembly_to_genome_name = get_query_assembly_to_genome_name (top_obj, clients)


    # get species reps by query id
    (all_sp_reps, sp_reps_by_query) = get_sp_rep_hits(summary_tables, query_assembly_to_genome_name)

    
    # copy over genome objs
    sp_rep_dst_upa = dict()
    for sp_rep_id in sorted(all_sp_reps.keys()):
        sp_rep_src_upa = genome_id_to_upa_map[sp_rep_id]
        (src_wsid_str, src_objid_str, src_ver) = sp_rep_src_upa.split('/')

        dst_obj_name = 'GTDB_SP_REP-'+sp_rep_id+'.Genome'
        found_free_name = False
        try:
            existing_obj = clients.ws().get_object_info3({'objects':[{'wsid':primary_wsid,'name':dst_obj_name}]})
            for extra_char in range(10):
                dst_obj_name = 'GTDB_SP_REP-'+sp_rep_id+'-'+str(extra_char)+'.Genome'
                try:
                    existing_obj = clients.ws().get_object_info3({'objects':[{'wsid':primary_wsid,'name':dst_obj_name}]})
                except:
                    found_free_name = True
                    break
        except:
            found_free_name = True

        if not found_free_name:
            raise ValueError ("unable to find available object name for GTDB Species Rep {}".format(sp_rep_id))

        # copy the object, using newest version (or does copy_object barf on assembly handle too?)
        dst_genome_obj_info =  _copy_genome_obj(sp_rep_src_upa,
                                                primary_wsid,
                                                dst_obj_name,
                                                clients)
        #dst_genome_obj_info = clients.ws().copy_object({'from':{'wsid':int(src_wsid_str), 'objid':int(src_objid_str)}, 'to':{'wsid':primary_wsid,'name':dst_obj_name}})

        sp_rep_dst_upa[sp_rep_id] = upa_from_info (dst_genome_obj_info)


    # create new genomesets with query and gtdb proximal genomes
    #
    all_genomeset_elements = dict()
    per_query_genomeset_elements = dict()
    for query_upa in upas:
        (this_wsid_str, this_objid_str, this_ver) = query_upa.split('/')
        genome_obj_info = clients.ws().get_object_info3({'objects': [{'wsid':int(this_wsid_str),'objid':int(this_objid_str)}]})['infos'][0]
        genome_name = genome_obj_info[NAME_I]
        genome_upa = upa_from_info (genome_obj_info)
        #query_names[genome_upa] = genome_name
        #query_upas[genome_name] = genome_upa
            
        all_genomeset_elements[genome_name] = {'ref':genome_upa}
        per_query_genomeset_elements[genome_name] = dict()
        per_query_genomeset_elements[genome_name][genome_name] = {'ref':genome_upa}

        if genome_name not in sp_reps_by_query:
            raise ValueError ("missing sp_rep_hts for query {}".format(genome_name))
            
        #print ("LOOKING FOR SP REP LIST for query {}".format(genome_name))  # DEBUG

        for sp_rep_id in sp_reps_by_query[genome_name]:
            all_genomeset_elements[sp_rep_id] = {'ref':sp_rep_dst_upa[sp_rep_id]}
            per_query_genomeset_elements[genome_name][sp_rep_id] = {'ref':sp_rep_dst_upa[sp_rep_id]}
            
    # build genomesets and save them, adding to created objects
    for query_name in sorted(per_query_genomeset_elements.keys()):
        genomeset_desc = 'Proximal GTDB species reps for '+query_name
        genomeset_name = 'GTDB_SP_REPS.'+query_name+'.GenomeSet'
        genomeset_type = 'KBaseSearch.GenomeSet'
        genomeset_data = {'description': genomeset_desc,
                          'elements': per_query_genomeset_elements[query_name]
                         }
        updated_obj_info = clients.dfu().save_objects(
            { 'id': primary_wsid,
              'objects': [{ 'type': genomeset_type,
                            'name': genomeset_name,
                            'data': genomeset_data
              }]})[0]
        new_ref = upa_from_info(updated_obj_info)
        obj_desc = genomeset_desc
        new_objects_created.append({'ref': new_ref, 'description': obj_desc})

    # make genomeset with all genome objs
    genomeset_desc = 'Proximal GTDB species reps for ALL query genomes in '+top_name
    genomeset_name = 'GTDB_SP_REPS-ALL.'+top_name+'.GenomeSet'
    genomeset_type = 'KBaseSearch.GenomeSet'
    genomeset_data = {'description': genomeset_desc,
                      'elements': all_genomeset_elements
                     }
    updated_obj_info = clients.dfu().save_objects(
        { 'id': primary_wsid,
          'objects': [{ 'type': genomeset_type,
                        'name': genomeset_name,
                        'data': genomeset_data
          }]})[0]
    new_ref = upa_from_info(updated_obj_info)
    obj_desc = genomeset_desc
    new_objects_created.append({'ref': new_ref, 'description': obj_desc})            

    return new_objects_created


# fix_unowned_shock_handles ()
#
def fix_unowned_shock_handles (obj_type, obj_data, clients):
    handle_fields = { 'genome': ['genbank_handle_ref', 'gff_handle_ref'],
                      'assembly': ['reads_handle_ref', 'fasta_handle_ref']
                    }

    for h_field in handle_fields[obj_type]:
        if obj_data.get(h_field):
            h_id = obj_data[h_field]
            s_id = clients.hs().hids_to_handles ([h_id])[0]['id']
            own_node_output = clients.dfu().own_shock_node({'shock_id': s_id, 'make_handle': 1})
            new_s_id = own_node_output['shock_id']
            new_h_id = own_node_output['handle']['hid']
            obj_data[h_field] = new_h_id

    return obj_data


# upa_from_info ()
#
def upa_from_info (obj_info):
    return '/'.join([str(obj_info[WSID_I]),
                     str(obj_info[OBJID_I]),
                     str(obj_info[VERSION_I])])
