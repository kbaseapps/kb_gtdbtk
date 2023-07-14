#!/opt/conda3/bin/python3
'''
Copyright 2023 Dylan Chivian

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
'''

import sys
import os
import argparse
import gzip
import re
import ete3


# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="trim tree to target leaves")

    parser.add_argument("-i", "--intree", help="tree to trim in newick format")
    parser.add_argument("-o", "--outtree", help="output tree in newick format")
    parser.add_argument("-l", "--leaflist", help="file with list of leaves to retain")
    parser.add_argument("-t", "--targetleafoutfile", help="output file with list of leaves with their new names")
    parser.add_argument("-g", "--gtdblineageoutfile", help="output file with lineage for leaves")
    parser.add_argument("-s", "--sisters", help="retain one leaf per sister lineage from target branches", action='store_true')
    parser.add_argument("-a", "--archaea_metadata_file", help="gtdb metadta for archaea (def: /data/ar53_metadata_r207.tsv)")
    parser.add_argument("-b", "--bacteria_metadata_file", help="gtdb metadta for bacteria (def: /data/bac210_metadata_r207.tsv)")

    args = parser.parse_args()

    args_pass = True

    if len(sys.argv) < 7:
        parser.print_help()
        sys.exit (-1)

    if args.intree is None:
        print ("must specify --{}\n".format('intree'))
        args_pass = False
    elif not os.path.exists(args.intree) or \
         not os.path.isfile(args.intree) or \
         not os.path.getsize(args.intree) > 0:
        print ("--{} {} must exist and not be empty\n".format('intree', args.intree))
        args_pass = False
    if args.outtree is None:
        print ("must specify --{}\n".format('outtree'))
        args_pass = False
    if args.leaflist is None:
        print ("must specify --{}\n".format('leaflist'))
        args_pass = False
    elif not os.path.exists(args.leaflist) or \
         not os.path.isfile(args.leaflist) or \
         not os.path.getsize(args.leaflist) > 0:
        print ("--{} {} must exist and not be empty\n".format('leaflist', args.leaflist))
        args_pass = False

    # defaults
    if args.archaea_metadata_file is None:
        args.archaea_metadata_file = os.path.join ("/data", "ar53_metadata_r207.tsv")
    if args.bacteria_metadata_file is None:
        args.bacteria_metadata_file = os.path.join ("/data", "bac120_metadata_r207.tsv")

        
    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# get_target_leaves ()
#
def get_target_leaves (leaflist_file):

    print ("reading target leaves from file {} ...".format(leaflist_file))

    target_leaves = dict()
    target_lineages = dict()
    
    if leaflist_file.lower().endswith('.gz'):
        f = gzip.open(leaflist_file, 'rt')
    else:
        f = open(leaflist_file, 'r')

    for line in f:
        line = line.rstrip()

        leaf_info = line.split("\t")
        if len(leaf_info) == 2:
            (leaf_id, genome_name) = leaf_info
        else:
            (leaf_id, genome_name, this_lineage) = leaf_info
            target_lineages[genome_name] = this_lineage

        target_leaves[leaf_id] = genome_name
    f.close()

    return (target_leaves, target_lineages)


# write_new_target_leaves ()
#
def write_new_target_leaves (new_target_leaves, leaflist_outfile):

    print ("writing new target leaves to file {} ...".format(leaflist_outfile))

    outbuf = []
    for target_id in sorted(new_target_leaves.keys()):
        #print ("TARGET: {}".format(target_id))  # DEBUG
        outbuf.append("\t".join([target_id, new_target_leaves[target_id]]))
    
    with open (leaflist_outfile, 'w') as leaflist_h:
        leaflist_h.write("\n".join(outbuf)+"\n")

    return leaflist_outfile


# write_gtdb_lineage_outfile ()
#
def write_gtdb_lineage_outfile (tree, taxa, target_lineages, gtdblineage_outfile):

    print ("writing new target lineages to file {} ...".format(gtdblineage_outfile))

    gtdb_lineage_outbuf = []
    for leaf in tree.get_leaves():
        leaf_id = leaf.name.strip('"')
        leaf_id = re.sub(' .*$', '', leaf_id)
        if leaf_id in taxa:
            gtdb_lineage_outbuf.append("\t".join([leaf_id, ";".join(taxa[leaf_id])]))
        elif leaf_id in target_lineages:
            gtdb_lineage_outbuf.append("\t".join([leaf_id, target_lineages[leaf_id]]))
            
    with open (gtdblineage_outfile, 'w') as lineage_h:
        lineage_h.write("\n".join(gtdb_lineage_outbuf)+"\n")
        
    return gtdblineage_outfile


# get_tree_obj_from_file ()
#
def get_tree_obj_from_file (tree_file):
    
    print ("reading tree from file {} ...".format(tree_file))

    """
    newick_buf = ''
    if tree_file.lower().endswith('.gz'):
        f = gzip.open(tree_file, 'rt')
    else:
        f = open(tree_file, 'r')
    for line in f:
        line = line.rstrip()
        newick_buf += line
    f.close()

    tree = ete3.Tree(newick_buf, quoted_node_names=True, format=1)
    """

    tree = ete3.Tree(tree_file, quoted_node_names=True, format=1)  # BAD: makes support the node name
    #tree = ete3.Tree(tree_file, quoted_node_names=True, format=2)  # fails

    # GTDB support vals are being read as names.  Fix right now!
    node_i = 0
    for n in tree.traverse():
        if not n.is_leaf():
            n.name = 'node_'+str(node_i)
            node_i += 1

    # if reading in output from this program, remove taxon and extra RS_ or GB_
    for n in tree.get_leaves():
        n.name = re.sub (' - .*$','',n.name)
        n.name = re.sub ('^RS_','',n.name)
        n.name = re.sub ('^GB_','',n.name)
        
    return tree
    

# write_tree_to_file()
#
def write_tree_to_file (tree, out_file):

    print ("writing tree to file {} ...".format(out_file))

    tree.ladderize()
    output_newick_buf = tree.write(features=[])
    if not output_newick_buf.endswith(';'):
        output_newick_buf += ';'

    with open (out_file, 'w') as newick_h:
        newick_h.write(output_newick_buf)

    return out_file


# read_gtdb_metadata()
#
def read_gtdb_metadata (arc_metadata_file, bac_metadata_file):
    print ("reading GTDB metadata files ...")

    gtdb_metadata = dict()
    taxa = dict()

    # sp reps
    read_sp_reps_only = True
    gtdb_sp_rep_genome_I = 15 
    gtdb_genome_rep_id_I = 14 

    # genome id and gtdb tax
    genome_id_I = 0
    gtdb_taxonomy_I = 16 
    
    # gtdb genome genome rep selection scores
    Type_species_I = 17 
    NCBI_type_strain_I = 84 
    NCBI_rep_I = 64 
    NCBI_assembly_level_I = 45 
    CheckM_completeness_I = 2 
    CheckM_contamination_I = 3 
    Genome_type_I = 55 
    Contig_count_I = 10 
    Undetermined_bases_I = 1 
    length_16S_I = 98  # using ssu_length >= 1400 instead of full_length
    #Many_frameshifted_proteins_I = NA  # can't find in metadata

    # extra score terms we're adding for global tree
    genus_rep_I = 19 # + 1,000,000
    #mimag_HQ_I = 40 # not really needed since have 16S and checkm quals

    for metadata_file in [arc_metadata_file, bac_metadata_file]:
        with open (metadata_file, 'r') as mh:
            for line in mh:
                if line.startswith('accession'):
                    continue
                line = line.rstrip()
                row = line.split("\t")

                # store sp reps only
                if read_sp_reps_only and row[gtdb_sp_rep_genome_I] != 't':
                    continue

                genome_id            = row[genome_id_I]
                genome_id            = re.sub('^RS_','',genome_id)
                genome_id            = re.sub('^GB_','',genome_id)
                gtdb_tax             = row[gtdb_taxonomy_I]
                
                genus_rep            = row[genus_rep_I]
            
                type_species         = row[Type_species_I]
                NCBI_type_strain     = row[NCBI_type_strain_I]
                NCBI_rep             = row[NCBI_rep_I]
                NCBI_assembly_level  = row[NCBI_assembly_level_I]
                CheckM_completeness  = row[CheckM_completeness_I]
                CheckM_contamination = row[CheckM_contamination_I]
                genome_type          = row[Genome_type_I]
                contig_count         = row[Contig_count_I]
                undetermined_bases   = row[Undetermined_bases_I]
                length_16S           = row[length_16S_I]
                

                # transform some fields to save space
                if type_species == 'type strain of species':
                    type_species = 't'
                else:
                    type_species = 'f'
                if NCBI_type_strain == 'assembly designated as reftype' or 'assembly from type material':
                    NCBI_type_strain = 't'
                else:
                    NCBI_type_strain = 'f'
                if NCBI_rep == 'representative genome':
                    NCBI_rep = 't'
                else:
                    NCBI_rep = 'f'
                if NCBI_assembly_level == 'Complete Genome':
                    complete_genome = 't'
                else:
                    complete_genome = 'f'
                if genome_type == 'derived from environmental sample' or \
                   genome_type == 'derived from metagenome':
                    genome_type = 'M'
                elif genome_type == 'derived from single cell':
                    genome_type = 'S'
                else:
                    genome_type = 'NA'
                if length_16S == 'none':
                    length_16S = '0'
                
                # store tax
                taxa[genome_id] = gtdb_tax.split(';')

                # store metadata of interest as tuple
                gtdb_metadata[genome_id] = (genus_rep,
                                            type_species,
                                            NCBI_type_strain,
                                            NCBI_rep,
                                            complete_genome,
                                            float(CheckM_completeness),
                                            float(CheckM_contamination),
                                            genome_type,
                                            int(contig_count),
                                            int(undetermined_bases),
                                            int(length_16S)
                                           )

    return (taxa, gtdb_metadata)
            

# get_gtdb_rep_scores()
#
def get_gtdb_rep_scores (genome_ids, gtdb_metadata):
    #print ("scoring genomes with GTDB preferences ...")

    gtdb_rep_scores = dict()

    # GTDB pref score (adding genus rep +1,000,000
    #  from: https://gtdb.ecogenomic.org/methods#updating-gtdb-species-representatives
    """
    Type species of genome	                       100,000
    Effective type strain of species according to NCBI	10,000
    NCBI representative of species	                 1,000
    Complete genome	                                   100
    CheckM quality estimate	completeness - 5*contamination
    MAG or SAG	                                          -100
    Contig count	                -5 * (no. contigs/100)
    Undetermined bases	  -5 * (no. undetermined bases/10,000)
    Full length 16S rRNA gene	                            10
    Many frameshifted proteins according to NCBI	   -25
    """

    for genome_id in genome_ids:
        score = 0.0

        (genus_rep,
         type_species,
         NCBI_type_strain,
         NCBI_rep,
         complete_genome,
         CheckM_completeness,
         CheckM_contamination,
         genome_type,
         contig_count,
         undetermined_bases,
         length_16S) = gtdb_metadata[genome_id]

        if genus_rep == 't':
            score += 1000000
        if type_species == 't':
            score += 100000
        if NCBI_type_strain == 't':
            score += 10000
        if NCBI_rep == 't':
            score += 1000
        if complete_genome == 't':
            score += 100
        score += CheckM_completeness - 5 * CheckM_contamination
        if genome_type == 'M' or genome_type == 'S':
            score -= 100
        score -= 5 * float(undetermined_bases) / 10000.0
        if length_16S >= 1400:
            score += 10
        elif length_16S > 0:
            score += 2
        # frameshifts not found in metadata.  maybe I missed it
            
        gtdb_rep_scores[genome_id] = score

    return gtdb_rep_scores


# get_trunk_nodes ()
#
def get_trunk_nodes (tree, target_leaves):
    trunk_nodes = dict()

    trunk_node_cnt = 0
    for n in tree.traverse():
        if n.is_leaf():
            if n.name in target_leaves:
                for ancestor_node in reversed(n.get_ancestors()):
                    if ancestor_node.name is None or not ancestor_node.name.startswith('TRUNK'):
                        # DEBUG
                        #if ancestor_node.is_root():
                        #    print ("found root")
                        if ancestor_node.is_leaf():
                            #print ("how can a leaf be an ancestor?")
                            continue
                        new_name = 'TRUNK_'+str(trunk_node_cnt)
                        # DEBUG
                        #if ancestor_node.name is not None:
                        #    print ("changing trunk node name from {} to {}".format(ancestor_node.name, new_name))
                        ancestor_node.name = new_name
                        trunk_nodes[new_name] = new_name
                        trunk_node_cnt += 1
            
    return trunk_nodes


# get_sister_leaves ()
#
def get_sister_leaves (tree, target_leaves, trunk_nodes, gtdb_metadata):
    these_sister_leaves = dict()

    for child in tree.get_children():
        if child.is_leaf():
            if child.name not in target_leaves:
                leaf_name = child.name
                print ("adding child leaf as sister {}".format(leaf_name))
                these_sister_leaves[leaf_name] = leaf_name
        elif child.name is not None and child.name in trunk_nodes:
            new_sister_leaves = get_sister_leaves (child,
                                                   target_leaves,
                                                   trunk_nodes,
                                                   gtdb_metadata)
            for leaf_name in new_sister_leaves:
                these_sister_leaves[leaf_name] = new_sister_leaves[leaf_name]
        else:
            leaf_name = find_sister_branch_rep (child, gtdb_metadata)
            print ("adding sister rep leaf {}".format(leaf_name))
            these_sister_leaves[leaf_name] = leaf_name

    return these_sister_leaves


# find_sister_branch_rep ()
#
def find_sister_branch_rep (branch_node, gtdb_metadata):
    high_score_leaf_name = None

    leaf_names = branch_node.get_leaf_names()
    leaf_scores = get_gtdb_rep_scores (leaf_names, gtdb_metadata)

    high_score = -10000000
    for leaf_name in sorted(leaf_names):
        if leaf_scores[leaf_name] > high_score:
            high_score = leaf_scores[leaf_name]
            high_score_leaf_name = leaf_name

    return high_score_leaf_name
    

# trim_tree_to_target_and_sister_leaves()
#
def trim_tree_to_target_leaves (tree, target_leaves, sister_leaves=None):

    print ("trimming tree ...")
    
    retain_list = []
    for leaf_name in tree.get_leaf_names():
        if leaf_name in target_leaves:
            retain_list.append(leaf_name)
        if sister_leaves is not None and leaf_name in sister_leaves:
            retain_list.append(leaf_name)

    if sister_leaves is not None:
        preserve_len = True
    else:
        preserve_len = False
                
    tree.prune (retain_list,preserve_branch_length=preserve_len)

    return tree


# replace_leaf_id_with_name()
#
def replace_leaf_id_with_name (tree, taxa, target_leaves, sister_leaves=None):
    print ("adjusting leaf names ...")

    new_target_leaves = dict()

    # DEBUG
    #for target_id in sorted(target_leaves.keys()):
    #    print ("TARGET_ID: '{}'".format(target_id))
    
    for leaf in tree.get_leaves():
        target_leaf_orig_name = None
        #print ("LEAF_NAME: '{}'".format(leaf.name))  # DEBUG
        if leaf.name in target_leaves:
            target_leaf_orig_name = leaf.name
            leaf.name = target_leaves[leaf.name]
        if sister_leaves is not None and leaf.name in sister_leaves:
            leaf.name = sister_leaves[leaf.name]

        if leaf.name in taxa:
            leaf.name = '"'+leaf.name+' - '+taxa[leaf.name][-1]+'"'
            
        if target_leaf_orig_name:
            new_target_leaves[target_leaf_orig_name] = leaf.name
                
    return (tree, new_target_leaves)


# add_internal_node_names ()
#
def add_internal_node_names (tree, taxa, target_lineages):

    nodes_to_unname = dict()

    # augment taxa with GTDB-Tk Classify defined lineages for queries
    for short_leaf_name in list(target_lineages.keys()):
        if short_leaf_name not in taxa:
            taxa[short_leaf_name] = target_lineages[short_leaf_name].split(';')

    # DEBUG
    """
    for short_leaf_name in sorted(taxa.keys()):
        print ("TAXA FOR {}".format(short_leaf_name))
        for tax_level_i,taxon in enumerate(taxa[short_leaf_name]):
            indent = ''
            for indent_i in range(tax_level_i):
                indent += "\t"
            print (indent+taxon)
    """
    
    # set lowest common taxon for each internal node
    for n in tree.traverse():
        if n.is_leaf():
            continue

        # don't name nodes without taxonomic lineage
        #  (note: GTDB-Tk has already given query a lineage)
        if n.name:
            unname_this = False
            for child in n.get_children():
                if child.is_leaf():
                    short_leaf_name = re.sub (' .*$', '', child.name)
                    short_leaf_name = short_leaf_name.lstrip('"')
                    if short_leaf_name not in taxa:
                        nodes_to_unname[n.name] = True
                        unname_this = True
                        break
            if unname_this:
                continue

        # trim leaf names to ids in taxa
        gtdb_defined_leaf_names = []
        for leaf_name in n.get_leaf_names():
            short_leaf_name = re.sub (' .*$', '', leaf_name)
            short_leaf_name = short_leaf_name.lstrip('"')
            if short_leaf_name in taxa:
                gtdb_defined_leaf_names.append(short_leaf_name)
                #print ("ADDING {}".format(short_leaf_name))  # DEBUG
                
        # determine common taxon level
        if len(gtdb_defined_leaf_names) > 0:
            common_taxon = None
            for tax_i,taxon_0 in enumerate(taxa[gtdb_defined_leaf_names[0]]):
                #print ("NODE {} TAXON_0: {} from base leaf {}".format(n.name, taxon_0, gtdb_defined_leaf_names[0])) # DEBUG
                if len(taxon_0) == 3:
                    continue
                level_match = True
                for leaf_name in gtdb_defined_leaf_names:
                    #print ("NODE {} TAXON: {} from leaf {}".format(n.name, taxa[leaf_name][tax_i], leaf_name)) # DEBUG
                    if len(taxa[leaf_name][tax_i]) == 3:  # e.g. 'g__'
                        level_match = False
                        break
                    elif taxa[leaf_name][tax_i] != taxon_0:
                        level_match = False
                        break
                if level_match:
                    common_taxon = taxon_0
                    #print ("NODE {} SETTING COMMON TAXON to {}".format(n.name, common_taxon))  # DEBUG
                else:
                    break
            n.name += ' taxon: '+common_taxon
        else:
            #print ("UNNAMING {}".format(n.name))
            nodes_to_unname[n.name] = True
                    
    # remove internal nodes names with same taxon as a parent node
    for leaf in tree.get_leaves():
        for n in leaf.get_ancestors():
            this_taxon = re.sub('^.* taxon: ','',n.name)
            if not n.is_root():
                up_taxon = re.sub('^.* taxon: ','',n.up.name)
                if this_taxon == up_taxon:
                    nodes_to_unname[n.name] = True
    for n in tree.traverse():
        if not n.is_leaf():
            if n.name in nodes_to_unname:
                #print ("NODE TO UNNAME: '{}'".format(n.name))  # DEBUG
                n.name = ''
            else:
                #print ("SETTING TAXON FOR NODE '{}'".format(n.name))  # DEBUG
                n.name = re.sub('^.* taxon: ','',n.name)

    return tree


# main()
#
def main() -> int:
    args = getargs()

    # read target leaves
    (target_leaves, target_lineages) = get_target_leaves (args.leaflist)

    # read full tree input
    tree = get_tree_obj_from_file (args.intree)

    # read gtdb metadata and taxa
    (taxa, gtdb_metadata) = read_gtdb_metadata (args.archaea_metadata_file,
                                                args.bacteria_metadata_file) 

    # determine extra leaves to act as sister reps
    sister_leaves = None
    if (args.sisters):
        trunk_nodes = get_trunk_nodes (tree, target_leaves)
        sister_leaves = get_sister_leaves (tree,
                                           target_leaves,
                                           trunk_nodes,
                                           gtdb_metadata)
        
    # trim to target leaves
    tree = trim_tree_to_target_leaves (tree,
                                       target_leaves,
                                       sister_leaves)

    # fix id names to genome names
    (tree, new_target_leaves) = replace_leaf_id_with_name (tree,
                                                           taxa,
                                                           target_leaves,
                                                           sister_leaves)

    # write new target leaves
    write_new_target_leaves (new_target_leaves,
                             args.targetleafoutfile)
    
    # write gtdb taxa file
    if args.gtdblineageoutfile:
        gtdb_lineage_path = write_gtdb_lineage_outfile (tree,
                                                        taxa,
                                                        target_lineages,
                                                        args.gtdblineageoutfile)
    
    # add internal node names
    tree = add_internal_node_names (tree, taxa, target_lineages)

    # write tree to file
    write_tree_to_file (tree, args.outtree)    
    
    
    print ("DONE")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())
