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
import copy
import subprocess

import ete3


# these need to be global for layout fxn
leaf_taxa = dict()
taxon_column = dict()
taxon_color = dict()


# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="trim tree to target leaves")

    parser.add_argument("-i", "--intree", help="tree to trim in newick format")
    parser.add_argument("-t", "--title", help="title for tree images")
    parser.add_argument("-o", "--outimgbase", help="image file base (will make circle and rectangle PNG and PDF)")
    parser.add_argument("-q", "--queryleaflist", help="query file with list of query leaves")
    parser.add_argument("-l", "--leaflist", help="file with list of target leaves )not all)")
    parser.add_argument("-g", "--gtdblineagefile", help="file with lineage for all leaves")
    
    args = parser.parse_args()

    args_pass = True

    if len(sys.argv) < 9:
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
    if args.outimgbase is None:
        print ("must specify --{}\n".format('outimgbase'))
        args_pass = False
    if args.queryleaflist is None:
        print ("must specify --{}\n".format('queryleaflist'))
        args_pass = False
    elif not os.path.exists(args.queryleaflist) or \
         not os.path.isfile(args.queryleaflist) or \
         not os.path.getsize(args.queryleaflist) > 0:
        print ("--{} {} must exist and not be empty\n".format('queryleaflist', args.queryleaflist))
        args_pass = False
    if args.leaflist is None:
        print ("must specify --{}\n".format('leaflist'))
        args_pass = False
    elif not os.path.exists(args.leaflist) or \
         not os.path.isfile(args.leaflist) or \
         not os.path.getsize(args.leaflist) > 0:
        print ("--{} {} must exist and not be empty\n".format('leaflist', args.leaflist))
        args_pass = False
    if args.gtdblineagefile is None:
        print ("must specify --{}\n".format('gtdblineagefile'))
        args_pass = False
    elif not os.path.exists(args.gtdblineagefile) or \
         not os.path.isfile(args.gtdblineagefile) or \
         not os.path.getsize(args.gtdblineagefile) > 0:
        print ("--{} {} must exist and not be empty\n".format('gtdblineagefile', args.gtdblineagefile))
        args_pass = False

    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# get_target_leaves ()
#
def get_target_leaves (leaflist_file):

    print ("reading target leaves from file {} ...".format(leaflist_file))

    target_leaves = dict()
    
    if leaflist_file.lower().endswith('.gz'):
        f = gzip.open(leaflist_file, 'rt')
    else:
        f = open(leaflist_file, 'r')

    for line in f:
        line = line.rstrip()

        leaf_info = line.split("\t")
        if len(leaf_info) == 2:
            (leaf_id, leaf_name) = leaf_info
        else:
            (leaf_id, leaf_name, this_lineage) = leaf_info
            
        base_leaf_id = leaf_name.split(' ')[0].lstrip('"')
        target_leaves[base_leaf_id] = leaf_name  # note: id has been replaced by trim program
    f.close()

    return target_leaves


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


# get_tree_obj_from_file ()
#
def get_tree_obj_from_file (tree_file):
    
    print ("reading tree from file {} ...".format(tree_file))

    tree = ete3.Tree(tree_file, quoted_node_names=True, format=1)

    return tree


# get_taxon_colors()
#
def get_taxon_colors (in_tree_path, lineage_file, outimgbase):

    # taxon colors bin
    get_taxon_colors_bin = os.path.join ('/kb', 'module', 'bin', 'get_taxon_colors.py')

    # configure taxon box color
    out_color_path = outimgbase + '-'+'taxon_colors.map'
    get_taxon_colors_cmd = [get_taxon_colors_bin,
                            '--intree', in_tree_path,
                            '--gtdblineagefile', lineage_file,
                            '--outcolorfile', out_color_path
                        ]
    env = dict(os.environ)
    subprocess.run(get_taxon_colors_cmd, check=True, env=env)
    with open (out_color_path,'r') as taxon_colors_h:
        for line in taxon_colors_h:
            line = line.rstrip()
            (taxon, hex_color) = line.split("\t")
            taxon_color[taxon] = hex_color


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

            
# Creates my own layout function.
#
def tax_box_layout(node):

    if node.is_leaf():
        leaf_name = node.name

        # add taxon boxes
        if leaf_name in leaf_taxa:
            for taxon in leaf_taxa[leaf_name]:

                this_column = taxon_column[taxon]
                this_color = '#'+taxon_color[taxon]
                size = 35
                this_box_face = ete3.faces.RectFace(size, size, this_color, this_color)
                this_box_face.margin_top = 2
                this_box_face.margin_bottom = 2
                this_box_face.margin_left = 2
                this_box_face.margin_right = 2
                this_box_face.border.margin = 1
                ete3.faces.add_face_to_node(this_box_face, node, column=this_column, aligned=True)
                

# write_tree_img_to_files()
#
def write_tree_img_to_files (t,
                             title_disp,
                             outimgbase,
                             query_target_leaves,
                             target_leaves,
                             all_leaf_lineages):
    out_files = dict()

    t.ladderize()
    ts = ete3.TreeStyle()

    # customize
    ts.show_leaf_name = True
    ts.show_branch_length = False
    ts.show_branch_support = False  # trimming killed support vals
    #ts.scale = 50 # 50 pixels per branch length unit
    ts.branch_vertical_margin = 5  # pixels between adjacent branches
    #title_disp = intree_name
    if title_disp is None:
        title_disp = os.path.basename (outimgbase)
    ts.title.add_face(ete3.TextFace(title_disp, fsize=30, fgcolor="#606060"), column=0)
    
    # color targets
    leaf_colors = dict()
    leaf_colors['white']     = "#FFFFFF"
    leaf_colors['mustard']   = "#FEE787"
    leaf_colors['violet']    = "#DFCFFC"
    leaf_colors['lightblue'] = "#A8EAFC"
    leaf_colors['darkblue']  = "#A8C1FC"
    default_user_genome_color = leaf_colors['mustard']
    default_reference_genome_color = leaf_colors['violet']

    node_style = ete3.NodeStyle()
    node_style["fgcolor"] = "#606060"  # for node balls
    #node_style["size"] = 10  # for node balls (gets reset based on support)
    node_style["size"] = 0  # support vals are meaningless after trimming and all 1.0
    node_style["vt_line_color"] = "#606060"
    node_style["hz_line_color"] = "#606060"
    node_style["vt_line_width"] = 2
    node_style["hz_line_width"] = 2
    node_style["vt_line_type"] = 0  # 0 solid, 1 dashed, 2 dotted
    node_style["hz_line_type"] = 0
    
    leaf_style = ete3.NodeStyle()
    leaf_style["fgcolor"] = "#ffffff"  # for node balls
    leaf_style["size"] = 3  # for node balls (we're using it to add space)
    leaf_style["vt_line_color"] = "#606060"  # unecessary
    leaf_style["hz_line_color"] = "#606060"
    leaf_style["vt_line_width"] = 2
    leaf_style["hz_line_width"] = 2
    leaf_style["vt_line_type"] = 0  # 0 solid, 1 dashed, 2 dotted
    leaf_style["hz_line_type"] = 0
    
    # set style
    for node in t.traverse():
        if not node.is_leaf():
            style = copy.copy(node_style)

        else:
            leaf_name = node.name

            style = copy.copy(leaf_style)
            node_id = leaf_name.split(' ')[0].lstrip('"')
            if node_id in query_target_leaves:
                style["bgcolor"] = default_user_genome_color
                style["fgcolor"] = default_user_genome_color  # for spacer ball
            elif node_id in target_leaves:
                style["bgcolor"] = default_reference_genome_color
                style["fgcolor"] = default_reference_genome_color  # for spacer ball
                
        node.set_style(style)


    # get complete lineages for taxa of interest
    full_parent_lineages = get_parents (all_leaf_lineages)

    
    # set taxon box to leaves of this node
    tax_column_I = { 'g': 2, 'f': 3, 'o': 4, 'c': 5, 'p': 6}
    taxon_added = dict()
    for n in t.traverse():
        if not n.is_leaf():
            if n.name:
                taxon = n.name
                taxon_added[taxon] = True

                # add lineage up to phylum
                parent_taxa = []
                if taxon in full_parent_lineages:
                    parent_taxon = full_parent_lineages[taxon]
                    if parent_taxon not in taxon_added:
                        parent_taxa.append(parent_taxon)
                        taxon_added[parent_taxon] = True
                    while parent_taxon in full_parent_lineages:
                        parent_taxon = full_parent_lineages[parent_taxon]
                        if parent_taxon not in taxon_added:
                            parent_taxa.append(parent_taxon)
                            taxon_added[parent_taxon] = True

                # decorate leaves below this node
                for this_taxon in [taxon] + parent_taxa:
                    print ("adding tax ring for: {}".format(this_taxon))
                    tax_level = this_taxon[0]
                    if tax_level not in tax_column_I:  # mostly 'd' but maybe also 's' someday
                        continue
                    taxon_column[this_taxon] = tax_column_I[tax_level]

                    for leaf_name in n.get_leaf_names():
                        #print ("ADDING TAX PLOT TO LEAF {}".format(leaf_name))  # DEBUG
                        if leaf_name not in leaf_taxa:
                            leaf_taxa[leaf_name] = []
                        # may have 2 leaves with same name
                        if this_taxon not in leaf_taxa[leaf_name]:
                            leaf_taxa[leaf_name].append(this_taxon)

    # add taxa to non-target leaves that match expanded rings from lower tax levels
    for leaf_name in t.get_leaf_names():
        leaf_id = re.sub(' .*$', '', leaf_name.strip('"'))
        if leaf_id not in all_leaf_lineages:
            raise ValueError ("missing lineage for leaf id {}".format(leaf_id))
        #print ("LEAF_ID: {}".format(leaf_id))  # DEBUG
        for this_taxon in all_leaf_lineages[leaf_id].split(';'):
            tax_level = this_taxon[0]
            if tax_level not in tax_column_I:  # mostly 'd' but maybe also 's' someday
                continue
            if this_taxon in taxon_added:
                if leaf_name not in leaf_taxa:
                    leaf_taxa[leaf_name] = []
                    print ("ADDING ADDITIONAL TAXA TO LEAF_NAME: {}".format(leaf_name))  # DEBUG
                # may have 2 leaves with same name                                       
                if this_taxon not in leaf_taxa[leaf_name]:
                    leaf_taxa[leaf_name].append(this_taxon)
                    print ("ADDING TAXON {} TO LEAF_NAME: {}".format(this_taxon, leaf_name))  # DEBUG

    # key step to make the taxon ring boxes get plotted
    ts.layout_fn = tax_box_layout

                                           
    # save rectangle tree image as PNG and PDF
    shape = 'rectangle'
    um_term = ''
    out_png_file_path = outimgbase + '-'+shape+um_term+'.PNG'
    out_pdf_file_path = outimgbase + '-'+shape+um_term+'.PDF'
    out_files[shape+um_term] = { 'png': out_png_file_path,
                                 'pdf': out_pdf_file_path
                               }
    dpi = 600
    img_units = "in"
    img_pix_width = 2400
    img_in_width = round(float(img_pix_width) / float(dpi), 1)
    img_html_width = img_pix_width // 2
    print ("writing tree img to file {} ...".format(out_png_file_path))
    t.render(out_png_file_path, w=img_in_width, units=img_units, dpi=dpi, tree_style=ts)
    print ("writing tree img to file {} ...".format(out_pdf_file_path))
    t.render(out_pdf_file_path, w=img_in_width, units=img_units, tree_style=ts)  # dpi irrelevant


    # save circle tree image as PNG and PDF
    ts.mode = "c"  # circular tree graph <-- THIS IS THE KEY
    shape = 'circle'
    um_term = ''
    out_png_file_path = outimgbase + '-'+shape+um_term+'.PNG'
    out_pdf_file_path = outimgbase + '-'+shape+um_term+'.PDF'
    out_files[shape+um_term] = { 'png': out_png_file_path,
                                 'pdf': out_pdf_file_path
                               }
    #ts.arc_start = -180 # 0 degrees = 3 o'clock
    #ts.arc_span = 180
    print ("writing tree img to file {} ...".format(out_png_file_path))
    t.render(out_png_file_path, w=img_in_width, units=img_units, dpi=dpi, tree_style=ts)
    print ("writing tree img to file {} ...".format(out_pdf_file_path))
    t.render(out_pdf_file_path, w=img_in_width, units=img_units, tree_style=ts)  # dpi irrelevant

    
    # save ultrametric circle tree image as PNG and PDF
    ts.mode = "c"  # circular tree graph <-- THIS IS THE KEY
    t.convert_to_ultrametric()  # make a dendrogram <-- THIS IS THE KEY 
    shape = 'circle'
    um_term = '-ultrametric'
    out_png_file_path = outimgbase + '-'+shape+um_term+'.PNG'
    out_pdf_file_path = outimgbase + '-'+shape+um_term+'.PDF'
    out_files[shape+um_term] = { 'png': out_png_file_path,
                                 'pdf': out_pdf_file_path
                               }
    #ts.arc_start = -180 # 0 degrees = 3 o'clock
    #ts.arc_span = 180
    print ("writing tree img to file {} ...".format(out_png_file_path))
    t.render(out_png_file_path, w=img_in_width, units=img_units, dpi=dpi, tree_style=ts)
    print ("writing tree img to file {} ...".format(out_pdf_file_path))
    t.render(out_pdf_file_path, w=img_in_width, units=img_units, tree_style=ts)  # dpi irrelevant
    

    return out_files


# main()
#
def main() -> int:
    args = getargs()

    # read query target leaves
    query_target_leaves = get_target_leaves (args.queryleaflist)

    # read all target leaves
    target_leaves = get_target_leaves (args.leaflist)

    # read all target lineages
    all_leaf_lineages = get_all_leaf_lineages (args.gtdblineagefile)

    # read full tree input
    tree = get_tree_obj_from_file (args.intree)

    # get taxon colors for nodes in tree
    get_taxon_colors (args.intree, args.gtdblineagefile, args.outimgbase)

    # write tree img to files
    out_files = write_tree_img_to_files (tree,
                                         args.title,
                                         args.outimgbase,
                                         query_target_leaves,
                                         target_leaves,
                                         all_leaf_lineages)
    
    print ("DONE")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())
