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
    parser.add_argument("-l", "--leaflist", help="file with list of leaves to retain")
    
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

        (old_leaf_id, leaf_name) = line.split("\t")

        base_leaf_id = leaf_name.split(' ')[0].lstrip('"')
        target_leaves[base_leaf_id] = leaf_name  # note: id has been replaced by trim program
    f.close()

    return target_leaves


# get_tree_obj_from_file ()
#
def get_tree_obj_from_file (tree_file):
    
    print ("reading tree from file {} ...".format(tree_file))

    tree = ete3.Tree(tree_file, quoted_node_names=True, format=1)

    return tree


# get_taxon_colors()
#
def get_taxon_colors (in_tree_path, outimgbase):

    # taxon colors bin
    get_taxon_colors_bin = os.path.join ('/kb', 'module', 'bin', 'get_taxon_colors.py')

    # configure taxon box color
    out_color_path = outimgbase + '-'+'taxon_colors.map'
    get_taxon_colors_cmd = [get_taxon_colors_bin,
                           '--intree', in_tree_path,
                           '--outcolorfile', out_color_path
                        ]
    env = dict(os.environ)
    subprocess.run(get_taxon_colors_cmd, check=True, env=env)
    with open (out_color_path,'r') as taxon_colors_h:
        for line in taxon_colors_h:
            line = line.rstrip()
            (taxon, hex_color) = line.split("\t")
            taxon_color[taxon] = hex_color

            
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
def write_tree_img_to_files (t, title_disp, outimgbase, query_target_leaves, target_leaves):
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

        
    # set taxon box to leaves of this node
    tax_column_I = { 'g': 2, 'f': 3, 'o': 4, 'c': 5, 'p': 6}
    for n in t.traverse():
        if not n.is_leaf():
            if n.name:
                taxon = n.name
                print ("adding tax ring for: {}".format(taxon))
                tax_level = taxon[0]
                if tax_level not in tax_column_I:  # mostly 'd' but may also 's' someday
                    continue
                taxon_column[taxon] = tax_column_I[tax_level]

                for leaf_node in n.get_leaves():
                    leaf_name = leaf_node.name
                    #print ("ADDING TAX PLOT TO LEAF {}".format(leaf_name))  # DEBUG
                    if leaf_name not in leaf_taxa:
                        leaf_taxa[leaf_name] = []
                    # may have 2 leaves with same name
                    if taxon not in leaf_taxa[leaf_name]:
                        leaf_taxa[leaf_name].append(taxon)

    ts.layout_fn = tax_box_layout

                                           
    # save rectangle tree image as PNG and PDF
    shape = 'rectangle'
    out_png_file_path = outimgbase + '-'+shape+'.PNG'
    out_pdf_file_path = outimgbase + '-'+shape+'.PDF'
    out_files[shape] = { 'png': out_png_file_path,
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
    shape = 'circle'
    out_png_file_path = outimgbase + '-'+shape+'.PNG'
    out_pdf_file_path = outimgbase + '-'+shape+'.PDF'
    out_files[shape] = { 'png': out_png_file_path,
                         'pdf': out_pdf_file_path
                         }
    ts.mode = "c"  # circular tree graph <-- THIS IS THE KEY
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

    # read full tree input
    tree = get_tree_obj_from_file (args.intree)

    # get taxon colors for nodes in tree
    get_taxon_colors (args.intree, args.outimgbase)

    # write tree img to files
    out_files = write_tree_img_to_files (tree, args.title, args.outimgbase, query_target_leaves, target_leaves)
    
    print ("DONE")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())
