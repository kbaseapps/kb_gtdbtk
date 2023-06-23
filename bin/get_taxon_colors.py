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

import ete3


# these need to be global for layout fxn
leaf_taxa = dict()
taxon_column = dict()
taxon_color = dict()


# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="get taxon colors, with higher levels controlling colors of lower levels")

    parser.add_argument("-i", "--intree", help="tree to trim in newick format")
    parser.add_argument("-o", "--outcolorfile", help="output taxon color file")
    
    args = parser.parse_args()

    args_pass = True

    if len(sys.argv) < 2:
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
        
    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# get_tree_obj_from_file ()
#
def get_tree_obj_from_file (tree_file):
    
    #print ("reading tree from file {} ...".format(tree_file))

    tree = ete3.Tree(tree_file, quoted_node_names=True, format=1)

    return tree


# get_taxon_colors()
#
def get_taxon_colors (t, out_file):
    out_files = dict()

    t.ladderize()
        
    # configure taxa boxes
    base_colors = [ ('OrangeRed', 'FF4500'),
                    ('DodgerBlue', '1E90FF'),
                    ('DarkViolet', '9400D3'),
                    ('Crimson', 'DC143C'),
                    ('MediumBlue', '0000CD'),
                    ('SeaGreen', '2E8B57'),
                    ('LightCoral', 'F08080'),
                    ('Magenta', 'FF00FF'),
                    ('DarkOrange', 'FF8C00'),
                    ('DeepSkyBlue', '00BFFF'),
                    ('MediumAquaMarine', '66CDAA'),
                    ('MediumVioletRed', 'C71585'),
                    ('Violet', 'EE82EE'),
                    ('Coral', 'FF7F50'),
                    ('Navy', '000080'),
                    ('Teal', '008080'),
                    ('PaleVioletRed', 'DB7093'),
                    ('Gold', 'FFD700'),
                    ('Green', '008000'),
                    ('LightSalmon', 'FFA07A'),
                    ('Purple', '800080'),
                    ('MediumSlateBlue', '7B68EE'),
                    ('Tomato', 'FF6347'),
                    ('CornflowerBlue', '6495ED'),
                    ('Lime', '00FF00'),
                    ('FireBrick', 'B22222'),
                    ('Indigo', '4B0082'),
                    ('Moccasin', 'FFE4B5'),
                    ('DarkGreen', '006400'),
                    ('DarkRed', '8B0000'),
                    ('Plum', 'DDA0DD'),
                    ('Yellow', 'FFFF00'),
                    ('DeepPink', 'FF1493'),
                    ('PaleVioletRed', 'DB8093'),
                    ('Salmon', 'FA8072')
                  ]

    tax_levels_supported = [ 'p', 'c', 'o', 'f', 'g']  # skipping domain and species

    # get tax hierarchies
    tax_lineages = dict()
    parent = dict()
    for leaf in t.get_leaves():
        last_taxon = None

        for ancestor in leaf.get_ancestors():  # must be g -> p so child is last_taxon
            if ancestor.name:
                taxon = ancestor.name
                tax_level = taxon[0]
                if tax_level in tax_levels_supported:  
                    if tax_level not in tax_lineages:
                        tax_lineages[tax_level] = dict()
                    if taxon not in tax_lineages[tax_level]:
                        tax_lineages[tax_level][taxon] = dict()
                    if last_taxon:
                        tax_lineages[tax_level][taxon][last_taxon] = True
                    last_taxon = taxon
                    
    # get colors
    out_buf = []
    num_levels = 5
    taxon_color = dict()
    top_taxon_color = dict()
    base_color_i = -1
    hex_to_val = {'A':10, 'B':11, 'C':12, 'D':13, 'E':14, 'F':15}
    val_to_hex = {10:'A', 11:'B', 12:'C', 13:'D', 14:'E', 15:'F'}
    for i in range(10):
        hex_to_val[str(i)] = i
        val_to_hex[i] = str(i)
    
    for tax_level in tax_levels_supported:  # order matters!
        if tax_level not in tax_lineages:
            continue
        for taxon in sorted(tax_lineages[tax_level].keys()):
            if taxon not in taxon_color:
                base_color_i += 1
                if not base_color_i < len(base_colors):
                    raise ValueError ("Ran out of colors for taxa")
                taxon_color[taxon] = base_colors[base_color_i][1]  # use the hex string
                top_taxon_color[taxon] = taxon_color[taxon]
                row = "\t".join([taxon, taxon_color[taxon]])
                print (row)
                out_buf.append(row)

            this_top_taxon_color = top_taxon_color[taxon]
            parent_color = taxon_color[taxon]
            for child_i,child_taxon in enumerate(sorted(tax_lineages[tax_level][taxon].keys())):
                top_taxon_color[child_taxon] = top_taxon_color[taxon]
                
                val_i = [0, 0, 0, 0, 0, 0]
                ci = child_i % 3
                if child_i < 3:
                    val_i[2*ci] = 1
                elif child_i < 6:
                    val_i[2*ci] = 1
                    if ci == 0 or ci == 1:
                        j = 2
                    else:
                        j = 0
                    val_i[2*j] = 1
                elif child_i < 9:
                    val_i[2*ci] = 1
                    val_i[2*ci+1] = 1
                elif child_i < 12:
                    val_i[2*ci] = 1
                    val_i[2*ci+1] = 1
                    if ci == 0 or ci == 1:
                        j = 2
                    else:
                        j = 0
                    val_i[2*j] = 1
                    val_i[2*j+1] = 1
                elif child_i < 15:
                    val_i[2*ci+1] = 1
                elif child_i < 18:
                    val_i[2*ci+1] = 1
                    if ci == 0 or ci == 1:
                        j = 2
                    else:
                        j = 0
                    val_i[2*j+1] = 1
                else:
                    cj = child_i % 6
                    val_i[cj] = 1
                    for k in range(cj):
                        val_i[k] = 1

                # make new color
                child_color = ''
                for k in range(6):
                    color_hex = parent_color[k]
                    top_taxon_color_hex = this_top_taxon_color[k]
                    if val_i[k]:
                        val = hex_to_val[color_hex]
                        top_taxon_val = hex_to_val[top_taxon_color_hex]
                        if top_taxon_val < 8:
                            color_step = int((16-top_taxon_val)/float(num_levels-1))
                            color_step = 1  if color_step < 1 else color_step
                            val += color_step
                            val = 15  if val > 15 else val
                        else:
                            color_step = int((top_taxon_val+1)/float(num_levels-1))
                            color_step = 1  if color_step < 1 else color_step
                            val -= color_step
                            val = 0  if val < 0 else val
                        color_hex = val_to_hex[val]
                    child_color += color_hex

                # save color for taxon
                taxon_color[child_taxon] = child_color
                row = "\t".join([child_taxon, child_color])
                print (row)
                out_buf.append(row)


    # write color map
    if out_file:
        with open (out_file, 'w') as out_h:
            out_h.write("\n".join(out_buf)+"\n")

    return out_file


# main()
#
def main() -> int:
    args = getargs()

    # read full tree input
    tree = get_tree_obj_from_file (args.intree)

    # write tree img to files
    out_file = get_taxon_colors (tree, args.outcolorfile)
    
    #print ("DONE")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())
