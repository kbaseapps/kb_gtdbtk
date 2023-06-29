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
    parser.add_argument("-g", "--gtdblineagefile", help="file with lineage for leaves")
    parser.add_argument("-o", "--outcolorfile", help="output taxon color file")
    
    args = parser.parse_args()

    args_pass = True

    if len(sys.argv) < 4:
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


# get_tree_obj_from_file ()
#
def get_tree_obj_from_file (tree_file):
    
    #print ("reading tree from file {} ...".format(tree_file))

    tree = ete3.Tree(tree_file, quoted_node_names=True, format=1)

    return tree


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


# get_complete_tax_lineages ()
#
def get_complete_tax_lineages (tax_lineages, full_parent_lineages, tax_levels_supported):
    complete_tax_lineages = dict()
    
    for tax_level in reversed(tax_levels_supported):  # reversed order matters!  want g -> p
        if tax_level not in tax_lineages:
            continue
        for child_taxon in tax_lineages[tax_level].keys():
            child_tax_level = child_taxon[0]
            if child_tax_level not in complete_tax_lineages:
                complete_tax_lineages[child_tax_level] = dict()
            if child_taxon not in complete_tax_lineages[child_tax_level]:
                complete_tax_lineages[child_tax_level][child_taxon] = dict()

            # climb
            this_child_taxon = child_taxon
            while this_child_taxon in full_parent_lineages:
                parent_taxon = full_parent_lineages[this_child_taxon]
                parent_tax_level = parent_taxon[0]
                if parent_tax_level not in complete_tax_lineages:
                    complete_tax_lineages[parent_tax_level] = dict()
                if parent_taxon not in complete_tax_lineages[parent_tax_level]:
                    complete_tax_lineages[parent_tax_level][parent_taxon] = dict()
                complete_tax_lineages[parent_tax_level][parent_taxon][this_child_taxon] = True

                this_child_taxon = parent_taxon
                
    return complete_tax_lineages


# get_all_taxon_colors()
#
def get_all_taxon_colors (t, all_leaf_lineages, out_file=None):
    out_files = dict()

    t.ladderize()
        
    # color options
    base_colors = [
        [   # blues
            ('DodgerBlue', '1E90FF'),
            ('MediumBlue', '0000CD'),
            ('DeepSkyBlue', '00BFFF'),
            ('MediumSlateBlue', '7B68EE'),
            ('CornflowerBlue', '6495ED'),
        ],
        [   # reds
            ('Crimson', 'DC143C'),
            ('PaleVioletRed', 'DB7093'),
            ('LightSalmon', 'FFA07A'),
            ('Tomato', 'FF6347'),
            ('FireBrick', 'B22222'),
            ('DarkRed', '8B0000'),
            ('DeepPink', 'FF1493'),
            ('PaleVioletRed', 'DB8093'),
            ('Salmon', 'FA8072')
        ],
        [   # purples / violets
            ('DarkViolet', '9400D3'),
            ('Magenta', 'FF00FF'),
            ('MediumVioletRed', 'C71585'),
            ('Violet', 'EE82EE'),
            ('Purple', '800080'),
            ('Indigo', '4B0082'),
            ('Plum', 'DDA0DD'),
        ],
        [   # oranges / yellows
            ('OrangeRed', 'FF4500'),
            ('LightCoral', 'F08080'),
            ('DarkOrange', 'FF8C00'),
            ('Coral', 'FF7F50'),
            ('Gold', 'FFD700'),
            ('Moccasin', 'FFE4B5'),
            ('Yellow', 'FFFF00'),
        ]
        #[   # greens  (skip to support red/green color blindness
        #    ('SeaGreen', '2E8B57'),
        #    ('MediumAquaMarine', '66CDAA'),
        #    ('Teal', '008080'),
        #    ('Green', '008000'),
        #    ('Lime', '00FF00'),
        #    ('DarkGreen', '006400'),
        #]
    ]

    tax_levels_supported = [ 'p', 'c', 'o', 'f', 'g']  # skipping domain and species
    
    # get tax hierarchies we care about from tree
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
                    
    # fill in any parent gaps in lineage
    full_parent_lineages = get_parents (all_leaf_lineages)
    complete_tax_lineages = get_complete_tax_lineages (tax_lineages,
                                                       full_parent_lineages,
                                                       tax_levels_supported)
    tax_lineages = complete_tax_lineages


    # get colors
    out_buf = []
    num_levels = 5
    taxon_color = dict()
    top_taxon_color = dict()
    base_colors_used = []
    for base_color_i in range(len(base_colors)):
        base_colors_used.append ({})
    hex_to_val = {'A':10, 'B':11, 'C':12, 'D':13, 'E':14, 'F':15}
    val_to_hex = {10:'A', 11:'B', 12:'C', 13:'D', 14:'E', 15:'F'}
    for i in range(10):
        hex_to_val[str(i)] = i
        val_to_hex[i] = str(i)

    base_color_i = 0
    color_shift_map = dict()
    for tax_level in tax_levels_supported:  # order matters!
        if tax_level not in tax_lineages:
            continue
        for taxon in sorted(tax_lineages[tax_level].keys()):
            if taxon not in taxon_color:

                (taxon_color[taxon], base_colors_used) = get_taxon_color (taxon,
                                                                          base_color_i,
                                                                          base_colors,
                                                                          base_colors_used)
                base_color_i += 1
                base_color_i = base_color_i % len(base_colors)
                top_taxon_color[taxon] = taxon_color[taxon]
                row = "\t".join([taxon, taxon_color[taxon]])
                print (row)
                out_buf.append(row)

            this_top_taxon_color = top_taxon_color[taxon]
            parent_color = taxon_color[taxon]
            
            for child_i,child_taxon in enumerate(sorted(tax_lineages[tax_level][taxon].keys())):
                top_taxon_color[child_taxon] = top_taxon_color[taxon]

                if child_i == 0 and taxon in color_shift_map:
                    val_i = color_shift_map[taxon]
                else:
                    val_i = [0, 0, 0, 0, 0, 0]
                    ci = child_i % 3
                    if child_i < 3:
                        #val_i[2*ci] = 1
                        if child_i == 0:
                            val_i[4] = 1
                        elif child_i == 1:
                            val_i[0] = 1
                        else:
                            val_i[2] = 1
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

                # save val_i vector
                color_shift_map[child_taxon] = val_i
                        
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


# get_taxon_color():
#
def get_taxon_color(taxon, this_base_color_i, base_colors, base_colors_used):
    this_taxon_color = None
    base_colors_unused = []
    for base_color_i in range(len(base_colors)):
        base_colors_unused.append([])
    for base_color_i in range(len(base_colors)):
        for color_tuple in base_colors[base_color_i]:
            hex_color = color_tuple[1]
            if hex_color not in base_colors_used[base_color_i]:
                base_colors_unused[base_color_i].append(hex_color)

    all_base_colors_used = True
    for base_color_i in range(len(base_colors)):
        if len(base_colors_unused[base_color_i]) > 0:
            all_base_colors_used = False
            break
    if all_base_colors_used:
        for base_color_i in range(len(base_colors)):
            base_colors_unused.append([])
            base_colors_used.append({})
        for base_color_i in range(len(base_colors)):
            for base_color_j in range(len(base_colors[base_color_i])):
                base_colors_unused[base_color_i].append(base_colors[base_color_i][base_color_j][1])

    # pick an unused color
    hash_taxon = hash(taxon)
    base_color_i = this_base_color_i
    if len(base_colors_unused[base_color_i]) == 0:
        for alt_i in range(len(base_colors_unused)):
            if len(base_colors_unused[alt_i]) != 0:
                base_color_i = alt_i
                break
    if len(base_colors_unused[base_color_i]) == 1:
        base_color_j = 0
    else:
        base_color_j = hash_taxon % len(base_colors_unused[base_color_i])

    # set the pick
    this_taxon_color = base_colors_unused[base_color_i][base_color_j]
    base_colors_used[base_color_i][this_taxon_color] = True

    return (this_taxon_color, base_colors_used)


# main()
#
def main() -> int:
    args = getargs()

    # read full tree input
    tree = get_tree_obj_from_file (args.intree)

    # read target lineages
    all_leaf_lineages = get_all_leaf_lineages (args.gtdblineagefile)
    
    # write tree img to files
    out_file = get_all_taxon_colors (tree, all_leaf_lineages, args.outcolorfile)
    
    #print ("DONE")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())
