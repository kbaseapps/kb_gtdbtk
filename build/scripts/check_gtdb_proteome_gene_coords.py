#!/usr/bin/python3
'''
Copyright 2022 Dylan Chivian

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
'''

import sys
import os
import errno
import argparse
import gzip
import re
import glob

# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="filter out short reads from fastq")
    parser.add_argument("-i", "--idlistfile", help="genome id list file", required=True)
    parser.add_argument("-a", "--assemblydir", help="assembly directory", required=True)
    parser.add_argument("-p", "--proteomedir", help="proteome directory", required=True)
    parser.add_argument("-n", "--ncbidir", help="ncbi gff+fna directory", required=True)
    parser.add_argument("-v", "--verbose", help="verbose mode", action="store_true")
    args = parser.parse_args()

    return args


# get_genome_ids()
#
def get_genome_ids (idslistfile):
    genome_ids = []
    print ("reading input file {} ...".format(idslistfile))
    if idslistfile.lower().endswith('.gz'):
        f = gzip.open(idslistfile, 'rt')
    else:
        f = open(idslistfile, 'r')

    for line in f:
        this_genome_id = line.rstrip()
        genome_ids.append(this_genome_id)
    f.close()
    print ("genomes found: {}".format(len(genome_ids)))
    return genome_ids


# process_genome()
#
def process_genome (genome_id, assemblydir, proteomedir, ncbidir, verbose=False):
    gtdb_contigs = get_contigs ('gtdb', genome_id, assemblydir, verbose)
    ncbi_contigs = get_contigs ('ncbi', genome_id, ncbidir, verbose)

    merged_contigs = match_contigs (genome_id, ncbi_contigs, gtdb_contigs)

    # adjust contig lens for circular contigs, but not the seqs
    corrected_contigs = revise_circular_contigs (genome_id, merged_contigs)
    
    genes = get_gene_aa_seqs (genome_id, proteomedir, verbose)

    find_overrun_genes (genome_id, genes, corrected_contigs)
    
    return
    

# get_contigs()
#
def get_contigs (source, genome_id, assemblydir, verbose=False):
    contigs = {'lens': dict(), 'seqs':dict()}
    contig_lens = dict()
    contig_seqs = dict()

    prefix = 'GCF'
    if genome_id.startswith('GCA'):
        prefix = 'GCA'
    genome_id_num = genome_id.replace(prefix+'_','')
    id_3 = genome_id_num[0:3]
    id_6 = genome_id_num[3:6]
    id_9 = genome_id_num[6:9]

    if source == 'gtdb':
        assembly_file = os.path.join(assemblydir, prefix, id_3, id_6, id_9, genome_id+'_genomic.fna.gz')
    elif source == 'ncbi':
        ncbi_assembly_base = os.path.join(assemblydir, genome_id)
        try:
            assembly_file = glob.glob(ncbi_assembly_base+'*_genomic.fna.gz')[0]
        except:
            raise ValueError ("no FNA file found for genome_id: {}".format(genome_id))
    else:
        raise ValueError ("bad source for assembly")

    if assembly_file.lower().endswith('.gz'):
        f = gzip.open(assembly_file, 'rt')
    else:
        f = open(assembly_file, 'r')

    last_id = None
    last_seq = ''
    contig_cnt = 0
    assembly_len = 0
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            if last_id is not None:
                last_id = last_id.split('.')[0]
                contig_seqs[last_id] = last_seq
                contig_lens[last_id] = len(last_seq)
                assembly_len += len(last_seq)
            last_seq = ''
            last_id = line.split(' ')[0].lstrip('>')
            contig_cnt += 1
        else:
            last_seq += line
    if last_id is not None:
        last_id = last_id.split('.')[0]
        contig_seqs[last_id] = last_seq
        contig_lens[last_id] = len(last_seq)
        assembly_len += len(last_seq)
        last_seq = ''
    f.close()

    if verbose:
        print ("{} contigs read {} total length {}".format(source, contig_cnt, assembly_len))

    contigs['lens'] = contig_lens
    contigs['seqs'] = contig_seqs
    return contigs
    

# match_contigs()
#
def match_contigs (genome_id, ncbi_contigs, gtdb_contigs):
    contigs = {'lens': dict(), 'seqs':dict()}
    contig_lens = dict()
    contig_seqs = dict()
    overlap_len = 100

    ncbi_contig_ids = sorted(ncbi_contigs['lens'].keys())
    gtdb_contig_ids = sorted(gtdb_contigs['lens'].keys())
    if ncbi_contig_ids != gtdb_contig_ids:
        raise ValueError ("mismatch in contig_ids for genome {}".format(genome_id))

    for contig_id in ncbi_contig_ids:
        if ncbi_contigs['lens'][contig_id] != gtdb_contigs['lens'][contig_id]:
            raise ValueError ("ncbi and gtdb contigs differ for genome {} contig_id {}".format(genome_id, contig_id))

        this_overlap_len = overlap_len
        if gtdb_contigs['lens'][contig_id] < overlap_len:
            this_overlap_len = gtdb_contigs['lens'][contig_id]
        if ncbi_contigs['seqs'][contig_id][0:this_overlap_len-1] != \
           gtdb_contigs['seqs'][contig_id][0:this_overlap_len-1]:
            raise ValueError ("ncbi and gtdb contigs differ at 5prime for genome {} contig_id {}".format(genome_id, contig_id))
        if ncbi_contigs['seqs'][contig_id][-(this_overlap_len+1):-1] != \
           gtdb_contigs['seqs'][contig_id][-(this_overlap_len+1):-1]:
            raise ValueError ("ncbi and gtdb contigs differ at 3prime for genome {} contig_id {}".format(genome_id, contig_id))

    # if no changes made to contigs, just return gtdb
    contigs['lens'] = gtdb_contigs['lens']
    contigs['seqs'] = gtdb_contigs['seqs']
    return contigs


# get_gene_aa_seqs()
#
def get_gene_aa_seqs (genome_id, proteomedir, verbose=False):
    genes = {'coords': dict(), 'seqs': dict()}

    gene_coords = dict()
    gene_aa_seqs = dict()

    db_prefix = 'RS'
    prefix = 'GCF'
    if genome_id.startswith('GCA'):
        db_prefix = 'GB'
        prefix = 'GCA'
    genome_id_num = genome_id.replace(prefix+'_','')
    id_3 = genome_id_num[0:3]
    id_6 = genome_id_num[3:6]
    id_9 = genome_id_num[6:9]

    domain = 'archaea'
    proteome_file = os.path.join(proteomedir, domain, db_prefix+'_'+prefix+'_'+genome_id_num+'_protein.faa.gz')
    if not os.path.exists(proteome_file):
        domain = 'bacteria'
        proteome_file = os.path.join(proteomedir, domain, db_prefix+'_'+prefix+'_'+genome_id_num+'_protein.faa.gz')
        if not os.path.exists(proteome_file):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), proteome_file)

        
    if proteome_file.lower().endswith('.gz'):
        f = gzip.open(proteome_file, 'rt')
    else:
        f = open(proteome_file, 'r')

    last_header = None
    last_seq = ''
    gene_cnt = 0
    total_gene_len = 0
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            if last_header is not None:
                header = last_header.lstrip('>').split(' ')
                last_id = header[0]
                last_contig_id = last_id.split('.')[0]
                if last_contig_id not in gene_coords:
                    gene_coords[last_contig_id] = []
                this_gene_coords = {'gene_id': last_id,
                                    'contig_id': last_contig_id,
                                    'beg': int(header[2]),
                                    'end': int(header[4]),
                                    'strand': header[6]}
                gene_coords[last_contig_id].append(this_gene_coords)
                gene_aa_seqs[last_id] = last_seq.rstrip('*')
                total_gene_len += len(last_seq)
            last_seq = ''
            last_header = line
            gene_cnt += 1
        else:
            last_seq += line
    if last_header is not None:
        header = last_header.lstrip('>').split(' ')
        last_id = header[0]
        last_contig_id = last_id.split('.')[0]
        if last_contig_id not in gene_coords:
            gene_coords[last_contig_id] = []
        this_gene_coords = {'gene_id': last_id,
                            'contig_id': last_contig_id,
                            'beg': int(header[2]),
                            'end': int(header[4]),
                            'strand': header[6]}
        gene_coords[last_contig_id].append(this_gene_coords)
        gene_aa_seqs[last_id] = last_seq.rstrip('*')
        total_gene_len += len(last_seq)
        last_seq = ''
    f.close()

    if verbose:
        print ("genes read {} total length {}".format(gene_cnt, total_gene_len))

    genes['coords'] = gene_coords
    genes['aa_seqs'] = gene_aa_seqs
    return genes


# revise_circular_contigs()
#
def revise_circular_contigs (genome_id,
                             contigs,
                             length_thresh=30000,
                             min_len=21,
                             mismatch_thresh=0,
                             local_mismatch_thresh=0,
                             local_window=100):

    ret_contigs = {'lens': dict(), 'seqs': dict()}
    contig_lens = contigs['lens']
    contig_seqs = contigs['seqs']
    
    for contig_id in sorted(contig_lens.keys()):
        contig_len = contig_lens[contig_id]
        contig_seq = contig_seqs[contig_id]

        for i in range(min_len,length_thresh):
            if i > contig_len / 2:
                break;
            if mismatch_thresh == 0:
                #print ("{} <> {}".format(contig_seq[0:i], contig_seq[-(i+1):-1]))
                if contig_seq[0:i+1] == contig_seq[contig_len-(i+1):contig_len]:
                    print ("FOUND ONE! genome_id: {} contig_id: {}  overlap: {}".format(genome_id, contig_id, contig_seq[0:i+1]))
                    contig_lens[contig_id] = contig_len - (i+1)

                    # don't adjust seq
                    #contig_seqs[contig_id] = contig_seq[0:-(i+1)]
                    #if len(contig_seqs[contig_id]) != contig_lens[contig_id]:
                    #    print ("BAD MATH {} != {}".format(len(contig_seqs[contig_id]),contig_lens[contig_id]))
                    #    sys.exit(-3)

                    break
            #else:
                # add logic for mismatch != 0 here

    ret_contigs['lens'] = contig_lens
    ret_contigs['seqs'] = contig_seqs
    return ret_contigs
    

# find_overrun_genes()
#
def find_overrun_genes (genome_id, genes, contigs):
    contig_lens = contigs['lens']
    gene_coords = genes['coords']
    gene_aa_seqs = genes['aa_seqs']
    discard_gene_ids = dict()

    for contig_id in sorted(gene_coords.keys()):

        overrun_i = -1
        for gene_i,gene in enumerate(gene_coords[contig_id]):
            gene_id = gene['gene_id']
            if gene['end'] > contig_lens[contig_id]:
                overrun_i += 1

                print ("CONSIDER genome_id: {} gene_id: {} beg: {} end: {}: strand: {} contig_id: {} contig_len: {}".format(
                    genome_id,
                    gene_id,
                    gene['beg'],
                    gene['end'],
                    gene['strand'],
                    contig_id,
                    contig_lens[contig_id]))

                num_genes = len(gene_coords[contig_id])
                if num_genes > 1:
                    #for pre_i in range(overrun_i,overrun_i+3):
                    for pre_i in range(0,overrun_i+3):
                        if pre_i > gene_i - 1:
                            break
                        if gene_coords[contig_id][pre_i]['strand'] == \
                           gene_coords[contig_id][gene_i]['strand']:

                            pre_id = gene_coords[contig_id][pre_i]['gene_id']

                            pre_len = gene_coords[contig_id][pre_i]['end'] - gene_coords[contig_id][pre_i]['beg'] + 1
                            overrun_len = gene_coords[contig_id][gene_i]['end'] - gene_coords[contig_id][gene_i]['beg'] + 1 
                            if overrun_len < pre_len:
                                compare_len = overrun_len
                            else:
                                compare_len = pre_len
                            compare_len = int(compare_len/3)

                            found_match = False
                            for i in range(0,int(compare_len/2)):
                                if gene_aa_seqs[gene_id][i:compare_len+i] == gene_aa_seqs[pre_id][0:compare_len]:
                                    print ("Found fit")
                                    if overrun_len <= pre_len:
                                        discard_gene_ids[gene_id] = True
                                    else:
                                        discard_gene_ids[pre_id] = True
                                    found_match = True
                                    break
                            if not found_match:
                                print ("ALERT: Sequence mismatch for genes pre:{} overrun:{} pre seq: {} overrun seq: {}".format(pre_id, gene_id, gene_aa_seqs[pre_id][0:compare_len], gene_aa_seqs[gene_id][0:compare_len]))

    # check discards
    for contig_id in sorted(gene_coords.keys()):
        for gene_i,gene in enumerate(gene_coords[contig_id]):
            gene_id = gene['gene_id']
            if gene_id in discard_gene_ids:
                print ("DISCARD genome_id: {} gene_id: {} beg: {} end: {}: strand: {} contig_id: {} contig_len: {}".format(
                    genome_id,
                    gene_id,
                    gene['beg'],
                    gene['end'],
                    gene['strand'],
                    contig_id,
                    contig_lens[contig_id]))

    return discard_gene_ids

    
# main()
#
def main() -> int:
    args = getargs()

    genome_cnt = 0
    for genome_id in get_genome_ids (args.idlistfile):
        genome_cnt += 1
        if genome_cnt % 100 == 0:
            print ("processed {}".format(genome_cnt))
        if args.verbose:
            print ("processing {}".format(genome_id))
        process_genome (genome_id,
                        args.assemblydir,
                        args.proteomedir,
                        args.ncbidir,
                        args.verbose)
                   
    print ("DONE")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())
