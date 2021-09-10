#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_table',
                        required=True,
                        help='CSV table with header and at least 2 columns: gene,locus_tag')
    parser.add_argument('-g', '--gbk_dir',
                        required=True,
                        help='Folder with genbanks.')
    parser.add_argument('-m', '--OG_matrix',
                        required=True,
                        help='matrix with OGs in rows, assemblies in columns and space-separated locus_tags at intersection.')
    parser.add_argument('-o', '--output_dir',
                        required=True,
                        help='output dir with fasta seqs of each OG')
    return parser.parse_args()


def main():
    # parse input
    targets_df = pd.read_csv(input_file, header=0, index_col=None)
    if 'gene' not in targets_df:
        sys.exit('Input table does not have a "gene" column.')
    if 'locus_tag' not in targets_df:
        sys.exit('Imput table does not have a "locus_tag" column.')
    targets_df = targets_df[['gene', 'locus_tag']]
    targets_df = targets_df.set_index('gene', drop=True)

    # parse matrix input
    og_df = pd.read_csv(OG_matrix_file, sep='\t', index_col=0)
    if 'function' in og_df:
        og_df = og_df.drop(columns='function')

    # melt matrix input
    og_df = pd.melt(og_df, var_name='assembly', value_name='locus_tags', ignore_index=False).dropna(how='any')

    # iterate rows and make new df with one locus_tag per row
    all_ogs = []
    all_assemblies = []
    all_locus_tags = []
    for cur_og, row in og_df.iterrows():
        cur_assembly = row.assembly
        cur_locus_tags = row.locus_tags
        for cur_tag in cur_locus_tags.rstrip().split():
            all_ogs.append(cur_og)
            all_assemblies.append(cur_assembly)
            all_locus_tags.append(cur_tag)
    
    og_df = pd.DataFrame({'OG': all_ogs, 'assembly': all_assemblies}, index=all_locus_tags)

    # iterate through input and get dictionary with key=gene, values=[locus_tags in OG]
    res = dict()
    for cur_gene, row in targets_df.iterrows():
        cur_locus_tag = row.locus_tag
        if cur_locus_tag in og_df.index:
            cur_og = og_df.OG.loc[cur_locus_tag]
            mask_og = og_df.OG == cur_og
            locus_tags_in_og = og_df.index[mask_og].to_list()
            if cur_gene in res:                                 # to avoid duplicate output filenames, modify genenames if there are duplicates
                n = 1
                new_name = cur_gene + '_DUP_' + str(n)
                while new_name in res:
                    n += 1
                    new_name = cur_gene + '_DUP_' + str(n)
                res[new_name] = locus_tags_in_og
            else:
                res[cur_gene] = locus_tags_in_og
    
    # iterate through results and prepare BioPython records lists
    res_records = dict()
    gbk_mem = dict()
    for cur_gene, cur_locus_tags in res.items():
        for cur_tag in cur_locus_tags:
            cur_assembly = og_df.assembly.loc[cur_tag]
            
            # load biopython assembly data into gbk_mem if it hasnt been done previously to avoid re-reading files.
            if cur_assembly not in gbk_mem:
                cur_assembly_file = os.path.join(gbk_dir, cur_assembly + '.gbk')
                if not os.path.exists(cur_assembly_file):
                    cur_assembly_file = os.path.join(gbk_dir, cur_assembly + '.gbff')
                    if not os.path.exists(cur_assembly_file):
                        print(cur_assembly + ' not found.')
                        continue
                
                gbk_mem[cur_assembly] = []        
                for record in SeqIO.parse(cur_assembly_file, 'genbank'):
                    gbk_mem[cur_assembly].append(record)
            
            # get locus_tag features from biopython
            for record in gbk_mem[cur_assembly]:
                for cur_feature in record.features:
                    if cur_feature.type == 'CDS':
                        if 'locus_tag' in cur_feature.qualifiers:
                            if cur_tag in cur_feature.qualifiers['locus_tag']:
                                # feature found, now init res_records if needed and save
                                if cur_gene not in res_records:
                                    res_records[cur_gene] = []
                                fasta_name = cur_tag
                                fasta_seq = Seq(cur_feature.qualifiers['translation'][0])
                                fasta_desc = []
                                if 'product' in cur_feature.qualifiers:
                                    fasta_desc.extend(cur_feature.qualifiers['product'])
                                if 'protein_id' in cur_feature.qualifiers:
                                    fasta_desc.extend(cur_feature.qualifiers['protein_id'])
                                if 'note' in cur_feature.qualifiers:
                                    fasta_desc.extend(cur_feature.qualifiers['note'])
                                fasta_desc = '; '.join(fasta_desc)

                                res_records[cur_gene].append(SeqRecord(fasta_seq, id=fasta_name, description=fasta_desc))

    # output
    for cur_res in res_records:
        cur_output = os.path.join(output_dir, cur_res + '.faa')
        SeqIO.write(res_records[cur_res], cur_output, 'fasta')
    
    return


if __name__ == '__main__':
    args = get_args()

    input_file = os.path.abspath(args.input_table)
    gbk_dir = os.path.abspath(args.gbk_dir)
    OG_matrix_file = os.path.abspath(args.OG_matrix)
    output_dir = os.path.abspath(args.output_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    main()
