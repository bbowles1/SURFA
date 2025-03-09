#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Goal: Make a python wrapper for calls to Bedtools using a
minimal set of dependencies, remove reliance on pybedtools.

@author: bbowles
"""

import subprocess
import pandas as pd
import os
import argparse

def produce_seqid_dict(seqid_map, key, value):
    """Generate a seqid dict from input .tsv of Refseq/Genbank identifiers for human chromosomes.
    Use with seqid_dict in get_seq.

    :param seqid_map: Dataframe containing Refseq/Genbank identifier mapping
    :type seqid_map: pandas.Dataframe
    :param key: Column to use as key in the seqid mapping.
    :type key: str
    :param value: Column to use as value in the seqid mapping. Should match FASTA seqid.
    :type value: str
    """

    return(seqid_map.set_index(key)[value].to_dict())

def fasta_from_stdout(fasta):
    """Unpack FASTA lines from stdout

    :param fasta: stdout stream from subprocess Bedtools call
    :type fasta: subprocess.CompletedProcess.stdout
    :yield: FASTA sequence return from subprocess .stdout
    :rtype: str
    """

    for line in fasta.splitlines():
        if line[0] == '>':
            # line is FASTA ID row
            continue
        else:
            yield line

def get_seq(BED_df, FASTA_path, working_dir, seqid_dict=None): 
    """Get FASTA sequence from a BED file using bedtools getfasta

    :param BED_df: Pandas dataframe with named chrom, chromStart, chromEnd columns
    :type BED_df: pandas.Dataframe
    :param FASTA_path: Path to FASTA file
    :type FASTA_path: str
    :param working_dir: Working dir for temporary files.
    :type working_dir: str
    :param seqid_dict: Dictionary to remap BED_df seq IDs to match input FASTA, defaults to None
    :type seqid_dict: dict, optional
    :return: list of FASTA strings generated from Bedtools getfasta
    :rtype: list
    """

    if seqid_dict:
        # remap identifiers
        BED_df.loc[:, 'chrom'] = BED_df.chrom.map(seqid_dict)

    # write temporary bedfile
    tmp_bed = os.path.join(working_dir, 'tmp.bed')
    BED_df.to_csv(tmp_bed, sep='\t', index=False, header=None)

    # format bedtools call
    cmd = ['bedtools', 'getfasta', '-fi', FASTA_path, '-bed', tmp_bed]

    # call bedtools, capture stdout
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
    # get seq from stdout
    FASTA_list = list(fasta_from_stdout(result.stdout))

    # remove tmp file
    if os.path.exists(tmp_bed):
        os.remove(tmp_bed)
    
    return FASTA_list

if __name__ == "__main__":

    # main is mostly just for testing
    # because it just runs the Bedtools CLI call with extra steps :p
    """Example Calls

    If seqids map between BED and FASTA:
    python get_seq.py --bed ./test_bed.bed \
        --fasta /Users/bbowles/Documents/Code/refdata/FASTA/GRCh37/release-113/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa

    If you need to remap seqids:
    python get_seq.py --bed ./test_bed.bed \
        --fasta /Users/bbowles/Documents/Code/refdata/FASTA/GRCh37/release-113/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa \
        --working_dir . \
        --seqid_map seqid_map.csv \
        --seqid_value 'genbank_sequence' \
        --seqid_key 'number'
    """

    parser = argparse.ArgumentParser(prog='get_seq.py')
    parser.add_argument('--bed', nargs='?', help='Path to BED file containing target regions (chrom/chromStart/chromEnd cols only).')
    parser.add_argument('--fasta', nargs='?', help='Path to FASTA file.')
    parser.add_argument('--working_dir', nargs='?', help='Working dir for intermediate files.', default=os.getcwd())
    parser.add_argument('--seqid_map', nargs='?', help='Dictionary of Genbank/Refseq identifiers for mapping seq ids.', default=None)
    parser.add_argument('--seqid_key', nargs='?', help='BED chromosome identifiers to remap using seqid_map.', default=None)
    parser.add_argument('--seqid_value', nargs='?', help='Output value for BED chromosome identifiers, must match a target col in seqid_map.', default=None)

    # parse args
    args = parser.parse_args()
    BED_df = args.bed
    FASTA_path = args.fasta
    working_dir = args.working_dir
    seqid_path = args.seqid_map

    # import BED file
    BED_df = pd.read_csv(BED_df, sep='\t', header=None, names=['chrom','chromStart','chromEnd'])

    if seqid_path:
        # import seqid map
        seqid_map = pd.read_csv(seqid_path)
        # convert to dict
        print("Remapping Seqid value in input BED file.")
        seqid_dict = produce_seqid_dict(seqid_map, args.seqid_key, args.seqid_value)
        BED_df['FASTA'] = get_seq(BED_df, FASTA_path, working_dir, seqid_dict = seqid_dict)

    else:
        try:
            # map FASTA return to input BED df
            FASTA_list = get_seq(BED_df, FASTA_path, working_dir)
            BED_df['FASTA'] = FASTA_list

            # print return
            print(BED_df.to_string(index=False))

        except:
            # handle cases where seqid does not match between FASTA and BED
            if len(FASTA_list) == 0:
                raise Exception('Return from BEDtools query is empty!')
            else:
                raise Exception(f'Return from BEDtools query is malformed! {BED_df.shape[0]} BED length does not match {len(FASTA_list)} values in FASTA return.')



