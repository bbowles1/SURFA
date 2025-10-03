#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 18 11:09:30 2025

I want to take my existing JSON structuring solution and turn it into a genome-wide annotation.

Example usage:
    python get_seq.py --bed ./test_bed.bed \
    --fasta /Users/bbowles/Documents/Code/refdata/FASTA/GRCh37/release-113/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa \
    --working_dir . \
    --seqid_map seqid_map.csv \
    --seqid_value 'genbank_sequence' \
    --seqid_key 'number'

@author: bbowles
"""



import pandas as pd
import numpy as np
import os
import math
import warnings
import json
import subprocess
import argparse
import sqlite3


#############
# FUNCTIONS #
#############


def chunker(seq, size):
# chunk data into n sizes
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))


class NpEncoder(json.JSONEncoder):
    # encoder used to write json
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)


def export_uorfs(uorfs):
    """Write the output uorf JSON

    :param uorfs: JSON structure containing uORF data
    :type uorfs: JSON
    """
    # Ensure the directory exists
    os.makedirs('./data', exist_ok=True)
    
    # Write to JSON file
    with open('./data/uorfs.json', 'w') as f:
        #json.dump(uorfs, f)
        json.dump(uorfs, f, cls=NpEncoder)
    

def complement_function(input_FASTA):  # This function translates negative strand nucleotides into their complements, but does not reverse the reading frame - must do this manually
    """FASTA string

    :param input_FASTA: FASTA nucleotide sequence
    :type input_FASTA: str
    :return: Output nucleotide string converted to reverse compliment seq
    :rtype: str
    """
    nucleotide_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
        
    input_FASTA = [nucleotide_dict[k.upper()] for k in input_FASTA]
       
    new_codon = ''.join(input_FASTA)
    
    return new_codon        # output new codons


    
def get_transcript_FASTA(ensg_df):
    """Convert FASTA sequence to reverse complement for negative strands.
    Many to one result where individual exons are concatenated together into transcript sequences.

    :param ensg_df: Ensembl GTF converted to pandas dataframe
    :type ensg_df: pandas.DataFrame
    :raises Exception: Strand Parsing Exception
    :return: Strand corrected ensg-df, where FASTA seq has been converted to reverse comp for negative strands
    :rtype: pandas.DataFrame
    """
    if len(ensg_df.groupby("transcript").strand.nunique().unique()) > 1:
        # there was an error inferring strand identity
        raise Exception("Multiple strand values associated with input transcripts!")
    
    # concatenate FASTA sequences
    ensg_df['transcript_FASTA'] = ensg_df.groupby('transcript')['FASTA'].transform(lambda x: ''.join(x))
    ensg_df = ensg_df[['transcript','strand', 'chrom', 'gene_id', 
                       'transcript_FASTA','gene_name']].drop_duplicates().copy()
    
    # subset based on strand
    pos_df = ensg_df.loc[ensg_df.strand=="+"].copy()
    neg_df = ensg_df.loc[ensg_df.strand=="-"].copy()
    
    # apply negative complement for negative strands
    neg_df["transcript_FASTA"] = neg_df.transcript_FASTA.str[::-1].apply(complement_function)
    
    # append sequences together
    ensg_df = pd.concat([pos_df, neg_df]).sort_index()
    
    return ensg_df



def score_kozak(codon):
    """
    Score binding strength of each start codon. For now, returns ordinal values
    but is a place holder for more advanced functionality in the future.
    """
    
    codon_dict = {'ATG':1, 'CTG':2, 'TTG':3, 'GTG':4, 'TGG':5, 
                  'TCG':6, 'UTU':7, 'TTT':8, 'TTC':9}
    
    return 1


def get_codons(seq, frame):
    # split sequence into codons for a given frame
    rna_seq = seq[frame:]
    return [rna_seq[i:i+3] for i in range(0, len(rna_seq)-2, 3)]


def fasta_codon_search(RNA, frame):
    # Begin RNA reading frame at Nth position, then iterate over in chunks of 3
    return ( list(map(''.join, zip(*[iter(RNA[frame:] )]*3))) )


def filter_codons(codons, targets):
    target_set = set(targets)  # O(1) lookups
    return [(codon, idx) for idx, codon in enumerate(codons) if codon in target_set]


def match_codons(start_codons, stop_codons):
    """
    Find next downstream stop codon for each start codon using numpy
    Returns: list of tuples (start_codon_info, stop_codon_info or None)
    """
    if not stop_codons:
        return [(start, None) for start in start_codons]
    if not start_codons:
        return []
    
    # extract indices from input codons
    start_indices = np.array([idx for _, idx in start_codons])
    stop_indices = np.array([idx for _, idx in stop_codons])
    
    # use searchsorted to find next stop codon for each start
    next_stop_positions = np.searchsorted(stop_indices, start_indices, side='right')
    
    result = []
    for i, start in enumerate(start_codons):
        stop_number = next_stop_positions[i]
        if stop_number < len(stop_codons):
            result.append((start, stop_codons[stop_number]))
        else:
            result.append((start, None))  # no downstream stop codon

    return result


def return_FASTA(FASTA, codon_tuple):
    # check if stop is None
    if codon_tuple[1] == None:
        return FASTA[codon_tuple[0][1]:]
    else:
        return FASTA[codon_tuple[0][1]:codon_tuple[1][1]+3]


def return_FASTA_optimized(FASTA, codon_tuple):
    start = codon_tuple[0][1]
    end = codon_tuple[1]
    
    # Direct slice assignment is faster than conditional
    return FASTA[start:] if end is None else FASTA[start:end[1]+3]


def get_uorfs(input_df):
    """    
    Analyze a table of transcript information and return a table of uORFs

    :param FASTA_df: Table containing `transcript_FASTA` olumn
    :type FASTA_df: pandas Data.Frame
    :return: Pandas dataframe containing uORF information
    :rtype: pandas Data.Frame
    """    

    # convert input to upper
    FASTA_df = input_df.copy()
    FASTA_df["transcript_FASTA"] = FASTA_df.transcript_FASTA.str.upper()
        
    # set start and stop codons
    stop_codons = {'UAA', 'UAG', 'UGA'}
    start_codons = {
        'AUG',  # Canonical499
        'CUG', 'UUG', 'GUG',  # Common non-canonical
        'UGG', 'UCG',  # Less common non-canonical
        'UUU', 'UUC', 'UUA'   # Rare non-canonical
    }
    all_codons = stop_codons.union(start_codons)
    
    # split data into three frame states

    n_rows = FASTA_df.shape[0]
    FASTA_df = FASTA_df.loc[FASTA_df.index.repeat(3)].reset_index(drop=True) # efficient way to copy dataframe 3x
    FASTA_df['frame'] = np.tile([0,1,2], n_rows)
    
    #########
    # STEPS #
    #########
    
    # find uORFs in frame:
        # 1. Convert to RNA seq
        # 1. get codons based on frame
        # 2. enumerate over to get start and indices

    
    # convert FASTA to RNA
    FASTA_df.loc[:, 'transcript_FASTA'] = FASTA_df.transcript_FASTA.str.replace("T","U")
    
    # apply fast codon search
    # faster than vectorized approaches for a chunk size of 1k
    FASTA_df["codons"] = FASTA_df.apply(lambda row: fasta_codon_search(row.transcript_FASTA, row.frame), axis=1)

    # retrieve start/stop codons from the list
    FASTA_df["start_codons"] = FASTA_df.codons.apply(lambda x: filter_codons(x, start_codons))
    FASTA_df["stop_codons"] = FASTA_df.codons.apply(lambda x: filter_codons(x, stop_codons))
    
    # determine which codon pairs are in frame with each other
    FASTA_df["codon_pairs"] = FASTA_df.apply(lambda row: match_codons(row.start_codons, row.stop_codons), axis=1 )
    FASTA_df.drop(columns=["start_codons","stop_codons"], inplace=True)
    
    # drop rows with no codon pairs
    FASTA_df = FASTA_df.loc[FASTA_df["codon_pairs"].apply(bool)]
    
    # explode rows --> 1 uORF per row
    uorf_df = FASTA_df.explode("codon_pairs")
    
    # unpack codon tuples -> need end, start POS
    uorf_df["start_codon"] = uorf_df.codon_pairs.str[0].str[0]
    uorf_df["start_codon_pos"] = uorf_df.codon_pairs.str[0].str[1]
    uorf_df["stop_codon"] = uorf_df.codon_pairs.str[1].str[0]
    uorf_df["stop_codon_pos"] = uorf_df.codon_pairs.str[1].str[1]
    
    # recode missing stop codons - None entries indicate uORF reads into main CDS.
    uorf_df.loc[uorf_df.stop_codon.isna(), "stop_codon"] = "NO_UTR_STOP"
    
    # get uORF sequence from transcript FASTA
    # entries where stop codon is missing are read from the uORF to the CDS start
    sequences = []
    adjusted_bp_start = []
    adjusted_bp_stop = []
    for fasta, start, end in zip(uorf_df['transcript_FASTA'], 
                                 uorf_df['start_codon_pos'], 
                                 uorf_df['stop_codon_pos']):
        bp_start = start*3
        if pd.isna(end):
            bp_end = np.nan
            sequences.append(fasta[bp_start:])
        else:
            bp_end = (end+1)*3
            sequences.append( fasta[bp_start:int(bp_end)] )
        adjusted_bp_start.append(bp_start)
        adjusted_bp_stop.append(bp_end)
            
    # append values back to uorf df
    uorf_df['uORF_FASTA'] = sequences
    uorf_df['start_bp_pos'] = adjusted_bp_start
    uorf_df['stop_bp_pos'] = adjusted_bp_stop
        
    # get total uorf length
    uorf_df['uORF_length'] = uorf_df.stop_bp_pos - uorf_df.start_bp_pos
    
    # generate anotations
    # 1. is entry confined to uORF or does it extend into CDS? -> "uorf_location"
    # 2. is entry in-frame with upstream CDS or out of frame? -> "CDS_frame_status"
    uorf_df["uorf_location"] = np.where(
        uorf_df.stop_bp_pos.isna(),
        "CDS_overlap",
        "UTR_only")
    
    uorf_df["CDS_frame_status"] = np.where(
        (uorf_df.utr_len - uorf_df.frame)%3 == 0,
        "in-frame",
        "out-frame")
    
    # drop transcript-level info: strand, FASTA, total_length (these are retrievable from the transcript table)
    # drop intermmediate info no longer needed: codon_pairs
    output_cols = ['transcript', 'uORF_length', 'start_bp_pos', 'stop_bp_pos', 
                   'uORF_FASTA','start_codon', 'stop_codon', 'frame',
                   'uorf_location', 'CDS_frame_status']
    uorf_df = uorf_df[output_cols]
    
    # rename columns to be more explicit
    # start/stop pos are relative to UTR start and not genomic coords
    uorf_df.rename(columns={
        'start_bp_pos':'rel_start_pos',
        'stop_bp_pos':'rel_stop_pos'}, inplace=True)
        
    return uorf_df


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

def fasta_from_stdout_old_old(fasta):
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

        
def fasta_from_stdout(fasta):
    """Yield (index, sequence) pairs - each sequence isolated"""
    current_index = None
    current_sequence = ""
    
    for line in fasta.split('\n'):
        if line.startswith('>'):
            # Yield the previous complete sequence if it exists
            if current_index is not None:
                yield (current_index, current_sequence)
            
            # Start a new sequence
            current_index = line.split("::")[0][1:]
            current_sequence = ""  # Reset for new sequence
        elif line.strip():  # Skip empty lines
            current_sequence += line  # This builds up the current sequence only
    
    # Yield the final sequence
    if current_index is not None:
        yield (current_index, current_sequence)
     
            
            

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
    :return: Input BED file with 1) unmapped chrom values dropped and 2) a new FASTA sequence column
    :rtype: Pandas Dataframe
    """

    if seqid_dict:
        
        # raise warning if there are empty chrom values
        empty_chrom_fields = BED_df.chrom.isna().sum()
        if empty_chrom_fields > 0:
            warnings.warn(f"GTF file used for FASTA retrieval includes {empty_chrom_fields} empty seqids.")
        
        # remap identifiers
        original_chroms = BED_df.chrom.copy()
        BED_df.loc[:, 'chrom'] = BED_df.chrom.map(seqid_dict)
        dropped_chroms = original_chroms.loc[BED_df.chrom.isna()].drop_duplicates()
        if not dropped_chroms.empty:
            print("\tThe following contigs could not be mapped to their FASTA equivalents and were removed:\n")
            for i in dropped_chroms:
                print(f"\t\t{i}")
            print("\n\tPlease use the seqid-map / seqid-key / seqid-value args to map between the GTF seqid \n\tand your FASTA contigs.")
        
        # drop unmapped cols
        BED_df = BED_df.loc[BED_df.chrom.notna()]

    # add identifier to bedfile
    BED_df["index"] = BED_df.index.astype(str)

    # write temporary bedfile
    tmp_bed = os.path.join(working_dir, 'tmp.bed')
    BED_df.to_csv(tmp_bed, sep='\t', index=False, header=None)

    # format bedtools call
    cmd = ['bedtools', 'getfasta', '-name', '-fi', FASTA_path, '-bed', tmp_bed]

    # call bedtools, capture stdout
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
    # get seq from stdout
    FASTA_list = pd.DataFrame(
        list(fasta_from_stdout(result.stdout)),
        columns=['index', 'FASTA'])
    FASTA_list.loc[:, "index"] = FASTA_list["index"].astype(str)
    
    # merge FASTA into original BED df
    BED_df = BED_df.merge(FASTA_list, on='index', how='outer')

    # remove tmp file
    if os.path.exists(tmp_bed):
        os.remove(tmp_bed)
    
    return BED_df



def make_bed(ensg_df):
    
    # make BED-compatable dataframe
    BED_df = ensg_df[['seqname','start','end']].rename(
        columns={'seqname':'chrom',
                 'start':'chromStart',
                 'end':'chromEnd'})
    BED_df.loc[:, 'chromStart'] = BED_df.chromStart.astype(int) - 1
    BED_df.loc[:, 'chromStart'] = BED_df.chromStart.astype(int).astype(str)
    BED_df.loc[:, 'chromEnd'] = BED_df.chromEnd.astype(int).astype(str)

    return BED_df    


def df_to_sequence(input_df, FASTA_path, output_dir, seqid_path, seqid_key, seqid_value):
    
    # import BED file
    BED_df = make_bed(input_df)
    
    if seqid_path:
        
        # import seqid map
        seqid_map = pd.read_csv(seqid_path)
        # convert to dict
        print("Remapping Seqid value in input BED file.")
        seqid_dict = produce_seqid_dict(seqid_map, seqid_key, seqid_value)
        BED_df = get_seq(BED_df, FASTA_path, output_dir, seqid_dict = seqid_dict)
        
        if BED_df.FASTA.isna().all():
            raise Exception("No sequences were retrieved from the FASTA input! Is your GTF correctly formatted?")
    
    else:
    
        # map FASTA return to input BED df
        BED_df = get_seq(BED_df, FASTA_path, output_dir)
        
    # join FASTA back to input data
    BED_df.loc[:, "index"] = BED_df["index"].astype(int)
    BED_df = BED_df.set_index("index").FASTA
    output_df = input_df.merge(BED_df, left_index=True, right_index=True)
        
    return output_df


########
# MAIN #
########

def main():

    ##########
    # INPUTS #
    ##########

    # Set up argument parser
    parser = argparse.ArgumentParser(description='Parse uORF information from a GTF into a structured JSON.')
    parser.add_argument('--gtf', required=True, help='Path to input GTF build.')
    parser.add_argument('--fasta', required=True, help='Path to input FASTA file.')
    parser.add_argument('--output-dir', required=True, help='Output Directory for files.')
    parser.add_argument('--ensembl-source', required=False, default='ensembl_havana',
                        help='Which GTF Data Source (ie Ensembl, Havana) to use.')
    parser.add_argument('--seqid-map', nargs='?', 
                        help='Dictionary of Genbank/Refseq identifiers for mapping seq ids.', 
                        default=None)
    parser.add_argument('--seqid-key', nargs='?', 
                        help='Column key in the input seqid_map to use for mapping GTF chrom identifiers to their FASTA equivalents.', 
                        default=None)
    parser.add_argument('--seqid-value', nargs='?', 
                        help='Column from seqid_map containing output values to remap GTF chrom identifiers to.', 
                        default=None)



    # Parse arguments
    args = parser.parse_args()
    gtf_path = args.gtf
    output_dir = args.output_dir
    FASTA_path = args.fasta
    source = args.ensembl_source
    seqid_path = args.seqid_map
    output_dir = args.output_dir
    seqid_key = args.seqid_key
    seqid_value = args.seqid_value
    
    # check that all paths exist
    for path in [gtf_path, output_dir, FASTA_path]:
        if not os.path.exists(path):
            raise Exception(f"Path {path} does not exist!")

    #######
    # TMP #
    #######
    
    # set params for testing
    if False:
        gtf_path = "/Users/bbowles/Documents/Code/refdata/MANE/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz"
        FASTA_path = '/Users/bbowles/Documents/Code/refdata/FASTA/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
        output_dir = "/Users/bbowles/Documents/Code/tmp"
        source = "ensembl_havana"
        seqid_path = "/Users/bbowles/Documents/Code/GitHub/d3-uORF-Viewer/deprecating-pybedtools/seqid_map.csv"
        seqid_value='number'
        seqid_key='chr_abbreviation'
        working_dir = output_dir
        
    # set params for testing
    if False:
        # required
        gtf_path = "/Users/bbowles/Documents/Code/GitHub/d3-uORF-Viewer/tests/mini_db/mini.gtf.gz"
        FASTA_path = '/Users/bbowles/Documents/Code/GitHub/d3-uORF-Viewer/tests/mini_db/minifasta.fa'
        output_dir = "/Users/bbowles/Documents/Code/tmp"
        source = "ensembl_havana"

        # optional
        seqid_path = "/Users/bbowles/Documents/Code/GitHub/d3-uORF-Viewer/deprecating-pybedtools/seqid_map.csv"
        seqid_value='number'
        seqid_key='chr_abbreviation'


    ###############
    # DATA IMPORT #
    ###############

    # Annotate gene_start: the start of coding at the canonical transcript ATG. This is the CDS start position in the Ensembl .gff3.
    ensg_df = pd.read_csv(gtf_path, 
                        header=None, 
                        sep='\t', comment='#', names=['seqname', 'source', 'feature', 
                                                        'start', 'end', 'score', 'strand', 'frame', 'attribute'])

    # subset to source
    ensg_df = ensg_df.loc[ensg_df.source==source]

    ####################
    # EXTRACT CDS INFO #
    ####################

    # get length of first CDS bound
    cds = ensg_df.loc[ensg_df.feature=="CDS"].copy()
    cds["transcript"] = cds.attribute.str.split(";").str[1].str.split(' ').str[2].str.strip('"')
    cds["exon"] = cds.attribute.str.split(';').str[6].str.split(' ').str[2].astype(int)
    cds["length"] = cds.end+1 - cds.start
    cds.sort_values(by=["transcript","exon"], ascending=True, inplace=True)

    # subset to first exon of CDS (which is split between 5' UTR and beginning of CDS)
    first_cds = cds.groupby("transcript").first().reset_index().drop_duplicates()
    first_cds = first_cds[['transcript','exon', 'length', 'start']]
    first_cds["type"] = "CDS"


    #################
    # GENE FEATURES #
    #################


    # Create exons table
    exons = ensg_df.loc[ensg_df.feature=='exon'].copy()
    exons['length'] = exons.end+1 - exons.start
    exons['exon'] = exons.attribute.str.split(';').str[6].str.split(' ').str[2]
    exons['transcript'] = exons.attribute.str.split(';').str[1].str.split(' ').str[2].str.strip('"')
    exons['rel_end'] = exons.groupby('transcript')['length'].cumsum()
    exons['rel_start'] = exons['rel_end'] - exons['length']

        

    ###################
    # 5' UTR FEATURES #
    ###################
    # Annotate features including transcript, exon, CDS start site, strand
    # collect 5'UTR regions specifically

    # subset to UTRs
    utr_df = ensg_df.loc[ensg_df.feature=='UTR'].copy()

    # unpack transcript and exon #
    utr_df['transcript'] = utr_df.attribute.str.split(';').str[1].str.split(' ').str[2].str.strip('"')
    utr_df['exon'] = utr_df.attribute.str.split(';').str[6].str.split(' ').str[2]
    utr_df.sort_values(by=["transcript","exon"], ascending=True, inplace=True)

    # merge in CDS start
    utr_df = pd.merge(utr_df,
                       first_cds[['transcript','start']].rename(columns={"start":"cds_start"}), 
                       on=['transcript'], how='left')
    
    # drop rows where CDS start could not be mapped
    print("Dropping rows where CDS start is not defined.")
    utr_df = utr_df.loc[utr_df.cds_start.notna()]

    # determine UTR status
    def check_identity(region_start, strand, CDS_start):
        if strand == '+':
            if region_start > CDS_start:
                return(3)
            elif region_start < CDS_start:
                return(5)
        if strand == '-':
            if region_start > CDS_start:
                return(5)
            elif region_start < CDS_start:
                return(3)
    utr_df['utr_type'] = utr_df.apply(lambda row: check_identity(row.start, row.strand, row.cds_start), axis=1)
    utr_df = utr_df.loc[utr_df.utr_type==5]
    

    #############
    # CDS FRAME #
    #############
    # now that we have UTR info, calculate CDS frame state
    
    utr_df['length'] = utr_df.end+1 - utr_df.start
    
    frame_state = (utr_df.groupby('transcript')["length"].sum() % 3).reset_index()
    frame_state.rename(columns={"length":"frame_state"}, inplace=True)
    utr_df = pd.merge(utr_df, 
             frame_state )
    
    
    ###################
    # UTR ANNOTATIONS #
    ###################
    
    # calculate relative start and stop of each UTR exon
    utr_df['rel_stop'] = utr_df.groupby("transcript").length.cumsum()
    utr_df['rel_start'] = utr_df['rel_stop'] - utr_df.length
    


    #########
    # FASTA #
    #########

    # get FASTA sequences for entire 5'UTR region
    if False:
        # import BED file
        BED_df = make_bed(utr_df)
    
        if seqid_path:
    
            seqid_key = args.seqid_key
            seqid_value = args.seqid_value
    
            # import seqid map
            seqid_map = pd.read_csv(seqid_path)
            # convert to dict
            print("Remapping Seqid value in input BED file.")
            seqid_dict = produce_seqid_dict(seqid_map, seqid_key, seqid_value)
            BED_df = get_seq(BED_df, FASTA_path, output_dir, seqid_dict = seqid_dict)
            
            if BED_df.FASTA.isna().all():
                raise Exception("No sequences were retrieved from the FASTA input! Is your GTF correctly formatted?")
    
        else:
    
            # map FASTA return to input BED df
            BED_df = get_seq(BED_df, FASTA_path, output_dir)
    
    
        # join FASTA back to input data
        BED_df.loc[:, "index"] = BED_df["index"].astype(int)
        BED_df = BED_df.set_index("index").FASTA
        utr_df = utr_df.merge(BED_df, left_index=True, right_index=True)
        
    # new method for retrieving FASTA seq
    utr_df = df_to_sequence(
        utr_df, FASTA_path, output_dir, seqid_path, seqid_key, seqid_value)
    

    # drop columns where FASTA sequence could not be mapped
    drop_rows = utr_df.loc[utr_df.FASTA.isna()].index
    if not drop_rows.empty:
        warnings.warn(f"{len(drop_rows)} FASTA sequences could not be mapped to the input GTF.")
    utr_df = utr_df.loc[utr_df.FASTA.notna()]


    # convert exon FASTA to transcript FASTA
    # 1. concatenate FASTA sequences for all exons of the same transcript
    # 2. Apply reverse complement for negative strands
    utr_df['gene_id'] = utr_df.attribute.str.split(' ').str[1].str.strip('"; ,')
    utr_df['gene_name'] = utr_df.attribute.str.split(';').str[3].str.split(' ').str[2].str.strip('"; ,')
    transcript_df = get_transcript_FASTA(utr_df.rename(columns={"seqname":"chrom"}))
    # determine transcript length
    transcript_df["length"] = transcript_df.transcript_FASTA.str.len()


    ##############
    # UTR LENGTH #
    ##############
    # calc UTR length and append to transcript table
    utr_len = utr_df.groupby("transcript")["length"].sum().reset_index()
    transcript_df = pd.merge(transcript_df,
             utr_len.rename(columns={'length':'utr_len'}))
    

    #############
    # TRANSFORM #
    #############
    # transform data into SQLite
    # 1. group by transcript
    # 2. generate transcript exon table
    # 3. generate uORF sequence table
    # 4. Write to sqlite
    # 5. Make function to unpack sqlite -> JSON

    
    # iterate over transcripts in chunks
    # extrat uORF sequenes
    uorf_list = []
    size = 1000
    count = 0
    n_chunks = math.ceil(transcript_df.shape[0]/size)
    for chunk in chunker(transcript_df, size):
                        
        progress = count / n_chunks * 100
        print(f"\rSearching for uORFs: {progress:.1f}%", end='', flush=True)
        
        tmp_table = get_uorfs(chunk)
        uorf_list.append(tmp_table)
        
        count+=1
    print("\rSearching for uORFs: 100%\n", end='', flush=True)
    uorf_table = pd.concat(uorf_list)

    # sort, add a unique uORF identifier
    uorf_table = uorf_table.sort_values(by=["transcript","start_codon",
                        "rel_start_pos","rel_stop_pos","uORF_length"]).rename(
                            columns={"uORF_length":"uorf_length"})
    
    # add a unique uORF identifier: ENST + codon + number
    uorf_table["uorf_count"] = uorf_table.groupby(["transcript","start_codon"]).cumcount()+1
    uorf_table["uorf_id"] = uorf_table.transcript + "." + uorf_table.start_codon + "." + uorf_table.uorf_count.astype(str)
    uorf_table.drop(columns=["uorf_count"], inplace=True)
    
    # count number of uORFs per transcript
    orfs_per_enst = uorf_table.groupby('transcript').size().reset_index()
    orfs_per_enst.columns = ['transcript','uorf_count']
    transcript_df = pd.merge(transcript_df, orfs_per_enst)
    


    ##########
    # SQLite #
    ##########    
    # export the following tables:
    # 1. transcript table
    # 2. UTR exon table
    # 3. uORF table
    
    # transcript table
    transcript_df
    
    # exon table
    exons.drop(columns=["attribute"], inplace=True)
    
    # uorf table
    uorf_table
    
    # cds_df
    first_cds
    
    # rename columns in utr_df
    utr_df.rename(columns={"seqname":"chrom"}, inplace=True)
    
    # force dtypes in UTR_df
    utr_df.loc[:, 'exon'] = utr_df.exon.astype(int)

    # set explicit output cols for SQL db
    utr_cols = ["transcript","exon","length","rel_start", "rel_stop", "start", "end", "chrom" ,
                     "frame_state", "FASTA"]

    # create database
    print("Writing output to SQL database.")
    db_path = os.path.join(output_dir, "uorfs.db")
    conn = sqlite3.connect(db_path)
    transcript_df.to_sql('transcripts', conn, if_exists='replace', index=False)
    utr_df[utr_cols].to_sql('utr', conn, if_exists='replace', index=False)
    uorf_table.to_sql('uorfs', conn, if_exists='replace', index=False)
    first_cds.to_sql('cds', conn, if_exists='replace', index=False)
    
    conn.close()
    print(f"Database saved to {db_path}")



    # locate example data for MEF2C
    transcript_df.loc[transcript_df.transcript == "ENST00000504921.7"].iloc[0]
    exons.loc[exons.transcript == "ENST00000504921.7"].iloc[0]
    uorf_table.loc[uorf_table.transcript == "ENST00000504921.7"].iloc[0]



    def query_uorf_db(database_path, table, transcript):
        """
        Safely query using context manager for automatic cleanup.
        """
        if table not in ['transcripts','utr','uorfs','cds']:
            raise Exception("SQL table not found in database.")
        try:
            with sqlite3.connect(database_path) as conn:
                query=f"SELECT * FROM {table} WHERE transcript = ?"
                df = pd.read_sql_query(query, conn, params=(transcript,))
                return df
        except sqlite3.Error as e:
            print(f"Database error: {e}")
            return []


    target_transcript = 'ENST00000504921.7' # MEF2C
    enst_query = query_uorf_db(db_path, "transcripts", target_transcript) # transcript
    utr_query = query_uorf_db(db_path, "utr", target_transcript) # utr
    cds_query = query_uorf_db(db_path, "cds", target_transcript) # cds
    uorf_query = query_uorf_db(db_path, "uorfs", target_transcript) # cds
    
    # check for multiple transcript returns
    enst = enst_query.iloc[0]
    utr = utr_query.iloc[0]
    cds = cds_query.iloc[0]
    
    def format_utr(utr_tuple):
        
        out_json = {
            "exon": utr_tuple.exon,
            "type": "utr",
            "start": utr_tuple.rel_start,
            "end": utr_tuple.rel_stop
            }
        return out_json
    
    
    def format_cds(cds_tuple, cds_start):
        out_json = {
            "exon": cds_tuple.exon,
            "type": "cds",
            "start": cds_start,
            "end": cds_start + cds_tuple.length
            }
        return out_json
    
        
    transcript_exons = [format_utr(row) for row in utr_query.itertuples()] + [
        format_cds(cds_query.iloc[0], utr_query.rel_stop.iloc[-1])]
        
    # format transcript block of JSON (transcript meta + exons)
    transcript_block = {
        "id": enst.transcript,
        "type": "transcript",
        "start": 0,
        "end": enst.length,
        "exons": transcript_exons
        }

    # format uORF blocks    
    def format_uorfs(uorf_tuple):
        uorf_start = uorf_tuple.rel_start_pos
        uorf_stop = uorf_tuple.rel_stop_pos
        out_json = {
            "id": uorf_tuple.uorf_id,
            "type": "uorf",
            "frame": uorf_tuple.frame,
            "start": uorf_start,
            "end": uorf_stop,
            "exons": [
                {
                    "exon": "1",
                    "type": "uorf",
                    "start_codon": uorf_tuple.start_codon,
                    "start": uorf_start,
                    "end": uorf_stop
                }
            ],
            "start_codon": uorf_tuple.start_codon,
            "stop_codon": uorf_tuple.stop_codon,
            "length": uorf_tuple.uorf_length,
            "sequence": uorf_tuple.uORF_FASTA
            }
        return out_json

    uorf_block = [format_uorfs(row) for row in uorf_query.itertuples()]
    
    uorf_json = {
        "gene":enst.gene_id,
        "start":0,
        "end":enst.length,
        "regions": [
            [transcript_block] + uorf_block
            ]
        }

    # print first element of uORF block for debugging
    print(json.dumps(uorf_block[0], indent=4))


    ########
    # MAIN #
    ########
            
        
if __name__ == "__main__":
    main()

