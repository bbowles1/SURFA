#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 13:35:37 2024

Goal: Extract...
1. region length
2. Codon frame state
3. Kozak strength (and start codon)
... then functionalize this script.

TSPEAR ENSG = ENSG00000175894
MEF2c ENSG = ENSG00000081189

Non-canonical (near cognate) start codons:
    CTG, TTG, GTG, UGG, UCG, UTU, UTT, UTC
    

@author: bbowles
"""

import pandas as pd
import numpy as np
import pybedtools
import os
import json

# hard code input params
gtf_path = '/Users/bbowles/Documents/Code/refdata/MANE/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz'
source = 'ensembl_havana'
geneid = 'ENSG00000081189' # MEF2C
#geneid = 'ENSG00000175894' # TSPEAR

FASTA_path = '/Users/bbowles/Documents/Code/refdata/FASTA/GRCh37/GRCh37_latest_genomic.fna'
FASTA_dict = '/Users/bbowles/Documents/Code/refdata/FASTA/GRCh37/FASTA_chrom_identifiers.txt'


def export_uorfs_for_d3(uorfs, export_type='file'):
    """
    Export uORFs data in format suitable for D3.js
    
    Parameters:
    uorfs (dict): The uORFs dictionary
    export_type (str): 'file', 'embedded', or 'both'
    
    Returns:
    str: Path to output file or HTML string
    """
    
    # Flatten and prepare data
    flattened_data = []
    for uorf_id, info in uorfs.items():
        row = {
            'uorf_id': uorf_id,
            'frame': info['frame'],
            'start_index': info['start_index'],
            'stop_index': info['stop_index'],
            'start_codon': info['start_codon'],
            'stop_codon': info['stop_codon'],
            'length': info['length'],
            'kozak_score': info['kozak_score'],
            'sequence': info['sequence'],
            'total_sequence_length': info['total_sequence_length']
        }
        flattened_data.append(row)
    
    if export_type in ['file', 'both']:
        # Export as separate JSON file
        with open('uorfs.json', 'w') as f:
            json.dump(flattened_data, f, indent=2)
    
    if export_type in ['embedded', 'both']:
        # Create HTML with embedded data
        html_template = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>uORF Visualization</title>
            <script src="https://d3js.org/d3.v7.min.js"></script>
            <style>
                .label {{ font-size: 12px; font-family: sans-serif; }}
            </style>
        </head>
        <body>
            <div id="visualization"></div>
            <script>
                const uorfData = {json.dumps(flattened_data, indent=2)};
                
                // Your D3.js visualization code here
                document.addEventListener('DOMContentLoaded', function() {{
                    // Create visualization using uorfData
                }});
            </script>
        </body>
        </html>
        """
        
        with open('uorfs_visualization.html', 'w') as f:
            f.write(html_template)
    
    return 'uorfs.json' if export_type == 'file' else 'uorfs_visualization.html'


def score_kozak(codon):
    
    """
    Score binding strength of each start codon. For now, returns ordinal values
    but is a place holder for more advanced functionality in the future.
    """
    
    codon_dict = {'ATG':1, 'CTG':2, 'TTG':3, 'GTG':4, 'TGG':5, 
                  'TCG':6, 'UTU':7, 'TTT':8, 'TTC':9}
    
    return 1


def get_uorfs(sequence, transcript_id):
    """
    Analyze a DNA sequence for potential uORFs in all three reading frames,
    including non-canonical start codons and Kozak sequence scoring.
    
    Parameters:
    sequence (str): DNA sequence string containing only G, C, A, T
    transcript_id (str): Identifier for the transcript (e.g., "ENST00000504921")
    score_kozak (callable): Function that takes a sequence context and returns binding strength score
    
    Returns:
    dict: Nested dictionary containing uORF information
    """
    # Step 1: Clean and validate input
    sequence = sequence.upper().replace(',', '')
    valid_nucleotides = set('GCTA')
    if not all(nuc in valid_nucleotides for nuc in sequence):
        raise ValueError("Sequence contains invalid nucleotides")
    
    # Store total sequence length
    total_length = len(sequence)
    
    # Step 2: Define start and stop codons (in RNA format)
    stop_codons = {'UAA', 'UAG', 'UGA'}
    start_codons = {
        'AUG',  # Canonical
        'CUG', 'UUG', 'GUG',  # Common non-canonical
        'UGG', 'UCG',  # Less common non-canonical
        'UUU', 'UUC', 'UUA'   # Rare non-canonical
    }
    
    # Step 3: Function to convert DNA to RNA
    def dna_to_rna(dna_seq):
        return dna_seq.replace('T', 'U')
    
    # Step 4: Function to convert RNA to DNA
    def rna_to_dna(rna_seq):
        return rna_seq.replace('U', 'T')
    
    # Step 5: Function to get sequence context around start codon
    def get_start_context(seq, start_pos):
        # Get sequence context while handling edge cases
        context_start = max(0, start_pos)
        context_end = min(len(seq), start_pos + 3)
        return seq[context_start:context_end]
    
    # Step 6: Function to split sequence into codons for a given frame
    def get_codons(seq, frame):
        rna_seq = dna_to_rna(seq[frame:])
        return [rna_seq[i:i+3] for i in range(0, len(rna_seq)-2, 3)]
    
    # Step 7: Function to find all start and stop positions in a frame
    def find_orfs_in_frame(seq, frame):
        codons = get_codons(seq, frame)
        starts = []  # Will now store tuples of (index, start_codon)
        stops = []   # Will continue to store tuples of (index, stop_codon)
        
        for i, codon in enumerate(codons):
            if codon in start_codons:
                starts.append((i, codon))
            elif codon in stop_codons:
                stops.append((i, codon))
        
        return starts, stops
    
    # Step 8: Create counters for each start codon type
    start_codon_counters = {}
    
    # Step 9: Analyze all three frames and build uORF dictionary
    uorfs = {}
    
    for frame in range(3):
        starts, stops = find_orfs_in_frame(sequence, frame)
        
        # Match starts with their corresponding stops
        for start_idx, start_codon in starts:
            start_pos = (start_idx * 3) + frame
            
            # Find the first stop codon that comes after this start
            for stop_idx, stop_codon in stops:
                if stop_idx > start_idx:
                    stop_pos = (stop_idx * 3) + frame
                    
                    # Calculate the sequence of this ORF
                    orf_sequence = sequence[start_pos:stop_pos+3]
                    
                    # Get context for Kozak sequence scoring
                    start_context = get_start_context(sequence, start_pos)
                    kozak_score = score_kozak(start_context)
                    
                    # Convert start codon to DNA format
                    dna_start_codon = rna_to_dna(start_codon)
                    
                    # Update counter for this start codon type
                    if dna_start_codon not in start_codon_counters:
                        start_codon_counters[dna_start_codon] = 1
                    else:
                        start_codon_counters[dna_start_codon] += 1
                    
                    # Generate uORF ID using new format
                    uorf_id = f"{transcript_id}.{dna_start_codon}.{start_codon_counters[dna_start_codon]}"
                    
                    uorfs[uorf_id] = {
                        'frame': frame + 1,
                        'start_index': start_pos,
                        'stop_index': stop_pos,
                        'start_codon': dna_start_codon,
                        'stop_codon': stop_codon,
                        'length': stop_pos - start_pos + 3,
                        'sequence': orf_sequence,
                        'total_sequence_length': total_length,
                        'kozak_score': kozak_score,
                        'start_context': start_context
                    }
                    break
    
    return uorfs




def get_MANE_FASTA(ensg_df, FASTA_path, FASTA_dict):

    # Goal: get UTR sequence and UTR exon bounds for each ENST in the ensembl GRCh37 release
    # provided input is a MANE GTF for a gene ID of interest

    # get FASTA sequence for each 5'UTR exon
    BED_df = ensg_df[['seqname','start','end']]
    BED_df.rename(columns = {'seqname':'chrom', 'start':'chromStart', 'end':'chromEnd'}, inplace=True)
    BED_df.loc[:, 'chrom'] = 'chr'+ BED_df.chrom
    BED_df.loc[:, 'chromStart'] = BED_df.chromStart - 1 # adjusts to zero-based BEDtools nucleotide assignment

    # get FASTA sequences
    print("Obtaining FASTA sequences for all 5'UTR regions.")
    FASTA_list = get_seq(BED_df, FASTA_path, FASTA_dict)

    ensg_df['FASTA'] = np.asarray(FASTA_list) # read FASTA elements back to the pos_df column


    # cut out duplicate entries in df, preferentially keeping entries in ensembl_havana, then havana, then ensembl
    ID_list = pd.Series(ensg_df.loc[(ensg_df.source == 'ensembl_havana')].transcript.unique())
    ensg_df.drop(ensg_df.loc[(ensg_df.source == 'havana') & (ensg_df.transcript.isin(ID_list))].index, axis=0, inplace=True)
    ID_list = ID_list.append((pd.Series(ensg_df.loc[(ensg_df.source == 'havana')].transcript.unique())))
    ensg_df.drop(ensg_df.loc[(ensg_df.source == 'ensembl') & (ensg_df.transcript.isin(ID_list))].index, axis=0, inplace=True)
    ID_list = ID_list.append((pd.Series(ensg_df.loc[(ensg_df.source == 'ensembl')].transcript.unique())))

    # rank values by start site in ascending order
    ensg_df = ensg_df.sort_values(by=['start'], ascending=True)

    # append FASTA sequences to the ensg_df - also appends range, not sure this is necessary
    ENST_FASTA_dict = {}
    range_dict = {}
    for ID in ensg_df.transcript.unique():
        FASTA_list = []
        range_list = []
        for row in ensg_df.loc[ensg_df.transcript == ID].itertuples():
            FASTA_list.append(row.FASTA)
            range_list.append(range((row.start), (row.stop + 1)))
        ENST_FASTA_dict.update({row.transcript:FASTA_list})
        range_dict.update({row.transcript:range_list})

    # create a range_df with transcript, start, stop, stranf, FASTA, and exon junctions ('FASTA_range' column)
    range_df = ensg_df[['transcript', 'seqid','start', 'stop', 'strand']].drop_duplicates(subset = 'transcript', keep = 'first')
    range_df.rename(columns={"seqid":"CHROM"}, inplace=True)
    range_df.reset_index(drop=True, inplace=True)
    range_df['FASTA'] = range_df.transcript.map(ENST_FASTA_dict)
    range_df.loc[:, 'FASTA'] = range_df.apply(lambda row : ''.join(row.FASTA), axis = 1)
    range_df['FASTA_range'] = range_df.transcript.map(range_dict)



def get_seq(BED_df, FASTA_path, FASTA_id_path):
    
    # get seq using pybedtools

    # create a BED object and attempt to get sequence
    Bed_obj = pybedtools.BedTool.from_dataframe(BED_df) # create a BedTools object for pybedtools
    fasta = pybedtools.example_filename(FASTA_path)
    Bed_obj = Bed_obj.sequence(fi=fasta)
    FASTA_list = ((open(Bed_obj.seqfn).read())).split('\n') # opens FASTA text file and splits it to list
    FASTA_list = [ x for x in FASTA_list if ">" not in x ] # remove FASTA line header
    FASTA_list = [i for i in FASTA_list if i] # remove empty elements from list

    # if pybedtools return is empty, search FASTA for RefSeq identifiers instead
    if not bool(FASTA_list):
        
        # map FASTA names to correct identifiers
        FASTA_ids = pd.read_csv(FASTA_id_path, sep='\t')[['Abbreviation', 'RefSeq sequence']].set_index(
            'Abbreviation').to_dict()['RefSeq sequence']
        
        # map FASTA IDs
        BED_df.loc[:, 'chrom'] = BED_df.chrom.map(FASTA_ids)

        # retry search
        Bed_obj = pybedtools.BedTool.from_dataframe(BED_df)  # create a BedTools object for pybedtools
        fasta = pybedtools.example_filename(FASTA_path)
        Bed_obj = Bed_obj.sequence(fi=fasta)
        FASTA_list = ((open(Bed_obj.seqfn).read())).split('\n')  # opens FASTA text file and splits it to list
        FASTA_list = [x for x in FASTA_list if ">" not in x]  # remove FASTA line header
        FASTA_list = [i for i in FASTA_list if i]  # remove empty elements from list
        FASTA_list = list(map(lambda x: x.upper(), FASTA_list)) # convert all elements to upper case

    return FASTA_list

def make_bed(ensg_df):
    
    # make BED-compatable dataframe
    BED_df = ensg_df[['seqname','start','end']].rename(
        columns={'seqname':'chrom',
                 'start':'chromStart',
                 'end':'chromEnd'})
    BED_df.loc[:, 'chromStart'] = BED_df.chromStart.astype(int) - 1

    return BED_df    
    

# Annotate gene_start: the start of coding at the canonical transcript ATG. This is the CDS start position in the Ensembl .gff3.
ensg_df = pd.read_csv(gtf_path, 
                      header=0, 
                      sep='\t', comment='#', names=['seqname', 'source', 'feature', 
                                                    'start', 'end', 'score', 'strand', 'frame', 'attribute'])

# subset to source, UTRs
ensg_df = ensg_df.loc[ensg_df.source==source]

# Extract gene
ensg_df['gene_id_w_version'] = ensg_df.attribute.str.split(';').str[0].str.split(' ').str[1].str.strip('"')

if '.' not in geneid:
    # assume most recent version
    ensg_df['gene_id'] = ensg_df.gene_id_w_version.str.split('.').str[0]
    ensg_df['id_version'] = ensg_df.gene_id_w_version.str.split('.').str[1].astype(int)
    
    ensg_df = ensg_df.loc[ensg_df.gene_id == geneid]
    genid_w_version = f"{geneid}.{str(ensg_df.id_version.max())}"
    ensg_df = ensg_df.loc[ensg_df.gene_id_w_version == genid_w_version]
     
else:
    ensg_df = ensg_df.loc[ensg_df.gene_id_w_version == geneid]
    
if ensg_df.empty:
    raise Exception(f'{geneid} not found in source {source}!')


# get CDS start site
CDS_start = ensg_df.loc[ensg_df.feature=='CDS'].start.min()
    
# subset to source, UTRs
ensg_df = ensg_df.loc[ensg_df.feature=='UTR']

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
ensg_df['utr_type'] = ensg_df.apply(lambda row: check_identity(row.start, row.strand, CDS_start), axis=1)
ensg_df = ensg_df.loc[ensg_df.utr_type==5]
    
# unpack transcript and exon #
ensg_df['transcript'] = ensg_df.attribute.str.split(';').str[1].str.split(' ').str[2].str.strip('"')
ensg_df['exon'] = ensg_df.attribute.str.split(';').str[6].str.split(' ').str[2]


# get FASTA sequences
BED_df = make_bed(ensg_df)
FASTA_list = get_seq(BED_df, FASTA_path, FASTA_dict)
ensg_df['FASTA'] = np.asarray(FASTA_list) # read FASTA elements back to the pos_df column


for row in ensg_df.itertuples():
    
    uorfs = get_uorfs(row.FASTA, row.transcript)


# 1. Get total UTR length
utr_len = (ensg_df.start-ensg_df.end).abs().sum()






