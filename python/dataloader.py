#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 13:35:37 2024

STATUS 12/29: This is the testbed for code, the parser+cli script is parse_gtf.py

Dataloader for uORF regions

Extract gene features:
1. Length of UTR


Extract UTR features:
1. region length
2. Codon frame state
3. Kozak strength (and start codon)
4. Mouseover info: uORF start/stop regions

Test Sets:
TSPEAR ENSG = ENSG00000175894
MEF2c ENSG = ENSG00000081189

Non-canonical (near cognate) start codons:
    CTG, TTG, GTG, UGG, UCG, UTU, UTT, UTC
    
Notes:
    Negative strands

@author: bbowles
"""

import pandas as pd
import numpy as np
import pybedtools
import os
import json


##########
# INPUTS #
##########

# set wd
os.chdir('/Users/bbowles/Documents/Code/GitHub/d3-uORF-Viewer/')

# hard code input params
gtf_path = '/Users/bbowles/Documents/Code/refdata/MANE/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz'
source = 'ensembl_havana'
geneid = 'ENSG00000081189' # MEF2C
#geneid = 'ENSG00000175894' # TSPEAR

FASTA_path = '/Users/bbowles/Documents/Code/refdata/FASTA/GRCh37/GRCh37_latest_genomic.fna'
FASTA_dict = '/Users/bbowles/Documents/Code/refdata/FASTA/GRCh37/FASTA_chrom_identifiers.txt'


#############
# FUNCTIONS #
#############


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
    
    # Ensure the directory exists
    os.makedirs('./data', exist_ok=True)
    
    # Write to JSON file
    with open('./data/uorfs.json', 'w') as f:
        #json.dump(uorfs, f)
        json.dump(uorfs, f, cls=NpEncoder)
    

def complement_function(input_FASTA):  # This function translates negative strand nucleotides into their complements, but does not reverse the reading frame - must do this manually
    
    nucleotide_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
        
    input_FASTA = [nucleotide_dict[k.upper()] for k in input_FASTA]
       
    new_codon = ''.join(input_FASTA)
    
    return new_codon        # output new codons

    
def get_transcript_FASTA(ensg_df):
    
    if ensg_df.strand.nunique() != 1:
        raise Exception(f"Multiple strand values associated with transcripts: {','.join(ensg_df.transcript.unique())}")
    strand = ensg_df.strand.iloc[0]
    
    # concatenate FASTA sequences
    ensg_df['transcript_FASTA'] = ensg_df.groupby('transcript')['FASTA'].transform(lambda x: ''.join(x))
    ensg_df = ensg_df[['transcript','transcript_FASTA']].drop_duplicates()
    
    if strand == '+':
        return ensg_df
    elif strand == '-':
        # apply reverse complement
        reverse_FASTA = ensg_df.transcript_FASTA.str[::-1].apply(complement_function)
        ensg_df['transcript_FASTA'] = reverse_FASTA
        return ensg_df



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
    
    # Clean and validate input
    sequence = sequence.upper().replace(',', '')
    valid_nucleotides = set('GCTAN')
    if not all(nuc in valid_nucleotides for nuc in sequence):
        raise ValueError("Sequence contains invalid nucleotides")
    
    # Store total sequence length
    total_length = len(sequence)
    
    # set start and stop codons
    stop_codons = {'UAA', 'UAG', 'UGA'}
    start_codons = {
        'AUG',  # Canonical
        'CUG', 'UUG', 'GUG',  # Common non-canonical
        'UGG', 'UCG',  # Less common non-canonical
        'UUU', 'UUC', 'UUA'   # Rare non-canonical
    }
    
    # convert DNA to RNA nucleotides
    def dna_to_rna(dna_seq):
        return dna_seq.replace('T', 'U')
    
    def rna_to_dna(rna_seq):
        return rna_seq.replace('U', 'T')
    
    # get sequence context around start codon
    def get_start_context(seq, start_pos):
        context_start = max(0, start_pos)
        context_end = min(len(seq), start_pos + 3)
        return seq[context_start:context_end]
    
    # split sequence into codons for a given frame
    def get_codons(seq, frame):
        rna_seq = dna_to_rna(seq[frame:])
        return [rna_seq[i:i+3] for i in range(0, len(rna_seq)-2, 3)]
    
    # find all start and stop positions in a frame
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
    
    # create counters for each start codon type
    start_codon_counters = {}
    
    # analyze all three frames and build uORF dictionary
    uorfs = pd.DataFrame(columns=[
        'uorf_id','frame',
        'start_pos','stop_pos',
        'dna_start_codon','stop_codon',
        'orf_sequence','total_length',
        'kozak_score','start_context'])
    
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
                    
                    # create a pandas dataframe row
                    uorf_row = pd.DataFrame([[uorf_id,
                                 frame,
                                 start_pos,
                                 stop_pos,
                                 dna_start_codon,
                                 stop_codon,
                                 orf_sequence,
                                 total_length,
                                 kozak_score,
                                 start_context
                                 ]], columns=[
                                     'uorf_id','frame',
                                     'start_pos','stop_pos',
                                     'dna_start_codon','stop_codon',
                                     'orf_sequence','total_length',
                                     'kozak_score','start_context'])
                    
                    # append row to pandas dataframe
                    if uorfs.empty:
                        uorfs = uorf_row
                    else:
                        uorfs = pd.concat([uorfs, uorf_row], axis=0)
                        uorfs.reset_index(drop=True)
                    
                    break
    
    # sort, remove duplicates 
    uorfs.drop_duplicates(inplace=True)
    uorfs.sort_values(by=['start_pos','stop_pos'], ascending=False, inplace=True)
    
    # convert to json
    uorf_list = []
    
    for row in uorfs.itertuples():
        uorf_row = {
            'id': row.uorf_id,
            'type': 'uorf',
            'frame': row.frame + 1,
            'start': row.start_pos,
            'end': row.stop_pos,
            'exons': [{
                'exon':'1',
                'type':'uorf',
                'start': int(row.start_pos),
                'end': int(row.stop_pos)
                }],
            'start_codon': row.dna_start_codon,
            'stop_codon': row.stop_codon,
            'length': row.stop_pos - row.start_pos + 3,
            'sequence': row.orf_sequence,
            'total_sequence_length': row.total_length,
            'kozak_score': row.kozak_score,
            'start_context': row.start_context
        }
        
        uorf_list.append( uorf_row )
    
    return uorf_list


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
    range_df = ensg_df[['transcript', 'seqid','start', 'end', 'strand']].drop_duplicates(subset = 'transcript', keep = 'first')
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



########
# MAIN #
########


# Annotate gene_start: the start of coding at the canonical transcript ATG. This is the CDS start position in the Ensembl .gff3.
ensg_df = pd.read_csv(gtf_path, 
                      header=0, 
                      sep='\t', comment='#', names=['seqname', 'source', 'feature', 
                                                    'start', 'end', 'score', 'strand', 'frame', 'attribute'])

# subset to source
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


#################
# GENE FEATURES #
#################

# get CDS start site - used to identify strand
CDS_start = ensg_df.loc[ensg_df.feature=='CDS'].start.min()

# Get UTR start site
exons = ensg_df.loc[ensg_df.feature=='exon']
exons['length'] = (exons.start - exons.end).abs()
exons['exon'] = exons.attribute.str.split(';').str[6].str.split(' ').str[2]
exons['transcript'] = exons.attribute.str.split(';').str[1].str.split(' ').str[2].str.strip('"')

    

###################
# 5' UTR FEATURES #
###################

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


#################
# EXON FEATURES #
#################

# 1. exon #
# 2. exon length
# 3. UTR (color 1) or CDS (color 2)

# create objects for exon(uorf) and exon(CDS) to create shading
CDS_exons = exons.loc[exons.exon.isin(ensg_df.exon)][['start','end','exon','transcript']]
CDS_exons['length'] = (CDS_exons.start - CDS_exons.end).abs()

uorf_exons = ensg_df[['start','end','exon', 'transcript']]
uorf_exons['type'] = 'utr'
uorf_exons['length'] = (uorf_exons.start - uorf_exons.end).abs()

# determine "gap" between end of uorf exons and beginning of CDS
gap_exons = pd.merge(CDS_exons[['exon','length','transcript']].rename(columns={'length':'cds_length'}),
         uorf_exons[['exon','length','transcript']].rename(columns={'length':'uorf_length'}))

gap_exons['gap'] = gap_exons.cds_length - gap_exons.uorf_length
gap_exons = gap_exons.loc[gap_exons.gap > 0].rename(columns={'gap':'length'})
gap_exons['type'] = 'cds'
gap_exons = gap_exons[['exon','length','type','transcript']]

# concat gap and regular exon bounds
gap_exons = pd.concat([uorf_exons[['exon','length','type','transcript']], gap_exons])

# sort
gap_exons['type'] = pd.Categorical(gap_exons['type'], ["utr", "cds"])
gap_exons.sort_values(by=['exon','type'], inplace=True)

# determine start and stop coords
sorted_exons = []
for transcript in gap_exons.transcript.unique():
    coord_exons = gap_exons.loc[gap_exons.transcript == transcript]
    coord_exons['start'] = np.nan
    coord_exons['end'] = np.nan
    start = 0
    for row in coord_exons.itertuples():
        start = start
        end = start+row.length
        coord_exons.loc[row.Index, 'start'] = start
        coord_exons.loc[row.Index, 'end'] = end
        start = start + row.length
    sorted_exons.append(coord_exons)
sorted_exons = pd.concat(sorted_exons)



#########
# FASTA #
#########

# get FASTA sequences
BED_df = make_bed(ensg_df)
FASTA_list = get_seq(BED_df, FASTA_path, FASTA_dict)
ensg_df['FASTA'] = np.asarray(FASTA_list) # read FASTA elements back to the pos_df column


# convert exon FASTA to transcript FASTA
# 1. append FASTA sequences
# 2. Apply reverse complement for negative strands
transcript_df = get_transcript_FASTA(ensg_df)


# iterate over exons to identify uORFs
for row in transcript_df.itertuples():
            
    # identify all exons for transcript
    exon_array = sorted_exons.loc[sorted_exons.transcript == row.transcript]
    
    # get exon information into JSON format
    exon_array.rename(columns={'transcript':'id'}, inplace=True)
    transcript_dict = {}
    for transcript in exon_array.id.unique():
        dict_entry = exon_array.loc[exon_array.id == transcript]        
        transcript_len = dict_entry.length.sum()
        dict_entry = dict_entry.drop(columns=['length','id']).to_dict(orient='records')
        transcript_dict.update(
            {'id' : transcript,
            'type' : 'transcript',
            'start' : 0,
            'end' : transcript_len,
            'exons' : dict_entry})
        
    # identify uORF sites in UTR-FASTA
    uorf_array = get_uorfs(row.transcript_FASTA, row.transcript)
    
    # structure UTR array
    utr_array = {
        'gene' : geneid,
        'start' : 0,
        'end' : exon_array.end.max(),
        'regions' : [transcript_dict] + uorf_array
        }
        
    export_uorfs( utr_array )

