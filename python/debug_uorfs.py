#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 24 20:39:26 2025

@author: bbowles
"""
import pandas as pd
import numpy as np



input_df = pd.read_csv("/Users/bbowles/Documents/Code/GitHub/d3-uORF-Viewer/vectorizing_get_uorfs/chunk_input.csv")


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
    
    # Store total sequence length
    FASTA_df['total_length'] = FASTA_df.transcript_FASTA.str.len()
    
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
    FASTA_df.to_csv("/Users/bbowles/Documents/Code/GitHub/d3-uORF-Viewer/debug-uorf-iterator/codon_search_input.csv", index=False)
    
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
    uorf_df.to_csv("/Users/bbowles/Documents/Code/GitHub/d3-uORF-Viewer/debug-uorf-iterator/uorf_df.csv", index=False)

    
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
            
    uorf_df['uORF_FASTA'] = sequences
    uorf_df['start_bp_pos'] = adjusted_bp_start
    uorf_df['stop_bp_pos'] = adjusted_bp_stop
    uorf_df['uORF_length'] = uorf_df.stop_bp_pos - uorf_df.start_bp_pos
        
    
    # generate anotations
    # 1. is entry confined to uORF or does it extend into CDS? -> "uorf_location"
    # 2. is entry in-frame with upstream CDS or out of frame? -> "CDS_frame_status"
    uorf_df["uorf_location"] = np.where(
        uorf_df.stop_pos.isna(),
        "CDS_overlap",
        "UTR_only")
    
    uorf_df["CDS_frame_status"] = np.where(
        (uorf_df.total_length - uorf_df.frame)%3 == 0,
        "in-frame",
        "out-frame")
    
    
    # drop transcript-level info: strand, FASTA, total_length (these are retrievable from the transcript table)
    # drop intermmediate info no longer needed: codon_pairs
    output_cols = ['transcript', 'uORF_length', 'start_bp_pos', 'stop_bp_pos', 
                   'uORF_FASTA','start_codon', 'stop_codon', 'frame']
    uorf_df = uorf_df[output_cols]
    
    # rename columns to be more explicit
    # start/stop pos are relative to UTR start and not genomic coords
    uorf_df.rename(columns={
        'start_bp_pos':'rel_start_pos',
        'stop_bp_pos':'rel_stop_pos'}, inplace=True)
    
    # convert RNA to DNA
    #uorf_df["uORF_FASTA"] = uorf_df.uORF_FASTA.str.replace("U","T")
    
    # get total uorf length
    uorf_df["uorf_length"] = (uorf_df.rel_stop_pos - uorf_df.rel_start_pos).astype(int)
    
    return uorf_df
