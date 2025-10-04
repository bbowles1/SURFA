#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 18 11:09:30 2025

Generate a SQL database of uORF given an input Ensembl GTF file.

Example usage:
gtf_to_json.py --gtf "/Users/bbowles/Documents/Code/refdata/MANE/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz" \
    --fasta  '/Users/bbowles/Documents/Code/refdata/FASTA/GRCh37/release-113/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa' \
    --output-dir "/Users/bbowles/Documents/Code/tmp" \
    --ensembl-source "ensembl_havana" \
    --seqid-map  "/Users/bbowles/Documents/Code/GitHub/d3-uORF-Viewer/deprecating-pybedtools/seqid_map.csv" \
    --seqid-key "chr_abbreviation" \
    --seqid-value "number"

@author: bbowles
"""

import pandas as pd
import os
import math
import warnings
import json
import argparse
import sqlite3
import logging

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
    
    # setup logging
    logger = logging.getLogger(__name__)
    logging.basicConfig(
        format='%(asctime)s - %(name)s - %(funcName)s - %(levelname)s - %(message)s',
        level=logging.INFO,
        handlers=[
            logging.FileHandler('gtf_to_json.log'),
            logging.StreamHandler()
        ])
    logger.info('Beginning database build.')

    
    # check that all paths exist
    for path in [gtf_path, output_dir, FASTA_path]:
        if not os.path.exists(path):
            logger.error(f"Input path {path} does not exist!")
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
        
    # set params for testing
    if False:
        # required
        gtf_path = "/Users/bbowles/Documents/Code/GitHub/d3-uORF-Viewer/tests/mini.gtf.gz"
        FASTA_path = '/Users/bbowles/Documents/Code/GitHub/d3-uORF-Viewer/tests/minifasta.fa'
        output_dir = "/Users/bbowles/Documents/Code/tmp"
        source = "ensembl_havana"

        # optional
        seqid_path = None
        seqid_value=None
        seqid_key=None


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
    logger.info(f'Subsetting to dataframe source {source}.')


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
    first_cds = first_cds[['seqname','transcript','exon', 'length', 'start','end']]
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
    no_cds_rows = utr_df.loc[utr_df.cds_start.notna()].shape[0]
    if no_cds_rows > 0:
        logger.warning(f"Removed {no_cds_rows} rows where CDS start is not defined in the input GTF.")
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
    
    logger.info("Matching GTF to FASTA sequences...")    
    
    # new method for retrieving FASTA seq
    utr_df = gtf_to_sequence(
        utr_df, FASTA_path, output_dir, seqid_path, seqid_key, seqid_value)
    first_cds = gtf_to_sequence(
        first_cds, FASTA_path, output_dir, seqid_path, seqid_key, seqid_value)
    
    logger.info("FASTA sequences retrieved.")    


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

    logger.info("Searching transcript sequences for upstream open reading frames...")    
    
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
    
    logger.info("Upstream open reading frames identified.")    

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
    first_cds.drop(columns=["seqname","end"], inplace=True)
    
    # rename columns in utr_df
    utr_df.rename(columns={"seqname":"chrom"}, inplace=True)
    
    # force dtypes in UTR_df
    utr_df.loc[:, 'exon'] = utr_df.exon.astype(int)

    # set explicit output cols for SQL db
    utr_cols = ["transcript","exon","length","rel_start", "rel_stop", "start", "end", "chrom" ,
                     "frame_state", "FASTA"]

    # create database
    logger.info("Writing output to SQL database.")
    
    db_path = os.path.join(output_dir, "uorfs.db")
    conn = sqlite3.connect(db_path)
    transcript_df.to_sql('transcripts', conn, if_exists='replace', index=False)
    utr_df[utr_cols].to_sql('utr', conn, if_exists='replace', index=False)
    uorf_table.to_sql('uorfs', conn, if_exists='replace', index=False)
    first_cds.to_sql('cds', conn, if_exists='replace', index=False)
    
    conn.close()
    logger.info(f"Database saved to {db_path}.")
    print(f"Database saved to {db_path}")


    if False:
    
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

