import json
import numpy as np
import os
import sqlite3
import pandas as pd

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

def export_uorfs(uorfs, outpath):
    """Write the output uorf JSON

    :param uorfs: JSON structure containing uORF data
    :type uorfs: JSON
    :param outpath: full export path, including JSON extension
    :type outpath: str
    """
    
    # Write to JSON file
    with open(outpath, 'w') as f:
        json.dump(uorfs, f, cls=NpEncoder)
    

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


def assemble_json_from_transcript(database_path, target_transcript, outpath):

    enst_query = query_uorf_db(database_path, "transcripts", target_transcript) # transcript
    utr_query = query_uorf_db(database_path, "utr", target_transcript) # utr
    cds_query = query_uorf_db(database_path, "cds", target_transcript) # cds
    uorf_query = query_uorf_db(database_path, "uorfs", target_transcript) # cds
    
    # check for multiple transcript returns
    if enst_query.shape[0] != 1:
        raise Exception("Multiple transcripts returned!")
    enst = enst_query.iloc[0]
        
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

    # save to outpath
    export_uorfs(uorf_json, outpath)


__all__ = ['NpEncoder','export_uorfs','query_uorf_db']