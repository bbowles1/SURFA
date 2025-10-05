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
