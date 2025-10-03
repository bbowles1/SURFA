#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 20 14:45:57 2025

run tests

@author: bbowles
"""

import sqlite3
import pandas as pd
import unittest

db_path = '/app/tests/uorfs.db'

def import_reference(path):
    # convert FASTA entry to sequence
    sequence = []
    with open(path) as f:
        lines = f.readlines()
        for line in lines:
            if (">" in line) or ("#" in line):
                pass
            else:
                sequence.append(line.strip(" \n"))
    return "".join(sequence)


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


class TestSequences(unittest.TestCase):
    
    # NEW TESTS THAT ARE INFORMATIVE
    def test_mef2c_utr_sequence_identity(self):
        # reference file = mef2c_noncoding_exons.fa
        mef2c_ref_utr = import_reference("/app/tests/mef2c_utr_mrna_sequence.fa")
        
        # extract sequences from db
        utr_query = query_uorf_db(db_path, "utr", 'ENST00000504921.7')
        mef2c_test_exon_1_dna = utr_query.loc[utr_query.exon==1].FASTA.iloc[0]
        mef2c_test_exon_2_dna = utr_query.loc[utr_query.exon==2].FASTA.iloc[0]
        
        # get reverse complement
        mef2c_test_exon_1_rna = complement_function(mef2c_test_exon_1_dna)[::-1]
        mef2c_test_exon_2_rna = complement_function(mef2c_test_exon_2_dna)[::-1]
        
        # assemble complete UTR
        mef2c_test_utr = mef2c_test_exon_1_rna + mef2c_test_exon_2_rna
        
        self.assertEqual(mef2c_ref_utr, mef2c_test_utr)    
        
        
    def test_nrg1_utr_sequence_identity(self):
        # reference file = mef2c_noncoding_exons.fa
        nrg1_ref_utr = import_reference("/app/tests/nrg1_utr_mrna_sequence.fa")
        
        # extract sequences from db, no processing necessary
        # because A) NRG1 UTR is entirely within exon 1 and 2) NRG1 is positive standed
        nrg1_test_utr = ''.join(query_uorf_db(db_path, "utr", 'ENST00000405005.8').FASTA)
        
        self.assertEqual(nrg1_ref_utr, nrg1_test_utr)    


    def test_mef2c_utr_seq_length(self):
        test_transcript = 'ENST00000504921.7'
        # len(UTR_exons) + len(CDS) = len(reference)
        mef2c_ref_noncoding_exons_len = len(
            import_reference("/app/tests/mef2c_noncoding_exons.fa")
            )
        
        # retrieve UTR len
        utr_test_len = query_uorf_db(db_path, "utr", test_transcript).length.sum()
        
        # retrieve CDS len
        cds_test_len = query_uorf_db(db_path, "cds", test_transcript).length.sum()
        
        # determine overall seq length
        mef2c_test_noncoding_exons_len = utr_test_len + cds_test_len
        
        self.assertEqual(mef2c_ref_noncoding_exons_len, mef2c_test_noncoding_exons_len)
     

    def test_nrg1_utr_seq_length(self):
        test_transcript = 'ENST00000405005.8'
        # len(UTR_exons) + len(CDS) = len(reference)
        nrg1_ref_noncoding_exons_len = len(
            import_reference("/app/tests/nrg1_noncoding_exons.fa")
            )
        
        # retrieve UTR len
        utr_test_len = query_uorf_db(db_path, "utr", test_transcript).length.sum()
        
        # retrieve CDS len
        cds_test_len = query_uorf_db(db_path, "cds", test_transcript).length.sum()
        
        # determine overall seq length
        nrg1_test_noncoding_exons_len = utr_test_len + cds_test_len
        
        self.assertEqual(nrg1_ref_noncoding_exons_len, nrg1_test_noncoding_exons_len)
 
# run checks
if __name__ == '__main__':
    unittest.main()


