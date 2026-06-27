#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 20 14:45:57 2025

run tests

@author: bbowles
"""

import unittest

# custom imports
from surfa.fasta_utils import complement_function
from surfa.uorf_utils import import_reference
from surfa.json_converter import query_uorf_db

db_path = "/app/tests/uorfs.db"


class TestSequences(unittest.TestCase):
    # NEW TESTS THAT ARE INFORMATIVE
    def test_mef2c_utr_sequence_identity(self):
        # reference file = mef2c_noncoding_exons.fa
        mef2c_ref_utr = import_reference("/app/tests/data/mef2c_utr_mrna_sequence.fa")

        # extract sequences from db
        utr_query = query_uorf_db(db_path, "utr", "ENST00000504921.7")
        mef2c_test_exon_1_dna = utr_query.loc[utr_query.exon == 1].FASTA.iloc[0]
        mef2c_test_exon_2_dna = utr_query.loc[utr_query.exon == 2].FASTA.iloc[0]

        # get reverse complement
        mef2c_test_exon_1_rna = complement_function(mef2c_test_exon_1_dna)[::-1]
        mef2c_test_exon_2_rna = complement_function(mef2c_test_exon_2_dna)[::-1]

        # assemble complete UTR
        mef2c_test_utr = mef2c_test_exon_1_rna + mef2c_test_exon_2_rna

        self.assertEqual(mef2c_ref_utr, mef2c_test_utr)

    def test_nrg1_utr_sequence_identity(self):
        # reference file = mef2c_noncoding_exons.fa
        nrg1_ref_utr = import_reference("/app/tests/data/nrg1_utr_mrna_sequence.fa")

        # extract sequences from db, no processing necessary
        # because A) NRG1 UTR is entirely within exon 1 and 2) NRG1 is positive standed
        nrg1_test_utr = "".join(
            query_uorf_db(db_path, "utr", "ENST00000405005.8").FASTA
        )

        self.assertEqual(nrg1_ref_utr, nrg1_test_utr)

    def test_mef2c_utr_seq_length(self):
        test_transcript = "ENST00000504921.7"
        # len(UTR_exons) + len(CDS) = len(reference)
        mef2c_ref_noncoding_exons_len = len(
            import_reference("/app/tests/data/mef2c_noncoding_exons.fa")
        )

        # retrieve UTR len
        utr_test_len = query_uorf_db(db_path, "utr", test_transcript).length.sum()

        # retrieve CDS len
        cds_test_len = query_uorf_db(db_path, "cds", test_transcript).length.sum()

        # determine overall seq length
        mef2c_test_noncoding_exons_len = utr_test_len + cds_test_len

        self.assertEqual(mef2c_ref_noncoding_exons_len, mef2c_test_noncoding_exons_len)

    def test_nrg1_utr_seq_length(self):
        test_transcript = "ENST00000405005.8"
        # len(UTR_exons) + len(CDS) = len(reference)
        nrg1_ref_noncoding_exons_len = len(
            import_reference("/app/tests/data/nrg1_noncoding_exons.fa")
        )

        # retrieve UTR len
        utr_test_len = query_uorf_db(db_path, "utr", test_transcript).length.sum()

        # retrieve CDS len
        cds_test_len = query_uorf_db(db_path, "cds", test_transcript).length.sum()

        # determine overall seq length
        nrg1_test_noncoding_exons_len = utr_test_len + cds_test_len

        self.assertEqual(nrg1_ref_noncoding_exons_len, nrg1_test_noncoding_exons_len)


# run checks
if __name__ == "__main__":
    unittest.main()
