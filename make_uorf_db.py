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

import argparse

# custom imports
from gtf_to_db.uorf_utils import gtf_to_uorf_db

########
# MAIN #
########

if __name__ == "__main__":

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
    FASTA_path = args.fasta
    output_dir = args.output_dir
    source = args.ensembl_source
    seqid_path = args.seqid_map
    seqid_key = args.seqid_key
    seqid_value = args.seqid_value

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

    # call main function
    gtf_to_uorf_db(gtf_path,
                   FASTA_path,
                   output_dir,
                   source,
                   seqid_path,
                   seqid_key,
                   seqid_value)
