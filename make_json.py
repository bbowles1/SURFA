#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 18 11:09:30 2025

Take the uorf.db and export a JSON of the target sequence.

@author: bbowles
"""

import argparse
import logging

# custom imports
from gtf_to_db.json_converter import assemble_json_from_transcript

########
# MAIN #
########

if __name__ == "__main__":

    # Set up argument parser
    parser = argparse.ArgumentParser(description='Export a target uORF as a JSON.')
    parser.add_argument('--db', required=True, help='Path to uorfs.db.')
    parser.add_argument('--transcript', required=True, help='Target transcript (incl version number).')
    parser.add_argument('--output', required=True, help='Output file name (incl JSON extension).')
    parser.add_argument('--log-level', default='INFO', 
                    choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                    help='Set the logging level.')

    # Parse arguments
    args = parser.parse_args()
    database_path = args.db
    target_transcript = args.transcript
    outpath = args.output

    # setup logging
    logger = logging.getLogger(__name__)
    logging.basicConfig(
        format='%(asctime)s - %(name)s - %(funcName)s - %(levelname)s - %(message)s',
        level=logging.getLevelName(args.log_level),
        handlers=[
            logging.FileHandler('make_json.log'),
            logging.StreamHandler()
        ])

    # call main function
    assemble_json_from_transcript(database_path, 
                                  target_transcript, 
                                  outpath)
