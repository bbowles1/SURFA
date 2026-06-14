#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 18 11:09:30 2025

Query: Take the uorf.db and export a JSON of the target sequence.

@author: bbowles
"""

import argparse
import logging

# custom imports
from surfa.json_converter import assemble_json_from_transcript

logger = logging.getLogger(__name__)

#########
# QUERY #
#########

def register(subparsers: argparse._SubParsersAction) -> None:  # noqa: SLF001
    """Attach the `query` subcommand to *subparsers*.

    :param subparsers: _description_
    :type subparsers: argparse._SubParsersAction
    """

    parser: argparse.ArgumentParser = subparsers.add_parser(
        "query",
        help="Export a JSON given a target uORF sequence.",
        description=("Export a target uORF as a JSON."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        '--db', 
        required=True, 
        help='Path to uorfs.db (created using surfa build command).')
    
    parser.add_argument(
        '--transcript', 
        required=True, 
        help='Target transcript (including version number).')
    
    parser.add_argument(
        '--output', 
        default="query.json",
        help='Output file name (incl JSON extension).')
    
    parser.set_defaults(func=run)


def run(args: argparse.Namespace) -> int:
    """Execute the query command.

    :param args: Parsed command line args.
    :type args: argparse.Namespace
    :return: Exit status
    :rtype: int
    """

    logger.info("Starting query: output=%r", args.output)
    logger.debug("Full args: %s", args)

    # Parse arguments
    args = parser.parse_args()
    database_path = args.db
    target_transcript = args.transcript
    outpath = args.output

    # call main function
    assemble_json_from_transcript(database_path, 
                                  target_transcript, 
                                  outpath)

    return 0