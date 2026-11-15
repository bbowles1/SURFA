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
from surfa.db_utils import export_db

logger = logging.getLogger(__name__)

#########
# QUERY #
#########


def register(subparsers, parents=None) -> None:
    """Attach the `export` subcommand to *subparsers*.

    :param subparsers: _description_
    :type subparsers: argparse._SubParsersAction
    """

    parser = subparsers.add_parser(
        "export",
        parents=parents or [],
        help="Export a JSON given a target uORF sequence.",
        description=("Export a target uORF as a JSON."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--db",
        required=True,
        help="Path to uorfs.db (created using surfa build command).",
    )

    parser.add_argument("--output", required=True, help="Output file name.")

    parser.add_argument(
        "--format",
        default="bed",
        required=False,
        help="Output format (csv or bed).",
    )

    parser.set_defaults(func=run)


def run(args: argparse.Namespace) -> int:
    """Export the UTR table from a SURFA uorfs database.

    :param args: Parsed command line args.
    :type args: argparse.Namespace
    :return: Exit status
    :rtype: int
    """

    logger.info("Starting export: output=%r", args.output)
    logger.debug("Full args: %s", args)

    # Parse arguments
    database_path = args.db
    target_format = args.format
    outpath = args.output

    # call main function
    export_db(target_format, database_path, outpath)

    return 0
