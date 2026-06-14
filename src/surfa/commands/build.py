#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 18 11:09:30 2025

build: Generate a SQL database of uORF given an input Ensembl GTF file and FASTA sequence.

Example usage:
make_uorf_db.py --gtf "/Users/bbowles/Documents/Code/refdata/MANE/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz" \
    --fasta  '/Users/bbowles/Documents/Code/refdata/FASTA/GRCh37/release-113/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa' \
    --output-dir "/Users/bbowles/Documents/Code/tmp" \
    --ensembl-source "ensembl_havana" \
    --seqid-map  "/Users/bbowles/Documents/Code/GitHub/d3-uORF-Viewer/deprecating-pybedtools/seqid_map.csv" \
    --seqid-key "chr_abbreviation" \
    --seqid-value "number"

@author: bbowles
"""

from __future__ import annotations

import argparse
import logging

# import surfa functions
from surfa.uorf_utils import gtf_to_uorf_db

logger = logging.getLogger(__name__)

#########
# BUILD #
#########


def register(subparsers: argparse._SubParsersAction) -> None:  # noqa: SLF001
    """Attach the `build` subcommand to *subparsers*.

    :param subparsers: _description_
    :type subparsers: argparse._SubParsersAction
    """

    parser: argparse.ArgumentParser = subparsers.add_parser(
        "build",
        help="Build a sqilite database of uORF calls from an input FASTA and GTF sequence.",
        description=(
            """
            Provide an input Ensembl GTF and FASTA sequence file
            The output .db file will contain a table of uORF calls.
            """
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--gtf", required=True, help="Path to Ensembl-format GTF file.")

    parser.add_argument("--fasta", required=True, help="Path to input FASTA file.")

    parser.add_argument(
        "--output-dir", required=True, help="Output Directory for files."
    )

    parser.add_argument(
        "--ensembl-source",
        required=False,
        default="ensembl_havana",
        help="Which Ensembl GTF Data Source (ie Ensembl, Havana) to use.",
    )

    parser.add_argument(
        "--seqid-map",
        nargs="?",
        help="Dictionary of Genbank/Refseq identifiers for mapping seq ids.",
        default=None,
    )

    parser.add_argument(
        "--seqid-key",
        nargs="?",
        help="Column key in the input seqid_map to use for mapping GTF chrom identifiers to their FASTA equivalents.",
        default=None,
    )

    parser.add_argument(
        "--seqid-value",
        nargs="?",
        help="Column from seqid_map containing output values to remap GTF chrom identifiers to.",
        default=None,
    )

    parser.set_defaults(func=run)


def run(args: argparse.Namespace) -> int:
    """Execute the build command.

    :param args: Parsed command line args.
    :type args: argparse.Namespace
    :return: Exit status
    :rtype: int
    """

    logger.info("Starting build: output_dir=%r", args.output_dir)
    logger.debug("Full args: %s", args)

    # Parse arguments
    gtf_path = args.gtf
    FASTA_path = args.fasta
    output_dir = args.output_dir
    source = args.ensembl_source
    seqid_path = args.seqid_map
    seqid_key = args.seqid_key
    seqid_value = args.seqid_value

    # exec build
    gtf_to_uorf_db(
        gtf_path, FASTA_path, output_dir, source, seqid_path, seqid_key, seqid_value
    )

    return 0
