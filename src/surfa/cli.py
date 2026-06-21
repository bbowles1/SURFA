#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
surfa — command-line interface.

Sub-parsers must contain register() and run() definitions.

"""

from __future__ import annotations

import argparse
import sys
import logging
import os

from surfa.commands import build, query


def make_parser() -> argparse.ArgumentParser:
    """Construct and return the top-level argument parser.

    :return: Main argument parser to delegate downstream function calls.
    :rtype: argparse.ArgumentParser
    """

    # parent carries shared flags (log-level, log-file)
    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument(
        "--log-level",
        "-l",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level.",
    )
    global_parser.add_argument(
        "--log-file",
        metavar="FILE",
        default="surfa.log",
        help="Write logs to this file (default: surfa.log in the current directory).",
    )

    # root parser owns the user-facing identity and --version
    parser = argparse.ArgumentParser(
        prog="surfa",
        parents=[global_parser],
        description="Surfa Library: identify upstream open reading frames in sequence data.",
        epilog="Run `surfa <command> --help` for per-command options.",
    )
    parser.add_argument(
        "--version",
        "-V",
        action="version",
        version="%(prog)s 0.1.0",
    )

    # subcommands
    subparsers = parser.add_subparsers(
        title="commands",
        metavar="<command>",
        dest="command",
    )

    build.register(subparsers, parents=[global_parser])
    query.register(subparsers, parents=[global_parser])

    return parser


def _configure_logging(log_level: str, log_file: str | None) -> None:
    """Set up logging

    :param log_level: "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"
    :type log_level: str
    :param log_file: Output path for log file.
    :type log_file: str | None
    """

    log_dir = os.path.dirname(log_file)
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)

    logging.basicConfig(
        format="%(asctime)s - %(name)s - %(funcName)s - %(levelname)s - %(message)s",
        level=logging.getLevelName(log_level),
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_file),
        ],
    )


def main(argv: list[str] | None = None) -> int:
    """Parse *argv* and dispatch to the appropriate subcommand.

    :param argv: surfa command to call (ie build or export), defaults to None
    :type argv: list[str] | None, optional
    :return: Sub-parser to use for downstream commands.
    :rtype: args.func
    """

    parser = make_parser()
    args = parser.parse_args(argv)

    # no subcommand supplied -> print help and exit cleanly.
    if args.command is None:
        parser.print_help()
        return 0

    _configure_logging(args.log_level, args.log_file)

    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
