# Copyright (c) 2026 Danny van Dyk
#
# This file is part of the EOS project. EOS is free software;
# you can redistribute it and/or modify it under the terms of the GNU General
# Public License version 2, as published by the Free Software Foundation.
#
# EOS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place, Suite 330, Boston, MA  02111-1307  USA

from __future__ import annotations

import argparse
from collections.abc import Sequence
from contextlib import redirect_stderr, redirect_stdout
import logging
import sys
from typing import TextIO


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Download and manage public EOS datasets')
    common_subparser = argparse.ArgumentParser(add_help=False)
    common_subparser.add_argument(
        '-v', '--verbose',
        help='Increase the verbosity of the script',
        dest='verbose', action='count', default=0,
    )
    subparsers = parser.add_subparsers(title='commands', dest='command')

    parser_download = subparsers.add_parser(
        'download',
        parents=[common_subparser],
        description='Download a single public EOS dataset from Github or Zenodo.',
        help='Download a public EOS dataset.',
    )
    parser_download.add_argument('id', metavar='ID', help='The unique id of the public EOS dataset.')
    parser_download.set_defaults(cmd=cmd_download)

    parser_list = subparsers.add_parser(
        'list',
        parents=[common_subparser],
        description='List the ids of all public EOS datasets.',
        help='List the ids of all public EOS datasets.',
    )
    parser_list.set_defaults(cmd=cmd_list)

    parser_update = subparsers.add_parser(
        'update',
        parents=[common_subparser],
        description='Update the list of public EOS datasets from Github.',
        help='Update the list of public EOS datasets from Github.',
    )
    parser_update.set_defaults(cmd=cmd_update)

    return parser


def _configure_logging(verbosity: int) -> None:
    levels = {
        0: logging.ERROR,
        1: logging.WARNING,
        2: logging.INFO,
        3: logging.DEBUG,
    }
    logging.basicConfig(level=levels[min(verbosity, 3)])


def _datasets():
    from ..datasets import DataSets

    return DataSets()


def cmd_download(args: argparse.Namespace, **_kwargs) -> int:
    _datasets().download(args.id)
    return 0


def cmd_list(args: argparse.Namespace, *, stdout: TextIO = sys.stdout, **_kwargs) -> int:
    for dataset_id, dataset in _datasets().datasets():
        print(f'{dataset_id:<9}  -  {dataset.title}', file=stdout)
        print(f'{"":<9}     {", ".join(dataset.authors)}', file=stdout)
    return 0


def cmd_update(args: argparse.Namespace, **_kwargs) -> int:
    _datasets().update()
    return 0


def main(
    argv: Sequence[str] | None = None,
    *,
    stdout: TextIO = sys.stdout,
    stderr: TextIO = sys.stderr,
) -> int:
    parser = _parser()
    try:
        with redirect_stdout(stdout), redirect_stderr(stderr):
            args = parser.parse_args(argv)
    except SystemExit as error:
        return int(error.code)

    if not hasattr(args, 'cmd') or not callable(args.cmd):
        parser.print_help(file=stdout)
        return 2

    _configure_logging(args.verbose)
    return args.cmd(
        args,
        stdout=stdout,
        stderr=stderr,
    )
