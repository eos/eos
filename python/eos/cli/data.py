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
from contextlib import contextmanager, redirect_stderr, redirect_stdout
import logging
from pathlib import Path
import sys
from typing import TextIO

from .data_checks import (
    CheckContext,
    CheckFactory,
    CheckScope,
    PlainTextRenderer,
    register_basic_checks,
    run_checks,
)
from .data_checks_dataset import register_dataset_checks


def create_check_factory() -> CheckFactory:
    """Return the built-in check registry."""
    factory = CheckFactory()
    register_basic_checks(factory)
    register_dataset_checks(factory)
    return factory


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

    parser_check = subparsers.add_parser(
        'check',
        parents=[common_subparser],
        description='Check an EOS dataset candidate.',
        help='Check an EOS dataset candidate.',
    )
    parser_check.add_argument(
        '--analysis-file',
        action='append',
        default=None,
        metavar='PATH',
        help='Select an analysis file relative to DIRECTORY; repeat for multiple files.',
    )
    parser_check.add_argument(
        '--main-analysis-file',
        metavar='PATH',
        help='Select the main analysis file from the analysis files.',
    )
    parser_check.add_argument(
        '-i', '--interactive',
        action='store_true',
        help='Use interactive checking (not yet available).',
    )
    parser_check.add_argument(
        'directory',
        nargs='?',
        default='.',
        metavar='DIRECTORY',
        help='Dataset directory to check (default: current directory).',
    )
    parser_check.set_defaults(cmd=cmd_check)

    return parser


@contextmanager
def _configured_logging(verbosity: int, stderr: TextIO):
    import eos

    levels = {
        0: logging.ERROR,
        1: logging.WARNING,
        2: logging.INFO,
        3: logging.DEBUG,
    }
    handler = eos.stderr_handler
    previous_level = handler.level
    previous_stream = handler.stream
    handler.setLevel(levels[min(verbosity, 3)])
    handler.setStream(stderr)
    try:
        yield
    finally:
        handler.setStream(previous_stream)
        handler.setLevel(previous_level)


def _datasets():
    from ..datasets import DataSets

    return DataSets()


def cmd_download(args: argparse.Namespace, **_kwargs) -> int:
    _datasets().download(args.id)
    return 0


def cmd_list(args: argparse.Namespace, *, stdout: TextIO | None = None, **_kwargs) -> int:
    stdout = sys.stdout if stdout is None else stdout
    for dataset_id, dataset in _datasets().datasets():
        print(f'{dataset_id:<9}  -  {dataset.title}', file=stdout)
        print(f'{"":<9}     {", ".join(dataset.authors)}', file=stdout)
    return 0


def cmd_update(args: argparse.Namespace, **_kwargs) -> int:
    _datasets().update()
    return 0


def cmd_check(
    args: argparse.Namespace,
    *,
    stdout: TextIO | None = None,
    stderr: TextIO | None = None,
    check_factory: CheckFactory | None = None,
) -> int:
    stdout = sys.stdout if stdout is None else stdout
    stderr = sys.stderr if stderr is None else stderr
    if args.interactive:
        print('eos-data check: interactive operation is not yet available', file=stderr)
        return 2

    context = CheckContext(
        dataset_root=Path(args.directory),
        analysis_paths=tuple(Path(path) for path in (args.analysis_file or ())),
        main_analysis_path=(
            Path(args.main_analysis_file) if args.main_analysis_file is not None else None
        ),
    )
    factory = check_factory if check_factory is not None else create_check_factory()
    result = run_checks(factory, context, scope=CheckScope.COMPLETE)
    PlainTextRenderer().write(result, stdout)
    return result.exit_status


def main(
    argv: Sequence[str] | None = None,
    *,
    stdout: TextIO | None = None,
    stderr: TextIO | None = None,
    check_factory: CheckFactory | None = None,
) -> int:
    stdout = sys.stdout if stdout is None else stdout
    stderr = sys.stderr if stderr is None else stderr
    parser = _parser()
    try:
        with redirect_stdout(stdout), redirect_stderr(stderr):
            args = parser.parse_args(argv)
    except SystemExit as error:
        return int(error.code)

    if not hasattr(args, 'cmd') or not callable(args.cmd):
        parser.print_help(file=stdout)
        return 0

    with _configured_logging(args.verbose, stderr):
        return args.cmd(
            args,
            stdout=stdout,
            stderr=stderr,
            check_factory=check_factory,
        )
