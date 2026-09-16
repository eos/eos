# Copyright (c) 2016-2026 Danny van Dyk
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
import sys
import traceback
from typing import TextIO

import eos

try:
    from termcolor import colored
except ImportError:
    colored = lambda s, *args, **kwargs: s


def _parser():
    parser = argparse.ArgumentParser(description='List references used in EOS')
    parser.add_argument('refs', metavar='REFS', type=str, nargs='*', help='Key for a reference')

    return parser


def main(argv: Sequence[str] | None = None, *, stdout: TextIO | None = None, stderr: TextIO | None = None) -> int:
    stdout = sys.stdout if stdout is None else stdout
    stderr = sys.stderr if stderr is None else stderr
    parser = _parser()
    try:
        with redirect_stdout(stdout), redirect_stderr(stderr):
            args = parser.parse_args(argv)
    except SystemExit as e:
        return int(e.code)

    with redirect_stdout(stdout), redirect_stderr(stderr):
        try:
            list_references(args)
        except Exception as e:
            print(colored('✖ Encountered an unrecoverable error:\n', 'red', attrs=['bold']), f'{e}', file=stdout)
            if not isinstance(e, ValueError):
                traceback.print_exception(e, e, e.__traceback__)
            return 1

    return 0


def list_references(args):
    """List the references used in EOS, optionally filtered by their keys."""
    refs = eos.References()

    for rn, r in refs:
        if args.refs and rn not in args.refs:
            continue

        eprint = None
        if r.eprint_archive() == 'arXiv':
            eprint = 'arXiv:' + r.eprint_id().rsplit(':', 1)[-1]

        print(f'[{rn}] : {r.authors()}')
        print(f'    {r.title()}')
        if eprint:
            print(f'    {eprint}')
        print()
