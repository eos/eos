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
    parser = argparse.ArgumentParser(description='List the signal PDFs implemented in EOS')
    parser.add_argument('-n', '--filter-by-name', metavar='NAME',
        help = 'Add a filter for the name part of a signal PDF; repeat for several names.',
        dest = 'names', action = 'append', default = []
    )
    parser.add_argument('-p', '--filter-by-prefix', metavar='PREFIX',
        help = 'Add a filter for the prefix part of a signal PDF; repeat for several prefixes.',
        dest = 'prefixes', action = 'append', default = []
    )

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
            list_signal_pdfs(args)
        except Exception as e:
            print(colored('✖ Encountered an unrecoverable error:\n', 'red', attrs=['bold']), f'{e}', file=stdout)
            if not isinstance(e, ValueError):
                traceback.print_exception(e, e, e.__traceback__)
            return 1

    return 0


def print_signal_pdf(entry):
    """Print the name, the description, and the kinematic variables of a single signal PDF."""
    numerator   = ', '.join(entry.numerator_kinematic_variables())
    denominator = ', '.join(entry.denominator_kinematic_variables())

    print(entry.name())
    print(f'    {entry.description()}')
    print()
    print(f'    Numerator kinematic variables:{numerator}')
    print(f'    Denominator kinematic variables:{denominator}')
    print()


def list_signal_pdfs(args):
    """List the signal PDFs implemented in EOS, optionally filtered by name or prefix."""
    signal_pdfs = eos.SignalPDFs()

    names    = set(args.names)
    prefixes = set(args.prefixes)

    for qn, entry in signal_pdfs:
        if names or prefixes:
            if str(qn.prefix_part()) not in prefixes and str(qn.name_part()) not in names:
                continue

        print_signal_pdf(entry)
