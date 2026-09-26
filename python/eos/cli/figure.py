# Copyright (c) 2025-2026 Danny van Dyk
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

import argcomplete
import argparse
from collections.abc import Sequence
from contextlib import contextmanager, redirect_stderr, redirect_stdout
import logging
import os
import sys
import traceback
from typing import TextIO
import yaml

import eos

try:
    from termcolor import colored
except ImportError:
    colored = lambda s, *args, **kwargs: s

# return the value of the environment variable, or a default value if the variable is unset.
def get_from_env(envvar, default):
    if not envvar in os.environ:
        return default
    if envvar == "EOS_VERBOSITY":
        return int(os.environ[envvar])

    return os.environ[envvar]


def _parser():
    parser = argparse.ArgumentParser(description='Create figures using EOS.')
    # 'parent' parser for common arguments
    common_subparser = argparse.ArgumentParser(add_help=False)
    # add verbosity arg to all commands
    common_subparser.add_argument('-v', '--verbose',
        help = 'Increases the verbosity of the script. Can also be set via the EOS_VERBOSITY environment variable.',
        dest = 'verbose', action = 'count', default = None
    )
    subparsers = parser.add_subparsers(title = 'commands')

    ## begin of commands

    # draw
    parser_draw = subparsers.add_parser('draw',
        parents = [common_subparser],
        description = '''
Draws the figure based on the YAML input provided.
''',
        help = 'Draws the figure.'
    )
    parser_draw.add_argument('input_file', metavar='INPUT_FILE',
        help = 'The YAML input file that specifies the figure to be drawn.'
    )
    parser_draw.add_argument('output_file', metavar='OUTPUT_FILE',
        help = 'The output file where the figure shall be stored.'
    )
    parser_draw.set_defaults(cmd = cmd_draw)

    ## end of commands

    return parser


class CustomLogFormatter(logging.Formatter):
    _MAP_LEVEL_TO_COLOR = {
        logging.ERROR:      ('✖', 'red'),
        logging.WARNING:    ('⚠', 'yellow'),
        logging.SUCCESS:    ('🗸', 'green'),
        logging.COMPLETED:  ('🗸', 'green'),
        logging.INPROGRESS: ('…', 'green'),
        logging.INFO:       ('ℹ', 'blue'),
        logging.DEBUG:      ('🤖', None),
    }
    def __init__(self):
        super().__init__(fmt='%(levelname)s %(message)s', datefmt=None, style='%')

    def format(self, record):
        levelno = record.levelno
        if record.levelno not in self._MAP_LEVEL_TO_COLOR:
            levelno = logging.ERROR

        symbol, color = self._MAP_LEVEL_TO_COLOR[record.levelno]
        record.levelname = colored(symbol, color, attrs=['bold'])
        record.msg = colored(record.msg, color)

        return super().format(record)

_LOG_LEVELS = {
    0: logging.ERROR,
    1: logging.WARNING,
    2: logging.SUCCESS,
    3: logging.INPROGRESS,
    4: logging.INFO,
    5: logging.DEBUG
}


@contextmanager
def _configured_logging(verbosity, stderr):
    handler            = eos.default_log_handler
    previous_level     = handler.level
    previous_stream    = handler.stream
    previous_formatter = handler.formatter
    eos.set_log_level(_LOG_LEVELS[verbosity])
    handler.setStream(stderr)
    handler.setFormatter(CustomLogFormatter())
    try:
        yield
    finally:
        handler.setFormatter(previous_formatter)
        handler.setStream(previous_stream)
        eos.set_log_level(previous_level)


def main(argv: Sequence[str] | None = None, *, stdout: TextIO | None = None, stderr: TextIO | None = None) -> int:
    stdout = sys.stdout if stdout is None else stdout
    stderr = sys.stderr if stderr is None else stderr
    parser = _parser()
    argcomplete.autocomplete(parser)
    try:
        with redirect_stdout(stdout), redirect_stderr(stderr):
            args = parser.parse_args(argv)
    except SystemExit as e:
        return int(e.code)

    if not hasattr(args, 'cmd') or not callable(args.cmd):
        parser.print_help(file=stdout)
        return 0

    if not args.verbose:
        args.verbose = get_from_env('EOS_VERBOSITY', 0)
    if args.verbose > 5:
        args.verbose = 5

    with _configured_logging(args.verbose, stderr), redirect_stdout(stdout), redirect_stderr(stderr):
        try:
            args.cmd(args)
        except Exception as e:
            print(colored('✖ Encountered an unrecoverable error:\n', 'red', attrs=['bold']), f'{e}', file=stdout)
            if not isinstance(e, ValueError):
                traceback.print_exception(e, e, e.__traceback__)
            return 1

    return 0


# Draw command
def cmd_draw(args):
    """Draw the figure based on the YAML input provided."""
    from eos.figure import FigureFactory

    with open(args.input_file) as f:
        yaml_input = yaml.safe_load(f)

    figure = FigureFactory.from_dict(**yaml_input)
    figure.draw()
    figure.save(args.output_file)
