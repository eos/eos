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
import functools
import math
import sys
import traceback
from typing import TextIO

import eos

try:
    from termcolor import colored
except ImportError:
    colored = lambda s, *args, **kwargs: s

_MIN_NAME_LENGTH = 20
# how big to choose the range around the mode
_NUMBER_OF_SIGMAS = 2.0
_MAX_PRIOR_LENGTH = len('log-gamma')


class _KinematicsAction(argparse.Action):
    """Collect a kinematic variable and its value for the next observable."""

    def __call__(self, parser, namespace, values, option_string=None):
        name, value = values
        try:
            namespace.kinematics[name] = float(value)
        except ValueError:
            raise argparse.ArgumentError(self, f"value '{value}' of kinematic variable '{name}' is not a number")


class _ObservableAction(argparse.Action):
    """Instantiate an observable with the kinematics collected so far, then start afresh."""

    def __call__(self, parser, namespace, values, option_string=None):
        namespace.observables.append((values, dict(namespace.kinematics)))
        namespace.kinematics.clear()


def _parser():
    parser = argparse.ArgumentParser(
        description='List the parameters used by EOS, or those a given observable depends on'
    )
    parser.add_argument('-k', '--kinematics', metavar=('NAME', 'VALUE'), nargs=2,
        help = 'Set a kinematic variable for the next observable; repeat for several variables.',
        dest = 'kinematics', action = _KinematicsAction, default = {}
    )
    parser.add_argument('-o', '--observable', metavar='NAME',
        help = 'List only the parameters this observable depends on; repeat for several observables.',
        dest = 'observables', action = _ObservableAction, default = []
    )
    parser.add_argument('-s', '--scan-format',
        help = 'Print the parameters as a list of scan ranges and priors.',
        dest = 'scan_format', action = 'store_true', default = False
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
            list_parameters(args)
        except Exception as e:
            print(colored('✖ Encountered an unrecoverable error:\n', 'red', attrs=['bold']), f'{e}', file=stdout)
            if not isinstance(e, ValueError):
                traceback.print_exception(e, e, e.__traceback__)
            return 1

    return 0


def _compare(x, y):
    """Order two parameter names, qualified names first and by their parts."""
    x_is_qualified_name = False
    try:
        qnx = eos.QualifiedName(x)
        x_is_qualified_name = True
        qny = eos.QualifiedName(y)

        for part in ('prefix_part', 'suffix_part', 'name_part'):
            px, py = str(getattr(qnx, part)()), str(getattr(qny, part)())
            if px != py:
                return -1 if px < py else 1

        return 0
    except RuntimeError:
        if x_is_qualified_name:
            return 1

        return -1 if x < y else (0 if x == y else 1)


def used_parameter_names(parameters, observables):
    """Return the names of the parameters the given observables depend on."""
    names = set()

    for name, kinematics in observables:
        observable = eos.Observable.make(name, parameters, eos.Kinematics(**kinematics), eos.Options())
        if observable is None:
            raise ValueError(f"Unknown observable '{name}'")

        names.update(str(parameters.by_id(i).name()) for i in observable.used_parameter_ids())

    return names


def print_scan_format(parameters, max_name_length):
    """Print every parameter as a scan range with a matching prior."""
    for p in parameters:
        delta_down = p.evaluate() - p.min()
        delta_up   = p.max() - p.evaluate()

        if delta_down == 0 and delta_up == 0:
            prior = 'flat'
            minimum, maximum = 'MIN\t', 'MAX\t'
        else:
            # for asymmetries larger than 5%, use log-gamma
            ratio = delta_up / delta_down if delta_down != 0.0 else math.inf
            prior = 'gaussian' if math.fabs(ratio - 1.0) < 0.05 else 'log-gamma'
            minimum = f'{p.evaluate() - _NUMBER_OF_SIGMAS * delta_down:.4g}'
            maximum = f'{p.evaluate() + _NUMBER_OF_SIGMAS * delta_up:.4g}'

        name = f'"{p.name()}"'
        line = f'    --scan\t{name:<{max_name_length}}\t{minimum}\t{maximum}\t--prior\t{prior:<{_MAX_PRIOR_LENGTH}}'

        if prior != 'flat':
            line += f'\t{"":<6}{p.min():+.4g}\t{p.evaluate():+.4g}\t{p.max():+.4g}'

        print(f'{line} \\')


def sections(parameters):
    """Return the parameter sections by name; the order EOS reads their files in is arbitrary."""
    return sorted(parameters.sections(), key=lambda section: section.name())


def print_sections(parameters, names, max_name_length):
    """Print the selected parameters, grouped by section and group."""
    for section in sections(parameters):
        title = section.name()
        print('=' * len(title))
        print(title)
        print('=' * len(title))
        print()

        for group in section:
            title = group.name()
            print(title)
            print('-' * len(title))
            print()

            # nasty hack to sort
            # TODO: remove entirely once all parameter names are QualifiedName compatible
            for p in sorted(group, key=functools.cmp_to_key(lambda x, y: _compare(x.name(), y.name()))):
                if p.name() not in names:
                    continue

                print(f'{p.name():>{max_name_length}}\t{p.min():+.4e}\t{p.evaluate():+.4e}\t{p.max():+.4e}')

            print()


def list_parameters(args):
    """List the parameters known to EOS, optionally restricted to those of given observables."""
    parameters = eos.Parameters.Defaults()

    if args.observables:
        names = used_parameter_names(parameters, args.observables)
    else:
        names = {str(p.name()) for section in sections(parameters) for group in section for p in group}

    max_name_length = max([_MIN_NAME_LENGTH, *(len(name) for name in names)])

    if args.scan_format:
        print_scan_format((p for section in sections(parameters) for group in section for p in group), max_name_length)
        return

    print_sections(parameters, names, max_name_length)
