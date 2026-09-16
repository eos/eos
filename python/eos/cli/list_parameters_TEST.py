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

import io
import os
from pathlib import Path
import re
import unittest

import eos
from eos.cli import list_parameters


_FIXTURES = Path(os.environ.get('SOURCE_DIR', Path(__file__).parents[2])) / 'eos/cli/list_parameters_TEST.d'

KNOWN_OBSERVABLE = 'B->Dlnu::dBR/dq2'
KNOWN_PARAMETER  = 'mass::B_d'

_PARAMETER_LINE = re.compile(r'^\s*(?P<name>\S+)\t(?P<min>[-+]\S+)\t(?P<value>[-+]\S+)\t(?P<max>[-+]\S+)$')
_SCAN_LINE      = re.compile(r'^    --scan\t"(?P<name>[^"]+)"\s*\t.*--prior\t(?P<prior>\S+)')


def run(*argv):
    """Runs one command line and returns its exit status together with the captured output."""
    stdout, stderr = io.StringIO(), io.StringIO()
    status = list_parameters.main([*argv], stdout=stdout, stderr=stderr)
    return status, stdout.getvalue(), stderr.getvalue()


def fixture(name):
    """Returns the lines of the expected output stored in the fixture of the given name."""
    return (_FIXTURES / name).read_text().rstrip('\n').splitlines()


def parameter_names(stdout):
    """Returns the names on the parameter lines of a listing."""
    return [m.group('name') for m in map(_PARAMETER_LINE.match, stdout.splitlines()) if m]


def section_titles(stdout):
    """Returns the section titles of a listing, each framed by a rule of equals signs."""
    lines = stdout.splitlines()
    return [line for i, line in enumerate(lines)
            if i > 0 and lines[i - 1] == '=' * len(line) and len(line) > 0]


def all_parameter_names():
    """Returns the names of every parameter EOS knows about."""
    parameters = eos.Parameters.Defaults()
    return {str(p.name()) for section in parameters.sections() for group in section for p in group}


class ParserTests(unittest.TestCase):

    def test_help(self):
        "'--help' exits successfully and prints the usage message."

        status, stdout, _ = run('--help')

        self.assertEqual(status, 0)
        self.assertIn('usage:', stdout)

    def test_unknown_option(self):
        "An unknown option is rejected by the parser."

        status, _, stderr = run('--no-such-option')

        self.assertEqual(status, 2)
        self.assertIn('usage:', stderr)

    def test_kinematics_needs_two_values(self):
        "'--kinematics' takes a name and a value."

        status, _, stderr = run('--kinematics', 'q2')

        self.assertEqual(status, 2)
        self.assertIn('usage:', stderr)

    def test_kinematics_value_must_be_a_number(self):
        "A non-numeric kinematic value is rejected by the parser."

        status, _, stderr = run('--kinematics', 'q2', 'not-a-number')

        self.assertEqual(status, 2)
        self.assertIn('is not a number', stderr)


class ListParametersTests(unittest.TestCase):

    def test_for_observable(self):
        "Restricting to one observable reproduces the recorded listing."

        status, stdout, _ = run('--kinematics', 'q2', '1.0', '--observable', KNOWN_OBSERVABLE)

        self.assertEqual(status, 0)
        self.assertEqual(stdout.rstrip('\n').splitlines(), fixture('b-to-dlnu.txt'))

    def test_for_observable_is_a_subset(self):
        "An observable lists fewer parameters than the complete listing, and only known ones."

        _, complete, _   = run()
        _, restricted, _ = run('--kinematics', 'q2', '1.0', '--observable', KNOWN_OBSERVABLE)

        names = parameter_names(restricted)
        self.assertGreater(len(names), 0)
        self.assertLess(len(names), len(parameter_names(complete)))
        self.assertTrue(set(names).issubset(all_parameter_names()))

    def test_without_observable(self):
        "Without an observable every parameter is listed exactly once."

        status, stdout, _ = run()

        self.assertEqual(status, 0)
        names = parameter_names(stdout)
        self.assertEqual(sorted(names), sorted(all_parameter_names()))
        self.assertIn(KNOWN_PARAMETER, names)

    def test_several_observables(self):
        "Several observables select the union of the parameters they depend on."

        _, one, _  = run('--kinematics', 'q2', '1.0', '--observable', KNOWN_OBSERVABLE)
        _, two, _  = run('--kinematics', 'q2', '1.0', '--observable', KNOWN_OBSERVABLE,
                         '--kinematics', 'q2', '1.0', '--observable', 'B->pilnu::dBR/dq2')

        self.assertTrue(set(parameter_names(one)).issubset(set(parameter_names(two))))
        self.assertGreater(len(parameter_names(two)), len(parameter_names(one)))

    def test_unknown_observable(self):
        "An unknown observable is reported and yields a non-zero exit status."

        status, stdout, _ = run('--observable', 'B->Dlnu::no-such-observable')

        self.assertEqual(status, 1)
        self.assertIn('unrecoverable error', stdout)

    def test_sections_are_sorted(self):
        "The sections are listed by name, since the order in which EOS reads their files is arbitrary."

        status, stdout, _ = run()

        self.assertEqual(status, 0)
        titles = section_titles(stdout)
        self.assertGreater(len(titles), 1)
        self.assertEqual(titles, sorted(titles))

    def test_sections_and_groups(self):
        "Every section and group of the parameter file heads its own block."

        status, stdout, _ = run()

        self.assertEqual(status, 0)
        parameters = eos.Parameters.Defaults()
        for section in parameters.sections():
            self.assertIn(f'{"=" * len(section.name())}\n{section.name()}\n', stdout)
            for group in section:
                self.assertIn(f'{group.name()}\n{"-" * len(group.name())}\n', stdout)


class ScanFormatTests(unittest.TestCase):

    def test_one_line_per_parameter(self):
        "The scan format prints one continued line per parameter, the filters do not apply."

        status, stdout, _ = run('--scan-format')

        self.assertEqual(status, 0)
        lines = stdout.splitlines()
        self.assertEqual(len(lines), len(all_parameter_names()))
        for line in lines:
            self.assertTrue(line.endswith(' \\'))

    def test_names_and_priors(self):
        "Every scan line names a known parameter and one of the three priors."

        status, stdout, _ = run('--scan-format')

        self.assertEqual(status, 0)
        names = all_parameter_names()
        matches = [_SCAN_LINE.match(line) for line in stdout.splitlines()]
        self.assertTrue(all(matches))
        for m in matches:
            self.assertIn(m.group('name'), names)
            self.assertIn(m.group('prior'), ('flat', 'gaussian', 'log-gamma'))

    def test_flat_parameters_carry_placeholders(self):
        "A parameter without a range is given a flat prior and MIN/MAX placeholders."

        status, stdout, _ = run('--scan-format')

        self.assertEqual(status, 0)
        flat = [line for line in stdout.splitlines() if '--prior\tflat' in line]
        self.assertGreater(len(flat), 0)
        for line in flat:
            self.assertIn('MIN', line)
            self.assertIn('MAX', line)


if __name__ == '__main__':
    unittest.main(verbosity=5)
