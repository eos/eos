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
import unittest

from eos.cli import list_observables


KNOWN_OBSERVABLE = 'B->Dlnu::dBR/dq2'


def run(*argv):
    """Runs one command line and returns its exit status together with the captured output."""
    stdout, stderr = io.StringIO(), io.StringIO()
    status = list_observables.main([*argv], stdout=stdout, stderr=stderr)
    return status, stdout.getvalue(), stderr.getvalue()


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


class ListObservablesTests(unittest.TestCase):

    def test_known_observable(self):
        "Filtering by a single qualified name lists exactly that observable."

        status, stdout, _ = run(KNOWN_OBSERVABLE)

        self.assertEqual(status, 0)
        self.assertEqual(stdout.splitlines(), [KNOWN_OBSERVABLE])

    def test_without_filter(self):
        "Without a filter every observable is listed, among them the known one."

        status, stdout, _ = run()

        self.assertEqual(status, 0)
        lines = stdout.splitlines()
        self.assertGreater(len(lines), 1)
        self.assertIn(KNOWN_OBSERVABLE, lines)

    def test_unknown_observable(self):
        "A qualified name that no observable carries lists nothing."

        status, stdout, _ = run('B->Dlnu::no-such-observable')

        self.assertEqual(status, 0)
        self.assertEqual(stdout, '')

    def test_several_names(self):
        "Several qualified names may be given at once."

        status, stdout, _ = run(KNOWN_OBSERVABLE, 'B->pilnu::dBR/dq2')

        self.assertEqual(status, 0)
        self.assertEqual(sorted(stdout.splitlines()), sorted([KNOWN_OBSERVABLE, 'B->pilnu::dBR/dq2']))


if __name__ == '__main__':
    unittest.main(verbosity=5)
