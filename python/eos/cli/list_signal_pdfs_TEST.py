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
import unittest

import eos
from eos.cli import list_signal_pdfs


_FIXTURES = Path(os.environ.get('SOURCE_DIR', Path(__file__).parents[2])) / 'eos/cli/list_signal_pdfs_TEST.d'

KNOWN_PREFIX = 'B->Dlnu'
KNOWN_NAME   = 'P(q2,cos(theta_l))'
KNOWN_PDF    = 'B->Dlnu::P(q2)'


def run(*argv):
    """Runs one command line and returns its exit status together with the captured output."""
    stdout, stderr = io.StringIO(), io.StringIO()
    status = list_signal_pdfs.main([*argv], stdout=stdout, stderr=stderr)
    return status, stdout.getvalue(), stderr.getvalue()


def fixture(name):
    """Returns the lines of the expected output stored in the fixture of the given name."""
    return (_FIXTURES / name).read_text().rstrip('\n').splitlines()


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


class ListSignalPDFsTests(unittest.TestCase):

    def test_filter_by_prefix(self):
        "Filtering by a prefix part reproduces the recorded output for those signal PDFs."

        status, stdout, _ = run('-p', KNOWN_PREFIX)

        self.assertEqual(status, 0)
        self.assertEqual(stdout.rstrip('\n').splitlines(), fixture('b-to-dlnu.txt'))

    def test_filter_by_name(self):
        "Filtering by a name part lists only signal PDFs carrying that name."

        status, stdout, _ = run('-n', KNOWN_NAME)

        self.assertEqual(status, 0)
        names = [line for line in stdout.splitlines() if line and not line.startswith(' ')]
        self.assertGreater(len(names), 0)
        for name in names:
            self.assertTrue(name.endswith(f'::{KNOWN_NAME}'))

    def test_filters_are_combined(self):
        "Name and prefix filters select the union of their matches."

        _, by_name, _   = run('-n', KNOWN_NAME)
        _, by_prefix, _ = run('-p', KNOWN_PREFIX)
        _, by_both, _   = run('-n', KNOWN_NAME, '-p', KNOWN_PREFIX)

        names = {line for line in by_both.splitlines() if line and not line.startswith(' ')}
        self.assertEqual(names, {line for line in (by_name + by_prefix).splitlines() if line and not line.startswith(' ')})

    def test_without_filter(self):
        "Without a filter every signal PDF is listed."

        status, stdout, _ = run()

        self.assertEqual(status, 0)
        names = [line for line in stdout.splitlines() if line and not line.startswith(' ')]
        self.assertEqual(len(names), sum(1 for _ in eos.SignalPDFs()))
        self.assertIn(KNOWN_PDF, names)

    def test_unknown_filter(self):
        "A filter that matches no signal PDF lists nothing."

        status, stdout, _ = run('-n', 'no-such-signal-pdf-name')

        self.assertEqual(status, 0)
        self.assertEqual(stdout, '')

    def test_kinematic_variables(self):
        "The kinematic variables of numerator and denominator are listed separately."

        status, stdout, _ = run('-p', KNOWN_PREFIX)

        self.assertEqual(status, 0)
        self.assertIn('    Numerator kinematic variables:q2, cos(theta_l)', stdout)
        self.assertIn('    Denominator kinematic variables:q2_min, q2_max', stdout)


if __name__ == '__main__':
    unittest.main(verbosity=5)
