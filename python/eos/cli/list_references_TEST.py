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

from eos.cli import list_references


_FIXTURES = Path(os.environ.get('SOURCE_DIR', Path(__file__).parents[2])) / 'eos/cli/list_references_TEST.d'

KNOWN_REFERENCE  = 'BGJvD:2019A'
LEGACY_REFERENCE = 'BFS:2001A'


def run(*argv):
    """Runs one command line and returns its exit status together with the captured output."""
    stdout, stderr = io.StringIO(), io.StringIO()
    status = list_references.main([*argv], stdout=stdout, stderr=stderr)
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


class ListReferencesTests(unittest.TestCase):

    def test_known_reference(self):
        "Filtering by a single key reproduces the recorded output for that reference."

        status, stdout, _ = run(KNOWN_REFERENCE)

        self.assertEqual(status, 0)
        self.assertEqual(stdout.rstrip('\n').splitlines(), fixture('bgjvd-2019a.txt'))

    def test_legacy_eprint_id(self):
        "A pre-2007 eprint id keeps its archive prefix in the arXiv identifier."

        status, stdout, _ = run(LEGACY_REFERENCE)

        self.assertEqual(status, 0)
        self.assertEqual(stdout.rstrip('\n').splitlines(), fixture('bfs-2001a.txt'))

    def test_without_filter(self):
        "Without a filter every reference is listed, among them the known one."

        status, stdout, _ = run()

        self.assertEqual(status, 0)
        lines = stdout.splitlines()
        self.assertGreater(len(lines), 4)
        self.assertIn(f'[{KNOWN_REFERENCE}] : Bordone, M. and Gubernari, N. and van Dyk, D. and Jung, M.', lines)

    def test_filter_is_exclusive(self):
        "Filtering by one key suppresses every other reference."

        _, filtered, _ = run(KNOWN_REFERENCE)
        _, complete, _ = run()

        self.assertLess(len(filtered.splitlines()), len(complete.splitlines()))
        self.assertEqual(len([line for line in filtered.splitlines() if line.startswith('[')]), 1)

    def test_malformed_key(self):
        "A key that is not a valid reference name is reported and yields a non-zero exit status."

        status, stdout, _ = run('Asatryan:2001zw')

        self.assertEqual(status, 1)
        self.assertIn('unrecoverable error', stdout)


if __name__ == '__main__':
    unittest.main(verbosity=5)
