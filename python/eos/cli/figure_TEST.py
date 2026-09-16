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
import tempfile
import unittest

import matplotlib
matplotlib.use('Agg')

from eos.cli import figure


_FIXTURES = Path(os.environ.get('SOURCE_DIR', Path(__file__).parents[2])) / 'eos/cli/figure_TEST.d'

SINGLE_FILE            = str(_FIXTURES / 'single.yaml')
UNKNOWN_TYPE_FILE      = str(_FIXTURES / 'unknown-type.yaml')
INVALID_STRUCTURE_FILE = str(_FIXTURES / 'invalid-structure.yaml')


def run(*argv):
    """Runs one command line and returns its exit status together with the captured output."""
    stdout, stderr = io.StringIO(), io.StringIO()
    status = figure.main([*argv], stdout=stdout, stderr=stderr)
    return status, stdout.getvalue(), stderr.getvalue()


class ParserTests(unittest.TestCase):

    def test_without_command(self):
        "Without a command the usage message is printed on stdout."

        status, stdout, _ = run()

        self.assertEqual(status, 0)
        self.assertIn('usage:', stdout)
        self.assertIn('draw', stdout)

    def test_help(self):
        "'--help' exits successfully and prints the usage message."

        status, stdout, _ = run('--help')

        self.assertEqual(status, 0)
        self.assertIn('usage:', stdout)

    def test_unknown_command(self):
        "An unknown command is rejected by the parser."

        status, _, stderr = run('no-such-command')

        self.assertEqual(status, 2)
        self.assertIn('usage:', stderr)

    def test_no_base_directory(self):
        "'draw' does not accept a base directory; that is analysis-file territory."

        status, _, stderr = run('draw', '-b', './', SINGLE_FILE, 'figure.pdf')

        self.assertEqual(status, 2)
        self.assertIn('unrecognized arguments', stderr)


class DrawTests(unittest.TestCase):

    def test_output_file_is_written(self):
        "'draw' renders the figure described by the input file into the output file."

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'figure.pdf'

            status, _, _ = run('draw', SINGLE_FILE, str(output))

            self.assertEqual(status, 0)
            self.assertTrue(output.is_file())
            self.assertGreater(output.stat().st_size, 0)

    def test_verbosity_is_accepted(self):
        "'-v' is accepted and does not change the outcome."

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'figure.pdf'

            status, _, _ = run('draw', '-v', SINGLE_FILE, str(output))

            self.assertEqual(status, 0)
            self.assertTrue(output.is_file())

    def test_missing_input_file(self):
        "A missing input file is reported and yields a non-zero exit status."

        with tempfile.TemporaryDirectory() as directory:
            missing = str(Path(directory) / 'no-such-file.yaml')
            output  = Path(directory) / 'figure.pdf'

            status, stdout, _ = run('draw', missing, str(output))

            self.assertEqual(status, 1)
            self.assertIn('unrecoverable error', stdout)
            self.assertFalse(output.exists())

    def test_malformed_input_file(self):
        "An input file that is not valid YAML is reported and yields a non-zero exit status."

        with tempfile.TemporaryDirectory() as directory:
            malformed = Path(directory) / 'malformed.yaml'
            malformed.write_text("type: 'single'\nplot:\n  items: [\n")
            output = Path(directory) / 'figure.pdf'

            status, stdout, _ = run('draw', str(malformed), str(output))

            self.assertEqual(status, 1)
            self.assertIn('unrecoverable error', stdout)
            self.assertFalse(output.exists())

    def test_unknown_figure_type(self):
        "An unknown figure type is reported and yields a non-zero exit status."

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'figure.pdf'

            status, stdout, _ = run('draw', UNKNOWN_TYPE_FILE, str(output))

            self.assertEqual(status, 1)
            self.assertIn('Unknown figure type', stdout)
            self.assertFalse(output.exists())

    def test_invalid_structure(self):
        "An input file with an unknown key is reported and yields a non-zero exit status."

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'figure.pdf'

            status, stdout, _ = run('draw', INVALID_STRUCTURE_FILE, str(output))

            self.assertEqual(status, 1)
            self.assertIn('unrecoverable error', stdout)
            self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main(verbosity=5)
