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

from eos.cli import analysis


_FIXTURES = Path(os.environ.get('SOURCE_DIR', Path(__file__).parents[2])) / 'eos/cli/analysis_TEST.d'

ANALYSIS_FILE           = str(_FIXTURES / 'analysis.yaml')
LIKELIHOOD_DETAILS_FILE = str(_FIXTURES / 'likelihood-details.yaml')


def run(*argv, analysis_file=ANALYSIS_FILE):
    """Runs one command line and returns its exit status together with the captured output."""
    stdout, stderr = io.StringIO(), io.StringIO()
    status = analysis.main([*argv, '-f', analysis_file], stdout=stdout, stderr=stderr)
    return status, stdout.getvalue(), stderr.getvalue()


class ListPriorsTests(unittest.TestCase):

    def test_names(self):
        "The named priors are listed in the order in which the analysis file declares them."

        status, stdout, _ = run('list-priors')

        self.assertEqual(status, 0)
        self.assertEqual(stdout.splitlines(), ['CKM', 'FF-norm', 'FF-shape'])


class ListLikelihoodsTests(unittest.TestCase):

    def test_names(self):
        "Without '--display-details' only the names are listed."

        status, stdout, _ = run('list-likelihoods')

        self.assertEqual(status, 0)
        self.assertEqual(stdout.splitlines(), ['lattice', 'lcsr', 'BaBar', 'Belle'])

    def test_details(self):
        "With '--display-details' each built-in constraint is listed below its likelihood."

        status, stdout, _ = run('list-likelihoods', '-d')

        self.assertEqual(status, 0)
        self.assertEqual(stdout.splitlines(), [
            'lattice',
            ' - B->pi::f_++f_0@RBC+UKQCD:2015A',
            'lcsr',
            ' - B->pi::f_+@KMO:2006A',
            'BaBar',
            ' - B^0->pi^+lnu::BR@BaBar:2010B',
            'Belle',
            ' - B^0->pi^+lnu::BR@Belle:2010A',
        ])

    def test_details_of_every_contribution(self):
        "Manual constraints and pyhf workspaces are listed alongside the built-in constraints."

        status, stdout, _ = run('list-likelihoods', '-d', analysis_file=LIKELIHOOD_DETAILS_FILE)

        self.assertEqual(status, 0)
        self.assertEqual(stdout.splitlines(), [
            'builtin',
            ' - B->pi::f_+@KMO:2006A',
            'manual',
            ' - test::manual-constraint [manual]',
            'with-pyhf',
            ' - B->pi::f_+@KMO:2006A',
            ' - workspace.json [pyhf]',
        ])


class ListPosteriorsTests(unittest.TestCase):

    def test_names(self):
        "The named posteriors are listed in the order in which the analysis file declares them."

        status, stdout, _ = run('list-posteriors')

        self.assertEqual(status, 0)
        self.assertEqual(stdout.splitlines(), ['FF', 'CKM+FF'])


class ListPredictionsTests(unittest.TestCase):

    def test_names(self):
        "The named prediction sets are listed in the order in which the analysis file declares them."

        status, stdout, _ = run('list-predictions')

        self.assertEqual(status, 0)
        self.assertEqual(stdout.splitlines(), ['PDF', 'dBRdq2', 'BR'])


class ListFiguresTests(unittest.TestCase):

    def test_names(self):
        "The named figures are listed in the order in which the analysis file declares them."

        status, stdout, _ = run('list-figures')

        self.assertEqual(status, 0)
        self.assertEqual(stdout.splitlines(), ['CKM-Vub', 'CKM-Vub-v-FF'])


class ListStepsTests(unittest.TestCase):

    def test_ids(self):
        "The steps are listed by id, in the order in which the analysis file declares them."

        status, stdout, _ = run('list-steps')

        self.assertEqual(status, 0)
        self.assertEqual(stdout.splitlines(), ['CKM+FF.find-mode', 'CKM+FF.sample', 'both.sample'])


class ListStepDependenciesTests(unittest.TestCase):

    def test_dependencies(self):
        "The ids a step depends on are listed."

        status, stdout, _ = run('list-step-dependencies', 'CKM+FF.sample')

        self.assertEqual(status, 0)
        self.assertEqual(stdout.splitlines(), ['CKM+FF.find-mode'])

    def test_without_dependencies(self):
        "A step without dependencies produces no output."

        status, stdout, _ = run('list-step-dependencies', 'CKM+FF.find-mode')

        self.assertEqual(status, 0)
        self.assertEqual(stdout, '')

    def test_unknown_id(self):
        "An unknown step id is reported as an error."

        status, stdout, _ = run('list-step-dependencies', 'nosuch.step')

        self.assertEqual(status, 1)
        self.assertIn("Step with id 'nosuch.step' not found", stdout)


class CommonOptionTests(unittest.TestCase):

    def test_default_analysis_file(self):
        "The analysis file defaults to 'analysis.yaml'."

        args = analysis._parser().parse_args(['list-priors'])

        self.assertEqual(args.analysis_file, 'analysis.yaml')
        self.assertIsNone(args.verbose)

    def test_missing_analysis_file(self):
        "A missing analysis file is reported as an error."

        status, stdout, _ = run('list-priors', analysis_file='/nonexistent/analysis.yaml')

        self.assertEqual(status, 1)
        self.assertIn('does not exist', stdout)

    def test_unknown_command(self):
        "An unknown command is rejected by the argument parser."

        stdout, stderr = io.StringIO(), io.StringIO()

        self.assertEqual(analysis.main(['nosuchcommand'], stdout=stdout, stderr=stderr), 2)
        self.assertIn('invalid choice', stderr.getvalue())

    def test_no_command(self):
        "Without a command the help is printed."

        stdout, stderr = io.StringIO(), io.StringIO()

        self.assertEqual(analysis.main([], stdout=stdout, stderr=stderr), 0)
        self.assertIn('usage: ', stdout.getvalue())


if __name__ == '__main__':
    unittest.main(verbosity=5)
