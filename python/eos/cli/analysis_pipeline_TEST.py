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
import math
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np

import eos

from eos.cli import analysis


_FIXTURES = Path(os.environ.get('SOURCE_DIR', Path(__file__).parents[2])) / 'eos/cli/analysis_TEST.d'

ANALYSIS_FILE = str(_FIXTURES / 'analysis.yaml')

# The posterior used throughout, and the parameters it varies in the order in which the
# analysis file's priors declare them.
POSTERIOR  = 'FF'
PARAMETERS = [
    'B->pi::f_+(0)@BCL2008',
    'B->pi::b_+^1@BCL2008',
    'B->pi::b_+^2@BCL2008',
    'B->pi::b_0^1@BCL2008',
    'B->pi::b_0^2@BCL2008',
]

# Nested-sampling settings that keep each run near a second while still yielding a few
# hundred samples.
NLIVE   = 40
MAXITER = 300
SEED    = 10

_BASE = None


def run(*argv, base=None, analysis_file=ANALYSIS_FILE):
    """Runs one command line and returns its exit status together with the captured output."""
    stdout, stderr = io.StringIO(), io.StringIO()
    argv = [*argv, '-f', analysis_file]
    if base is not None:
        argv += ['-b', base]
    status = analysis.main(argv, stdout=stdout, stderr=stderr)
    return status, stdout.getvalue(), stderr.getvalue()


def sample_nested(base, *argv):
    """Draws nested samples into the given base directory."""
    return run('sample-nested', POSTERIOR, '-n', str(NLIVE), '-m', str(MAXITER), '-s', str(SEED),
               *argv, base=base)


def setUpModule():
    """Draws one set of nested samples for the tests that consume it."""
    global _BASE
    _BASE = tempfile.mkdtemp()
    status, _, _ = sample_nested(_BASE)
    assert status == 0, 'failed to draw the nested samples the pipeline tests build on'


def tearDownModule():
    shutil.rmtree(_BASE, ignore_errors=True)


def data_path(*parts, base=None):
    return os.path.join(_BASE if base is None else base, 'data', POSTERIOR, *parts)


class FindModeTests(unittest.TestCase):

    def mode(self, label):
        return eos.data.Mode(data_path(f'mode-{label}'))

    def test_single_optimization(self):
        "A single call to the optimization algorithm is enough to record a mode."

        status, stdout, _ = run('find-mode', POSTERIOR, '-S', '-o', '1', '-L', 'single', base=_BASE)

        self.assertEqual(status, 0, stdout)
        self.assertEqual(len(self.mode('single').mode), len(PARAMETERS))

    def test_more_optimizations_do_not_worsen_the_mode(self):
        "Iterating the optimization does not record a worse chi^2 than a single call."

        self.assertEqual(run('find-mode', POSTERIOR, '-S', '-o', '1', '-L', 'one', base=_BASE)[0], 0)
        self.assertEqual(run('find-mode', POSTERIOR, '-S', '-o', '3', '-L', 'three', base=_BASE)[0], 0)

        self.assertLessEqual(self.mode('three').global_chi2, self.mode('one').global_chi2 + 1.0e-6)

    def test_non_finite_chi_square(self):
        "An optimization whose chi^2 is not finite still records the point it found."

        real = eos.GoodnessOfFit

        class _NonFiniteGoodnessOfFit:
            "Stands in for eos.GoodnessOfFit, reporting a NaN total chi^2."

            def __init__(self, log_posterior):
                self._gof = real(log_posterior)

            def total_chi_square(self):
                return float('nan')

            def total_degrees_of_freedom(self):
                return self._gof.total_degrees_of_freedom()

            def __iter__(self):
                return iter(self._gof)

        with tempfile.TemporaryDirectory() as base:
            with mock.patch.object(eos, 'GoodnessOfFit', _NonFiniteGoodnessOfFit):
                status, stdout, _ = run('find-mode', POSTERIOR, '-o', '2', '-L', 'nan', base=base)

            self.assertEqual(status, 0, stdout)
            mode = eos.data.Mode(data_path('mode-nan', base=base))
            self.assertEqual(len(mode.mode), len(PARAMETERS))
            self.assertTrue(np.all(np.isfinite(mode.mode)))
            self.assertTrue(math.isnan(mode.global_chi2))


if __name__ == '__main__':
    unittest.main(verbosity=5)
