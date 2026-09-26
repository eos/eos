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

# These tests cover the 'eos-analysis' commands that make up the pypmc-based pipeline:
# 'sample-mcmc', 'find-clusters', 'sample-pmc' and 'find-mode --from-mcmc'.

import io
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

POSTERIOR = 'CKM+FF'
CHAINS    = 2
MCMC_N    = 200
PMC_N     = 1000

_BASE = None


def run(*argv, base=None):
    """Runs one command line and returns its exit status together with the captured output."""
    stdout, stderr = io.StringIO(), io.StringIO()
    argv = [*argv, '-f', ANALYSIS_FILE]
    if base is not None:
        argv += ['-b', base]
    status = analysis.main(argv, stdout=stdout, stderr=stderr)
    return status, stdout.getvalue(), stderr.getvalue()


def setUpModule():
    """Runs the MCMC, clustering, and PMC commands once for the tests that inspect their results."""
    global _BASE
    _BASE = tempfile.mkdtemp()
    for chain in range(CHAINS):
        status, stdout, _ = run('sample-mcmc', POSTERIOR, str(chain), '-N', str(MCMC_N), base=_BASE)
        assert status == 0, f'sample-mcmc failed:\n{stdout}'
    status, stdout, _ = run('find-clusters', POSTERIOR, base=_BASE)
    assert status == 0, f'find-clusters failed:\n{stdout}'
    status, stdout, _ = run('sample-pmc', POSTERIOR, '-n', str(PMC_N), '-s', '3', '-N', str(PMC_N), base=_BASE)
    assert status == 0, f'sample-pmc failed:\n{stdout}'


def tearDownModule():
    shutil.rmtree(_BASE, ignore_errors=True)


def data_path(*parts):
    return os.path.join(_BASE, 'data', POSTERIOR, *parts)


class SampleMCMCTests(unittest.TestCase):

    def test_artifacts(self):
        "Each chain records the requested number of finite samples."

        for chain in range(CHAINS):
            mc = eos.data.MarkovChain(data_path(f'mcmc-{chain:04}'))
            self.assertEqual(mc.samples.shape, (MCMC_N, 6))
            self.assertTrue(np.all(np.isfinite(mc.samples)))

    def test_unknown_posterior(self):
        "An unknown posterior is reported as an error."

        with tempfile.TemporaryDirectory() as base:
            status, stdout, _ = run('sample-mcmc', 'NOSUCH', '0', '-N', '10', base=base)

        self.assertEqual(status, 1)
        self.assertIn('NOSUCH', stdout)

    def test_arguments_are_forwarded(self):
        "Every option reaches the task under its documented name."

        with mock.patch.object(eos, 'sample_mcmc') as task:
            status, _, _ = run('sample-mcmc', POSTERIOR, '3', '-N', '7', '-S', '2', '-p', '4', '-n', '9',
                               '-s', '0.1,0.2', '-c', '0.5', base='/base')

        self.assertEqual(status, 0)
        task.assert_called_once_with(
            analysis_file=ANALYSIS_FILE, posterior=POSTERIOR, chain=3, base_directory='/base',
            N=7, stride=2, preruns=4, pre_N=9, start_point=[0.1, 0.2], cov_scale=0.5,
        )


class FindClustersTests(unittest.TestCase):

    def test_artifacts(self):
        "The clusters form a normalized mixture density over the varied parameters."

        density = eos.data.MixtureDensity(data_path('clusters'))

        self.assertGreater(len(density.components), 0)
        self.assertAlmostEqual(sum(density.weights), 1.0)
        self.assertEqual(len(density.varied_parameters), 6)

    def test_missing_chains(self):
        "Clustering without MCMC samples is reported as an error."

        with tempfile.TemporaryDirectory() as base:
            status, _, _ = run('find-clusters', POSTERIOR, base=base)

        self.assertEqual(status, 1)

    def test_arguments_are_forwarded(self):
        "Every option reaches the task under its documented name."

        with mock.patch.object(eos, 'find_clusters') as task:
            status, _, _ = run('find-clusters', POSTERIOR, '-t', '1.5', '-c', '2', base='/base')

        self.assertEqual(status, 0)
        task.assert_called_once_with(
            analysis_file=ANALYSIS_FILE, posterior=POSTERIOR, base_directory='/base',
            threshold=1.5, K_g=2,
        )


class SamplePMCTests(unittest.TestCase):

    def test_artifacts(self):
        "The PMC run records its proposal and the final importance samples."

        sampler = eos.data.PMCSampler(data_path('pmc'))
        samples = eos.data.ImportanceSamples(data_path('samples'))

        self.assertGreater(len(sampler.components), 0)
        self.assertEqual(samples.samples.shape, (PMC_N, 6))
        self.assertEqual(samples.weights.shape, (PMC_N,))
        self.assertTrue(np.all(np.isfinite(samples.samples)))
        self.assertTrue(np.all(samples.weights >= 0.0))

    def test_missing_clusters(self):
        "Sampling without an initial proposal is reported as an error."

        with tempfile.TemporaryDirectory() as base:
            status, _, _ = run('sample-pmc', POSTERIOR, base=base)

        self.assertEqual(status, 1)

    def test_defaults(self):
        "The PMC arguments carry the documented defaults."

        args = analysis._parser().parse_args(['sample-pmc', POSTERIOR])

        self.assertEqual(args.step_N, 500)
        self.assertEqual(args.steps, 10)
        self.assertEqual(args.final_N, 5000)
        self.assertEqual(args.initial_proposal, 'clusters')
        self.assertIsNone(args.sigma_test_stat)


class FindModeFromMCMCTests(unittest.TestCase):

    def test_from_mcmc(self):
        "Starting from the most probable sample of a Markov chain yields a mode."

        status, stdout, _ = run('find-mode', POSTERIOR, '-c', '1', '-o', '1', '-L', 'from-mcmc', base=_BASE)

        self.assertEqual(status, 0, stdout)
        mc   = eos.data.MarkovChain(data_path('mcmc-0001'))
        mode = eos.data.Mode(data_path('mode-from-mcmc'))
        self.assertEqual([p['name'] for p in mode.varied_parameters], [p['name'] for p in mc.varied_parameters])
        self.assertTrue(np.all(np.isfinite(mode.mode)))

    def test_missing_chain(self):
        "Starting from a Markov chain that does not exist is reported as an error."

        status, _, _ = run('find-mode', POSTERIOR, '-c', str(CHAINS), base=_BASE)

        self.assertEqual(status, 1)

    def test_arguments_are_forwarded(self):
        "The chain index reaches the task as 'chain'."

        with mock.patch.object(eos.tasks, 'find_mode') as task:
            status, _, _ = run('find-mode', POSTERIOR, '-c', '1', base='/base')

        self.assertEqual(status, 0)
        self.assertEqual(task.call_args.kwargs['chain'], 1)


if __name__ == '__main__':
    unittest.main(verbosity=5)
