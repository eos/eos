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

# These tests cover the four commands that make up the typical 'eos-analysis' pipeline:
# 'validate', 'sample-nested', 'find-mode' and 'predict-observables'.

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
VALID_FILE    = str(_FIXTURES / 'valid.yaml')

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


class ValidateTests(unittest.TestCase):

    def test_valid_file(self):
        "A well-formed analysis file yields no diagnostics."

        status, stdout, _ = run('validate', analysis_file=VALID_FILE)

        self.assertEqual(status, 0)
        self.assertEqual(stdout, '')

    def test_valid_file_without_deep_phase(self):
        "'--no-deep' also yields no diagnostics for a well-formed analysis file."

        status, stdout, _ = run('validate', '--no-deep', analysis_file=VALID_FILE)

        self.assertEqual(status, 0)
        self.assertEqual(stdout, '')

    def test_structural_error(self):
        "A missing mandatory key is reported with the path at which it is missing."

        status, stdout, _ = run('validate', analysis_file=str(_FIXTURES / 'invalid-structure.yaml'))

        self.assertEqual(stdout.splitlines(), [
            "priors/FF/descriptions[0]/max: Missing mandatory key 'max'",
        ])

    def test_semantic_error(self):
        "A reference to an undefined prior is reported."

        status, stdout, _ = run('validate', analysis_file=str(_FIXTURES / 'invalid-semantics.yaml'))

        self.assertEqual(stdout.splitlines(), [
            "posteriors/FF/prior[0]: Posterior 'FF' references prior 'NONEXISTENT' which is not defined",
        ])

    def test_unused_entities(self):
        "Priors and likelihoods that no posterior uses are reported as warnings."

        status, stdout, _ = run('validate', analysis_file=str(_FIXTURES / 'unused-entities.yaml'))

        self.assertEqual(stdout.splitlines(), [
            "priors/unused-prior: Prior 'unused-prior' is unused",
            "likelihoods/unused-likelihood: Likelihood 'unused-likelihood' is unused",
        ])

    def test_no_deep_reports_the_same_errors(self):
        "The structural and semantic phases run irrespective of '--no-deep'."

        for fixture in ('invalid-structure.yaml', 'invalid-semantics.yaml', 'unused-entities.yaml'):
            with self.subTest(fixture=fixture):
                path = str(_FIXTURES / fixture)
                _, deep, _    = run('validate', analysis_file=path)
                _, shallow, _ = run('validate', '--no-deep', analysis_file=path)
                self.assertEqual(deep, shallow)

    def test_exit_status_does_not_report_validity(self):
        "The exit status is 0 even for an analysis file that cannot be loaded."

        for fixture in ('valid.yaml', 'invalid-structure.yaml', 'invalid-semantics.yaml'):
            with self.subTest(fixture=fixture):
                status, _, _ = run('validate', analysis_file=str(_FIXTURES / fixture))
                self.assertEqual(status, 0)

    def test_missing_file(self):
        "A missing analysis file is reported as an error."

        status, stdout, _ = run('validate', analysis_file='/nonexistent/analysis.yaml')

        self.assertEqual(status, 1)
        self.assertIn('does not exist', stdout)

    def test_defaults(self):
        "The deep phase is enabled unless '--no-deep' is given."

        self.assertTrue (analysis._parser().parse_args(['validate']).deep)
        self.assertFalse(analysis._parser().parse_args(['validate', '--no-deep']).deep)


class SampleNestedTests(unittest.TestCase):

    def test_artifacts(self):
        "Both the dynesty results and the exported importance samples are written."

        results = eos.data.DynestyResults(data_path('nested'))
        samples = eos.data.ImportanceSamples(data_path('samples'))

        self.assertEqual(results.type, 'DynestyResults')
        self.assertEqual([p['name'] for p in results.varied_parameters], PARAMETERS)
        self.assertEqual(samples.type, 'ImportanceSamples')
        self.assertEqual([p['name'] for p in samples.varied_parameters], PARAMETERS)
        self.assertEqual(samples.samples.shape[1], len(PARAMETERS))
        self.assertEqual(samples.weights.shape, (samples.samples.shape[0],))
        self.assertEqual(samples.posterior_values.shape, (samples.samples.shape[0],))

    def test_samples_are_finite_and_within_the_prior(self):
        "Every sample lies within the bounds the priors declare."

        samples = eos.data.ImportanceSamples(data_path('samples'))

        self.assertTrue(np.all(np.isfinite(samples.samples)))
        self.assertTrue(np.all(np.isfinite(samples.weights)))
        self.assertTrue(np.all(samples.samples[:, 0] >= 0.15))
        self.assertTrue(np.all(samples.samples[:, 0] <= 0.30))
        self.assertTrue(np.all(np.abs(samples.samples[:, 1:]) <= 4.0))

    def test_seed_is_reproducible(self):
        "Two runs with the same seed yield identical samples."

        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            self.assertEqual(sample_nested(first)[0], 0)
            self.assertEqual(sample_nested(second)[0], 0)

            np.testing.assert_array_equal(
                eos.data.ImportanceSamples(data_path('samples', base=first)).samples,
                eos.data.ImportanceSamples(data_path('samples', base=second)).samples,
            )

    def test_seed_changes_the_samples(self):
        "A different seed yields different samples."

        with tempfile.TemporaryDirectory() as base:
            self.assertEqual(sample_nested(base, '-s', str(SEED + 1))[0], 0)

            other = eos.data.ImportanceSamples(data_path('samples', base=base)).samples
            shared = eos.data.ImportanceSamples(data_path('samples')).samples
            self.assertFalse(other.shape == shared.shape and np.array_equal(other, shared))

    def test_maxiter_limits_the_run(self):
        "A smaller iteration limit yields fewer samples."

        with tempfile.TemporaryDirectory() as base:
            status, _, _ = run('sample-nested', POSTERIOR, '-n', str(NLIVE), '-m', '100',
                               '-s', str(SEED), base=base)

            self.assertEqual(status, 0)
            self.assertLess(
                len(eos.data.ImportanceSamples(data_path('samples', base=base)).samples),
                len(eos.data.ImportanceSamples(data_path('samples')).samples),
            )

    def test_unknown_posterior(self):
        "An unknown posterior is reported as an error."

        with tempfile.TemporaryDirectory() as base:
            status, stdout, _ = run('sample-nested', 'NOSUCH', '-n', str(NLIVE), base=base)

            self.assertEqual(status, 1)
            self.assertIn('NOSUCH', stdout)

    def test_defaults(self):
        "The sampler's arguments carry the documented defaults."

        args = analysis._parser().parse_args(['sample-nested', POSTERIOR])

        self.assertEqual(args.posterior, POSTERIOR)
        self.assertEqual(args.bound,   'multi')
        self.assertEqual(args.nlive,   250)
        self.assertEqual(args.dlogz,   1.0)
        self.assertIsNone(args.maxiter)
        self.assertEqual(args.miniter, 0)
        self.assertEqual(args.min_ess, 0)
        self.assertEqual(args.seed,    10)
        self.assertEqual(args.sample,  'auto')

    def test_arguments_are_forwarded(self):
        "Every option reaches the task under its documented name."

        with mock.patch.object(eos, 'sample_nested') as task:
            status, _, _ = run(
                'sample-nested', POSTERIOR, '-B', 'single', '-n', '17', '-d', '0.5',
                '-m', '99', '-l', '11', '-e', '13', '-s', '23', '-M', 'rwalk', base='/base',
            )

        self.assertEqual(status, 0)
        task.assert_called_once_with(
            analysis_file=ANALYSIS_FILE, posterior=POSTERIOR, base_directory='/base',
            bound='single', nlive=17, dlogz=0.5, maxiter=99, miniter=11, min_ess=13, seed=23,
            sample='rwalk',
        )


class FindModeTests(unittest.TestCase):

    def mode(self, label):
        return eos.data.Mode(data_path(f'mode-{label}'))

    def test_from_samples(self):
        "Starting from the importance samples yields a mode over the varied parameters."

        status, _, _ = run('find-mode', POSTERIOR, '-S', '-L', 'from-samples', base=_BASE)

        self.assertEqual(status, 0)
        mode = self.mode('from-samples')
        self.assertEqual(len(mode.mode), len(PARAMETERS))
        self.assertEqual([p['name'] for p in mode.varied_parameters], PARAMETERS)
        self.assertTrue(np.all(np.isfinite(mode.mode)))
        self.assertGreater(mode.dof, 0)
        self.assertGreaterEqual(mode.pvalue, 0.0)
        self.assertLessEqual(mode.pvalue, 1.0)

    def test_from_point(self):
        "Starting from an explicit point yields a mode."

        status, _, _ = run('find-mode', POSTERIOR, '-p', '0.25,0.0,0.0,0.0,0.0',
                           '-L', 'from-point', base=_BASE)

        self.assertEqual(status, 0)
        self.assertEqual(len(self.mode('from-point').mode), len(PARAMETERS))

    def test_start_point_is_parsed_as_floats(self):
        "The starting point is a comma-separated list of floats."

        args = analysis._parser().parse_args(['find-mode', POSTERIOR, '-p', '0.25,-1.5,2'])

        self.assertEqual(args.start_point, [0.25, -1.5, 2.0])

    def test_random_start_is_reproducible(self):
        "Two runs from the same random seed yield the same mode."

        self.assertEqual(run('find-mode', POSTERIOR, '-s', '17', '-L', 'seed-a', base=_BASE)[0], 0)
        self.assertEqual(run('find-mode', POSTERIOR, '-s', '17', '-L', 'seed-b', base=_BASE)[0], 0)

        np.testing.assert_allclose(self.mode('seed-a').mode, self.mode('seed-b').mode)

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

    def test_zero_optimizations(self):
        "Requesting no optimization at all is rejected."

        status, stdout, _ = run('find-mode', POSTERIOR, '-S', '-o', '0', base=_BASE)

        self.assertEqual(status, 1)
        self.assertIn('larger than zero', stdout)

    def test_mutually_exclusive_start_points(self):
        "At most one way of choosing the starting point may be given."

        cases = (
            ['-s', '17', '-p', '0.25,0,0,0,0'],
            ['-s', '17', '-c', '0'],
            ['-s', '17', '-S'],
            ['-p', '0.25,0,0,0,0', '-c', '0'],
            ['-p', '0.25,0,0,0,0', '-S'],
            ['-c', '0', '-S'],
        )
        for case in cases:
            with self.subTest(arguments=case):
                status, stdout, _ = run('find-mode', POSTERIOR, *case, base=_BASE)
                self.assertEqual(status, 1)
                self.assertIn('mutually exclusive', stdout)

    def test_mask_requires_importance_samples(self):
        "A mask may only be applied to the importance samples."

        status, stdout, _ = run('find-mode', POSTERIOR, '-M', 'nosuchmask', base=_BASE)

        self.assertEqual(status, 1)
        self.assertIn('mask-name can only be used with importance_samples', stdout)

    def test_unknown_posterior(self):
        "An unknown posterior is reported as an error."

        status, stdout, _ = run('find-mode', 'NOSUCH', base=_BASE)

        self.assertEqual(status, 1)
        self.assertIn('NOSUCH', stdout)

    def test_defaults(self):
        "The optimizer's arguments carry the documented defaults."

        args = analysis._parser().parse_args(['find-mode', POSTERIOR])

        self.assertEqual(args.optimizations, 3)
        self.assertEqual(args.label, 'default')
        self.assertIsNone(args.chain)
        self.assertIsNone(args.importance_samples)
        self.assertIsNone(args.start_point)
        self.assertIsNone(args.seed)
        self.assertIsNone(args.mask_name)

    def test_arguments_are_forwarded(self):
        "Every option reaches the task under its documented name."

        with mock.patch.object(eos.tasks, 'find_mode') as task:
            status, _, _ = run('find-mode', POSTERIOR, '-o', '5', '-S', '-L', 'label',
                               '-M', 'mask', base='/base')

        self.assertEqual(status, 0)
        task.assert_called_once_with(
            analysis_file=ANALYSIS_FILE, posterior=POSTERIOR, base_directory='/base',
            optimizations=5, chain=None, importance_samples=True, start_point=None,
            seed=None, label='label', mask_name='mask',
        )


class PredictObservablesTests(unittest.TestCase):

    def prediction(self, name):
        return eos.data.Prediction(data_path(f'pred-{name}'))

    def test_artifacts(self):
        "One prediction is recorded per sample and per observable."

        status, _, _ = run('predict-observables', POSTERIOR, 'dBRdq2', base=_BASE)

        self.assertEqual(status, 0)
        samples    = eos.data.ImportanceSamples(data_path('samples'))
        prediction = self.prediction('dBRdq2')

        self.assertEqual(prediction.type, 'Prediction')
        self.assertEqual(prediction.samples.shape, (len(samples.samples), 5))
        self.assertEqual(prediction.weights.shape, (len(samples.samples),))
        np.testing.assert_array_equal(prediction.weights, samples.weights)
        self.assertEqual(
            [str(o['name']) for o in prediction.varied_parameters],
            ['B->pilnu::dBR/dq2'] * 5,
        )
        self.assertTrue(np.all(np.isfinite(prediction.samples)))

    def test_fixed_parameters(self):
        "A prediction set may fix parameters of its own."

        status, _, _ = run('predict-observables', POSTERIOR, 'BR', base=_BASE)

        self.assertEqual(status, 0)
        prediction = self.prediction('BR')
        self.assertEqual(prediction.samples.shape[1], 1)
        self.assertTrue(np.all(prediction.samples > 0.0))

    def test_sample_range(self):
        "'--begin-index' and '--end-index' restrict the samples that are used."

        status, _, _ = run('predict-observables', POSTERIOR, 'dBRdq2', '-B', '10', '-E', '20',
                           base=_BASE)

        self.assertEqual(status, 0)
        self.assertEqual(self.prediction('dBRdq2').samples.shape[0], 10)

    def test_mask_excludes_a_sample_range(self):
        "A mask may not be combined with an explicit sample range."

        status, stdout, _ = run('predict-observables', POSTERIOR, 'dBRdq2', '-M', 'mask',
                                '-B', '5', base=_BASE)

        self.assertEqual(status, 1)
        self.assertIn('mutually exclusive', stdout)

    def test_unknown_prediction(self):
        "An unknown prediction set is reported as an error."

        status, stdout, _ = run('predict-observables', POSTERIOR, 'NOSUCH', base=_BASE)

        self.assertEqual(status, 1)
        self.assertIn('NOSUCH', stdout)

    def test_missing_samples(self):
        "Predicting without importance samples is reported as an error."

        with tempfile.TemporaryDirectory() as base:
            status, stdout, _ = run('predict-observables', POSTERIOR, 'dBRdq2', base=base)

            self.assertEqual(status, 1)
            self.assertIn('samples does not exist', stdout)

    def test_defaults(self):
        "The prediction's arguments carry the documented defaults."

        args = analysis._parser().parse_args(['predict-observables', POSTERIOR, 'dBRdq2'])

        self.assertEqual(args.prediction, 'dBRdq2')
        self.assertEqual(args.begin, 0)
        self.assertIsNone(args.end)
        self.assertIsNone(args.mask_name)

    def test_arguments_are_forwarded(self):
        "Every option reaches the task under its documented name."

        with mock.patch.object(eos, 'predict_observables') as task:
            status, _, _ = run('predict-observables', POSTERIOR, 'dBRdq2', '-B', '3', '-E', '7',
                               base='/base')

        self.assertEqual(status, 0)
        task.assert_called_once_with(
            analysis_file=ANALYSIS_FILE, posterior=POSTERIOR, prediction='dBRdq2',
            base_directory='/base', begin=3, end=7, mask_name=None,
        )


if __name__ == '__main__':
    unittest.main(verbosity=5)
