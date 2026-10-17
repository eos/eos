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

import unittest

import eos
import eos.data
from eos.diagnostic import Severity
import numpy as np
import os
import shutil
import tempfile
import yaml


def nabu_is_missing():
    try:
        import nabu
        return False
    except ImportError:
        return True


def _model_path():
    return os.path.join(os.environ['SOURCE_DIR'], 'eos/data/nabu_likelihood_TEST.d/nabu')


# These tests exercise the hardening of the description parsing, which happens before the (optional)
# nabu payload is loaded, so they do not require the nabu module to be installed.
class NabuLikelihoodTests(unittest.TestCase):

    @staticmethod
    def _write(directory, description):
        "Materialize a (possibly malformed) NabuLikelihood description in the given directory."
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, 'description.yaml'), 'w') as f:
            yaml.safe_dump(description, f, default_flow_style=False)

    def test_invalid_path(self):
        "Loading from a non-existent path raises an error."
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(RuntimeError):
                eos.data.NabuLikelihood(os.path.join(d, 'does-not-exist'))

    def test_wrong_type(self):
        "A description whose type is not 'NabuLikelihood' is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'MarkovChain', 'parameters': []})
            with self.assertRaises(ValueError):
                eos.data.NabuLikelihood(d)

    def test_unknown_key(self):
        "An unexpected key in the description is rejected rather than silently ignored."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'NabuLikelihood',
                'parameters': [{'name': 'a', 'min': 0.0, 'max': 1.0}], 'bogus': 42,
            })
            with self.assertRaises(ValueError):
                eos.data.NabuLikelihood(d)

    def test_missing_parameters(self):
        "A description missing the required 'parameters' key is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'NabuLikelihood'})
            with self.assertRaises(ValueError):
                eos.data.NabuLikelihood(d)

    def test_missing_payload(self):
        "A valid description without the likelihood payload is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'NabuLikelihood',
                'parameters': [{'name': 'a', 'min': 0.0, 'max': 1.0}],
            })
            with self.assertRaises(RuntimeError):
                eos.data.NabuLikelihood(d)

    def test_non_mapping_description(self):
        "A description.yaml that is not a top-level mapping is rejected."
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'description.yaml'), 'w') as f:
                yaml.safe_dump(['not', 'a', 'mapping'], f)
            with self.assertRaises(RuntimeError):
                eos.data.NabuLikelihood(d)


@unittest.skipIf(nabu_is_missing(), "Test requires the 'nabu' module")
class NabuLogPriorTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.affinity = os.sched_getaffinity(0) if hasattr(os, 'sched_getaffinity') else None
        cls.nabu = eos.data.NabuLikelihood(_model_path())
        cls.lh = cls.nabu.likelihood

    def _prior(self):
        parameters = eos.Parameters()
        prior = self.nabu.log_prior(parameters)
        return parameters, prior, list(prior.varied_parameters())

    def test_affinity_is_restored(self):
        "Loading the likelihood leaves the CPU affinity of every thread unchanged."
        if self.affinity is None or not os.path.isdir('/proc/self/task'):
            self.skipTest('platform lacks sched_getaffinity or /proc/self/task')
        for tid in os.listdir('/proc/self/task'):
            self.assertEqual(os.sched_getaffinity(int(tid)), self.affinity, f'thread {tid}')

    def test_varied_parameters(self):
        "The prior varies the parameters of description.yaml, in their order."
        _, prior, varied = self._prior()
        self.assertEqual([p.name() for p in varied], ['mass::b(MSbar)', 'mass::c'])

    def test_evaluate(self):
        "evaluate() matches nabu's log_prob at the current parameter values."
        _, prior, varied = self._prior()
        for x in self.lh.sample(5, random_seed=7):
            for p, v in zip(varied, x):
                p.set(v)
            self.assertAlmostEqual(prior.evaluate(), float(self.lh.log_prob(x)), places=12)

    def test_round_trips(self):
        "sample() and compute_cdf() invert each other and follow nabu's unit-cube maps."
        _, prior, varied = self._prior()
        for u in np.random.default_rng(3).uniform(size=(5, 2)):
            for p, uv in zip(varied, u):
                p.set_generator(uv)
            prior.sample()
            x = np.array([p.evaluate() for p in varied])
            np.testing.assert_allclose(x, self.lh.from_unit_cube(u), rtol=0, atol=1e-12)
            prior.compute_cdf()
            np.testing.assert_allclose([p.evaluate_generator() for p in varied], u, rtol=0, atol=1e-6)

    def test_log_posterior(self):
        "The prior combines with other priors inside a LogPosterior."
        from math import log, pi, sqrt
        parameters, prior, varied = self._prior()
        posterior = eos.LogPosterior(eos.LogLikelihood(parameters))
        self.assertTrue(posterior.add(prior, False))
        self.assertTrue(posterior.add(eos.LogPrior.Gaussian(parameters, 'CKM::abs(V_cb)', 0.04, 0.002), False))
        x = self.lh.sample(1, random_seed=11)[0]
        for p, v in zip(varied, x):
            p.set(v)
        parameters['CKM::abs(V_cb)'].set(0.04)
        expected = float(self.lh.log_prob(x)) - log(0.002 * sqrt(2.0 * pi))
        self.assertAlmostEqual(posterior.evaluate(), expected, places=10)


@unittest.skipIf(nabu_is_missing(), "Test requires the 'nabu' module")
class NabuPriorAnalysisFileTests(unittest.TestCase):

    _analysis_file = {
        'priors': [
            { 'name': 'flow', 'descriptions': [ { 'type': 'nabu', 'path': 'data/flow/nabu' } ] },
            { 'name': 'Vcb',  'descriptions': [ { 'type': 'uniform', 'parameter': 'CKM::abs(V_cb)', 'min': 0.03, 'max': 0.05 } ] },
            { 'name': 'mc',   'descriptions': [ { 'type': 'uniform', 'parameter': 'mass::c', 'min': 1.0, 'max': 1.5 } ] },
        ],
        'likelihoods': [
            { 'name': 'mc', 'manual_constraints': { 'test::mc': {
                'type': 'Gaussian', 'observable': 'mass::c', 'kinematics': {}, 'options': {},
                'mean': 1.28, 'sigma-stat': { 'hi': 0.01, 'lo': 0.01 }, 'sigma-sys': { 'hi': 0.0, 'lo': 0.0 },
                'dof': 1,
            } } },
            { 'name': 'Vcb', 'manual_constraints': { 'test::Vcb': {
                'type': 'Gaussian', 'observable': 'CKM::abs(V_cb)', 'kinematics': {}, 'options': {},
                'mean': 0.04, 'sigma-stat': { 'hi': 0.001, 'lo': 0.001 }, 'sigma-sys': { 'hi': 0.0, 'lo': 0.0 },
                'dof': 1,
            } } },
        ],
        'posteriors': [
            { 'name': 'prior-only', 'prior': [ 'flow' ], 'likelihood': [] },
            { 'name': 'combined',   'prior': [ 'flow', 'Vcb' ], 'likelihood': [ 'mc', 'Vcb' ] },
            { 'name': 'repeated',   'prior': [ 'flow', 'mc' ], 'likelihood': [] },
            { 'name': 'repeated-reversed', 'prior': [ 'mc', 'flow' ], 'likelihood': [] },
        ],
    }

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.base = self._tmp.name
        shutil.copytree(_model_path(), os.path.join(self.base, 'data', 'flow', 'nabu'))
        self.file = os.path.join(self.base, 'analysis.yaml')
        with open(self.file, 'w') as f:
            yaml.safe_dump(self._analysis_file, f)
        self.af = eos.AnalysisFile(self.file)
        self.lh = eos.data.NabuLikelihood(os.path.join(self.base, 'data', 'flow', 'nabu')).likelihood

    def tearDown(self):
        self._tmp.cleanup()

    def test_combined_posterior(self):
        "The hypercube maps and the log(posterior) include the nabu prior alongside other priors and likelihoods."
        from scipy.stats import norm
        analysis = self.af.analysis('combined', base_directory=self.base)
        self.assertEqual(analysis.varied_parameter_names, ['mass::b(MSbar)', 'mass::c', 'CKM::abs(V_cb)'])
        u = np.array([0.3, 0.6, 0.25])
        x = analysis._u_to_par(u)
        np.testing.assert_allclose(x[:2], self.lh.from_unit_cube(u[:2]), rtol=0, atol=1e-12)
        self.assertAlmostEqual(x[2], 0.035, places=12)
        np.testing.assert_allclose(analysis._par_to_u(x), u, rtol=0, atol=1e-6)
        expected = (float(self.lh.log_prob(x[:2])) - np.log(0.02)
                    + norm.logpdf(x[1], 1.28, 0.01) + norm.logpdf(x[2], 0.04, 0.001))
        self.assertAlmostEqual(analysis.log_pdf(u), expected, places=8)

    def test_repeated_parameter(self):
        "A parameter varied by the nabu prior must not appear in another prior."
        for posterior in ('repeated', 'repeated-reversed'):
            with self.assertRaisesRegex(ValueError, "'mass::c' is repeated"):
                self.af.analysis(posterior, base_directory=self.base)
        diagnostics = self.af.validate(base_directory=self.base)
        self.assertEqual(
            sorted(d.path for d in diagnostics if d.severity == Severity.ERROR),
            [('posteriors', 'repeated'), ('posteriors', 'repeated-reversed')])

    def test_repr_html(self):
        "The HTML representation of an analysis renders the nabu prior by its path."
        analysis = self.af.analysis('combined', base_directory=self.base)
        self.assertIn(f"<tt>{os.path.join(os.path.abspath(self.base), 'data/flow/nabu')}</tt>", analysis._repr_html_())

    def test_start_point_outside_support(self):
        "optimize() reports a start point at which the prior's CDF is undefined."
        analysis = self.af.analysis('prior-only', base_directory=self.base)
        with self.assertLogs('EOS', level='ERROR') as logs:
            analysis.optimize(start_point=[float('nan'), 1.27], options={'maxiter': 1})
        self.assertTrue(any('outside its defined range' in m for m in logs.output))

    def test_missing_path(self):
        "A path that does not hold a NabuLikelihood is reported when the posterior is created."
        with tempfile.TemporaryDirectory() as other:
            with self.assertRaises(RuntimeError):
                self.af.analysis('prior-only', base_directory=other)

    def test_sample_prior(self):
        "The sample-prior task reproduces the statistics of nabu's own samples."
        from scipy.stats import ks_2samp
        eos.tasks.sample_prior(self.af, 'prior-only', base_directory=self.base, N=2000, seed=1701)
        samples = eos.data.ImportanceSamples(os.path.join(self.base, 'data', 'prior-only', 'samples')).samples
        reference = self.lh.sample(2000, random_seed=1701)
        for i in range(2):
            self.assertGreater(ks_2samp(samples[:, i], reference[:, i]).pvalue, 0.01)
        self.assertAlmostEqual(np.corrcoef(samples.T)[0, 1], np.corrcoef(reference.T)[0, 1], delta=0.1)

    def test_sample_nested(self):
        "The sample-nested task runs on a posterior with a nabu prior."
        try:
            import dynesty
        except ImportError:
            self.skipTest("Test requires the 'dynesty' module")
        eos.tasks.sample_nested(self.af, 'combined', base_directory=self.base, nlive=50, dlogz=1.0, maxiter=1500, seed=10)
        samples = eos.data.ImportanceSamples(os.path.join(self.base, 'data', 'combined', 'samples'))
        self.assertEqual(samples.samples.shape[1], 3)
        mean = np.average(samples.samples, weights=samples.weights, axis=0)
        self.assertAlmostEqual(mean[2], 0.04, delta=0.001)


if __name__ == '__main__':
    unittest.main(verbosity=5)
