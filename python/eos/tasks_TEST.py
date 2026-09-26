import os
# The report task draws a corner figure; select a head-less backend before eos (and thus
# matplotlib) is imported.
os.environ.setdefault('MPLBACKEND', 'Agg')

import numpy as np
import shutil
import sys
import tempfile
import types
import unittest
import unittest.mock
import eos
import yaml
from pathlib import Path

class ClassMethodTests(unittest.TestCase):

    _analysis_file = str(Path(__file__).parent / "analysis_file_TEST.d/valid-analysis-file.yaml")

    def test_validate_task(self):
        eos.tasks.validate(self._analysis_file)

    def test_list_steps_task(self):
        steps = eos.tasks.list_steps(self._analysis_file)
        self.assertEqual(steps, {'CKM-all,WET-all.sample', 'CKM-all.draw-figure', 'WET-all.mode,draw-figure'})


class ReportTaskTests(unittest.TestCase):

    _fixture = Path(__file__).parent / "reporting_TEST.d"

    def test_report_task(self):
        "Render a report from the reporting_TEST.d fixtures (without generating a PDF)."
        # Copy the fixture into a temporary directory and render there, so that neither the
        # intermediate file, the drawn figures, nor the task's log pollute the source tree.
        workdir = tempfile.mkdtemp(prefix='eos-report-task-')
        self.addCleanup(shutil.rmtree, workdir, ignore_errors=True)
        base = os.path.join(workdir, 'base')
        shutil.copytree(self._fixture, base)

        cwd = os.getcwd()
        self.addCleanup(os.chdir, cwd)
        os.chdir(base)

        eos.tasks.report('analysis.yaml', 'report.md.jinja', base_directory='.', generate_pdf=False)

        with open(os.path.join(base, 'report.md')) as f:
            rendered = f.read()

        # Both posteriors are discovered from the fixture's data/ tree.
        self.assertIn('## Posterior `CKM`', rendered)
        self.assertIn('## Posterior `FF`', rendered)
        # 'CKM' has a recorded mode, hence a goodness-of-fit table with a local p-value and a total row.
        self.assertIn('`B->pilnu::BR`', rendered)
        self.assertIn('**total**', rendered)
        # 'FF' has no mode, hence the goodness-of-fit fallback.
        self.assertIn('No goodness-of-fit information has been recorded', rendered)
        # Parameter names are typeset as math; emitting the raw LaTeX as text breaks the PDF.
        self.assertIn('- `CKM::abs(V_ub)` ($|V_{ub}|$)', rendered)
        self.assertIn('- `B->pi::f_+(0)@BCL2008` ($f_+^{B\\to \\pi}(0)$)', rendered)
        # A corner figure is produced for 'CKM' and referenced from the report.
        self.assertIn('figures/corner-CKM.pdf', rendered)
        self.assertTrue(os.path.isfile(os.path.join(base, 'figures', 'corner-CKM.pdf')))


class _Parameter:
    "Minimal stand-in for eos.Parameter, exposing the name/min/max accessors that ImportanceSamples.create() uses."

    def __init__(self, name, minimum, maximum):
        self._name = name
        self._min  = minimum
        self._max  = maximum

    def name(self):
        return self._name

    def min(self):
        return self._min

    def max(self):
        return self._max


# Three real, registered EOS parameters (reused from other fixtures in this test suite), used here
# purely as labels for arbitrary points in an abstract 3D parameter space -- create_constraint never
# evaluates them against the physical model, so their generated values need not lie within any
# physically meaningful range.
_CREATE_CONSTRAINT_PARAMETER_NAMES = ['CKM::abs(V_ub)', 'B->pi::f_+(0)@BCL2008', 'decay-constant::B_u']

# A correlated (non-diagonal) covariance, shared by every generated fixture below, so that a bug
# that only fit the diagonal would be caught.
_CREATE_CONSTRAINT_COVARIANCE = np.array([
    [1.0, 0.5, 0.2],
    [0.5, 2.0, 0.3],
    [0.2, 0.3, 1.5],
])

_CREATE_CONSTRAINT_ANALYSIS_FILE_YAML = '''
likelihoods:
  - name: dummy-likelihood
    manual_constraints:
      test::create-constraint-dummy@Test:2026A:
        type: Gaussian
        observable: mass::b(MSbar)
        kinematics: {}
        options: {}
        mean: 4.18
        sigma-stat: {hi: 0.03, lo: -0.03}
        sigma-sys: {hi: 0.0, lo: 0.0}

priors:
  - name: uniform-prior
    descriptions:
      - { parameter: 'CKM::abs(V_ub)', min: -10.0, max: 10.0, type: uniform }
      - { parameter: 'B->pi::f_+(0)@BCL2008', min: -10.0, max: 10.0, type: uniform }
      - { parameter: 'decay-constant::B_u', min: -10.0, max: 10.0, type: uniform }

posteriors:
  - name: restricted
    prior: [ uniform-prior ]
    likelihood: [ dummy-likelihood ]
  - name: restricted-bad
    prior: [ uniform-prior ]
    likelihood: [ dummy-likelihood ]

steps:
  - title: 'Create constraint via a dict goodness-of-fit threshold (expected to pass)'
    id: 'restricted.create-constraint'
    tasks:
      - task: 'create-constraint'
        arguments:
          posterior: 'restricted'
          constraint_name: 'test-run-threshold-dict'
          goodness_of_fit_threshold:
            KS-test: 0.03

  - title: 'Create constraint via a dict goodness-of-fit threshold (expected to fail)'
    id: 'restricted-bad.create-constraint'
    tasks:
      - task: 'create-constraint'
        arguments:
          posterior: 'restricted-bad'
          constraint_name: 'test-run-threshold-dict-fail'
          goodness_of_fit_threshold:
            KS-test: 0.99
'''


class CreateConstraintTaskTests(unittest.TestCase):
    "Fixture-based tests of the create-constraint task (issue #1101)."

    def setUp(self):
        workdir = tempfile.mkdtemp(prefix='eos-create-constraint-task-')
        self.addCleanup(shutil.rmtree, workdir, ignore_errors=True)
        self.base = os.path.join(workdir, 'base')
        os.makedirs(self.base, exist_ok=True)
        self.analysis_file = os.path.join(self.base, 'analysis.yaml')
        with open(self.analysis_file, 'w') as f:
            f.write(_CREATE_CONSTRAINT_ANALYSIS_FILE_YAML)

    def _write_samples(self, posterior, samples, weights):
        parameters = [_Parameter(name, -10.0, 10.0) for name in _CREATE_CONSTRAINT_PARAMETER_NAMES]
        path = os.path.join(self.base, 'data', posterior, 'samples')
        eos.data.ImportanceSamples.create(path, parameters, samples, weights)

    @staticmethod
    def _good_fit_samples():
        "~6000 correlated, unequally-weighted samples genuinely drawn from a known Gaussian."
        rng = np.random.default_rng(1701)
        mean = np.array([0.2, -0.1, 0.05])
        samples = rng.multivariate_normal(mean, _CREATE_CONSTRAINT_COVARIANCE, size=6000)
        weights = rng.exponential(1.0, size=6000)
        return samples, weights, mean, _CREATE_CONSTRAINT_COVARIANCE

    @staticmethod
    def _bad_fit_samples():
        "A well-separated two-component mixture -- deliberately not a single Gaussian."
        rng = np.random.default_rng(2026)
        n = 3000
        half = n // 2
        component_a = rng.multivariate_normal([0.0, 0.0, 0.0], _CREATE_CONSTRAINT_COVARIANCE, size=half)
        component_b = rng.multivariate_normal([8.0, 0.0, 0.0], _CREATE_CONSTRAINT_COVARIANCE, size=n - half)
        samples = np.vstack([component_a, component_b])
        weights = np.ones(n)
        return samples, weights

    def _constraint_path(self, constraint_name):
        return os.path.join(self.base, 'constraints', constraint_name, 'constraint.yaml')

    def _load_constraint(self, constraint_name):
        with open(self._constraint_path(constraint_name)) as f:
            return yaml.safe_load(f)

    def test_basic_correctness(self):
        "The fitted mean/covariance are close to the known generating values, off-diagonal included."
        samples, weights, mean, covariance = self._good_fit_samples()
        self._write_samples('restricted', samples, weights)

        body = eos.create_constraint(self.analysis_file, 'restricted', 'test-basic-correctness', base_directory=self.base)

        written = self._load_constraint('test-basic-correctness')['test-basic-correctness']
        self.assertEqual(written['type'], 'MultivariateGaussian(Covariance)')
        self.assertEqual(written['observables'], _CREATE_CONSTRAINT_PARAMETER_NAMES)
        self.assertEqual(written['kinematics'], [{}, {}, {}])
        self.assertEqual(written['options'], [{}, {}, {}])
        np.testing.assert_allclose(written['means'], mean, atol=0.15)
        np.testing.assert_allclose(written['covariance'], covariance, atol=0.3)
        self.assertEqual(body, written)

    def test_strict_aborts_and_writes_nothing(self):
        "strict=True on a deliberately non-Gaussian fixture raises and leaves no constraint file behind."
        samples, weights = self._bad_fit_samples()
        self._write_samples('restricted', samples, weights)

        with self.assertRaises(RuntimeError):
            eos.create_constraint(self.analysis_file, 'restricted', 'test-strict-abort', base_directory=self.base, strict=True)

        self.assertFalse(os.path.exists(self._constraint_path('test-strict-abort')))

    def test_non_strict_warns_but_writes(self):
        "strict=False on the same non-Gaussian fixture does not raise, and the constraint is written."
        samples, weights = self._bad_fit_samples()
        self._write_samples('restricted', samples, weights)

        eos.create_constraint(self.analysis_file, 'restricted', 'test-non-strict', base_directory=self.base, strict=False)

        self.assertTrue(os.path.exists(self._constraint_path('test-non-strict')))

    def test_unknown_test_name_raises_before_any_fit(self):
        "An unrecognized entry in `tests` raises ValueError, naming the offending test, before any fit."
        samples, weights, _, _ = self._good_fit_samples()
        self._write_samples('restricted', samples, weights)

        with self.assertRaises(ValueError) as ctx:
            eos.create_constraint(self.analysis_file, 'restricted', 'test-unknown-test', base_directory=self.base,
                                   tests=['not-a-real-test'])
        self.assertIn('not-a-real-test', str(ctx.exception))
        self.assertFalse(os.path.exists(self._constraint_path('test-unknown-test')))

    def test_empty_tests_list_skips_checking(self):
        "tests=[] runs no goodness-of-fit check at all and writes the constraint unconditionally."
        samples, weights = self._bad_fit_samples()
        self._write_samples('restricted', samples, weights)

        eos.create_constraint(self.analysis_file, 'restricted', 'test-empty-tests', base_directory=self.base, tests=[])

        self.assertTrue(os.path.exists(self._constraint_path('test-empty-tests')))

    def test_parameters_subsetting(self):
        "Passing `parameters` restricts the fit to a single column."
        samples, weights, mean, covariance = self._good_fit_samples()
        self._write_samples('restricted', samples, weights)

        eos.create_constraint(self.analysis_file, 'restricted', 'test-subset', base_directory=self.base,
                               parameters=['B->pi::f_+(0)@BCL2008'])

        written = self._load_constraint('test-subset')['test-subset']
        self.assertEqual(written['observables'], ['B->pi::f_+(0)@BCL2008'])
        self.assertEqual(len(written['means']), 1)
        self.assertEqual(np.array(written['covariance']).shape, (1, 1))
        np.testing.assert_allclose(written['means'], [mean[1]], atol=0.15)
        np.testing.assert_allclose(written['covariance'], [[covariance[1, 1]]], atol=0.3)

    def test_pred_source(self):
        "source='pred-<name>' fits a Prediction fixture, carrying real kinematics/options through."
        params = eos.Parameters.Defaults()
        q2_values = (1.0, 2.0, 3.0)
        observables = [
            eos.Observable.make(
                'B->pilnu::dBR/dq2', params,
                eos.Kinematics({'q2': q2}),
                eos.Options({'P': 'pi', 'form-factors': 'BCL2008', 'model': 'CKM'}),
            )
            for q2 in q2_values
        ]

        samples, weights, mean, covariance = self._good_fit_samples()
        path = os.path.join(self.base, 'data', 'restricted', 'pred-fit-test')
        eos.data.Prediction.create(path, observables, samples, weights)

        eos.create_constraint(self.analysis_file, 'restricted', 'test-pred-source', base_directory=self.base,
                               source='pred-fit-test')

        written = self._load_constraint('test-pred-source')['test-pred-source']
        self.assertEqual(written['observables'], ['B->pilnu::dBR/dq2'] * 3)
        self.assertEqual(written['kinematics'], [{'q2': q2} for q2 in q2_values])
        for options in written['options']:
            self.assertEqual(options, {'P': 'pi', 'form-factors': 'BCL2008', 'model': 'CKM'})
        np.testing.assert_allclose(written['means'], mean, atol=0.15)
        np.testing.assert_allclose(written['covariance'], covariance, atol=0.3)

    def test_unsupported_source_raises_before_any_file_access(self):
        "A source that is neither 'samples' nor 'pred-'-prefixed raises ValueError; no data directory is needed."
        with self.assertRaises(ValueError):
            eos.create_constraint(self.analysis_file, 'restricted', 'test-bad-source', base_directory=self.base,
                                   source='mcmc-0001')
        self.assertFalse(os.path.exists(self._constraint_path('test-bad-source')))

    def test_ks_test_is_directly_callable_and_diagnostic(self):
        "eos.ks_test can be called directly on reloaded samples to diagnose a failed fit."
        samples, weights = self._bad_fit_samples()
        mean = np.average(samples, axis=0, weights=weights)
        covariance = np.atleast_2d(np.cov(samples.T, aweights=weights))

        self.assertIs(eos.goodness_of_fit_tests['KS-test'], eos.ks_test)
        result = eos.ks_test(samples, weights, mean, covariance)

        self.assertLess(result.pvalue, 0.03)
        self.assertGreater(result.statistic, 0.0)
        self.assertEqual(result.squared_distances.shape, (len(samples),))
        self.assertTrue(np.all(result.squared_distances >= 0.0))

    def test_dict_threshold_python_api(self):
        "A dict-valued goodness_of_fit_threshold behaves like the scalar form, and must cover every test."
        samples, weights, _, _ = self._good_fit_samples()
        self._write_samples('restricted', samples, weights)

        eos.create_constraint(self.analysis_file, 'restricted', 'test-dict-threshold', base_directory=self.base,
                               tests=['KS-test'], goodness_of_fit_threshold={'KS-test': 0.03})
        self.assertTrue(os.path.exists(self._constraint_path('test-dict-threshold')))

        with self.assertRaises(ValueError):
            eos.create_constraint(self.analysis_file, 'restricted', 'test-dict-threshold-missing', base_directory=self.base,
                                   tests=['KS-test'], goodness_of_fit_threshold={})
        self.assertFalse(os.path.exists(self._constraint_path('test-dict-threshold-missing')))

    def test_run_task_drives_dict_threshold_from_analysis_file(self):
        "A goodness_of_fit_threshold dict declared in an analysis file's steps: reaches create_constraint via run()."
        good_samples, good_weights, _, _ = self._good_fit_samples()
        self._write_samples('restricted', good_samples, good_weights)
        bad_samples, bad_weights = self._bad_fit_samples()
        self._write_samples('restricted-bad', bad_samples, bad_weights)

        eos.tasks.run(self.analysis_file, 'restricted.create-constraint', base_directory=self.base)
        self.assertTrue(os.path.exists(self._constraint_path('test-run-threshold-dict')))

        with self.assertRaises(RuntimeError):
            eos.tasks.run(self.analysis_file, 'restricted-bad.create-constraint', base_directory=self.base)
        self.assertFalse(os.path.exists(self._constraint_path('test-run-threshold-dict-fail')))


class _SuppressingOutput:
    "Stands in for ipywidgets.Output, which displays and then suppresses what is raised within it."

    def __init__(self, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, etype, evalue, tb):
        return True


class _Accordion:
    "Stands in for ipywidgets.Accordion."

    def __init__(self, children=None):
        self.selected_index = 0

    def set_title(self, index, title):
        pass


class TaskFailureTests(unittest.TestCase):

    def test_failure_escapes_the_output_widget(self):
        "A task that raises under IPython fails its caller, although the output widget suppresses."
        widgets         = types.ModuleType('ipywidgets')
        widgets.Output  = _SuppressingOutput
        widgets.Accordion = _Accordion

        ipython         = types.ModuleType('IPython')
        ipython.display = types.ModuleType('IPython.display')
        ipython.display.display = lambda *args, **kwargs: None

        modules = { 'ipywidgets': widgets, 'IPython': ipython, 'IPython.display': ipython.display }

        @eos.tasks.task('test-failure', '', logfile=False)
        def failing_task():
            raise RuntimeError('the task failed')

        self.addCleanup(eos.tasks._tasks.pop, 'test-failure', None)
        self.addCleanup(eos.tasks._task_outputs.pop, 'test-failure', None)

        with unittest.mock.patch.dict(sys.modules, modules), \
             unittest.mock.patch.object(eos.tasks, '__ipython__', True):
            with self.assertRaises(RuntimeError):
                failing_task()


_MODEL_COMPARISON_ANALYSIS = '''
likelihoods:
  - name: EXP
    constraints: [ 'B^+->tau^+nu::BR@Belle:2014A' ]
  - name: OTHER
    constraints: [ 'B^0->pi^-l^+nu::BR@HFLAV:2019A;form-factors=BCL2008-4' ]

priors:
  - name: CKM
    descriptions:
      - { parameter: 'CKM::abs(V_ub)', min: 3.0e-3, max: 4.5e-3, type: uniform }
  - name: CKM-wide
    descriptions:
      - { parameter: 'CKM::abs(V_ub)', min: 3.0e-3, max: 5e-3, type: uniform }
  - name: DC
    descriptions:
      - { parameter: 'decay-constant::B_u', central: 0.1894, sigma: 0.0014, type: gaussian }
  - name: DC-flat
    descriptions:
      - { parameter: 'decay-constant::B_u', min: 0.18, max: 0.20, type: uniform }

posteriors:
  - { name: A,     prior: [ CKM, DC ],      likelihood: [ EXP ] }
  - { name: B,     prior: [ CKM-wide, DC ], likelihood: [ EXP ] }
  - { name: C,     prior: [ CKM, DC ],      likelihood: [ EXP ] }
  - { name: FLAT,  prior: [ CKM, DC-flat ], likelihood: [ EXP ] }
  - { name: OTHER, prior: [ CKM, DC ],      likelihood: [ EXP, OTHER ] }
  - { name: TWICE, prior: [ CKM, DC ],      likelihood: [ EXP, EXP ] }
'''


class ModelComparisonTaskTests(unittest.TestCase):

    def setUp(self):
        self.base = tempfile.mkdtemp(prefix='eos-model-comparison-task-')
        self.addCleanup(shutil.rmtree, self.base, ignore_errors=True)
        self.analysis_file = os.path.join(self.base, 'analysis.yaml')
        with open(self.analysis_file, 'w') as f:
            f.write(_MODEL_COMPARISON_ANALYSIS)

    def _write_nested(self, posterior, logz, logzerr, vub_max=4.2e-3):
        "Record synthetic nested-sampling results with the given final evidence estimate."
        import dynesty
        N = 100
        rng = np.random.default_rng(1701)
        samples = np.stack([rng.uniform(3.2e-3, vub_max, N), rng.normal(0.1894, 0.0014, N)], axis=1)
        results = dynesty.results.Results(dict(
            samples=samples, samples_u=np.zeros((N, 2)), samples_it=np.arange(N), samples_id=np.arange(N),
            logwt=np.full(N, logz - np.log(N)), logl=np.zeros(N), logvol=np.zeros(N), information=np.zeros(N),
            logz=np.full(N, logz), logzerr=np.full(N, logzerr), ncall=np.ones(N, dtype=int), nlive=10, niter=N, eff=1.0,
        ))
        analysis = eos.AnalysisFile(self.analysis_file).analysis(posterior)
        eos.data.DynestyResults.create(os.path.join(self.base, 'data', posterior, 'nested'), analysis.varied_parameters, results)

    def _run(self, posteriors, **kwargs):
        return eos.tasks.model_comparison(self.analysis_file, posteriors, base_directory=self.base, **kwargs)

    def test_adjustment_and_bayes_factors(self):
        "Differing uniform ranges of a shared parameter are corrected to their common range."
        delta = np.log(2.0 / 1.5)
        self._write_nested('A', -1.0, 0.1)
        self._write_nested('B', -1.5 - delta, 0.1)
        self._write_nested('C', -5.3, 0.1)
        mc = self._run(['A', 'B', 'C'], group='grp')

        entries = { e['posterior']: e for e in mc.posteriors }
        self.assertEqual(mc.group, 'grp')
        self.assertEqual(mc.reference, 'A')
        self.assertAlmostEqual(entries['A']['log_prior_volume_adjustment'], 0.0)
        self.assertAlmostEqual(entries['B']['log_prior_volume_adjustment'], delta)
        self.assertAlmostEqual(entries['B']['adjusted_log_evidence'], -1.5)
        self.assertAlmostEqual(entries['C']['log_bayes_factor'], -4.3)
        self.assertAlmostEqual(entries['C']['log_bayes_factor_uncertainty'], np.hypot(0.1, 0.1))
        self.assertEqual(entries['A']['strength'], 'reference')
        self.assertEqual(entries['B']['strength'], 'barely worth mentioning')
        self.assertEqual(entries['C']['strength'], 'very strong')
        self.assertEqual(len(mc.comparisons), 3)
        self.assertEqual([c['name'] for c in mc.checks], ['uncertainty', 'prior-volume-adjustment'])
        self.assertTrue(mc.stable)

        from eos.reporting import AnalysisData
        ad = AnalysisData(base_directory=self.base)
        self.assertEqual(list(ad.model_comparisons), ['grp'])
        self.assertEqual([c['status'] for c in ad.model_comparisons['grp'].checks], ['passed', 'passed'])
        self.assertNotIn('model-comparison', ad)

    def test_uncertainty_check(self):
        "A pair whose strength changes within one standard deviation fails the uncertainty check."
        self._write_nested('A', -1.0, 0.1)
        self._write_nested('C', -2.2, 0.5)
        mc = self._run(['A', 'C'])
        checks = { c['name']: c for c in mc.checks }
        self.assertEqual(checks['uncertainty']['status'], 'failed')
        self.assertEqual(checks['uncertainty']['failed_pairs'], [['A', 'C']])
        self.assertEqual(checks['prior-volume-adjustment']['status'], 'passed')
        self.assertEqual(mc.comparisons[0]['failed_checks'], ['uncertainty'])

    def test_prior_volume_adjustment_check(self):
        "A pair whose strength changes without the prior-volume adjustment fails that check."
        self._write_nested('A', -1.0, 0.01)
        self._write_nested('B', -2.3, 0.01)
        mc = self._run(['A', 'B'])
        checks = { c['name']: c for c in mc.checks }
        self.assertEqual(checks['uncertainty']['status'], 'passed')
        self.assertEqual(checks['prior-volume-adjustment']['status'], 'failed')
        self.assertEqual(mc.group, 'default')
        self.assertEqual(mc.comparisons[0]['strength'], 'barely worth mentioning')
        self.assertAlmostEqual(mc.comparisons[0]['unadjusted_log_bayes_factor'], 1.3)

    def test_errors(self):
        "Invalid groups of posteriors are rejected."
        for p in ['A', 'FLAT', 'OTHER']:
            self._write_nested(p, -1.0, 0.1)
        self._write_nested('B', -1.0, 0.1, vub_max=4.8e-3)

        with self.assertRaisesRegex(ValueError, 'at least two distinct'):
            self._run(['A', 'A'])
        with self.assertRaisesRegex(ValueError, 'same likelihood'):
            self._run(['A', 'OTHER'])
        with self.assertRaisesRegex(ValueError, 'same likelihood'):
            self._run(['A', 'TWICE'])
        with self.assertRaisesRegex(ValueError, 'not all uniform'):
            self._run(['A', 'FLAT'])
        with self.assertRaisesRegex(ValueError, 'outside the common range'):
            self._run(['A', 'B'])
        with self.assertRaisesRegex(RuntimeError, 'sample-nested'):
            self._run(['A', 'C'])

    def test_report_table(self):
        "The example report template places the model comparison ahead of the posteriors."
        template = Path(__file__).parents[2] / 'examples' / 'inference.md.jinja'
        if not template.is_file():
            self.skipTest('example report template not available')
        self._write_nested('A', -1.0, 0.1)
        self._write_nested('C', -4.5, 0.1)
        self._run(['A', 'C'])

        cwd = os.getcwd()
        self.addCleanup(os.chdir, cwd)
        os.chdir(self.base)
        shutil.copy(template, 'inference.md.jinja')
        eos.tasks.report('analysis.yaml', 'inference.md.jinja', base_directory='.', generate_pdf=False)
        with open('inference.md') as f:
            rendered = f.read()

        self.assertNotIn('## Group', rendered)
        self.assertIn('| `C` | $-4.50 \\pm 0.10$ | $-3.50 \\pm 0.14$ | very strong |', rendered)
        self.assertLess(rendered.index('# Model comparison'), rendered.index('# Posteriors'))


if __name__ == '__main__':
    unittest.main(verbosity=5)
