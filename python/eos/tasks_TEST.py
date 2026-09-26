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



class TaskOutputTests(unittest.TestCase):
    "Tests of the output-directory handling of the task decorator (issue #1291)."

    def setUp(self):
        self.base = tempfile.mkdtemp(prefix='eos-task-output-')
        self.addCleanup(shutil.rmtree, self.base, ignore_errors=True)

        @eos.tasks.task('test-output', 'data/{posterior}/out-{label}', load_analysis_file=False)
        def output_task(posterior, label, base_directory, content, fail=False):
            path = os.path.join(base_directory, 'data', posterior, f'out-{label}')
            with open(os.path.join(path, 'result'), 'w') as f:
                f.write(content)
            if fail:
                raise RuntimeError('the task failed')

        self.addCleanup(eos.tasks._tasks.pop, 'test-output', None)
        self.addCleanup(eos.tasks._task_outputs.pop, 'test-output', None)
        self.task = output_task

    def test_invalid_names_are_rejected(self):
        for label in ('foo/bar', '../x', 'a b', '.', '..', ''):
            with self.subTest(label=label), self.assertRaises(ValueError):
                self.task(posterior='P', label=label, base_directory=self.base, content='x')
        for posterior in ('.', '..'):
            with self.subTest(posterior=posterior), self.assertRaises(ValueError):
                self.task(posterior=posterior, label='l', base_directory=self.base, content='x')
        self.assertEqual([], os.listdir(self.base))

if __name__ == '__main__':
    unittest.main(verbosity=5)
