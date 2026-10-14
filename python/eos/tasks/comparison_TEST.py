import os
# The report task draws a corner figure; select a head-less backend before eos (and thus
# matplotlib) is imported.
os.environ.setdefault('MPLBACKEND', 'Agg')

import numpy as np
import shutil
import tempfile
import unittest
import eos
from pathlib import Path


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
  - name: CKM-low
    descriptions:
      - { parameter: 'CKM::abs(V_ub)', min: 3.0e-3, max: 4.0e-3, type: uniform }
  - name: CKM-high
    descriptions:
      - { parameter: 'CKM::abs(V_ub)', min: 3.9e-3, max: 5.9e-3, type: uniform }
  - name: CKM-narrow
    descriptions:
      - { parameter: 'CKM::abs(V_ub)', min: 1.0e-3, max: 9.0e-3, type: uniform }
  - name: CKM-broad
    descriptions:
      - { parameter: 'CKM::abs(V_ub)', min: 1.0e-3, max: 41.0e-3, type: uniform }
  - name: DC
    descriptions:
      - { parameter: 'decay-constant::B_u', central: 0.1894, sigma: 0.0014, type: gaussian }
  - name: DC-flat
    descriptions:
      - { parameter: 'decay-constant::B_u', min: 0.18, max: 0.20, type: uniform }
  - name: TR
    descriptions:
      - { parameters: [ 'CKM::abs(V_ub)', 'decay-constant::B_u' ], shift: [ 3.0e-3, 0.18 ], transform: [ [ 1.0e-3, 0.0 ], [ 0.0, 0.01 ] ],
          min: [ 0.0, 0.0 ], max: [ 1.5, 2.0 ], type: transform }
  - name: TR-wide
    descriptions:
      - { parameters: [ 'CKM::abs(V_ub)', 'decay-constant::B_u' ], shift: [ 3.0e-3, 0.18 ], transform: [ [ 1.0e-3, 0.0 ], [ 0.0, 0.01 ] ],
          min: [ 0.0, 0.0 ], max: [ 1.5, 4.0 ], type: transform }
  - name: TR-skew
    descriptions:
      - { parameters: [ 'CKM::abs(V_ub)', 'decay-constant::B_u' ], shift: [ 3.0e-3, 0.18 ], transform: [ [ 1.0e-3, 0.0 ], [ 0.01, 0.01 ] ],
          min: [ 0.0, 0.0 ], max: [ 1.0, 1.0 ], type: transform }
  - name: TR-extra
    descriptions:
      - { parameters: [ 'CKM::abs(V_ub)', 'decay-constant::B_u', 'CKM::abs(V_cb)' ], shift: [ 3.0e-3, 0.18, 0.04 ],
          transform: [ [ 1.0e-3, 0.0, 0.0 ], [ 0.0, 0.01, 0.0 ], [ 0.0, 0.0, 0.001 ] ], min: [ 0.0, 0.0, 0.0 ], max: [ 1.5, 2.0, 1.0 ], type: transform }

posteriors:
  - { name: A,     prior: [ CKM, DC ],      likelihood: [ EXP ] }
  - { name: B,     prior: [ CKM-wide, DC ], likelihood: [ EXP ] }
  - { name: C,     prior: [ CKM, DC ],      likelihood: [ EXP ] }
  - { name: FLAT,  prior: [ CKM, DC-flat ], likelihood: [ EXP ] }
  - { name: OTHER, prior: [ CKM, DC ],      likelihood: [ EXP, OTHER ] }
  - { name: TWICE, prior: [ CKM, DC ],      likelihood: [ EXP, EXP ] }
  - { name: LOW,   prior: [ CKM-low, DC ],  likelihood: [ EXP ] }
  - { name: HIGH,  prior: [ CKM-high, DC ], likelihood: [ EXP ] }
  - { name: UNI,   prior: [ CKM, DC-flat ], likelihood: [ EXP ] }
  - { name: TR,    prior: [ TR ],          likelihood: [ EXP ] }
  - { name: TR-wide,  prior: [ TR-wide ],  likelihood: [ EXP ] }
  - { name: TR-skew,  prior: [ TR-skew ],  likelihood: [ EXP ] }
  - { name: TR-extra, prior: [ TR-extra ], likelihood: [ EXP ] }
  - { name: NARROW, prior: [ CKM-narrow ], likelihood: [ EXP ], global_options: { model: CKM } }
  - { name: WIDE,   prior: [ CKM-broad ],  likelihood: [ EXP ], global_options: { model: CKM } }
'''


class ModelComparisonTaskTests(unittest.TestCase):

    def setUp(self):
        self.base = tempfile.mkdtemp(prefix='eos-model-comparison-task-')
        self.addCleanup(shutil.rmtree, self.base, ignore_errors=True)
        self.analysis_file = os.path.join(self.base, 'analysis.yaml')
        with open(self.analysis_file, 'w') as f:
            f.write(_MODEL_COMPARISON_ANALYSIS)

    def _write_nested(self, posterior, logz, logzerr, shift=0.0):
        "Record synthetic nested-sampling results with the given final evidence estimate."
        import dynesty
        N = 100
        rng = np.random.default_rng(1701)
        analysis = eos.AnalysisFile(self.analysis_file).analysis(posterior)
        columns = {
            'CKM::abs(V_ub)':      lambda: rng.uniform(3.2e-3, 3.8e-3, N) + shift,
            'decay-constant::B_u': lambda: rng.normal(0.1894, 0.0014, N),
        }
        samples = np.stack([columns[p.name()]() if p.name() in columns else np.full(N, p.evaluate()) for p in analysis.varied_parameters], axis=1)
        D = samples.shape[1]
        results = dynesty.results.Results(dict(
            samples=samples, samples_u=np.zeros((N, D)), samples_it=np.arange(N), samples_id=np.arange(N),
            logwt=np.full(N, logz - np.log(N)), logl=np.zeros(N), logvol=np.zeros(N), information=np.zeros(N),
            logz=np.full(N, logz), logzerr=np.full(N, logzerr), ncall=np.ones(N, dtype=int), nlive=10, niter=N, eff=1.0,
        ))
        eos.data.DynestyResults.create(os.path.join(self.base, 'data', posterior, 'nested'), analysis.varied_parameters, results)

    def _run(self, posteriors, **kwargs):
        return eos.tasks.model_comparison(self.analysis_file, posteriors, base_directory=self.base, **kwargs)

    def test_adjustment_and_bayes_factors(self):
        "Differing uniform ranges of a shared parameter are corrected to the span of all ranges."
        delta = np.log(2.0 / 1.5)
        self._write_nested('A', -1.0, 0.1)
        self._write_nested('B', -1.5 - delta, 0.1)
        self._write_nested('C', -5.3, 0.1)
        mc = self._run(['A', 'B', 'C'], group='grp')

        entries = { e['posterior']: e for e in mc.posteriors }
        self.assertEqual(mc.group, 'grp')
        self.assertEqual(mc.reference, 'A')
        self.assertAlmostEqual(entries['A']['log_prior_volume_adjustment'], -delta)
        self.assertAlmostEqual(entries['B']['log_prior_volume_adjustment'], 0.0)
        self.assertAlmostEqual(entries['B']['adjusted_log_evidence'], -1.5 - delta)
        self.assertAlmostEqual(entries['C']['log_bayes_factor'], -4.3)
        self.assertAlmostEqual(entries['C']['log_bayes_factor_uncertainty'], np.hypot(0.1, 0.1))
        self.assertEqual(entries['A']['strength'], 'reference')
        self.assertEqual(entries['B']['strength'], 'barely worth mentioning')
        self.assertEqual(entries['C']['strength'], 'very strong')
        self.assertEqual(len(mc.comparisons), 3)
        self.assertEqual([c['name'] for c in mc.checks], ['uncertainty', 'prior-volume-adjustment'])
        self.assertEqual(len(mc.log_prior_volumes), 1)
        self.assertEqual(mc.log_prior_volumes[0]['parameters'], ['CKM::abs(V_ub)'])
        self.assertAlmostEqual(mc.log_prior_volumes[0]['log_reference_volume'], np.log(2.0e-3))
        self.assertAlmostEqual(mc.log_prior_volumes[0]['log_volumes']['A'], np.log(1.5e-3))
        self.assertTrue(mc.stable)

        # left behind by a killed task
        os.makedirs(os.path.join(self.base, 'model-comparison', '.grp.0123.tmp'))

        from eos.reporting import AnalysisData
        ad = AnalysisData(base_directory=self.base)
        self.assertEqual(list(ad.model_comparisons), ['grp'])
        self.assertEqual([c['status'] for c in ad.model_comparisons['grp'].checks], ['passed', 'passed'])
        self.assertNotIn('model-comparison', ad)

    def test_barely_overlapping_ranges(self):
        "Identical but shifted posteriors whose prior ranges barely overlap yield a vanishing Bayes factor."
        self._write_nested('LOW',  -1.0,             0.1)
        self._write_nested('HIGH', -1.0 - np.log(2), 0.1, shift=1.5e-3)
        mc = self._run(['LOW', 'HIGH'])
        entries = { e['posterior']: e for e in mc.posteriors }
        self.assertAlmostEqual(entries['LOW']['log_prior_volume_adjustment'],  np.log(1.0 / 2.9))
        self.assertAlmostEqual(entries['HIGH']['log_prior_volume_adjustment'], np.log(2.0 / 2.9))
        self.assertAlmostEqual(mc.comparisons[0]['log_bayes_factor'], 0.0)
        self.assertAlmostEqual(abs(mc.comparisons[0]['unadjusted_log_bayes_factor']), np.log(2))

    def test_transform_priors(self):
        "Transform priors are corrected by their volume within the box that contains all supports."
        cases = [
            # posteriors, expected adjustments, expected log reference volume
            (['UNI', 'TR'],      [0.0, 0.0],                 np.log(3.0e-5)),
            (['TR', 'TR-wide'],  [np.log(0.5), 0.0],         np.log(6.0e-5)),
            (['UNI', 'TR-skew'], [0.0, np.log(1.0 / 3.0)],   np.log(3.0e-5)),
        ]
        for p in ['UNI', 'TR', 'TR-wide', 'TR-skew']:
            self._write_nested(p, -1.0, 0.1)
        for posteriors, adjustments, log_reference_volume in cases:
            with self.subTest(posteriors=posteriors):
                mc = self._run(posteriors)
                entries = { e['posterior']: e for e in mc.posteriors }
                for p, adjustment in zip(posteriors, adjustments):
                    self.assertAlmostEqual(entries[p]['log_prior_volume_adjustment'], adjustment)
                self.assertEqual(len(mc.log_prior_volumes), 1)
                self.assertEqual(mc.log_prior_volumes[0]['parameters'], ['CKM::abs(V_ub)', 'decay-constant::B_u'])
                self.assertAlmostEqual(mc.log_prior_volumes[0]['log_reference_volume'], log_reference_volume)

    def test_prior_volume_with_nested_sampling(self):
        "Posteriors that differ only by the prior volume of a shared parameter are equally good after the adjustment."
        for p in ['NARROW', 'WIDE']:
            eos.tasks.sample_nested(self.analysis_file, p, base_directory=self.base, nlive=100, dlogz=0.5, seed=1701)
        mc = self._run(['NARROW', 'WIDE'])
        pair = mc.comparisons[0]
        sigma = pair['log_bayes_factor_uncertainty']
        self.assertAlmostEqual(abs(pair['unadjusted_log_bayes_factor']), np.log(5.0), delta=3 * sigma)
        self.assertEqual(eos.tasks._bayes_factor_strength(pair['unadjusted_log_bayes_factor']), 'substantial')
        self.assertLess(pair['log_bayes_factor'], 3 * sigma)
        self.assertEqual(pair['strength'], 'barely worth mentioning')
        self.assertEqual(pair['failed_checks'], ['prior-volume-adjustment'])

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
        for p in ['A', 'FLAT', 'OTHER', 'UNI', 'TR-extra']:
            self._write_nested(p, -1.0, 0.1)

        with self.assertRaisesRegex(ValueError, 'at least two distinct'):
            self._run(['A', 'A'])
        with self.assertRaisesRegex(ValueError, 'same likelihood'):
            self._run(['A', 'OTHER'])
        with self.assertRaisesRegex(ValueError, 'same likelihood'):
            self._run(['A', 'TWICE'])
        with self.assertRaisesRegex(ValueError, 'not all uniform'):
            self._run(['A', 'FLAT'])
        with self.assertRaisesRegex(ValueError, 'model-specific'):
            self._run(['UNI', 'TR-extra'])
        with self.assertRaisesRegex(RuntimeError, 'sample-nested'):
            self._run(['A', 'C'])
        for group in ['a/b', '..', 'a b']:
            with self.assertRaisesRegex(ValueError, 'single path component'):
                self._run(['A', 'UNI'], group=group)

        # failed comparisons leave no output behind
        directory = os.path.join(self.base, 'model-comparison')
        self.assertEqual([e for e in os.listdir(directory) if not e.startswith('.')], [])

    def test_report_table(self):
        "The example report template places the model comparison ahead of the posteriors."
        template = Path(__file__).parents[3] / 'examples' / 'inference.md.jinja'
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
