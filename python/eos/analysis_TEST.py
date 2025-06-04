import unittest

from numpy import random
import dynesty
import eos
import numpy as np
import yaml
import json
from pathlib import Path

class ClassMethodTests(unittest.TestCase):

    def test_priors(self):

        # Note: the cases are not physically meaningful
        prior_cases = [
            { 'parameter': 'B->D::alpha^f+_0@BSZ2015', 'min':  0.0, 'max': 1.0, 'type': 'uniform'},
            { 'parameter': 'B->D::alpha^f+_0@BSZ2015', 'min':  0.0, 'max': 1.0, 'central': 0.0, 'sigma': 1.0 , 'type': 'gaussian'},
            { 'parameter': 'B->D::alpha^f+_0@BSZ2015', 'min':  0.0, 'max': 1.0, 'central': 0.0, 'sigma': (0.3, 0.5) , 'type': 'gaussian'},
            { 'parameter': 'B->D::alpha^f+_0@BSZ2015', 'min':  0.0, 'max': 1.0, 'mu_0': 1.0, 'lambda' :2.0, 'type': 'scale'},
        ]

        for prior in prior_cases:
            analysis_args = {
                'priors': [ prior ],
                'likelihood': [ ]
            }

            try:
                analysis = eos.Analysis(**analysis_args)
            except Exception as e:
                self.fail(f"Can not construct analysis with prior '{ prior }'")


    def test_repeating_priors_parameters(self):

         # Note: the cases are not physically meaningful
         priors = [
             { 'parameter': 'B->D::alpha^f+_0@BSZ2015', 'min':  0.0, 'max': 1.0, 'type': 'uniform'},
             { 'parameter': 'B->D::alpha^f+_0@BSZ2015', 'min':  0.0, 'max': 1.0, 'type': 'uniform'}
         ]

         analysis_args = {
             'priors': priors ,
             'likelihood': [ ]
         }

         errorraised = False
         try:
             analysis = eos.Analysis(**analysis_args)
         except Exception as e:
             errorraised=True
         if errorraised==False:
             self.fail(f"Can construct analysis with prior '{ priors } eventhough prior is repeated'")

    def test_repeating_priors_constraints(self):

         # Note: the cases are not physically meaningful
         priors = [
             { 'constraint': 'D_s^*::decay-constants@PZ:2021A' },
             { 'constraint': 'D_s^*_T::decay-constants@PZ:2021A' }
         ]

         analysis_args = {
             'priors': priors ,
             'likelihood': [ ]
         }

         errorraised = False
         try:
             analysis = eos.Analysis(**analysis_args)
         except Exception as e:
             errorraised=True
         if errorraised==False:
             self.fail(f"Can construct analysis with prior '{ priors } eventhough prior is repeated'")

    def test_repeating_priors(self):

         # Note: the cases are not physically meaningful
         priors = [
             { 'parameter': 'decay-constant::D_s',      'min':  0.24781, 'max':  0.25197, 'type': 'uniform'},
             { 'constraint': 'D_s^*::decay-constants@PZ:2021A' }
         ]

         analysis_args = {
             'priors': priors ,
             'likelihood': [ ]
         }

         errorraised = False
         try:
             analysis = eos.Analysis(**analysis_args)
         except Exception as e:
             errorraised=True
         if errorraised==False:
             self.fail(f"Can construct analysis with prior '{ priors } eventhough prior is repeated'")

    def test_parameter_scaling(self):

        analysis_args = {
            'global_options': { 'form-factors': 'BSZ2015', 'model': 'CKM' },
            'priors': [
                { 'parameter': 'CKM::abs(V_cb)',           'min':  38e-3, 'max':  45e-3 , 'type': 'uniform'},
                { 'parameter': 'B->D::alpha^f+_0@BSZ2015', 'min':  0.0,   'max':  1.0   , 'type': 'uniform'},
                { 'parameter': 'B->D::alpha^f+_1@BSZ2015', 'min': -4.0,   'max': -1.0   , 'type': 'uniform'},
                { 'parameter': 'B->D::alpha^f+_2@BSZ2015', 'min': +4.0,   'max': +6.0   , 'type': 'uniform'},
                { 'parameter': 'B->D::alpha^f0_1@BSZ2015', 'min': -1.0,   'max': +2.0   , 'type': 'uniform'},
                { 'parameter': 'B->D::alpha^f0_2@BSZ2015', 'min': -2.0,   'max':  0.0   , 'type': 'uniform'}
            ],
            'likelihood': [
                'B->D::f_++f_0@HPQCD:2015A',
                'B->D::f_++f_0@FNAL+MILC:2015B',
                'B^0->D^+e^-nu::BRs@Belle:2015A',
                'B^0->D^+mu^-nu::BRs@Belle:2015A'
            ]
        }

        # Test analysis definition
        analysis = eos.Analysis(**analysis_args)

        # Test analysis optimization
        bfp = analysis.optimize()

        # Test analysis optimization from a random point
        bfp = analysis.optimize(start_point='random', rng=np.random.mtrand.RandomState(123))

        # Test parameter scaling
        point = bfp.point

        self.assertTrue(
            np.max(analysis._par_to_u(point)) < 1.0
            )
        self.assertTrue(
            np.min(analysis._par_to_u(point)) >= 0.0
            )
        self.assertTrue(
            np.max(analysis._u_to_par(analysis._par_to_u(point)) - point) < 1e-10
            )


    def test_sanitize_manual_input(self):

        types  = np.array(["Gaussian", "Gaussian"])
        values = np.array([4.123, 3.141])

        input_output_cases = [
            {'dirty': np.double(1.61), 'clean': float(1.61)},
            {'dirty': np.str_('Gaussian'), 'clean': "Gaussian"},
            {
             'dirty': {'type': types[0], 'mean': values[0]}, # Simplified manual contraint with numpy types
             'clean': {'type': 'Gaussian', 'mean': float(4.123)}
            },
            {
             'dirty': {'type': types, 'mean': values}, # Simplified manual contraint with numpy arrays
             'clean': {'type': ['Gaussian', 'Gaussian'], 'mean': [float(4.123), float(3.141)]}
            },
        ]

        for case in input_output_cases:
            # check value and type
            self.assertEqual(
                eos.Analysis._sanitize_manual_input(case['dirty']),
                case['clean']
            )
            self.assertEqual(
                type(eos.Analysis._sanitize_manual_input(case['dirty'])),
                type(case['clean'])
            )
            # check Yaml string
            self.assertNotEqual(
                yaml.dump(case['dirty']),
                yaml.dump(case['clean'])
            )
            self.assertEqual(
                yaml.dump(eos.Analysis._sanitize_manual_input(case['dirty'])),
                yaml.dump(case['clean'])
            )


    def test_univariate_priors_with_infinite_support(self):

        # Note: the cases are not physically meaningful
        prior_cases = [
            { 'parameter': 'B->D::alpha^f+_0@BSZ2015', 'central': 0.0, 'sigma': 1.0 , 'type': 'gaussian'},
            { 'parameter': 'B->D::alpha^f+_0@BSZ2015', 'k':       5.0,                'type': 'poisson'},
            { 'parameters': ['scnuee::Re{cVL}','scnuee::Re{cVR}'], 'shift': [0.0,0.0], 'transform': [[ 0.707106, 0.707106 ], [-0.707106, 0.707106]], 'min': [-2.0,-2.0], 'max': [2.0,2.0], 'type': 'transform'}
        ]

        for prior in prior_cases:
            analysis_args = {
                'priors': [ prior ],
                'likelihood': [ ]
            }

            try:
                analysis = eos.Analysis(**analysis_args)
            except Exception as e:
                self.fail(f"Can not construct analysis with prior '{ prior }'")


    def test_multivariate_priors_1(self):

        import yaml, numpy as np
        analysis_args = {
            'global_options': { },
            'manual_constraints': {
                'test::test': {
                    'type': 'MultivariateGaussian(Covariance)',
                    'observables': ['mass::c', 'mass::b(MSbar)'],
                    'kinematics': [{}, {}],
                    'options': [{}, {}],
                    'means': [1.27, 4.18],
                    'covariance': [[0.03**2, 0.0], [0.0, 0.02**2]],
                }
            },
            'priors': [
                { 'constraint': 'B->K::FormFactors[parametric,LCSR]@GKvD:2018A' },
                { 'parameter': 'mass::c', 'min': 1.18, 'max': 1.36, 'type': 'uniform' },
                { 'parameter': 'mass::b(MSbar)', 'min': 4.12, 'max': 4.24, 'type': 'uniform' },
            ],
            'likelihood': [
                # no entry; ``test::test`` is automatically selected as a manual constraint
            ]
        }

        # Test analysis definition
        analysis = eos.Analysis(**analysis_args)

        # Sample from the prior
        results = analysis.sample_nested(bound='multi', nlive=250, dlogz=0.01, miniter=10000, maxiter=100000, seed=10, print_progress=False)

        # Test samples against constraint
        entry = eos.Constraints()['B->K::FormFactors[parametric,LCSR]@GKvD:2018A']
        entry_data = yaml.load(entry.serialize(), Loader=yaml.SafeLoader)
        means  = np.array(entry_data['means'])
        cov    = np.array(entry_data['covariance'])
        sigmas = np.diag(cov)

        means  = np.append(means,  [1.27, 4.18])
        sigmas = np.append(sigmas, [0.03, 0.02])

        avg    = np.median(results.samples_equal(), axis=0)
        delta  = avg - means
        chi2_1 = np.dot(delta[0:8], np.dot(np.linalg.inv(cov), delta[0:8]))
        chi2_2 = (avg[8] - means[8])**2 / sigmas[8]**2
        chi2_3 = (avg[9] - means[9])**2 / sigmas[9]**2
        self.assertLess(chi2_1, 2.0325e-0, 'chi^2 for BSZ parameters exceeds 2% integrated probability for 8 degrees of freedom')
        self.assertLess(chi2_2, 1.5791e-2, 'chi^2 for mass::c exceeds 10% integrated probability for 1 degree of freedom')
        self.assertLess(chi2_3, 1.5791e-2, 'chi^2 for mass::b(MSbar) exceeds 10% integrated probability for 1 degree of freedom')

        # Result for log(Z) obtained by numerically integrating the posterior in Mathematica: 3.82966
        logz_analytic = 3.82966
        chi2_4 = (results['logz'][-1] - logz_analytic)**2 / results['logzerr'][-1]**2
        self.assertLess(chi2_2, 4.5494e-1, 'chi^2 for log(Z) exceeds 50% integrated probability for 1 degree of freedom')


    def test_multivariate_priors_2(self):

        import yaml, numpy as np
        # inject new constraint
        eos.Constraints().insert('test::test2', '''
        'type': 'MultivariateGaussian(Covariance)'
        'observables': ['mass::c', 'mass::b(MSbar)']
        'kinematics': [{}, {}]
        'options': [{}, {}]
        'means': [1.25, 4.19]
        'covariance': [[0.0009, 0.0003], [0.0003, 0.0004]]
        ''')
        analysis_args = {
            'global_options': { },
            'manual_constraints': {
                'test::test': {
                    'type': 'MultivariateGaussian(Covariance)',
                    'observables': ['mass::c', 'mass::b(MSbar)'],
                    'kinematics': [{}, {}],
                    'options': [{}, {}],
                    'means': [1.28, 4.17],
                    'covariance': [[0.03**2, 0.0], [0.0, 0.02**2]],
                }
            },
            'priors': [
                { 'constraint': 'test::test2' },
            ],
            'likelihood': [
                # no entry; ``test::test`` is automatically selected as a manual constraint
            ]
        }

        # Test analysis definition
        analysis = eos.Analysis(**analysis_args)

        # Sample from the prior
        results = analysis.sample_nested(bound='multi', nlive=250, dlogz=0.01, seed=10, print_progress=False)

        # Test samples against analytic results:
        means  = [1.26000, 4.18333]
        sigmas = [0.02049, 0.01366]
        avg    = np.median(results.samples_equal(), axis=0)

        chi2_0 = (avg[0] - means[0])**2 / sigmas[0]**2
        chi2_1 = (avg[1] - means[1])**2 / sigmas[1]**2
        self.assertLess(chi2_0, 3.9321e-3, 'chi^2 for mass::c exceeds 5% integrated probability for 1 degree of freedom')
        self.assertLess(chi2_1, 3.9921e-3, 'chi^2 for mass::b(MSbar) exceeds 5% integrated probability for 1 degree of freedom')

        # Result for log(Z) obtained by numerically integrating the posterior in Mathematica: 4.2532
        logz_analytic = 4.2532
        chi2_2 = (results['logz'][-1] - logz_analytic)**2 / results['logzerr'][-1]**2 # Assuming 2% error on the log(Z) value
        self.assertLess(chi2_2, 4.5494e-1, 'chi^2 for log(Z) exceeds 50% integrated probability for 1 degree of freedom')

    def test_pyhf_likelihood(self):

        try:
            import pyhf, os, contextlib
        except ModuleNotFoundError:
            raise RuntimeError('Attempting to use a PyHF likelihood without the `pyhf` module installed')

        @contextlib.contextmanager
        def cd(path):
            cwd = os.getcwd()
            os.chdir(path)
            try:
                yield
            finally:
                os.chdir(cwd)

        with cd(Path(__file__).parent / 'analysis_TEST.d'):
            with open('workspace.json') as workspace:
                spec = json.load(workspace)

            workspace = pyhf.Workspace(spec)
            model = workspace.model()
            data = np.array(workspace.data(model, include_auxdata=False))
            init_pars = np.array(model.config.suggested_init())

            af = eos.AnalysisFile('analysis_file_pyhf.yaml')
            analysis = af.analysis('EXP-pyhf')

            # Compare log likelihood values for the main term
            eos_likelihood = analysis.log_likelihood([])
            pyhf_likelihood = model.mainlogpdf(data, init_pars).item()
            self.assertAlmostEqual(eos_likelihood, pyhf_likelihood)

            # Check if the best fit point is close to pyhf fit results
            best_fit_point = analysis.optimize().point
            expected_best_fit_point = [2., 0., 1., 1., 1., 1.]
            eps = 5

            for expected, obtained in zip(expected_best_fit_point, best_fit_point):
                self.assertAlmostEqual(expected, obtained, eps)

if __name__ == '__main__':
    unittest.main(verbosity=5)
