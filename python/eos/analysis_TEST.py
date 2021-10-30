import unittest

from numpy import random
import eos
import numpy as np
import yaml

class ClassMethodTests(unittest.TestCase):

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
            np.max(analysis._par_to_x(point)) < 1.0
            )
        self.assertTrue(
            np.min(analysis._par_to_x(point)) > -1.0
            )
        self.assertTrue(
            np.max(analysis._x_to_par(analysis._par_to_x(point)) - point) < 1e-10
            )


    def test_sanitize_manual_input(self):

        types  = np.array(["Gaussian", "Gaussian"])
        values = np.array([4.123, 3.141])

        input_output_cases = [
            {'dirty': np.float64(1.61), 'clean': float(1.61)},
            {'dirty': np.str_('Gaussian'), 'clean': str("Gaussian")},
            {
             'dirty': {'type': types[0], 'mean': values[0]}, # Simplified manual contraint with numpy types
             'clean': {'type': str('Gaussian'), 'mean': float(4.123)}
            },
            {
             'dirty': {'type': types, 'mean': values}, # Simplified manual contraint with numpy arrays
             'clean': {'type': [str('Gaussian'), str('Gaussian')], 'mean': [float(4.123), float(3.141)]}
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


if __name__ == '__main__':
    unittest.main(verbosity=5)
