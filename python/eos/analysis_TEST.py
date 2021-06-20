import unittest
import eos
import numpy as np
import yaml

class ClassMethodTests(unittest.TestCase):

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
