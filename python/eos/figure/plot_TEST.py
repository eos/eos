import unittest

import eos
import eos.figure
import yaml

from matplotlib import pyplot as plt

class PlotTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            items:
              - type: 'observable'
                observable: 'B->Dlnu::dBR/dq2'
                options: { 'l': 'e' }
                variable: 'q2'
                x_values: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
              - type: 'observable'
                observable: 'B->Dlnu::dBR/dq2'
                options: { 'l': 'mu' }
                variable: 'q2'
                x_values: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
            """
            plot = eos.figure.PlotFactory.from_yaml(input)
            plot.prepare()
            fig, ax = plt.subplots()
            plot.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing plot of type '2D': {e}")

if __name__ == '__main__':
    unittest.main(verbosity=5)
