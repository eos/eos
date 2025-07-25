import unittest

import eos
import eos.figure
import yaml

from matplotlib import pyplot as plt

class ObservableItemTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: observable
            observable: 'B->Dlnu::dBR/dq2'
            variable: q2
            range: [0.1, 1.0]
            resolution: 100
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare()
            fig, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'observable': {e}")

if __name__ == '__main__':
    unittest.main(verbosity=5)
