import unittest

import eos
import os


class ImportanceSamplesTests(unittest.TestCase):

    def test_importing_samples(self):
        "Test the import of importance samples."

        file = eos.data.ImportanceSamples(os.path.join(os.environ['SOURCE_DIR'], "eos/data/importance_samples_TEST.d/samples"))


if __name__ == '__main__':
    unittest.main(verbosity=5)
