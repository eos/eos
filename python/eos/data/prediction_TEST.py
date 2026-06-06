import unittest

import eos
import os


class PredictionTests(unittest.TestCase):

    def test_importing_predictions(self):
        "Test the import of prediction samples."

        file = eos.data.Prediction(os.path.join(os.environ['SOURCE_DIR'], "eos/data/prediction_TEST.d/predictions-binned"))
        file = eos.data.Prediction(os.path.join(os.environ['SOURCE_DIR'], "eos/data/prediction_TEST.d/predictions"))


if __name__ == '__main__':
    unittest.main(verbosity=5)
