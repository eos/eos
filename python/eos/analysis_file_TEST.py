import unittest

import contextlib
import eos
import os
from pathlib import Path

class TestAnalysisFile(unittest.TestCase):

    def test_pyhf_likelihood(self):

        af = eos.AnalysisFile(Path(__file__).parent / 'analysis_file_TEST.d' / 'valid-analysis-file.yaml')
        af.validate()


if __name__ == '__main__':
    unittest.main(verbosity=5)
