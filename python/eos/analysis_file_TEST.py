import unittest

import contextlib
import eos
import os
from pathlib import Path

class TestAnalysisFile(unittest.TestCase):

    def test_analysis_file(self):

        af = eos.AnalysisFile(Path(__file__).parent / 'analysis_file_TEST.d' / 'valid-analysis-file.yaml')
        af.validate()

    def test_minimal_analysis_file(self):

        af = eos.AnalysisFile(Path(__file__).parent / 'analysis_file_TEST.d' / 'minimal-analysis-file.yaml')
        af.validate()


if __name__ == '__main__':
    unittest.main(verbosity=5)
