import unittest
import eos
import numpy as np
import yaml

class ClassMethodTests(unittest.TestCase):

    def test_importing_from_wcxf(self):

        from wilson import Wilson

        inputs = {
	    'cbtaunutau::cSL':  0.3 + 0.1j,
            'ubenue::cVL':     -0.1 + 0.2j,
        }
        w = Wilson(inputs, scale=4.2, eft="WET", basis="EOS")
        p = eos.Parameters.FromWCxf(w)


if __name__ == '__main__':
    unittest.main(verbosity=5)