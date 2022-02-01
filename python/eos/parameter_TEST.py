import unittest
import eos
import numpy as np
import yaml

def wilson_is_missing():
    try:
        import wilson
        return False
    except ModuleNotFoundError:
        return True

class ClassMethodTests(unittest.TestCase):

    @unittest.skipIf(wilson_is_missing(), "Test is missing the module 'wilson'")
    def test_importing_from_wcxf(self):

        from wilson import Wilson

        inputs = {
            'cbtaunutau::cSL':  0.3 + 0.1j,
            'ubenue::cVL':     -0.1 + 0.2j,
        }
        w = Wilson(inputs, scale=4.2, eft="WET", basis="EOS")
        p = eos.Parameters.FromWCxf(w)

        pdefaults = eos.Parameters.Defaults()
        outputs = {
            'cbtaunutau::Re{cSL}': +0.3,
            'cbtaunutau::Im{cSL}': +0.1,
            'ubenue::Re{cVL}':     -0.1,
            'ubenue::Im{cVL}':     +0.2,
        }
        outputs = { ok: ov + pdefaults[ok].evaluate() for ok, ov in outputs.items() }
        for ok, ov in outputs.items():
            self.assertAlmostEqual(
                p[ok].evaluate(),
                ov,
                delta = 1.0e-10
            )

if __name__ == '__main__':
    unittest.main(verbosity=5)
