# vim: set sts=4 tw=120 :

import unittest
import inspect
import os
import eos

class DocTests(unittest.TestCase):
    def test_000_Parameters(self):
        """Check the latex representation of all the parameters"""
        from matplotlib.texmanager import TexManager

        parameters = eos.Parameters.Defaults()
        texmanager = TexManager()

        for section in eos.Parameters.Defaults().sections():
            for group in section:
                for parameter in group:
                    try:
                        texmanager.get_text_width_height_descent(parameter.latex(), fontsize=12)
                    except Exception as e:
                        self.fail(f"Cannot compile latex representation of parameter {parameter.name()}, caucht exception of type {type(e).__name__}: {e}")


# Run new tests
if __name__ == '__main__':
    unittest.main()
