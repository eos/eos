# vim: set sts=4 tw=120 :

import unittest
import inspect
import os
import eos
import matplotlib

# EOS uses \text{...} throughout, which plain LaTeX does not provide
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}\usepackage{amssymb}'

class DocTests(unittest.TestCase):
    def test_000_Parameters(self):
        """Check the latex representation of all the parameters"""
        from matplotlib.texmanager import TexManager

        parameters = eos.Parameters.Defaults()
        texmanager = TexManager()

        for section in eos.Parameters.Defaults().sections():
            for group in section:
                for parameter in group:
                    latex_string = parameter.latex()
                    if latex_string:
                        try:
                            texmanager.get_text_width_height_descent('$' + latex_string + '$', fontsize=12)
                        except Exception as e:
                            self.fail(f"Cannot compile latex representation of parameter {parameter.name()}, caucht exception of type {type(e).__name__}: {e}")

    def test_001_Observables(self):
        """Check the latex representation of all the observables"""
        from matplotlib.texmanager import TexManager

        texmanager = TexManager()

        for section in eos.Observables().sections():
            for group in section:
                for qn, entry in group:
                    latex_string = entry.latex()
                    if latex_string:
                        try:
                            texmanager.get_text_width_height_descent('$' + latex_string + '$', fontsize=12)
                        except Exception as e:
                            self.fail(f"Cannot compile latex representation of observable {qn}, caught exception of type {type(e).__name__}: {e}")

    def test_002_References(self):
        """Check the latex representation of all reference titles"""
        from matplotlib.texmanager import TexManager

        texmanager = TexManager()

        for ref_id, ref in eos.References():
            latex_string = ref.title()
            if latex_string:
                try:
                    texmanager.get_text_width_height_descent(latex_string, fontsize=12)
                except Exception as e:
                    self.fail(f"Cannot compile latex representation of reference {ref_id} title, caught exception of type {type(e).__name__}: {e}")


# Run new tests
if __name__ == '__main__':
    unittest.main()
