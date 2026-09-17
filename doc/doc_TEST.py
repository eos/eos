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

    def test_003_math_conversion(self):
        """Check for any lingering ':math:' in the generated html documentation"""
        html_directory = os.path.join(os.environ.get("BUILDDIR"), 'html')
        unconverted_math = []
        checked_files = []
        for root, _, files in os.walk(html_directory):
            for filename in files:
                if not filename.endswith('.html'):
                    continue
                checked_files.append(filename)

                path = os.path.join(root, filename)
                with open(path, encoding='utf-8') as html_file:
                    for line_number, line in enumerate(html_file, start=1):
                        if ':math:' in line:
                            unconverted_math.append(f'{os.path.abspath(path)}:{line_number}')

        if not checked_files:
            self.fail('Found no HTML files in the documentation build')

        if unconverted_math:
            self.fail('Found unconverted :math: markup in the following HTML files:\n' + '\n'.join(unconverted_math))


# Run new tests
if __name__ == '__main__':
    unittest.main()
