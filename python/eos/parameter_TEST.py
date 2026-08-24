import os
import tempfile
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


class _StubParameter:
    """Minimal stand-in for :class:`eos.Parameter`, used to unit-test the sort key."""

    def __init__(self, name):
        self._name = name

    def name(self):
        return self._name


class _StubWilsonCoefficients:
    """Minimal stand-in for the ``wc`` attribute of a :class:`wilson.Wilson` object."""

    def __init__(self, values, basis='EOS', scale=4.2):
        self.basis = basis
        self.scale = scale
        self.dict  = values


class _StubWilson:
    """Minimal stand-in for :class:`wilson.Wilson`."""

    def __init__(self, *args, **kwargs):
        self.wc = _StubWilsonCoefficients(*args, **kwargs)


class FilterTests(unittest.TestCase):

    _QN = 'B->pi::f_+(0)@BCL2008'

    def test_no_filter(self):
        "Without a filter every entry is accepted."

        self.assertTrue(eos.Parameters().filter_entry(self._QN))

    def test_prefix(self):
        "The prefix filter matches a substring of the prefix part."

        self.assertTrue (eos.Parameters(prefix='B->pi').filter_entry(self._QN))
        self.assertFalse(eos.Parameters(prefix='B->D').filter_entry(self._QN))

    def test_name(self):
        "The name filter matches a substring of the name part."

        self.assertTrue (eos.Parameters(name='f_+').filter_entry(self._QN))
        self.assertFalse(eos.Parameters(name='f_0').filter_entry(self._QN))

    def test_suffix(self):
        "The suffix filter matches a substring of the suffix part."

        self.assertTrue (eos.Parameters(suffix='BCL2008').filter_entry(self._QN))
        self.assertFalse(eos.Parameters(suffix='BSZ2015').filter_entry(self._QN))

    def test_combined(self):
        "All provided filters must match."

        self.assertTrue (eos.Parameters(prefix='B->pi', name='f_+', suffix='BCL2008').filter_entry(self._QN))
        self.assertFalse(eos.Parameters(prefix='B->pi', name='f_+', suffix='BSZ2015').filter_entry(self._QN))


class ContainsTests(unittest.TestCase):

    def test_known_parameter(self):
        "A known parameter is found, by string as well as by qualified name."

        self.assertIn('mass::tau', eos.Parameters.Defaults())
        self.assertIn(eos.QualifiedName('mass::tau'), eos.Parameters.Defaults())

    def test_unknown_parameter(self):
        "An unknown parameter is not found."

        self.assertNotIn('mass::boing747', eos.Parameters.Defaults())


class KeyTests(unittest.TestCase):

    def test_qualified_name(self):
        "A qualified name is sorted by prefix, then suffix, then name."

        parameter = eos.Parameters.Defaults()['B->pi::f_+(0)@BCL2008']

        self.assertEqual(eos.Parameters._key(parameter), ('B->pi', 'BCL2008', 'f_+(0)'))

    def test_unqualified_name(self):
        "A name that is not a qualified name is sorted by the name alone."

        self.assertEqual(eos.Parameters._key(_StubParameter('not-a-qualified-name')), ('', '', 'not-a-qualified-name'))


class LatexRefineTests(unittest.TestCase):

    def test_gev(self):
        "The \\GeV macro is replaced by plain text."

        self.assertEqual(eos.Parameters()._latex_refine(r'\GeV'), r'\text{GeV}')

    def test_display_mode(self):
        "Inline math is promoted to display mode."

        self.assertEqual(eos.Parameters()._latex_refine(r'$m_b$'), r'$$m_b$$')

    def test_combined(self):
        "Both replacements are applied."

        self.assertEqual(eos.Parameters()._latex_refine(r'$m_b / \GeV$'), r'$$m_b / \text{GeV}$$')


class ReprHTMLTests(unittest.TestCase):

    def test_complete_list(self):
        "The complete list of parameters renders with one row per parameter."

        result = eos.Parameters()._repr_html_()

        self.assertIn('<th>qualified name</th>', result)
        self.assertIn('<tt style="color:grey">B->pi::f_+(0)@BCL2008</tt>', result)
        self.assertIn('<tt style="color:grey">mass::b(MSbar)</tt>',        result)

    def test_filtered_list(self):
        "A filter removes both the entries it rejects and the groups left empty."

        result = eos.Parameters(prefix='B->pi', suffix='BCL2008')._repr_html_()

        self.assertIn   ('<tt style="color:grey">B->pi::f_+(0)@BCL2008</tt>', result)
        self.assertNotIn('<tt style="color:grey">mass::b(MSbar)</tt>',        result)

    def test_empty_list(self):
        "A filter that rejects every parameter yields a table without rows."

        result = eos.Parameters(prefix='nosuchprefix')._repr_html_()

        self.assertNotIn('color:grey', result)

    def test_unit(self):
        "A parameter with a unit shows it in display mode; one without shows a dash."

        parameters = eos.Parameters.Defaults()
        latex      = parameters['mass::pi^0@HME'].unit().latex()

        self.assertIn(fr'$$\left[ {latex} \right]$$', eos.Parameters(name='pi^0')._repr_html_())
        self.assertIn('&mdash;', eos.Parameters(name='R^1')._repr_html_())

    def test_without_symbol(self):
        "A parameter without a LaTeX symbol shows a placeholder."

        parameters = eos.Parameters.Defaults()
        self.assertEqual(parameters['Lambda_b::polarisation@unpolarised'].latex(), '')

        result = eos.Parameters(name='polarisation')._repr_html_()
        self.assertIn('>---<', result)


class ToYAMLTests(unittest.TestCase):

    def test_selected_names(self):
        "Only the requested parameters are converted."

        parameters = eos.Parameters()
        contents   = yaml.safe_load(parameters.to_yaml(names=['mass::b(MSbar)']))

        self.assertEqual(list(contents.keys()), ['mass::b(MSbar)'])
        entry = contents['mass::b(MSbar)']
        self.assertEqual(entry['central'], parameters['mass::b(MSbar)'].evaluate())
        self.assertEqual(entry['min'],     parameters['mass::b(MSbar)'].min())
        self.assertEqual(entry['max'],     parameters['mass::b(MSbar)'].max())
        self.assertEqual(entry['latex'],   parameters['mass::b(MSbar)'].latex())

    def test_all_names(self):
        "Without a selection every parameter is converted."

        parameters = eos.Parameters()
        contents   = yaml.safe_load(parameters.to_yaml())

        self.assertEqual(len(contents), len([p for p in parameters]))
        self.assertIn('mass::b(MSbar)',        contents)
        self.assertIn('B->pi::f_+(0)@BCL2008', contents)

    def test_dump(self):
        "Dumping to a file yields the same contents as the conversion."

        parameters = eos.Parameters()
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'parameters.yaml')
            parameters.dump(path, names=['mass::b(MSbar)'])
            with open(path) as f:
                contents = f.read()

        self.assertEqual(contents, parameters.to_yaml(names=['mass::b(MSbar)']))


class FromWCxfStubTests(unittest.TestCase):

    def test_wrong_basis(self):
        "Wilson coefficients must be provided in the EOS basis."

        with self.assertRaisesRegex(ValueError, 'EOS basis'):
            eos.Parameters.FromWCxf(_StubWilson({}, basis='flavio'))

    def test_wrong_scale(self):
        "Wilson coefficients must be provided at the reference scale."

        with self.assertRaisesRegex(ValueError, '4.2 GeV'):
            eos.Parameters.FromWCxf(_StubWilson({}, scale=5.0))

    def test_real_coefficient(self):
        "A coefficient that EOS treats as real-valued is added to its default value."

        defaults = eos.Parameters.Defaults()
        result   = eos.Parameters.FromWCxf(_StubWilson({ 'b->s::c1': complex(0.25, 0.0) }))

        self.assertAlmostEqual(
            result['b->s::c1'].evaluate(),
            defaults['b->s::c1'].evaluate() + 0.25,
            delta = 1.0e-10
        )

    def test_real_coefficient_with_imaginary_part(self):
        "A coefficient that EOS treats as real-valued must not carry an imaginary part."

        with self.assertRaisesRegex(ValueError, 'imaginary part'):
            eos.Parameters.FromWCxf(_StubWilson({ 'b->s::c1': complex(0.25, 0.1) }))

    def test_complex_coefficient(self):
        "A complex-valued coefficient is split into its real and imaginary parts."

        defaults = eos.Parameters.Defaults()
        result   = eos.Parameters.FromWCxf(_StubWilson({ 'cbtaunutau::cSL': complex(0.3, 0.1) }))

        self.assertAlmostEqual(
            result['cbtaunutau::Re{cSL}'].evaluate(),
            defaults['cbtaunutau::Re{cSL}'].evaluate() + 0.3,
            delta = 1.0e-10
        )
        self.assertAlmostEqual(
            result['cbtaunutau::Im{cSL}'].evaluate(),
            defaults['cbtaunutau::Im{cSL}'].evaluate() + 0.1,
            delta = 1.0e-10
        )


if __name__ == '__main__':
    unittest.main(verbosity=5)
