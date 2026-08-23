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
