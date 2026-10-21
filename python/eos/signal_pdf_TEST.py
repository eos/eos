# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2026 Danny van Dyk
#
# This file is part of the EOS project. EOS is free software;
# you can redistribute it and/or modify it under the terms of the GNU General
# Public License version 2, as published by the Free Software Foundation.
#
# EOS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place, Suite 330, Boston, MA  02111-1307  USA

import unittest

import eos
import numpy as np
import scipy.special


def pypmc_is_missing():
    try:
        import pypmc
        return False
    except ModuleNotFoundError:
        return True


def _make_pdf():
    parameters = eos.Parameters.Defaults()
    kinematics = eos.Kinematics(**{ 'q2': 5.0, 'q2_min': 0.02, 'q2_max': 11.6 })
    options    = eos.Options(**{ 'l': 'mu' })

    return eos.SignalPDF.make('B->Dlnu::P(q2)', parameters, kinematics, options), kinematics


class MakeTests(unittest.TestCase):

    def test_variables_and_bounds(self):
        "The variables exclude the bounds, which are paired up separately."

        pdf, kinematics = _make_pdf()

        self.assertEqual([v.name() for v in pdf.variables], ['q2'])
        self.assertEqual([(lo.name(), hi.name()) for lo, hi in pdf.bounds], [('q2_min', 'q2_max')])
        self.assertEqual([lo.evaluate() for lo, hi in pdf.bounds], [0.02])
        self.assertEqual([hi.evaluate() for lo, hi in pdf.bounds], [11.6])

    def test_unknown_name(self):
        "An unknown name yields no PDF."

        parameters = eos.Parameters.Defaults()

        self.assertIsNone(eos.SignalPDF.make('nosuch::pdf', parameters, eos.Kinematics(), eos.Options()))


class LogPDFTests(unittest.TestCase):

    def test_sets_the_variables(self):
        "Evaluating the adapter sets the kinematic variables to the given phase space point."

        pdf, kinematics = _make_pdf()
        result = pdf.log_pdf([7.5])

        self.assertEqual(kinematics['q2'].evaluate(), 7.5)
        self.assertEqual(result, pdf.evaluate())

    def test_agrees_with_evaluate(self):
        "The adapter agrees with a direct evaluation at the same phase space point."

        pdf, kinematics = _make_pdf()
        kinematics['q2'].set(3.0)
        expected = pdf.evaluate()

        self.assertEqual(pdf.log_pdf([3.0]), expected)


@unittest.skipIf(pypmc_is_missing(), "Test is missing the module 'pypmc'")
class SampleMCMCTests(unittest.TestCase):

    def test_shapes(self):
        "The requested number of samples is returned, thinned by the stride."

        pdf, _ = _make_pdf()
        rng = np.random.RandomState(1701)
        samples, weights = pdf.sample_mcmc(N=20, stride=2, pre_N=20, preruns=1, rng=rng)

        self.assertEqual(samples.shape, (20, 1))
        self.assertEqual(weights.shape, (20,))

    def test_within_bounds(self):
        "All samples lie within the bounds of the phase space."

        pdf, _ = _make_pdf()
        np.random.seed(1702)
        samples, weights = pdf.sample_mcmc(N=20, stride=1, pre_N=20, preruns=1, start_point=[5.0])

        self.assertTrue(np.all(samples[:, 0] >= 0.02))
        self.assertTrue(np.all(samples[:, 0] <= 11.6))
        self.assertTrue(np.all(np.isfinite(weights)))


class DetectorLevelPDFPolesTests(unittest.TestCase):
    """
    For p(x) = |1 + R / (x - z)|^2 = 1 + Re[c / (x - z)] with c = 2 R + i |R|^2 / gamma, the remainder is constant, and the
    smeared density is 1 + Re[c S(x - z)], with S(zeta) = -i sqrt(pi / 2) / sigma w(zeta / (sqrt(2) sigma)) for a Gaussian resolution.
    """

    x_p, gamma, sigma = 0.3, 1.0e-3, 0.5
    R = 0.01 * np.exp(0.25j * np.pi)
    c = 2.0 * R + 1.0j * abs(R)**2 / gamma
    # off the grid nodes, within 3 sigma of the pole
    points = [-1.2, -0.13, 0.3013, 0.67, 1.7]

    @classmethod
    def setUpClass(cls):
        R, x_p, gamma = cls.R, cls.x_p, cls.gamma
        pole_re = f'(({R.real} * ({{x}} - {x_p}) + {R.imag} * {gamma}) / (({{x}} - {x_p})^2 + {gamma}^2))'
        pole_im = f'(({R.imag} * ({{x}} - {x_p}) - {R.real} * {gamma}) / (({{x}} - {x_p})^2 + {gamma}^2))'
        expressions = {
            'p(x)': f'(1 + {pole_re})^2 + {pole_im}^2', 'N': '{x_max} - {x_min}',
            'Re_z': f'{x_p}', 'Im_z': f'{-gamma}', 'Re_c': f'{cls.c.real}', 'Im_c': f'{cls.c.imag}', 'zero': '0',
        }
        for name, expression in expressions.items():
            eos.Observables().insert(eos.QualifiedName(f'DetectorLevelPDFPolesTest::{name}'), '', eos.Unit.Unity(), eos.Options(), expression)
        eos.SignalPDFs().insert(eos.QualifiedName('DetectorLevelPDFPolesTest::P(x)'), '', eos.Options(), eos.QualifiedName('DetectorLevelPDFPolesTest::p(x)'), ['x'],
                                eos.QualifiedName('DetectorLevelPDFPolesTest::N'), ['x_min', 'x_max'])

    def _pole(self, **coefficients):
        return { 'variable': 'x', 'position_re': 'DetectorLevelPDFPolesTest::Re_z', 'position_im': 'DetectorLevelPDFPolesTest::Im_z',
                 **{ key: f'DetectorLevelPDFPolesTest::{value}' for key, value in coefficients.items() } }

    def _make(self, poles):
        cache = eos.ObservableCache(eos.Parameters.Defaults())
        axes  = [{ 'variable': 'x', 'min': -20.0, 'max': 20.0, 'points': 1024 }]
        return eos.DetectorLevelPDF.make_from_expression(cache, 'DetectorLevelPDFPolesTest::P(x)', eos.Options(), axes,
                                                         f'Kernel::Gaussian({{x}}, 0, {self.sigma})', poles=poles), cache

    def _values(self, poles):
        pdf, cache = self._make(poles)
        cache.update()
        variable = pdf.kinematics()['x']
        values = []
        for x in self.points:
            variable.set(x)
            values.append(pdf.evaluate_linear())
        return np.array(values)

    def test_poles(self):
        "Check that a pole given as a dictionary reaches the PDF, against the smeared density in closed form."
        zeta     = np.array(self.points) - complex(self.x_p, -self.gamma)
        S        = -1.0j * np.sqrt(np.pi / 2.0) / self.sigma * scipy.special.wofz(zeta / (np.sqrt(2.0) * self.sigma))
        expected = 1.0 + (self.c * S).real
        np.testing.assert_allclose(self._values([self._pole(coefficient_re='Re_c', coefficient_im='Im_c')]), expected, rtol=1.0e-10)

    def test_no_coefficient(self):
        "Check that a pole without a subtraction coefficient has it extracted from the grid."
        np.testing.assert_allclose(self._values([self._pole()]), self._values([self._pole(coefficient_re='Re_c', coefficient_im='Im_c')]), rtol=1.0e-10)

    def test_zero_coefficient(self):
        "Check that a pole with a vanishing coefficient leaves the PDF without poles unchanged."
        np.testing.assert_allclose(self._values([self._pole(coefficient_re='zero', coefficient_im='zero')]), self._values([]), rtol=1.0e-14)

    def test_one_coefficient(self):
        "Check that a pole with only one part of the coefficient is rejected."
        with self.assertRaises(ValueError):
            self._make([self._pole(coefficient_re='Re_c')])

    def test_unknown_or_missing_key(self):
        "Check that a pole with a misspelled or a missing key is rejected."
        with self.assertRaisesRegex(ValueError, 'coefficent_re'):
            self._make([self._pole(coefficent_re='Re_c', coefficent_im='Im_c')])
        pole = self._pole()
        del pole['variable']
        with self.assertRaisesRegex(ValueError, 'variable'):
            self._make([pole])


class SignalPDFsContainsTests(unittest.TestCase):

    def test_known_signal_pdf(self):
        "A known signal PDF is found, by string as well as by qualified name."

        self.assertIn('B->Dlnu::P(q2)', eos.SignalPDFs())
        self.assertIn(eos.QualifiedName('B->Dlnu::P(q2)'), eos.SignalPDFs())

    def test_unknown_signal_pdf(self):
        "An unknown signal PDF is not found."

        self.assertNotIn('B->pi::NOSUCH', eos.SignalPDFs())


class SignalPDFsFilterTests(unittest.TestCase):

    _QN = eos.QualifiedName('B->Dlnu::P(q2)')

    def test_no_filter(self):
        "Without a filter every entry is accepted."

        self.assertTrue(eos.SignalPDFs().filter_entry(self._QN))

    def test_prefix(self):
        "The prefix filter matches a substring of the prefix part."

        self.assertTrue (eos.SignalPDFs(prefix='B->Dlnu').filter_entry(self._QN))
        self.assertFalse(eos.SignalPDFs(prefix='B->K^*ll').filter_entry(self._QN))

    def test_name(self):
        "The name filter matches a substring of the name part."

        self.assertTrue (eos.SignalPDFs(name='P(q2)').filter_entry(self._QN))
        self.assertFalse(eos.SignalPDFs(name='cos(theta_l)').filter_entry(self._QN))

    def test_suffix(self):
        "The suffix filter rejects an entry without a suffix part."

        self.assertFalse(eos.SignalPDFs(suffix='ABC:2020A').filter_entry(self._QN))


class SignalPDFsReprHTMLTests(unittest.TestCase):

    def test_complete_list(self):
        "The complete list of PDFs renders with one row per entry."

        result = eos.SignalPDFs()._repr_html_()

        self.assertIn('<tt>B->Dlnu::P(q2)</tt>',   result)
        self.assertIn('<tt>B->D^*lnu::P(q2)</tt>', result)

    def test_filtered_list(self):
        "A filter removes the entries it rejects from the table."

        result = eos.SignalPDFs(prefix='B->Dlnu')._repr_html_()

        self.assertIn   ('<tt>B->Dlnu::P(q2)</tt>',   result)
        self.assertNotIn('<tt>B->D^*lnu::P(q2)</tt>', result)

    def test_empty_list(self):
        "A filter that rejects every entry yields a table without rows."

        result = eos.SignalPDFs(prefix='nosuchprefix')._repr_html_()

        self.assertNotIn('<tt>', result)


if __name__ == '__main__':
    unittest.main(verbosity=5)
