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
