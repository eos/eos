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

import builtins
import importlib
import unittest

import eos
import eos.ipython

# Alias the formatters here, since an attribute reference to a name with two leading
# underscores is mangled within a class body.
from eos.ipython import __format_Parameter         as format_Parameter
from eos.ipython import __format_KinematicVariable as format_KinematicVariable
from eos.ipython import __format_Kinematics        as format_Kinematics
from eos.ipython import __format_Options           as format_Options
from eos.ipython import __format_ObservableEntry   as format_ObservableEntry
from eos.ipython import __format_Observable        as format_Observable
from eos.ipython import __format_GoodnessOfFit     as format_GoodnessOfFit
from eos.ipython import __format_Reference         as format_Reference


class ParameterTests(unittest.TestCase):

    def test_with_latex(self):
        "A parameter that carries a LaTeX symbol is titled by that symbol."

        parameters = eos.Parameters.Defaults()
        parameter  = parameters['mass::b(MSbar)']
        result     = format_Parameter(parameter)

        self.assertIn('(eos.Parameter)', result)
        self.assertIn(f'${parameter.latex()}$', result)
        self.assertIn(str(parameter.evaluate()), result)
        self.assertIn(str(parameter.central()),  result)

    def test_without_latex(self):
        "A parameter without a LaTeX symbol is titled by its qualified name."

        parameters = eos.Parameters.Defaults()
        parameter  = parameters['Lambda_b::polarisation@unpolarised']
        self.assertEqual(parameter.latex(), '')

        result = format_Parameter(parameter)
        self.assertIn('Lambda_b::polarisation@unpolarised', result)


class KinematicsTests(unittest.TestCase):

    def test_kinematic_variable(self):
        "A single kinematic variable is rendered with its name and value."

        kinematics = eos.Kinematics(**{ 'q2': 4.5 })
        result     = format_KinematicVariable(kinematics['q2'])

        self.assertIn('(eos.KinematicVariable)', result)
        self.assertIn('q2',  result)
        self.assertIn('4.5', result)

    def test_kinematics(self):
        "A set of kinematic variables is rendered as one row per variable."

        kinematics = eos.Kinematics(**{ 'q2_min': 0.02, 'q2_max': 11.6 })
        result     = format_Kinematics(kinematics)

        self.assertEqual(result.count('<tr>'), 2)
        self.assertIn('q2_min', result)
        self.assertIn('q2_max', result)
        self.assertIn('11.6',   result)

    def test_empty_kinematics(self):
        "An empty set of kinematic variables yields an empty table."

        result = format_Kinematics(eos.Kinematics())

        self.assertNotIn('<tr>', result)
        self.assertIn('<table>',  result)
        self.assertIn('</table>', result)


class OptionsTests(unittest.TestCase):

    def test_options(self):
        "A set of options is rendered as one row per key."

        result = format_Options(eos.Options(**{ 'l': 'mu', 'q': 'd' }))

        self.assertEqual(result.count('<tr>'), 2)
        self.assertIn('<td>mu</td>', result)
        self.assertIn('<td>d</td>',  result)

    def test_empty_options(self):
        "An empty set of options yields an empty table."

        result = format_Options(eos.Options())

        self.assertNotIn('<tr>', result)


class ObservableEntryTests(unittest.TestCase):

    def test_with_kinematic_variables(self):
        "The kinematic variables of an entry are spanned by a single header cell."

        entry  = eos.Observables()['B->Dlnu::dBR/dq2;l=e,q=d']
        result = format_ObservableEntry(entry)

        kvs = [kv for kv in entry.kinematic_variables()]
        self.assertGreater(len(kvs), 0)
        self.assertIn('B->Dlnu::dBR/dq2', result)
        self.assertIn(f'<th rowspan={len(kvs)}>Kinematic Variables</th>', result)
        for kv in kvs:
            self.assertIn(f'<td>{kv}</td>', result)

    def test_with_several_kinematic_variables(self):
        "Every kinematic variable beyond the first is given a row of its own."

        entry  = eos.Observables()['0->Kpi::Im{f_+}(Re{q2},Im{q2})']
        result = format_ObservableEntry(entry)

        kvs = [kv for kv in entry.kinematic_variables()]
        self.assertEqual(len(kvs), 2)
        self.assertIn(f'<th rowspan={len(kvs)}>Kinematic Variables</th><td>{kvs[0]}</td></tr>', result)
        self.assertIn(f'<tr><td>{kvs[1]}</td></tr>', result)

    def test_without_kinematic_variables(self):
        "An entry without kinematic variables omits the corresponding header cell."

        entry  = eos.Observables()['0->Kpi::Delta_CT@KSvD2025']
        result = format_ObservableEntry(entry)

        self.assertEqual([kv for kv in entry.kinematic_variables()], [])
        self.assertIn('0->Kpi::Delta_CT@KSvD2025', result)
        self.assertNotIn('Kinematic Variables', result)


class ObservableTests(unittest.TestCase):

    def test_with_kinematics_and_options(self):
        "Kinematics and options are each spanned by a single header cell."

        parameters = eos.Parameters.Defaults()
        kinematics = eos.Kinematics(**{ 'q2': 4.5 })
        options    = eos.Options(**{ 'l': 'mu', 'q': 'd' })
        observable = eos.Observable.make('B->Dlnu::dBR/dq2', parameters, kinematics, options)
        result     = format_Observable(observable)

        options = [ok for ok, ov in observable.options()]
        self.assertIn('(eos.Observable)', result)
        self.assertIn('<th rowspan="1">kinematics</th>',              result)
        self.assertIn('<th>q2</th><td>4.5</td>',                      result)
        self.assertIn(f'<th rowspan="{len(options)}">options</th>',   result)
        self.assertIn('<tr><th>q</th><td>d</td></tr>',                result)

    def test_without_kinematics_and_options(self):
        "An observable with neither kinematics nor options reports both as 'none'."

        parameters = eos.Parameters.Defaults()
        observable = eos.Observable.make('0->pipi::b_0@KKRvD2024', parameters, eos.Kinematics(), eos.Options())
        result     = format_Observable(observable)

        self.assertEqual(result.count('<td colspan=2>none</td>'), 2)


class GoodnessOfFitTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        analysis = eos.Analysis(
            priors     = [{ 'parameter': 'B->pi::f_+(0)@BCL2008', 'min': 0.0, 'max': 1.0, 'type': 'uniform' }],
            likelihood = ['B->pi::f_+@IKMvD:2014A'],
        )
        cls.gof = eos.GoodnessOfFit(analysis._log_posterior)

    def test_totals(self):
        "The totals of the goodness of fit are reported below the per-constraint table."

        result = format_GoodnessOfFit(self.gof)

        self.assertIn(f'<td>{self.gof.total_chi_square():6.4f}</td>',        result)
        self.assertIn(f'<td>{self.gof.total_degrees_of_freedom()}</td>',     result)
        self.assertIn('<tr><th>p-value</th>',                               result)

    def test_per_constraint(self):
        "Each constraint contributes one row; an undefined signed chi is shown as a dash."

        result = format_GoodnessOfFit(self.gof)

        self.assertIn('<tt>B->pi::f_+@IKMvD:2014A</tt>', result)
        self.assertIn('&mdash;', result)


class ReferenceTests(unittest.TestCase):

    def test_on_arxiv(self):
        "A reference with an arXiv eprint is linked to its arXiv abstract page."

        reference = eos.References()['IKMvD:2014A']
        result    = format_Reference(reference)

        self.assertEqual(reference.eprint_archive(), 'arXiv')
        arxiv_id = reference.eprint_id().split(':')[-1]
        self.assertIn(f'<a href="http://arxiv.org/abs/{arxiv_id}">', result)
        self.assertIn(reference.title(), result)

    def test_not_on_arxiv(self):
        "A reference without an arXiv eprint is rendered without a link."

        reference = eos.References()['ATLAS:2013A']
        result    = format_Reference(reference)

        self.assertNotEqual(reference.eprint_archive(), 'arXiv')
        self.assertNotIn('<a href', result)
        self.assertIn(reference.title(), result)


class DetectionTests(unittest.TestCase):

    def test_outside_ipython(self):
        "Outside of IPython the integration is reported as unavailable."

        self.assertFalse(eos.ipython.__ipython__)

    def test_inside_ipython(self):
        "Within IPython the integration is reported as available."

        builtins.__IPYTHON__ = True
        try:
            self.assertTrue(importlib.reload(eos.ipython).__ipython__)
        finally:
            del builtins.__IPYTHON__
            importlib.reload(eos.ipython)


if __name__ == '__main__':
    unittest.main(verbosity=5)
