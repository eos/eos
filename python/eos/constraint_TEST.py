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


class FilterTests(unittest.TestCase):

    _QN = eos.QualifiedName('B->pi::f_+@IKMvD:2014A')

    def test_no_filter(self):
        "Without a filter every entry is accepted."

        self.assertTrue(eos.Constraints().filter_entry(self._QN))

    def test_prefix(self):
        "The prefix filter matches a substring of the prefix part."

        self.assertTrue (eos.Constraints(prefix='B->pi').filter_entry(self._QN))
        self.assertFalse(eos.Constraints(prefix='B->D').filter_entry(self._QN))

    def test_name(self):
        "The name filter matches a substring of the name part."

        self.assertTrue (eos.Constraints(name='f_+').filter_entry(self._QN))
        self.assertFalse(eos.Constraints(name='f_0').filter_entry(self._QN))

    def test_suffix(self):
        "The suffix filter matches a substring of the suffix part."

        self.assertTrue (eos.Constraints(suffix='IKMvD').filter_entry(self._QN))
        self.assertFalse(eos.Constraints(suffix='KMO').filter_entry(self._QN))

    def test_combined(self):
        "All provided filters must match."

        self.assertTrue (eos.Constraints(prefix='B->pi', name='f_+', suffix='IKMvD').filter_entry(self._QN))
        self.assertFalse(eos.Constraints(prefix='B->pi', name='f_+', suffix='KMO').filter_entry(self._QN))


class ContainsTests(unittest.TestCase):

    _NAME = 'B^0->pi^+lnu::BR[0.0,4.0]@BaBar:2010A'

    def test_known_constraint(self):
        "A known constraint is found, by string as well as by qualified name."

        self.assertIn(self._NAME, eos.Constraints())
        self.assertIn(eos.QualifiedName(self._NAME), eos.Constraints())

    def test_unknown_constraint(self):
        "An unknown constraint is not found."

        self.assertNotIn('B->pi::NOSUCH@Nobody:2000A', eos.Constraints())

    def test_runtime_insertion(self):
        "A constraint inserted at run time is found by every instance created thereafter."

        name = 'B->pi::TEST@BaBar:2010A'
        entry = """
        type: Gaussian
        observable: B->pilnu::BR
        kinematics: {q2_min: 0.0, q2_max: 4.0}
        options: {q: d}
        mean: 1.0e-04
        sigma-stat: {hi: 1.0e-05, lo: -1.0e-05}
        sigma-sys: {hi: 0.0, lo: -0.0}
        references: []
        """

        constraints = eos.Constraints()
        constraints.insert(eos.QualifiedName(name), entry)

        # the inserting instance holds the entries as of its own creation
        self.assertNotIn(name, constraints)
        self.assertIn(name, eos.Constraints())


class ReprHTMLTests(unittest.TestCase):

    def test_complete_list(self):
        "The complete list of constraints renders, including entries without a reference."

        result = eos.Constraints()._repr_html_()

        self.assertIn('<tt>B->pi::f_+@IKMvD:2014A</tt>', result)
        # 'b->c::Bound[0^+]' has an empty suffix part and thus no reference
        self.assertIn('<tt>b->c::Bound[0^+]</tt>', result)
        self.assertIn('<a></a>', result)

    def test_reference_link(self):
        "A constraint taken from an arXiv eprint links to the arXiv abstract page."

        result    = eos.Constraints(name='f_+', suffix='IKMvD:2014A')._repr_html_()
        reference = eos.References()['IKMvD:2014A']
        arxiv_id  = reference.eprint_id().split(':')[-1]

        self.assertIn(f'<a href="https://arxiv.org/abs/{arxiv_id}">IKMvD:2014A</a>', result)

    def test_filtered_list(self):
        "A filter removes the entries it rejects from the table."

        result = eos.Constraints(prefix='B->pi', name='f_+')._repr_html_()

        self.assertIn   ('<tt>B->pi::f_+@IKMvD:2014A</tt>', result)
        self.assertNotIn('<tt>b->c::Bound[0^+]</tt>',       result)

    def test_observables(self):
        "The observables constrained by an entry are listed in the table."

        result = eos.Constraints(name='f_+', suffix='IKMvD:2014A')._repr_html_()
        entry  = eos.Constraints()['B->pi::f_+@IKMvD:2014A']

        for observable in entry.observables():
            self.assertIn(f'<tt>{observable}</tt>', result)


if __name__ == '__main__':
    unittest.main(verbosity=5)
