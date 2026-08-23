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

from _eos import ReferenceName


class FilterTests(unittest.TestCase):

    _RN = ReferenceName('IKMvD:2014A')

    def test_no_filter(self):
        "Without a filter every entry is accepted."

        self.assertTrue(eos.References().filter_entry(self._RN))

    def test_year(self):
        "The year filter matches a substring of the year part."

        self.assertTrue (eos.References(year='2014').filter_entry(self._RN))
        self.assertTrue (eos.References(year='201').filter_entry(self._RN))
        self.assertFalse(eos.References(year='2015').filter_entry(self._RN))

    def test_index(self):
        "The index filter matches a substring of the index part."

        self.assertTrue (eos.References(index='A').filter_entry(self._RN))
        self.assertFalse(eos.References(index='B').filter_entry(self._RN))

    def test_combined(self):
        "All provided filters must match."

        self.assertTrue (eos.References(year='2014', index='A').filter_entry(self._RN))
        self.assertFalse(eos.References(year='2014', index='B').filter_entry(self._RN))


class ReprHTMLTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.result = eos.References()._repr_html_()

    def test_eprint_link(self):
        "A reference with an arXiv eprint links to the arXiv abstract page."

        reference = eos.References()['IKMvD:2014A']
        arxiv_id  = reference.eprint_id().split(':')[-1]

        self.assertEqual(reference.eprint_archive(), 'arXiv')
        self.assertIn(f'<a  href="https://arxiv.org/abs/{arxiv_id}"><tt>IKMvD:2014A</tt></a>', self.result)
        self.assertIn(reference.title(), self.result)

    def test_without_eprint_link(self):
        "A reference without an arXiv eprint is rendered without a link target."

        reference = eos.References()['ATLAS:2013A']

        self.assertNotEqual(reference.eprint_archive(), 'arXiv')
        self.assertIn('<a ><tt>ATLAS:2013A</tt></a>', self.result)
        self.assertIn(reference.title(), self.result)

    def test_authors(self):
        "Authors are separated onto individual lines."

        reference = eos.References()['IKMvD:2014A']
        authors   = [a.strip() for a in reference.authors().split(' and ')]

        self.assertGreater(len(authors), 1)
        self.assertIn('<br/>'.join(authors), self.result)

    def test_filtered_list(self):
        "A filter removes the entries it rejects from the table."

        result = eos.References(year='2014', index='A')._repr_html_()

        self.assertIn   ('<tt>IKMvD:2014A</tt>', result)
        self.assertNotIn('<tt>ATLAS:2013A</tt>', result)


if __name__ == '__main__':
    unittest.main(verbosity=5)
