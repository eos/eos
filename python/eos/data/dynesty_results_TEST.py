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
import eos.data
import os


def dynesty_is_missing():
    try:
        import dynesty
        return False
    except ImportError:
        return True


class DynestyResultsTests(unittest.TestCase):

    _path = os.path.join(os.environ['SOURCE_DIR'], 'eos/data/dynesty_results_TEST.d')

    @unittest.skipIf(dynesty_is_missing(), "Test requires the 'dynesty' module")
    def test_loading(self):
        "Load a DynestyResults object from a prepared fixture and check its contents."
        dr = eos.data.DynestyResults(self._path)

        self.assertEqual(dr.type, 'DynestyResults')
        self.assertEqual(dr.lookup_table, {'CKM::abs(V_ub)': 0, 'B->pi::f_+(0)@BCL2008': 1})
        # the reconstructed samples have one column per varied parameter
        self.assertEqual(dr.samples.shape[1], 2)
        # there is exactly one weight per sample
        self.assertEqual(len(dr.weights), dr.samples.shape[0])

    @unittest.skipIf(dynesty_is_missing(), "Test requires the 'dynesty' module")
    def test_invalid_path(self):
        "Loading from a non-existent path raises an error."
        with self.assertRaises(RuntimeError):
            eos.data.DynestyResults(os.path.join(self._path, 'does-not-exist'))


if __name__ == '__main__':
    unittest.main(verbosity=5)
