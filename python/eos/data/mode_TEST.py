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


class ModeTests(unittest.TestCase):

    _path = os.path.join(os.environ['SOURCE_DIR'], 'eos/data/mode_TEST.d')

    def test_loading(self):
        "Load a Mode from a prepared fixture and check its contents."
        m = eos.data.Mode(self._path)

        self.assertEqual(m.type, 'Mode')
        self.assertEqual([p['name'] for p in m.varied_parameters], ['CKM::abs(V_ub)', 'B->pi::f_+(0)@BCL2008'])
        self.assertEqual(m.mode, [3.67e-3, 0.27])
        self.assertAlmostEqual(m.pvalue, 0.25)
        self.assertEqual(m.local_pvalues, {'B->pilnu::BR': 0.5})
        self.assertAlmostEqual(m.global_chi2, 12.0)
        self.assertAlmostEqual(m.dof, 5.0)

    def test_invalid_path(self):
        "Loading from a non-existent path raises an error."
        with self.assertRaises(RuntimeError):
            eos.data.Mode(os.path.join(self._path, 'does-not-exist'))


if __name__ == '__main__':
    unittest.main(verbosity=5)
