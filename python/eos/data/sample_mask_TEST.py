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


class SampleMaskTests(unittest.TestCase):

    _path = os.path.join(os.environ['SOURCE_DIR'], 'eos/data/sample_mask_TEST.d')

    def test_loading(self):
        "Load a SampleMask from a prepared fixture and check its contents."
        sm = eos.data.SampleMask(self._path)

        self.assertEqual(sm.type, 'Mask')
        self.assertEqual(sm.observables,
            ['B->Dlnu::BR', 'B->pilnu::BR', 'B->Dlnu::dBR/dq2', 'B->pilnu::dBR/dq2', 'B_u->lnu::BR'])
        self.assertEqual(sm.mask.dtype, bool)
        self.assertEqual(sm.mask.shape, (5,))
        self.assertEqual(list(sm.mask), [True, False, True, True, False])

    def test_invalid_path(self):
        "Loading from a non-existent path raises an error."
        with self.assertRaises(RuntimeError):
            eos.data.SampleMask(os.path.join(self._path, 'does-not-exist'))


if __name__ == '__main__':
    unittest.main(verbosity=5)
