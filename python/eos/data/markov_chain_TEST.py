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


class MarkovChainTests(unittest.TestCase):

    _path = os.path.join(os.environ['SOURCE_DIR'], 'eos/data/markov_chain_TEST.d')

    def test_loading(self):
        "Load a MarkovChain from a prepared fixture and check its contents."
        mc = eos.data.MarkovChain(self._path)

        self.assertEqual(mc.type, 'MarkovChain')
        self.assertEqual(mc.lookup_table, {'CKM::abs(V_ub)': 0, 'B->pi::f_+(0)@BCL2008': 1})
        self.assertEqual(mc.samples.shape, (5, 2))
        self.assertEqual(mc.usamples.shape, (5, 2))
        self.assertIsNotNone(mc.weights)
        self.assertEqual(mc.weights.shape, (5,))

    def test_invalid_path(self):
        "Loading from a non-existent path raises an error."
        with self.assertRaises(RuntimeError):
            eos.data.MarkovChain(os.path.join(self._path, 'does-not-exist'))


if __name__ == '__main__':
    unittest.main(verbosity=5)
