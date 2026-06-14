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


def pypmc_is_missing():
    try:
        import pypmc
        return False
    except ImportError:
        return True


class MixtureDensityTests(unittest.TestCase):

    _path = os.path.join(os.environ['SOURCE_DIR'], 'eos/data/mixture_density_TEST.d')

    def test_loading(self):
        "Load a MixtureDensity from a prepared fixture and check its contents."
        md = eos.data.MixtureDensity(self._path)

        self.assertEqual(md.type, 'MixtureDensity')
        self.assertEqual(len(md.components), 2)
        self.assertEqual(md.components[0]['type'], 'gauss')
        self.assertEqual(md.components[0]['mu'], [3.7e-3, 0.27])
        self.assertEqual(md.components[0]['sigma'], [[1.0e-7, 0.0], [0.0, 1.0e-3]])
        self.assertEqual(md.weights, [0.4, 0.6])
        self.assertEqual(md.varied_parameters, ['CKM::abs(V_ub)', 'B->pi::f_+(0)@BCL2008'])
        self.assertEqual(md.test_statistics, {'sigma': [], 'densities': []})

    @unittest.skipIf(pypmc_is_missing(), "Test requires the 'pypmc' module")
    def test_density(self):
        "Build a pypmc mixture density from the loaded fixture."
        import pypmc

        md      = eos.data.MixtureDensity(self._path)
        density = md.density()
        self.assertIsInstance(density, pypmc.density.mixture.MixtureDensity)
        self.assertEqual(len(density.components), 2)

    def test_invalid_path(self):
        "Loading from a non-existent path raises an error."
        with self.assertRaises(RuntimeError):
            eos.data.MixtureDensity(os.path.join(self._path, 'does-not-exist'))


if __name__ == '__main__':
    unittest.main(verbosity=5)
