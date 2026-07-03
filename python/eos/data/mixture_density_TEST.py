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
import numpy as np
import os
import tempfile
import yaml


def pypmc_is_missing():
    try:
        import pypmc
        return False
    except ImportError:
        return True


class MixtureDensityTests(unittest.TestCase):

    _path = os.path.join(os.environ['SOURCE_DIR'], 'eos/data/mixture_density_TEST.d')

    @staticmethod
    def _write(directory, description):
        "Materialize a (possibly malformed) MixtureDensity description in the given directory."
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, 'description.yaml'), 'w') as f:
            yaml.safe_dump(description, f, default_flow_style=False)

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

    def test_loading_from_file(self):
        "The description.yaml file itself is also accepted as the path."
        md = eos.data.MixtureDensity(os.path.join(self._path, 'description.yaml'))
        self.assertEqual(md.type, 'MixtureDensity')
        self.assertEqual(len(md.components), 2)

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

    @unittest.skipIf(pypmc_is_missing(), "Test requires the 'pypmc' module")
    def test_roundtrip(self):
        "A mixture written by create() reads back identically (write/read symmetry)."
        import pypmc

        means   = [np.array([3.7e-3, 0.27]), np.array([4.0e-3, 0.30])]
        covs    = [np.array([[1.0e-7, 0.0], [0.0, 1.0e-3]]), np.array([[2.0e-7, 0.0], [0.0, 2.0e-3]])]
        weights = np.array([0.4, 0.6])
        density = pypmc.density.mixture.create_gaussian_mixture(means, covs, weights)

        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'clusters')
            eos.data.MixtureDensity.create(path, density, varied_parameters=['a', 'b'])

            md = eos.data.MixtureDensity(path)
            self.assertEqual(md.type, 'MixtureDensity')
            self.assertEqual(len(md.components), 2)
            self.assertEqual(md.components[0]['type'], 'gauss')
            np.testing.assert_allclose(md.components[0]['mu'], means[0])
            np.testing.assert_allclose(md.components[0]['sigma'], covs[0])
            np.testing.assert_allclose(md.weights, [0.4, 0.6])
            self.assertEqual(md.varied_parameters, ['a', 'b'])
            self.assertEqual(md.test_statistics, {'sigma': [], 'densities': []})

    def test_wrong_type(self):
        "A description whose type is not 'MixtureDensity' is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'PMCSampler', 'components': [], 'weights': []})
            with self.assertRaises(ValueError):
                eos.data.MixtureDensity(d)

    def test_unknown_key(self):
        "An unexpected key in the description is rejected rather than silently ignored."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'MixtureDensity',
                'components': [{'type': 'gauss', 'mu': [0.0], 'sigma': [[1.0]]}],
                'weights': [1.0], 'bogus': 42,
            })
            with self.assertRaises(ValueError):
                eos.data.MixtureDensity(d)

    def test_unsupported_component_type(self):
        "A mixture component of an unsupported type is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'MixtureDensity',
                'components': [{'type': 'student-t', 'mu': [0.0], 'sigma': [[1.0]]}],
                'weights': [1.0],
            })
            with self.assertRaises(ValueError):
                eos.data.MixtureDensity(d)

    def test_missing_weights(self):
        "A description missing the required 'weights' key is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'MixtureDensity',
                'components': [{'type': 'gauss', 'mu': [0.0], 'sigma': [[1.0]]}],
            })
            with self.assertRaises(ValueError):
                eos.data.MixtureDensity(d)

    def test_non_mapping_description(self):
        "A description.yaml that is not a top-level mapping is rejected."
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'description.yaml'), 'w') as f:
                yaml.safe_dump(['not', 'a', 'mapping'], f)
            with self.assertRaises(RuntimeError):
                eos.data.MixtureDensity(d)


if __name__ == '__main__':
    unittest.main(verbosity=5)
