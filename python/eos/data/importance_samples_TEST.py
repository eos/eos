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


class _Parameter:
    "Minimal stand-in for eos.Parameter, exposing the name/min/max accessors that create() uses."

    def __init__(self, name, minimum, maximum):
        self._name = name
        self._min  = minimum
        self._max  = maximum

    def name(self):
        return self._name

    def min(self):
        return self._min

    def max(self):
        return self._max


class ImportanceSamplesTests(unittest.TestCase):

    _path = os.path.join(os.environ['SOURCE_DIR'], 'eos/data/importance_samples_TEST.d/samples')

    _parameters = [
        _Parameter('CKM::abs(V_ub)',        3.0e-3, 4.5e-3),
        _Parameter('B->pi::f_+(0)@BCL2008', 0.21,   0.32  ),
    ]

    @staticmethod
    def _write(directory, description, samples=None, weights=None, posterior_values=None):
        "Materialize a (possibly malformed) ImportanceSamples fixture in the given directory."
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, 'description.yaml'), 'w') as f:
            yaml.safe_dump(description, f, default_flow_style=False)
        if samples is not None:
            np.save(os.path.join(directory, 'samples.npy'), samples)
        if weights is not None:
            np.save(os.path.join(directory, 'weights.npy'), weights)
        if posterior_values is not None:
            np.save(os.path.join(directory, 'posterior_values.npy'), posterior_values)

    def test_loading(self):
        "Load an ImportanceSamples object from a prepared fixture and check its contents."
        f = eos.data.ImportanceSamples(self._path)

        self.assertEqual(f.type, 'ImportanceSamples')
        self.assertEqual(len(f.lookup_table), 6)
        self.assertEqual(f.lookup_table['CKM::abs(V_ub)'], 0)
        self.assertEqual(f.samples.shape[1], 6)
        self.assertEqual(f.weights.shape[0], f.samples.shape[0])
        self.assertIsNone(f.posterior_values)

    def test_varied_parameters_are_dicts(self):
        "The public varied_parameters remain subscriptable dicts for backward compatibility."
        f = eos.data.ImportanceSamples(self._path)
        self.assertEqual(f.varied_parameters[0]['name'], 'CKM::abs(V_ub)')
        self.assertIn('min', f.varied_parameters[0])
        self.assertIn('max', f.varied_parameters[0])

    def test_invalid_path(self):
        "Loading from a non-existent path raises an error."
        with self.assertRaises(RuntimeError):
            eos.data.ImportanceSamples(os.path.join(self._path, 'does-not-exist'))

    def test_roundtrip(self):
        "A data set written by create() reads back identically (write/read symmetry)."
        samples          = np.array([[3.5e-3, 0.22], [3.7e-3, 0.27], [4.1e-3, 0.30]])
        weights          = np.array([1.0, 2.0, 3.0])
        posterior_values = np.array([-1.0, -2.0, -3.0])

        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'samples')
            eos.data.ImportanceSamples.create(path, self._parameters, samples, weights, posterior_values)

            f = eos.data.ImportanceSamples(path)
            self.assertEqual(f.type, 'ImportanceSamples')
            self.assertEqual(f.lookup_table, {'CKM::abs(V_ub)': 0, 'B->pi::f_+(0)@BCL2008': 1})
            self.assertEqual(f.varied_parameters[0]['min'], 3.0e-3)
            np.testing.assert_allclose(f.samples,          samples)
            np.testing.assert_allclose(f.weights,          weights)
            np.testing.assert_allclose(f.posterior_values, posterior_values)

    def test_roundtrip_without_posterior_values(self):
        "A data set created without posterior_values reads back with the attribute set to None."
        samples = np.array([[3.5e-3, 0.22], [3.7e-3, 0.27]])
        weights = np.array([1.0, 2.0])

        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'samples')
            eos.data.ImportanceSamples.create(path, self._parameters, samples, weights)

            f = eos.data.ImportanceSamples(path)
            self.assertIsNone(f.posterior_values)

    def test_wrong_type(self):
        "A description whose type is not 'ImportanceSamples' is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'MarkovChain', 'parameters': []})
            with self.assertRaises(ValueError):
                eos.data.ImportanceSamples(d)

    def test_unknown_key(self):
        "An unexpected key in the description is rejected rather than silently ignored."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'ImportanceSamples',
                'parameters': [{'name': 'a', 'min': 0.0, 'max': 1.0}], 'bogus': 42,
            })
            with self.assertRaises(ValueError):
                eos.data.ImportanceSamples(d)

    def test_missing_parameters(self):
        "A description missing the required 'parameters' key is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'ImportanceSamples'})
            with self.assertRaises(ValueError):
                eos.data.ImportanceSamples(d)

    def test_non_mapping_description(self):
        "A description.yaml that is not a top-level mapping is rejected."
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'description.yaml'), 'w') as f:
                yaml.safe_dump(['not', 'a', 'mapping'], f)
            with self.assertRaises(RuntimeError):
                eos.data.ImportanceSamples(d)

    def test_shape_mismatch(self):
        "Samples whose column count disagrees with the parameters are rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'ImportanceSamples',
                'parameters': [{'name': 'a', 'min': 0.0, 'max': 1.0}, {'name': 'b', 'min': 0.0, 'max': 1.0}],
            }, samples=np.zeros((5, 3)), weights=np.zeros(5))
            with self.assertRaises(RuntimeError):
                eos.data.ImportanceSamples(d)


if __name__ == '__main__':
    unittest.main(verbosity=5)
