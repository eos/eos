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


class PredictionTests(unittest.TestCase):

    _path = os.path.join(os.environ['SOURCE_DIR'], 'eos/data/prediction_TEST.d')

    @staticmethod
    def _write(directory, description, samples=None, weights=None):
        "Materialize a (possibly malformed) Prediction fixture in the given directory."
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, 'description.yaml'), 'w') as f:
            yaml.safe_dump(description, f, default_flow_style=False)
        if samples is not None:
            np.save(os.path.join(directory, 'samples.npy'), samples)
        if weights is not None:
            np.save(os.path.join(directory, 'weights.npy'), weights)

    def test_loading(self):
        "Load the prepared (legacy) fixtures and check their contents."
        for name in ('predictions', 'predictions-binned'):
            p = eos.data.Prediction(os.path.join(self._path, name))
            self.assertEqual(p.type, 'Prediction')
            self.assertEqual(p.format, 1)
            self.assertEqual(p.samples.shape[0], p.weights.shape[0])
            self.assertEqual(p.samples.shape[1], len(p.varied_parameters))
            # the fixture columns are all genuine observables (they carry options and kinematics)
            for vp in p.varied_parameters:
                self.assertEqual(vp['kind'], 'observable')
            # the lookup keys keep the qualified-name form 'name;options[kinematics]'
            first = next(iter(p.lookup_table))
            self.assertTrue(first.startswith('B->pilnu::'))
            self.assertIn(';', first)
            self.assertTrue(first.endswith(']'))

    def test_invalid_path(self):
        "Loading from a non-existent path raises an error."
        with self.assertRaises(RuntimeError):
            eos.data.Prediction(os.path.join(self._path, 'does-not-exist'))

    def test_roundtrip(self):
        "A prediction of a genuine observable and a parameter-clothed column round-trips with kinds."
        params = eos.Parameters.Defaults()
        obs = eos.Observable.make(
            'B->pilnu::dBR/dq2', params,
            eos.Kinematics({'q2': 1.0}),
            eos.Options({'P': 'pi', 'form-factors': 'BCL2008', 'model': 'CKM'}),
        )
        par = eos.Observable.make('CKM::abs(V_ub)', params, eos.Kinematics({}), eos.Options({}))
        self.assertIsNotNone(obs)
        self.assertIsNotNone(par)

        samples = np.array([[1.0e-6, 3.7e-3], [2.0e-6, 3.6e-3], [3.0e-6, 3.8e-3]])
        weights = np.array([1.0, 2.0, 3.0])

        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'pred-test')
            eos.data.Prediction.create(path, [obs, par], samples, weights)

            p = eos.data.Prediction(path)
            self.assertEqual(p.type, 'Prediction')
            self.assertEqual(p.format, 2)

            # the genuine observable is classified as such, the parameter as a parameter
            self.assertEqual(p.varied_parameters[0]['kind'], 'observable')
            self.assertEqual(p.varied_parameters[1]['kind'], 'parameter')
            self.assertEqual(p.varied_parameters[1]['name'], 'CKM::abs(V_ub)')
            self.assertEqual(p.varied_parameters[1]['options'], {})
            self.assertEqual(p.varied_parameters[1]['kinematics'], {})

            # lookup keys: the parameter column has empty options/kinematics parts
            self.assertEqual(p.lookup_table['CKM::abs(V_ub);[]'], 1)
            obs_keys = [k for k in p.lookup_table if k != 'CKM::abs(V_ub);[]']
            self.assertEqual(len(obs_keys), 1)
            self.assertEqual(p.lookup_table[obs_keys[0]], 0)
            self.assertTrue(obs_keys[0].startswith('B->pilnu::dBR/dq2;'))

            np.testing.assert_allclose(p.samples, samples)
            np.testing.assert_allclose(p.weights, weights)

    def test_legacy_reclassification(self):
        "A legacy (format-1) file without kinds is reclassified on read."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'Prediction',
                'observables': [
                    {'name': 'B->pilnu::dBR/dq2', 'kinematics': {'q2': 1.0},
                     'options': {'P': 'pi', 'form-factors': 'BCL2008', 'model': 'CKM'}},
                    {'name': 'CKM::abs(V_ub)', 'kinematics': {}, 'options': {}},
                ],
            }, samples=np.zeros((3, 2)), weights=np.zeros(3))

            p = eos.data.Prediction(d)
            self.assertEqual(p.format, 1)
            self.assertEqual(p.varied_parameters[0]['kind'], 'observable')
            self.assertEqual(p.varied_parameters[1]['kind'], 'parameter')

    def test_wrong_type(self):
        "A description whose type is not 'Prediction' is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'MarkovChain', 'observables': []})
            with self.assertRaises(ValueError):
                eos.data.Prediction(d)

    def test_unknown_key(self):
        "An unexpected key in the description is rejected rather than silently ignored."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'Prediction', 'format': 2,
                'observables': [], 'bogus': 42,
            })
            with self.assertRaises(ValueError):
                eos.data.Prediction(d)

    def test_missing_observables(self):
        "A description missing the required 'observables' key is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'Prediction'})
            with self.assertRaises(ValueError):
                eos.data.Prediction(d)

    def test_invalid_kind(self):
        "A column with an unrecognized kind is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'Prediction', 'format': 2,
                'observables': [{'name': 'a', 'kind': 'nonsense', 'kinematics': {}, 'options': {}}],
            })
            with self.assertRaises(ValueError):
                eos.data.Prediction(d)

    def test_shape_mismatch(self):
        "Samples whose column count disagrees with the observables are rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'Prediction', 'format': 2,
                'observables': [
                    {'name': 'a', 'kind': 'observable', 'kinematics': {}, 'options': {}},
                    {'name': 'b', 'kind': 'observable', 'kinematics': {}, 'options': {}},
                ],
            }, samples=np.zeros((5, 3)), weights=np.zeros(5))
            with self.assertRaises(RuntimeError):
                eos.data.Prediction(d)

    def test_non_mapping_description(self):
        "A description.yaml that is not a top-level mapping is rejected."
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'description.yaml'), 'w') as f:
                yaml.safe_dump(['not', 'a', 'mapping'], f)
            with self.assertRaises(RuntimeError):
                eos.data.Prediction(d)


if __name__ == '__main__':
    unittest.main(verbosity=5)
