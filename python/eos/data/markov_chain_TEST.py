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


class MarkovChainTests(unittest.TestCase):

    _path = os.path.join(os.environ['SOURCE_DIR'], 'eos/data/markov_chain_TEST.d')

    _parameters = [
        _Parameter('CKM::abs(V_ub)',        3.0e-3, 4.5e-3),
        _Parameter('B->pi::f_+(0)@BCL2008', 0.21,   0.32  ),
    ]

    @staticmethod
    def _write(directory, description, samples=None, usamples=None, weights=None):
        "Materialize a (possibly malformed) MarkovChain fixture in the given directory."
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, 'description.yaml'), 'w') as f:
            yaml.safe_dump(description, f, default_flow_style=False)
        if samples is not None:
            np.save(os.path.join(directory, 'samples.npy'), samples)
        if usamples is not None:
            np.save(os.path.join(directory, 'usamples.npy'), usamples)
        if weights is not None:
            np.save(os.path.join(directory, 'weights.npy'), weights)

    def test_loading(self):
        "Load a MarkovChain from a prepared fixture and check its contents."
        mc = eos.data.MarkovChain(self._path)

        self.assertEqual(mc.type, 'MarkovChain')
        self.assertEqual(mc.lookup_table, {'CKM::abs(V_ub)': 0, 'B->pi::f_+(0)@BCL2008': 1})
        self.assertEqual(mc.samples.shape, (5, 2))
        self.assertEqual(mc.usamples.shape, (5, 2))
        self.assertIsNotNone(mc.weights)
        self.assertEqual(mc.weights.shape, (5,))

    def test_varied_parameters_are_dicts(self):
        "The public varied_parameters remain subscriptable dicts for backward compatibility."
        mc = eos.data.MarkovChain(self._path)
        self.assertEqual(mc.varied_parameters[0]['name'], 'CKM::abs(V_ub)')
        self.assertIn('min', mc.varied_parameters[0])
        self.assertIn('max', mc.varied_parameters[0])

    def test_invalid_path(self):
        "Loading from a non-existent path raises an error."
        with self.assertRaises(RuntimeError):
            eos.data.MarkovChain(os.path.join(self._path, 'does-not-exist'))

    def test_roundtrip(self):
        "A chain written by create() reads back identically (write/read symmetry)."
        samples  = np.array([[3.5e-3, 0.22], [3.7e-3, 0.27], [4.1e-3, 0.30]])
        usamples = np.array([[0.10, 0.20], [0.30, 0.40], [0.50, 0.60]])
        weights  = np.array([1.0, 2.0, 3.0])

        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'mcmc-0000')
            eos.data.MarkovChain.create(path, self._parameters, samples, usamples, weights)

            mc = eos.data.MarkovChain(path)
            self.assertEqual(mc.type, 'MarkovChain')
            self.assertEqual(mc.lookup_table, {'CKM::abs(V_ub)': 0, 'B->pi::f_+(0)@BCL2008': 1})
            self.assertEqual(mc.varied_parameters[0]['min'], 3.0e-3)
            self.assertEqual(mc.varied_parameters[1]['max'], 0.32)
            np.testing.assert_allclose(mc.samples,  samples)
            np.testing.assert_allclose(mc.usamples, usamples)
            np.testing.assert_allclose(mc.weights,  weights)

    def test_roundtrip_unweighted(self):
        "A chain created without weights reads back with weights set to None."
        samples  = np.array([[3.5e-3, 0.22], [3.7e-3, 0.27]])
        usamples = np.array([[0.10, 0.20], [0.30, 0.40]])

        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'mcmc-0000')
            eos.data.MarkovChain.create(path, self._parameters, samples, usamples)

            mc = eos.data.MarkovChain(path)
            self.assertIsNone(mc.weights)

    def test_wrong_type(self):
        "A description whose type is not 'MarkovChain' is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'Prediction', 'parameters': [], 'has-weights': False})
            with self.assertRaises(ValueError):
                eos.data.MarkovChain(d)

    def test_unknown_key(self):
        "An unexpected key in the description is rejected rather than silently ignored."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'MarkovChain',
                'parameters': [{'name': 'a', 'min': 0.0, 'max': 1.0}],
                'has-weights': False, 'bogus': 42,
            })
            with self.assertRaises(ValueError):
                eos.data.MarkovChain(d)

    def test_missing_parameters(self):
        "A description missing the required 'parameters' key is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'MarkovChain', 'has-weights': False})
            with self.assertRaises(ValueError):
                eos.data.MarkovChain(d)

    def test_non_mapping_description(self):
        "A description.yaml that is not a top-level mapping is rejected."
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'description.yaml'), 'w') as f:
                yaml.safe_dump(['not', 'a', 'mapping'], f)
            with self.assertRaises(RuntimeError):
                eos.data.MarkovChain(d)

    def test_shape_mismatch(self):
        "Samples whose column count disagrees with the parameters are rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'MarkovChain',
                'parameters': [{'name': 'a', 'min': 0.0, 'max': 1.0}, {'name': 'b', 'min': 0.0, 'max': 1.0}],
                'has-weights': False,
            }, samples=np.zeros((5, 3)), usamples=np.zeros((5, 2)))
            with self.assertRaises(RuntimeError):
                eos.data.MarkovChain(d)


if __name__ == '__main__':
    unittest.main(verbosity=5)
