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


def dynesty_is_missing():
    try:
        import dynesty
        return False
    except ImportError:
        return True


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


class DynestyResultsTests(unittest.TestCase):

    _path = os.path.join(os.environ['SOURCE_DIR'], 'eos/data/dynesty_results_TEST.d')

    _parameters = [
        _Parameter('CKM::abs(V_ub)',        3.0e-3, 4.5e-3),
        _Parameter('B->pi::f_+(0)@BCL2008', 0.21,   0.32  ),
    ]

    @staticmethod
    def _write(directory, description):
        "Materialize a (possibly malformed) DynestyResults description in the given directory."
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, 'description.yaml'), 'w') as f:
            yaml.safe_dump(description, f, default_flow_style=False)

    @unittest.skipIf(dynesty_is_missing(), "Test requires the 'dynesty' module")
    def test_loading(self):
        "Load a DynestyResults object from a prepared fixture and check its contents."
        dr = eos.data.DynestyResults(self._path)

        self.assertEqual(dr.type, 'DynestyResults')
        self.assertEqual(dr.lookup_table, {'CKM::abs(V_ub)': 0, 'B->pi::f_+(0)@BCL2008': 1})
        self.assertEqual(dr.varied_parameters[0]['name'], 'CKM::abs(V_ub)')
        # the reconstructed samples have one column per varied parameter
        self.assertEqual(dr.samples.shape[1], 2)
        # there is exactly one weight per sample
        self.assertEqual(len(dr.weights), dr.samples.shape[0])

    def test_invalid_path(self):
        "Loading from a non-existent path raises an error."
        with self.assertRaises(RuntimeError):
            eos.data.DynestyResults(os.path.join(self._path, 'does-not-exist'))

    @unittest.skipIf(dynesty_is_missing(), "Test requires the 'dynesty' module")
    def test_roundtrip(self):
        "Results written by create() read back identically (write/read symmetry)."
        src = eos.data.DynestyResults(self._path)

        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'nested')
            eos.data.DynestyResults.create(path, self._parameters, src.results)

            dr = eos.data.DynestyResults(path)
            self.assertEqual(dr.type, 'DynestyResults')
            self.assertEqual(dr.lookup_table, {'CKM::abs(V_ub)': 0, 'B->pi::f_+(0)@BCL2008': 1})
            np.testing.assert_allclose(dr.samples, src.samples)
            np.testing.assert_allclose(dr.weights, src.weights)

    def test_wrong_type(self):
        "A description whose type is not 'DynestyResults' is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'MarkovChain', 'parameters': []})
            with self.assertRaises(ValueError):
                eos.data.DynestyResults(d)

    def test_unknown_key(self):
        "An unexpected key in the description is rejected rather than silently ignored."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'DynestyResults',
                'parameters': [{'name': 'a', 'min': 0.0, 'max': 1.0}], 'bogus': 42,
            })
            with self.assertRaises(ValueError):
                eos.data.DynestyResults(d)

    def test_missing_parameters(self):
        "A description missing the required 'parameters' key is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'DynestyResults'})
            with self.assertRaises(ValueError):
                eos.data.DynestyResults(d)

    def test_missing_payload(self):
        "A valid description without the results payload is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'DynestyResults',
                'parameters': [{'name': 'a', 'min': 0.0, 'max': 1.0}],
            })
            with self.assertRaises(RuntimeError):
                eos.data.DynestyResults(d)

    def test_non_mapping_description(self):
        "A description.yaml that is not a top-level mapping is rejected."
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'description.yaml'), 'w') as f:
                yaml.safe_dump(['not', 'a', 'mapping'], f)
            with self.assertRaises(RuntimeError):
                eos.data.DynestyResults(d)


if __name__ == '__main__':
    unittest.main(verbosity=5)
