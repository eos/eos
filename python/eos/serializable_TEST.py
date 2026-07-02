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

import numpy as np
import os
import tempfile
import yaml

from dataclasses import dataclass
from eos.deserializable import Deserializable
from eos.serializable import Serializable, _SimpleDumper, _EXCLUDED_MODULES


@dataclass
class _Holder(Serializable, Deserializable):
    "Minimal dataclass used to exercise Serializable.to_yaml_file."
    value:object


class SerializableTests(unittest.TestCase):

    @staticmethod
    def _dump(obj):
        return yaml.dump(obj, Dumper=_SimpleDumper, default_flow_style=False)

    def test_numpy_is_excluded(self):
        "numpy is on the list of modules that must not be serialized."
        self.assertIn('numpy', _EXCLUDED_MODULES)

    def test_native_types_serialize(self):
        "Native Python structures serialize (and round-trip) unchanged."
        obj = {'a': [1.0, 2, 'x'], 'b': None, 'c': {'d': True}}
        self.assertEqual(yaml.safe_load(self._dump(obj)), obj)

    def test_rejects_numpy_array(self):
        "A numpy array is rejected with a TypeError."
        with self.assertRaises(TypeError):
            self._dump({'x': np.array([1.0, 2.0])})

    def test_rejects_numpy_scalars(self):
        "numpy scalar types are rejected, even those subclassing Python builtins."
        for scalar in (np.float64(1.5), np.int64(3), np.bool_(True), np.str_('abc')):
            with self.assertRaises(TypeError):
                self._dump({'x': scalar})

    def test_rejects_numpy_when_nested(self):
        "A numpy value nested at any depth is rejected."
        with self.assertRaises(TypeError):
            self._dump({'a': {'b': [1.0, np.float64(2.0)]}})

    def test_to_yaml_file_native(self):
        "Serializable.to_yaml_file writes a description of native types."
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'x.yaml')
            _Holder(value=[1.0, 2.0]).to_yaml_file(path)
            with open(path) as f:
                self.assertEqual(yaml.safe_load(f), {'value': [1.0, 2.0]})

    def test_to_yaml_file_rejects_numpy(self):
        "Serializable.to_yaml_file refuses to serialize a numpy value."
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'x.yaml')
            with self.assertRaises(TypeError):
                _Holder(value=np.arange(3)).to_yaml_file(path)


if __name__ == '__main__':
    unittest.main(verbosity=5)
