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


class SampleMaskTests(unittest.TestCase):

    _path = os.path.join(os.environ['SOURCE_DIR'], 'eos/data/sample_mask_TEST.d')

    _observables = ['B->Dlnu::BR', 'B->pilnu::BR', 'B->Dlnu::dBR/dq2', 'B->pilnu::dBR/dq2', 'B_u->lnu::BR']

    @staticmethod
    def _write(directory, description, mask=None):
        "Materialize a (possibly malformed) SampleMask fixture in the given directory."
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, 'description.yaml'), 'w') as f:
            yaml.safe_dump(description, f, default_flow_style=False)
        if mask is not None:
            np.save(os.path.join(directory, 'mask.npy'), mask)

    def test_loading_legacy(self):
        "Load a legacy (type 'Mask') fixture: it is accepted, warns, and normalizes the type."
        with self.assertLogs('EOS', level='WARNING'):
            sm = eos.data.SampleMask(self._path)

        self.assertEqual(sm.type, 'SampleMask')
        self.assertEqual(sm.observables, self._observables)
        self.assertEqual(sm.mask.dtype, bool)
        self.assertEqual(sm.mask.shape, (5,))
        self.assertEqual(list(sm.mask), [True, False, True, True, False])

    def test_invalid_path(self):
        "Loading from a non-existent path raises an error."
        with self.assertRaises(RuntimeError):
            eos.data.SampleMask(os.path.join(self._path, 'does-not-exist'))

    def test_roundtrip(self):
        "A mask written by create() is stored as 'SampleMask' and reads back identically."
        mask = np.array([True, False, True, True, False])

        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'mask-default')
            eos.data.SampleMask.create(path, mask, self._observables)

            # the discriminator is written as the current 'SampleMask'
            with open(os.path.join(path, 'description.yaml')) as f:
                self.assertEqual(yaml.safe_load(f)['type'], 'SampleMask')

            sm = eos.data.SampleMask(path)
            self.assertEqual(sm.type, 'SampleMask')
            self.assertEqual(sm.observables, self._observables)
            np.testing.assert_array_equal(sm.mask, mask)

    def test_wrong_type(self):
        "A description whose type is neither 'SampleMask' nor 'Mask' is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'Prediction', 'observables': []}, mask=np.array([True]))
            with self.assertRaises(ValueError):
                eos.data.SampleMask(d)

    def test_unknown_key(self):
        "An unexpected key in the description is rejected rather than silently ignored."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'SampleMask', 'observables': [], 'bogus': 42},
                        mask=np.array([True]))
            with self.assertRaises(ValueError):
                eos.data.SampleMask(d)

    def test_missing_observables(self):
        "A description missing the required 'observables' key is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'SampleMask'}, mask=np.array([True]))
            with self.assertRaises(ValueError):
                eos.data.SampleMask(d)

    def test_missing_payload(self):
        "A valid description without the mask payload is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'SampleMask', 'observables': self._observables})
            with self.assertRaises(RuntimeError):
                eos.data.SampleMask(d)

    def test_non_mapping_description(self):
        "A description.yaml that is not a top-level mapping is rejected."
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'description.yaml'), 'w') as f:
                yaml.safe_dump(['not', 'a', 'mapping'], f)
            with self.assertRaises(RuntimeError):
                eos.data.SampleMask(d)


if __name__ == '__main__':
    unittest.main(verbosity=5)
