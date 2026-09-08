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

from eos.data.unbinned_likelihood import (
    UnbinnedAxisDescription,
    UnbinnedExpressionResolutionDescription,
    UnbinnedOffsetAxisDescription,
    UnbinnedPDFResolutionDescription,
    UnbinnedSampledResolutionDescription,
)


class UnbinnedLikelihoodTests(unittest.TestCase):

    _path = os.path.join(os.environ['SOURCE_DIR'], 'eos/data/unbinned_likelihood_TEST.d')

    @staticmethod
    def _write(directory, description, observations=None, resolution=None):
        "Materialize a (possibly malformed) UnbinnedLikelihood fixture in the given directory."
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, 'description.yaml'), 'w') as f:
            yaml.safe_dump(description, f, default_flow_style=False)
        if observations is not None:
            np.save(os.path.join(directory, 'observations.npy'), observations)
        if resolution is not None:
            np.save(os.path.join(directory, 'resolution.npy'), resolution)

    def test_load_fixture(self):
        "The checked-in 1D fixture (a sampled resolution, declared wider than required) loads correctly."
        ul = eos.data.UnbinnedLikelihood(self._path)

        self.assertEqual(ul.type, 'UnbinnedLikelihood')
        self.assertEqual(ul.rank, 1)
        self.assertEqual(len(ul.axes), 1)
        self.assertEqual(ul.axes[0].variable, 'q2')
        self.assertEqual(ul.axes[0].min, 0.0)
        self.assertEqual(ul.axes[0].max, 7.0)
        self.assertEqual(ul.axes[0].points, 8)
        self.assertIsInstance(ul.resolution, UnbinnedSampledResolutionDescription)
        self.assertEqual(ul.resolution.offsets[0].points, 10)
        self.assertEqual(ul.observations.shape, (3, 1))
        self.assertEqual(list(ul.observations[:, 0]), [1.0, 3.5, 6.0])
        self.assertEqual(ul.resolution_values.shape, (10,))

    def test_invalid_path(self):
        "Loading from a non-existent path raises an error."
        with self.assertRaises(RuntimeError):
            eos.data.UnbinnedLikelihood(os.path.join(self._path, 'does-not-exist'))

    def test_wrong_type(self):
        "A description whose 'type' is not 'UnbinnedLikelihood' is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': '1.0',
                'type':    'SomeOtherType',
                'axes':    [{'variable': 'q2', 'min': 0.0, 'max': 7.0, 'points': 8}],
                'resolution': {'kind': 'expression', 'expression': '1.0'},
            })
            with self.assertRaises(ValueError):
                eos.data.UnbinnedLikelihood(d)

    def test_malformed_axes(self):
        "An axis description missing a mandatory key (here 'points') is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': '1.0',
                'type':    'UnbinnedLikelihood',
                'axes':    [{'variable': 'q2', 'min': 0.0, 'max': 7.0}],
                'resolution': {'kind': 'expression', 'expression': '1.0'},
            }, observations=np.zeros((1, 1)))
            with self.assertRaises(ValueError):
                eos.data.UnbinnedLikelihood(d)

    def test_roundtrip_1d_samples(self):
        "A 1D data object with a sampled resolution, written by create(), reads back identically."
        with tempfile.TemporaryDirectory() as d:
            axes = [UnbinnedAxisDescription(variable='q2', min=0.0, max=3.0, points=4)]
            resolution = UnbinnedSampledResolutionDescription(
                offsets=[UnbinnedOffsetAxisDescription(variable='q2', min=-2.0, spacing=1.0, points=4)],
            )
            observations = np.array([[0.5], [1.5], [2.5]])
            resolution_values = np.array([0.1, 0.2, 0.3, 0.4])

            eos.data.UnbinnedLikelihood.create(d, axes, resolution, observations, resolution_values=resolution_values)
            ul = eos.data.UnbinnedLikelihood(d)

            self.assertEqual(ul.rank, 1)
            self.assertEqual(ul.axes[0].points, 4)
            np.testing.assert_array_equal(ul.observations, observations)
            np.testing.assert_array_equal(ul.resolution_values, resolution_values)

    def test_roundtrip_2d_expression(self):
        "A 2D data object with an evaluatable resolution, written by create(), reads back identically."
        with tempfile.TemporaryDirectory() as d:
            axes = [
                UnbinnedAxisDescription(variable='q2', min=0.0, max=7.0, points=8),
                UnbinnedAxisDescription(variable='cos_theta', min=-1.0, max=1.0, points=4),
            ]
            resolution = UnbinnedExpressionResolutionDescription(expression='1.0')
            observations = np.array([[1.0, -0.5], [5.0, 0.5]])

            eos.data.UnbinnedLikelihood.create(d, axes, resolution, observations)
            ul = eos.data.UnbinnedLikelihood(d)

            self.assertEqual(ul.rank, 2)
            self.assertIsInstance(ul.resolution, UnbinnedExpressionResolutionDescription)
            self.assertIsNone(ul.resolution_values)
            np.testing.assert_array_equal(ul.observations, observations)

    def test_create_rejects_missing_resolution_values(self):
        "create() refuses a sampled resolution without resolution_values."
        with tempfile.TemporaryDirectory() as d:
            axes = [UnbinnedAxisDescription(variable='q2', min=0.0, max=3.0, points=4)]
            resolution = UnbinnedSampledResolutionDescription(
                offsets=[UnbinnedOffsetAxisDescription(variable='q2', min=-2.0, spacing=1.0, points=4)],
            )
            with self.assertRaises(RuntimeError):
                eos.data.UnbinnedLikelihood.create(d, axes, resolution, np.array([[0.5]]))

    def test_cropped_axes_native_range(self):
        "Omitting a crop uses the full native range."
        ul = eos.data.UnbinnedLikelihood(self._path)
        cropped = ul.cropped_axes()
        self.assertEqual(cropped, [{'variable': 'q2', 'min': 0.0, 'max': 7.0, 'points': 8}])

    def test_cropped_axes_snaps_and_crops(self):
        "A crop snaps onto native nodes and yields an even, reduced point count."
        ul = eos.data.UnbinnedLikelihood(self._path)
        cropped = ul.cropped_axes({'q2': (2.0, 5.0)})
        self.assertEqual(cropped, [{'variable': 'q2', 'min': 2.0, 'max': 5.0, 'points': 4}])

    def test_cropped_axes_rejects_odd_count(self):
        "A crop that yields an odd point count is rejected, not silently snapped."
        ul = eos.data.UnbinnedLikelihood(self._path)
        with self.assertRaises(ValueError):
            ul.cropped_axes({'q2': (2.0, 4.0)})

    def test_resolution_kernel_slices_declared_grid(self):
        "The declared (wider) resolution grid is sliced onto the required, narrower offset grid."
        ul = eos.data.UnbinnedLikelihood(self._path)
        cropped = ul.cropped_axes()
        kernel = ul.resolution_kernel(cropped)
        # resolution.npy holds arange(10); the required offsets [-4 .. 3] sit at declared indices [1 .. 8]
        np.testing.assert_array_equal(kernel, np.arange(1, 9, dtype=float))

    def test_resolution_kernel_slices_cropped_grid(self):
        "Cropping the signal axis crops the required offset window accordingly."
        ul = eos.data.UnbinnedLikelihood(self._path)
        cropped = ul.cropped_axes({'q2': (2.0, 5.0)})
        kernel = ul.resolution_kernel(cropped)
        # points=4, required offsets [-2 .. 1] sit at declared indices [3 .. 6]
        np.testing.assert_array_equal(kernel, np.arange(3, 7, dtype=float))

    def test_resolution_kernel_rejects_misaligned_offsets(self):
        "A declared offset grid that does not align with the required grid raises, naming the offset."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': '1.0',
                'type':    'UnbinnedLikelihood',
                'axes':    [{'variable': 'q2', 'min': 0.0, 'max': 7.0, 'points': 8}],
                'resolution': {
                    'kind':    'samples',
                    'offsets': [{'variable': 'q2', 'min': -4.3, 'spacing': 1.0, 'points': 8}],
                },
            }, observations=np.zeros((1, 1)), resolution=np.arange(8, dtype=float))
            ul = eos.data.UnbinnedLikelihood(d)
            cropped = ul.cropped_axes()
            with self.assertRaisesRegex(RuntimeError, '-4'):
                ul.resolution_kernel(cropped)

    def test_resolution_kernel_rejects_non_sampled_resolution(self):
        "resolution_kernel() refuses an evaluatable resolution."
        with tempfile.TemporaryDirectory() as d:
            axes = [UnbinnedAxisDescription(variable='q2', min=0.0, max=3.0, points=4)]
            resolution = UnbinnedExpressionResolutionDescription(expression='1.0')
            eos.data.UnbinnedLikelihood.create(d, axes, resolution, np.array([[0.5]]))
            ul = eos.data.UnbinnedLikelihood(d)
            with self.assertRaises(RuntimeError):
                ul.resolution_kernel(ul.cropped_axes())

    def test_sampled_resolution_kernel_expression(self):
        "An evaluatable (expression) resolution is sampled once onto the required grid."
        with tempfile.TemporaryDirectory() as d:
            axes = [
                UnbinnedAxisDescription(variable='q2', min=0.0, max=7.0, points=8),
                UnbinnedAxisDescription(variable='cos_theta', min=-1.0, max=1.0, points=4),
            ]
            resolution = UnbinnedExpressionResolutionDescription(expression='1.0')
            eos.data.UnbinnedLikelihood.create(d, axes, resolution, np.array([[1.0, 0.0]]))
            ul = eos.data.UnbinnedLikelihood(d)

            cropped = ul.cropped_axes()
            kernel = ul.sampled_resolution_kernel(cropped)
            self.assertEqual(kernel.shape, (8, 4))
            np.testing.assert_allclose(kernel, np.ones((8, 4)))

    def test_sampled_resolution_kernel_rejects_sampled_resolution(self):
        "sampled_resolution_kernel() refuses a sampled resolution."
        ul = eos.data.UnbinnedLikelihood(self._path)
        with self.assertRaises(RuntimeError):
            ul.sampled_resolution_kernel(ul.cropped_axes())


if __name__ == '__main__':
    unittest.main()
