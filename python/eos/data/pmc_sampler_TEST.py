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


class PMCSamplerTests(unittest.TestCase):

    _parameters = [
        _Parameter('CKM::abs(V_ub)',        3.0e-3, 4.5e-3),
        _Parameter('B->pi::f_+(0)@BCL2008', 0.21,   0.32  ),
    ]

    @staticmethod
    def _write(directory, description):
        "Materialize a (possibly malformed) PMCSampler description in the given directory."
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, 'description.yaml'), 'w') as f:
            yaml.safe_dump(description, f, default_flow_style=False)

    @staticmethod
    def _proposal():
        import pypmc
        means   = [np.array([5.0, 0.01]), np.array([-4.0, 1.0])]
        covs    = [np.array([[0.01, 0.003], [0.003, 0.0025]]), np.array([[0.1, 0.0], [0.0, 0.02]])]
        weights = np.array([0.3, 0.7])
        return pypmc.density.mixture.create_gaussian_mixture(means, covs, weights)

    def test_evaluate_mixture_pdf(self):
        "Test the evaluation of a mixture PDF, used in the computation of test statistics."

        try:
            import pypmc
        except ImportError:
            raise unittest.SkipTest("skipping 'test_evaluate_mixture_pdf' - pypmc is not installed")

        component_weights = np.array([0.3, 0.7])

        mean0       = np.array ([ 5.0  , 0.01  ])
        covariance0 = np.array([[ 0.01 , 0.003 ],
                                [ 0.003, 0.0025]])

        mean1       = np.array ([-4.0  , 1.0   ])
        covariance1 = np.array([[ 0.1  , 0.    ],
                                [ 0.   , 0.02  ]])

        component_means = [mean0, mean1]
        component_covariances = [covariance0, covariance1]

        target_mixture = pypmc.density.mixture.create_gaussian_mixture(component_means, component_covariances, component_weights)

        self.assertAlmostEqual(
            eos.data.PMCSampler._evaluate_mixture_pdf(target_mixture, np.array([-3.5, 0.8])),
            0.26256727,
            delta=1e-5
        )

    @unittest.skipIf(pypmc_is_missing(), "Test requires the 'pypmc' module")
    def test_roundtrip(self):
        "A PMC sampler written by create() reads back identically and writes the modern key."
        proposal = self._proposal()

        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'pmc')
            eos.data.PMCSampler.create(path, self._parameters, proposal)

            # the modern key is written; the legacy misspelling is not
            with open(os.path.join(path, 'description.yaml')) as f:
                on_disk = yaml.safe_load(f)
            self.assertIn('test_statistics', on_disk)
            self.assertNotIn('test statistics', on_disk)

            s = eos.data.PMCSampler(path)
            self.assertEqual(s.type, 'PMCSampler')
            self.assertEqual(s.lookup_table, {'CKM::abs(V_ub)': 0, 'B->pi::f_+(0)@BCL2008': 1})
            self.assertEqual(len(s.components), 2)
            np.testing.assert_allclose(s.component_weights, [0.3, 0.7])
            self.assertEqual(s.test_statistics, {'sigma': [], 'densities': []})

            density = s.density()
            self.assertEqual(len(density.components), 2)

    @unittest.skipIf(pypmc_is_missing(), "Test requires the 'pypmc' module")
    def test_legacy_test_statistics_warns(self):
        "A legacy 'test statistics' (space) key is accepted, warns, and is normalized."
        proposal = self._proposal()

        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'pmc')
            eos.data.PMCSampler.create(path, self._parameters, proposal)

            # rewrite the description with the legacy misspelled key
            desc = os.path.join(path, 'description.yaml')
            with open(desc) as f:
                on_disk = yaml.safe_load(f)
            on_disk['test statistics'] = on_disk.pop('test_statistics')
            with open(desc, 'w') as f:
                yaml.safe_dump(on_disk, f)

            with self.assertLogs('EOS', level='WARNING'):
                s = eos.data.PMCSampler(path)
            self.assertEqual(s.test_statistics, {'sigma': [], 'densities': []})

    def test_invalid_path(self):
        "Loading from a non-existent path raises an error."
        with self.assertRaises(RuntimeError):
            eos.data.PMCSampler(os.path.join(tempfile.gettempdir(), 'eos-does-not-exist-pmc'))

    def test_wrong_type(self):
        "A description whose type is not 'PMCSampler' is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'MixtureDensity',
                            'parameters': [], 'proposal': {'components': [], 'weights': []}})
            with self.assertRaises(ValueError):
                eos.data.PMCSampler(d)

    def test_unknown_key(self):
        "An unexpected key in the description is rejected rather than silently ignored."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'PMCSampler',
                            'parameters': [], 'proposal': {'components': [], 'weights': []}, 'bogus': 42})
            with self.assertRaises(ValueError):
                eos.data.PMCSampler(d)

    def test_missing_proposal(self):
        "A description missing the required 'proposal' key is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'PMCSampler', 'parameters': []})
            with self.assertRaises(ValueError):
                eos.data.PMCSampler(d)

    def test_non_mapping_description(self):
        "A description.yaml that is not a top-level mapping is rejected."
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'description.yaml'), 'w') as f:
                yaml.safe_dump(['not', 'a', 'mapping'], f)
            with self.assertRaises(RuntimeError):
                eos.data.PMCSampler(d)


if __name__ == '__main__':
    unittest.main(verbosity=5)
