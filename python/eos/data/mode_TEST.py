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
import scipy.stats
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


class _Entry:
    "Minimal stand-in for a GoodnessOfFit test-statistic entry."

    def __init__(self, chi2, dof):
        self.chi2 = chi2
        self.dof  = dof


class _GoodnessOfFit:
    "Minimal stand-in for eos.GoodnessOfFit, iterating as (name, entry) pairs."

    def __init__(self, total_chi2, total_dof, entries):
        self._total_chi2 = total_chi2
        self._total_dof  = total_dof
        self._entries    = entries

    def total_chi_square(self):
        return self._total_chi2

    def total_degrees_of_freedom(self):
        return self._total_dof

    def __iter__(self):
        return iter(self._entries)


class _Analysis:
    def __init__(self, varied_parameters):
        self.varied_parameters = varied_parameters


class _BestFitPoint:
    def __init__(self, analysis, point):
        self.analysis = analysis
        self.point    = point


class ModeTests(unittest.TestCase):

    _path = os.path.join(os.environ['SOURCE_DIR'], 'eos/data/mode_TEST.d')

    _parameters = [
        _Parameter('CKM::abs(V_ub)',        3.0e-3, 4.5e-3),
        _Parameter('B->pi::f_+(0)@BCL2008', 0.21,   0.32  ),
    ]

    @staticmethod
    def _write(directory, description):
        "Materialize a (possibly malformed) Mode fixture in the given directory."
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, 'description.yaml'), 'w') as f:
            yaml.safe_dump(description, f, default_flow_style=False)

    def test_loading_format1(self):
        "Load a legacy (format-v1) Mode fixture and check that it is upgraded in memory."
        m = eos.data.Mode(self._path)

        self.assertEqual(m.type, 'Mode')
        self.assertEqual(m.format_version, 1)
        self.assertEqual([p['name'] for p in m.varied_parameters], ['CKM::abs(V_ub)', 'B->pi::f_+(0)@BCL2008'])
        self.assertEqual(m.mode, [3.67e-3, 0.27])
        self.assertAlmostEqual(m.pvalue, 0.25)
        self.assertAlmostEqual(m.global_chi2, 12.0)
        self.assertAlmostEqual(m.dof, 5.0)
        # the flat local_pvalues map is upgraded into test statistics (with no recorded value)
        self.assertEqual(m.local_pvalues, {'B->pilnu::BR': 0.5})
        self.assertEqual(m.test_statistics, [
            {'name': 'B->pilnu::BR', 'type': 'chi^2', 'value': None, 'local_pvalue': 0.5},
        ])

    def test_invalid_path(self):
        "Loading from a non-existent path raises an error."
        with self.assertRaises(RuntimeError):
            eos.data.Mode(os.path.join(self._path, 'does-not-exist'))

    def test_roundtrip(self):
        "A mode written by create() reads back identically as format 2 (write/read symmetry)."
        test_statistics = [
            {'name': 'B->pilnu::BR',    'type': 'chi^2', 'value': 3.50, 'local_pvalue': 0.06},
            {'name': 'B->Dlnu::R_D',    'type': 'chi^2', 'value': 1.20, 'local_pvalue': 0.27},
        ]
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'mode-default')
            eos.data.Mode.create(path, self._parameters, np.array([3.67e-3, 0.27]), 0.25, test_statistics, 12.0, 5.0)

            m = eos.data.Mode(path)
            self.assertEqual(m.type, 'Mode')
            self.assertEqual(m.format_version, 2)
            self.assertEqual([p['name'] for p in m.varied_parameters], ['CKM::abs(V_ub)', 'B->pi::f_+(0)@BCL2008'])
            self.assertEqual(m.mode, [3.67e-3, 0.27])
            self.assertAlmostEqual(m.pvalue, 0.25)
            self.assertAlmostEqual(m.global_chi2, 12.0)
            self.assertAlmostEqual(m.dof, 5.0)
            self.assertEqual(m.test_statistics, test_statistics)
            self.assertEqual(m.local_pvalues, {'B->pilnu::BR': 0.06, 'B->Dlnu::R_D': 0.27})

    def test_roundtrip_none_pvalue(self):
        "A mode created with pvalue None and no test statistics reads back accordingly."
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'mode-default')
            eos.data.Mode.create(path, self._parameters, np.array([3.67e-3, 0.27]), None, [], None, None)

            m = eos.data.Mode(path)
            self.assertIsNone(m.pvalue)
            self.assertIsNone(m.global_chi2)
            self.assertIsNone(m.dof)
            self.assertEqual(m.test_statistics, [])
            self.assertEqual(m.local_pvalues, {})

    def test_from_bfp_and_gof(self):
        "from_bfp_and_gof derives the mode, chi2/dof and p-values, then stores them."
        bfp = _BestFitPoint(_Analysis(self._parameters), np.array([3.67e-3, 0.27]))
        gof = _GoodnessOfFit(2.0, 3, [
            ('B->pilnu::BR', _Entry(1.0, 2)),
            ('B->Dlnu::R_D', _Entry(1.0, 1)),
        ])

        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'mode-default')
            eos.data.Mode.from_bfp_and_gof(path, bfp, gof)

            m = eos.data.Mode(path)
            self.assertEqual(m.format_version, 2)
            self.assertEqual([p['name'] for p in m.varied_parameters], ['CKM::abs(V_ub)', 'B->pi::f_+(0)@BCL2008'])
            self.assertEqual(m.mode, [3.67e-3, 0.27])
            self.assertAlmostEqual(m.global_chi2, 2.0)
            self.assertEqual(m.dof, 3)
            self.assertAlmostEqual(m.pvalue, float(1.0 - scipy.stats.chi2(3).cdf(2.0)))

            self.assertEqual([t['name']  for t in m.test_statistics], ['B->pilnu::BR', 'B->Dlnu::R_D'])
            self.assertEqual([t['type']  for t in m.test_statistics], ['chi^2', 'chi^2'])
            self.assertEqual([t['value'] for t in m.test_statistics], [1.0, 1.0])
            self.assertAlmostEqual(m.local_pvalues['B->pilnu::BR'], float(1.0 - scipy.stats.chi2(2).cdf(1.0)))
            self.assertAlmostEqual(m.local_pvalues['B->Dlnu::R_D'], float(1.0 - scipy.stats.chi2(1).cdf(1.0)))

    def test_wrong_type(self):
        "A description whose type is not 'Mode' is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'Prediction', 'format_version': 2,
                'parameters': [], 'mode': [], 'pvalue': None, 'test_statistics': [],
            })
            with self.assertRaises(ValueError):
                eos.data.Mode(d)

    def test_unknown_key(self):
        "An unexpected key in the description is rejected rather than silently ignored."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'Mode', 'format_version': 2,
                'parameters': [{'name': 'a', 'min': 0.0, 'max': 1.0}],
                'mode': [0.5], 'pvalue': 0.25, 'test_statistics': [], 'bogus': 42,
            })
            with self.assertRaises(ValueError):
                eos.data.Mode(d)

    def test_missing_mode(self):
        "A format-2 description missing the required 'mode' key is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'Mode', 'format_version': 2,
                'parameters': [{'name': 'a', 'min': 0.0, 'max': 1.0}],
                'pvalue': 0.25, 'test_statistics': [],
            })
            with self.assertRaises(ValueError):
                eos.data.Mode(d)

    def test_format1_missing_local_pvalues(self):
        "A format-1 description without the required 'local_pvalues' map is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'Mode', 'format_version': 1,
                'parameters': [{'name': 'a', 'min': 0.0, 'max': 1.0}],
                'mode': [0.5], 'pvalue': 0.25,
            })
            with self.assertRaises(ValueError):
                eos.data.Mode(d)

    def test_non_mapping_description(self):
        "A description.yaml that is not a top-level mapping is rejected."
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'description.yaml'), 'w') as f:
                yaml.safe_dump(['not', 'a', 'mapping'], f)
            with self.assertRaises(RuntimeError):
                eos.data.Mode(d)


if __name__ == '__main__':
    unittest.main(verbosity=5)
