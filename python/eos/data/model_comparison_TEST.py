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

import os
import shutil
import tempfile
import unittest

import eos
import eos.data
import yaml


_ENTRIES = [
    { 'posterior': 'A', 'log_evidence': -1.0, 'log_evidence_uncertainty': 0.1, 'log_prior_volume_adjustment': 0.0,
      'adjusted_log_evidence': -1.0, 'log_bayes_factor': 0.0, 'log_bayes_factor_uncertainty': 0.0, 'strength': 'reference' },
    { 'posterior': 'B', 'log_evidence': -3.0, 'log_evidence_uncertainty': 0.1, 'log_prior_volume_adjustment': 0.5,
      'adjusted_log_evidence': -2.5, 'log_bayes_factor': -1.5, 'log_bayes_factor_uncertainty': 0.14, 'strength': 'positive' },
]
_COMPARISONS = [
    { 'first': 'A', 'second': 'B', 'log_bayes_factor': 1.5, 'log_bayes_factor_uncertainty': 0.14,
      'unadjusted_log_bayes_factor': 2.0, 'strength': 'positive', 'failed_checks': [] },
]
_CHECKS = [
    { 'name': 'uncertainty', 'status': 'passed', 'description': 'some check', 'failed_pairs': [] },
    { 'name': 'other',       'status': 'failed', 'description': 'another check', 'failed_pairs': [['A', 'B']] },
]


class ModelComparisonTests(unittest.TestCase):

    def setUp(self):
        self.path = tempfile.mkdtemp(prefix='eos-model-comparison-')
        self.addCleanup(shutil.rmtree, self.path, ignore_errors=True)

    def test_round_trip(self):
        "A model comparison written to disk reads back unchanged."
        eos.data.ModelComparison.create(self.path, 'grp', 'A', _ENTRIES, _COMPARISONS, _CHECKS)
        mc = eos.data.ModelComparison(self.path)
        self.assertEqual(mc.type, 'ModelComparison')
        self.assertEqual(mc.group, 'grp')
        self.assertEqual(mc.reference, 'A')
        self.assertEqual(mc.posteriors, _ENTRIES)
        self.assertEqual(mc.comparisons, _COMPARISONS)
        self.assertEqual(mc.checks, _CHECKS)
        self.assertFalse(mc.stable)
        self.assertIn('<tt>grp</tt>', mc._repr_html_())

    def test_wrong_type_raises(self):
        "A description of a different type is rejected."
        eos.data.ModelComparison.create(self.path, 'grp', 'A', _ENTRIES, _COMPARISONS, _CHECKS)
        description_file = os.path.join(self.path, 'description.yaml')
        with open(description_file) as f:
            description = yaml.safe_load(f)
        description['type'] = 'Mode'
        with open(description_file, 'w') as f:
            yaml.safe_dump(description, f)
        with self.assertRaises(ValueError):
            eos.data.ModelComparison(self.path)


if __name__ == '__main__':
    unittest.main(verbosity=5)
