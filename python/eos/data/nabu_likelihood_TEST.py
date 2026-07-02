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
import os
import tempfile
import yaml


# These tests exercise the hardening of the description parsing, which happens before the (optional)
# nabu payload is loaded, so they do not require the nabu module to be installed.
class NabuLikelihoodTests(unittest.TestCase):

    @staticmethod
    def _write(directory, description):
        "Materialize a (possibly malformed) NabuLikelihood description in the given directory."
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, 'description.yaml'), 'w') as f:
            yaml.safe_dump(description, f, default_flow_style=False)

    def test_invalid_path(self):
        "Loading from a non-existent path raises an error."
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(RuntimeError):
                eos.data.NabuLikelihood(os.path.join(d, 'does-not-exist'))

    def test_wrong_type(self):
        "A description whose type is not 'NabuLikelihood' is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'MarkovChain', 'parameters': []})
            with self.assertRaises(ValueError):
                eos.data.NabuLikelihood(d)

    def test_unknown_key(self):
        "An unexpected key in the description is rejected rather than silently ignored."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'NabuLikelihood',
                'parameters': [{'name': 'a', 'min': 0.0, 'max': 1.0}], 'bogus': 42,
            })
            with self.assertRaises(ValueError):
                eos.data.NabuLikelihood(d)

    def test_missing_parameters(self):
        "A description missing the required 'parameters' key is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {'version': 'test', 'type': 'NabuLikelihood'})
            with self.assertRaises(ValueError):
                eos.data.NabuLikelihood(d)

    def test_missing_payload(self):
        "A valid description without the likelihood payload is rejected."
        with tempfile.TemporaryDirectory() as d:
            self._write(d, {
                'version': 'test', 'type': 'NabuLikelihood',
                'parameters': [{'name': 'a', 'min': 0.0, 'max': 1.0}],
            })
            with self.assertRaises(RuntimeError):
                eos.data.NabuLikelihood(d)

    def test_non_mapping_description(self):
        "A description.yaml that is not a top-level mapping is rejected."
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'description.yaml'), 'w') as f:
                yaml.safe_dump(['not', 'a', 'mapping'], f)
            with self.assertRaises(RuntimeError):
                eos.data.NabuLikelihood(d)


if __name__ == '__main__':
    unittest.main(verbosity=5)
