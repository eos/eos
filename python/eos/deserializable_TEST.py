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

from dataclasses import dataclass, field

from eos.deserializable import Deserializable, InvalidComponent
from eos.diagnostic import Severity


@dataclass
class _Description(Deserializable):
    first: str
    second: int
    third: float
    optional: bool = False
    type: str = field(repr=False, init=False, default='description')


@dataclass
class _ConstructionObserver(Deserializable):
    value: str

    constructed = False

    def __post_init__(self):
        _ConstructionObserver.constructed = True
        raise RuntimeError('constructor called')


class DeserializableTests(unittest.TestCase):

    def test_reports_all_missing_mandatory_keys(self):
        result = Deserializable.make_with_diagnostics(_Description)

        self.assertIsInstance(result, InvalidComponent)
        self.assertEqual(
            {diagnostic.path for diagnostic in result.diagnostics},
            {('first',), ('second',), ('third',)}
        )
        self.assertTrue(all(diagnostic.severity is Severity.ERROR for diagnostic in result.diagnostics))

    def test_reports_unknown_and_missing_keys_together(self):
        result = Deserializable.make_with_diagnostics(_Description, first='value', unknown=None)

        self.assertIsInstance(result, InvalidComponent)
        self.assertEqual(
            {diagnostic.path for diagnostic in result.validate_structure()},
            {('unknown',), ('second',), ('third',)}
        )

    def test_init_false_field_is_not_required_and_is_rejected_if_supplied(self):
        valid = Deserializable.make_with_diagnostics(_Description, first='value', second=2, third=3.0)
        self.assertIsInstance(valid, _Description)

        invalid = Deserializable.make_with_diagnostics(
            _Description, first='value', second=2, third=3.0, type='other'
        )
        self.assertIsInstance(invalid, InvalidComponent)
        self.assertEqual([('type',)], [diagnostic.path for diagnostic in invalid.diagnostics])

    def test_well_formed_mapping_returns_instance(self):
        result = Deserializable.make_with_diagnostics(
            _Description, first='value', second=2, third=3.0, optional=True
        )

        self.assertEqual(
            _Description(first='value', second=2, third=3.0, optional=True),
            result
        )

    def test_check_keys_does_not_construct(self):
        _ConstructionObserver.constructed = False

        diagnostics = Deserializable.check_keys(_ConstructionObserver, {'value': 'valid'})

        self.assertEqual([], diagnostics)
        self.assertFalse(_ConstructionObserver.constructed)


if __name__ == '__main__':
    unittest.main(verbosity=5)
