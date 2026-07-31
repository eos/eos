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

from eos.analysis_file_description import AnalysisFileDescription
from eos.diagnostic import Severity
from eos.validation_context import ValidationContext


class ValidationContextTests(unittest.TestCase):

    @staticmethod
    def _description(tag):
        return AnalysisFileDescription.from_dict(
            priors=[{
                'name': f'prior-{tag}',
                'descriptions': [{
                    'parameter': f'test::{tag}',
                    'type': 'uniform',
                    'min': 0.0,
                    'max': 1.0,
                }],
            }],
            likelihoods=[{
                'name': f'likelihood-{tag}',
                'manual_constraints': {
                    f'test::{tag}-constraint': {'type': 'Gaussian'},
                },
            }],
            posteriors=[{
                'name': f'posterior-{tag}',
                'prior': [f'prior-{tag}'],
                'likelihood': [f'likelihood-{tag}'],
                'fixed_parameters': {
                    f'test::{tag}-alias': 0.5,
                },
            }],
            predictions=[{
                'name': f'prediction-{tag}',
                'observables': [{
                    'name': f'test::{tag}-observable',
                }],
            }],
            observables={
                f'test::{tag}-observable': {
                    'latex': tag,
                    'unit': '1',
                    'expression': '1',
                },
            },
            parameters={
                f'test::{tag}': {
                    'latex': tag,
                    'unit': '1',
                    'central': 0.5,
                    'min': 0.0,
                    'max': 1.0,
                    'alias_of': [f'test::{tag}-alias'],
                },
            },
            masks=[{
                'name': f'mask-{tag}',
                'description': [{
                    'name': f'test::{tag}-mask-observable',
                    'expression': '1',
                }],
            }],
        )

    def test_lookup_uses_shadow_and_real_registries(self):
        context = ValidationContext(self._description('one'))

        for kind, name in (
            ('observable', 'test::one-observable'),
            ('observable', 'test::one-mask-observable'),
            ('parameter', 'test::one'),
            ('parameter', 'test::one-alias'),
            ('constraint', 'test::one-constraint'),
            ('prior', 'prior-one'),
            ('likelihood', 'likelihood-one'),
            ('mask', 'mask-one'),
            ('posterior', 'posterior-one'),
        ):
            self.assertTrue(context.lookup(kind, name), (kind, name))

        self.assertTrue(context.lookup('parameter', eos.QualifiedName('mass::e')))
        self.assertTrue(context.lookup('observable', eos.QualifiedName('B->Dlnu::BR')))
        self.assertTrue(context.lookup(
            'constraint',
            eos.QualifiedName('B->D::f_++f_0@HPQCD:2015A'),
        ))
        self.assertFalse(context.lookup('parameter', 'test::not-declared'))
        self.assertFalse(context.lookup('prior', 'prior-not-declared'))

    def test_contexts_are_file_local(self):
        first = ValidationContext(self._description('first'))
        second = ValidationContext(self._description('second'))

        for kind, first_name, second_name in (
            ('observable', 'test::first-observable', 'test::second-observable'),
            ('parameter', 'test::first', 'test::second'),
            ('constraint', 'test::first-constraint', 'test::second-constraint'),
            ('prior', 'prior-first', 'prior-second'),
            ('likelihood', 'likelihood-first', 'likelihood-second'),
            ('mask', 'mask-first', 'mask-second'),
            ('posterior', 'posterior-first', 'posterior-second'),
        ):
            self.assertTrue(first.lookup(kind, first_name))
            self.assertFalse(first.lookup(kind, second_name))
            self.assertTrue(second.lookup(kind, second_name))
            self.assertFalse(second.lookup(kind, first_name))

    def test_phases_one_and_two_are_idempotent_and_do_not_change_registries(self):
        parameter_names_before = tuple(str(parameter.name()) for parameter in eos.Parameters())
        observable_names_before = tuple(str(name) for name, _ in eos.Observables())
        custom_parameter = eos.QualifiedName('test::idempotent-first')
        custom_observable = eos.QualifiedName('test::idempotent-first-observable')
        self.assertFalse(eos.Parameters().has(custom_parameter))
        self.assertNotIn(custom_observable, eos.Observables())

        first_description = self._description('idempotent-first')
        first_diagnostics = (
            list(first_description.validate_structure()),
            list(first_description.validate_semantics(ValidationContext(first_description))),
        )
        repeated_diagnostics = (
            list(first_description.validate_structure()),
            list(first_description.validate_semantics(ValidationContext(first_description))),
        )
        self.assertEqual(first_diagnostics, repeated_diagnostics)
        # this test is about idempotence and registry immutability, not about the exact set of
        # diagnostics: assert only that neither phase reports an error. The fixture legitimately
        # contains unused entities, which are reported as warnings, and further warning-level
        # checks may be added without invalidating what is tested here.
        for phase in first_diagnostics:
            self.assertEqual([d for d in phase if d.severity is Severity.ERROR], [])

        second_description = self._description('idempotent-second')
        list(second_description.validate_structure())
        list(second_description.validate_semantics(ValidationContext(second_description)))

        self.assertEqual(
            parameter_names_before,
            tuple(str(parameter.name()) for parameter in eos.Parameters()),
        )
        self.assertEqual(
            observable_names_before,
            tuple(str(name) for name, _ in eos.Observables()),
        )
        self.assertFalse(eos.Parameters().has(custom_parameter))
        self.assertNotIn(custom_observable, eos.Observables())

        first_context = ValidationContext(first_description)
        second_context = ValidationContext(second_description)
        self.assertFalse(first_context.lookup('parameter', 'test::idempotent-second'))
        self.assertFalse(second_context.lookup('parameter', 'test::idempotent-first'))

    def test_unknown_kind_is_rejected(self):
        context = ValidationContext(self._description('one'))
        with self.assertRaisesRegex(ValueError, 'Unknown validation entity kind'):
            context.lookup('prediction', 'anything')


if __name__ == '__main__':
    unittest.main(verbosity=5)
