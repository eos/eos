#!/usr/bin/python
# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2024-2026 Danny van Dyk
# Copyright (c) 2025 Matthew Kirk
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

import contextlib
import io
import eos
import os
from pathlib import Path

from eos.diagnostic import Diagnostic, Severity
from eos.validation_context import ValidationContext

_TESTD = Path(__file__).parent / 'analysis_file_TEST.d'

class TestAnalysisFile(unittest.TestCase):

    def test_analysis_file(self):

        af = eos.AnalysisFile(_TESTD / 'valid-analysis-file.yaml')
        af.validate()

    def test_minimal_analysis_file(self):

        af = eos.AnalysisFile(_TESTD / 'minimal-analysis-file.yaml')
        af.validate()

    def test_prior_only_analysis_file(self):

        # a file without a 'likelihoods' section is a legitimate prior-only analysis; it must load
        # and expose no likelihoods. Regression test: this previously raised a KeyError because the
        # (absent) 'likelihoods' entry was accessed unconditionally.
        with self.assertNoLogs('EOS', level='WARNING'):
            af = eos.AnalysisFile(_TESTD / 'no-likelihoods-analysis-file.yaml')
        self.assertEqual(dict(af.likelihoods), {})
        warnings = [
            diagnostic
            for diagnostic in af._description.validate_structure()
            if diagnostic.severity is Severity.WARNING
        ]
        self.assertEqual(
            [(('likelihoods',), 'No likelihood components found in analysis file')],
            [(diagnostic.path, diagnostic.message) for diagnostic in warnings],
        )
        af.validate()

    def test_deprecations_are_located_warning_diagnostics(self):
        with self.assertNoLogs('EOS', level='WARNING'):
            description = eos.analysis_file_description.AnalysisFileDescription.from_dict(
                priors=[{
                    'name': 'legacy-prior',
                    'parameters': [
                        {'constraint': 'B->D::f_++f_0@HPQCD:2015A'},
                        {
                            'parameter': 'mass::e',
                            'type': 'gaussian',
                            'central': 0.5,
                            'sigma': 0.1,
                            'min': 0.0,
                            'max': 1.0,
                        },
                    ],
                }],
                likelihoods=[{
                    'name': 'likelihood',
                    'constraints': ['B->D::f_++f_0@HPQCD:2015A'],
                }],
                posteriors=[{
                    'name': 'posterior',
                    'prior': ['legacy-prior'],
                    'likelihood': ['likelihood'],
                }],
            )
            diagnostics = list(description.validate_structure())

        warnings = {
            diagnostic.path: diagnostic
            for diagnostic in diagnostics
            if diagnostic.severity is Severity.WARNING
        }
        self.assertIn(('priors', 'legacy-prior', 'parameters'), warnings)
        self.assertIn(('priors', 'legacy-prior', 'descriptions', 0, 'type'), warnings)
        self.assertIn(('priors', 'legacy-prior', 'descriptions', 1), warnings)

    def test_both_prior_description_keys_produce_warning_without_logging(self):
        with self.assertNoLogs('EOS', level='WARNING'):
            prior = eos.analysis_file_description.PriorComponent.from_dict(
                name='legacy-prior',
                descriptions=[{
                    'parameter': 'mass::e',
                    'type': 'uniform',
                    'min': 0.0,
                    'max': 1.0,
                }],
                parameters=[{
                    'constraint': 'B->D::f_++f_0@HPQCD:2015A',
                    'type': 'constraint',
                }],
            )

        warnings = [
            diagnostic
            for diagnostic in prior.validate_structure()
            if diagnostic.severity is Severity.WARNING
        ]
        self.assertEqual(2, len(warnings))
        self.assertTrue(all(diagnostic.path == ('parameters',) for diagnostic in warnings))

    def test_task_defaults_are_warning_diagnostics_without_logging(self):
        with self.assertNoLogs('EOS', level='WARNING'):
            task = eos.analysis_file_description.TaskComponent.from_dict(
                task='corner-plot',
                arguments={'posterior': 'posterior'},
            )

        warnings = [
            diagnostic
            for diagnostic in task.validate_structure()
            if diagnostic.severity is Severity.WARNING
        ]
        self.assertTrue(warnings)
        self.assertTrue(all(diagnostic.path[:1] == ('arguments',) for diagnostic in warnings))

    def test_format_version(self):

        # the fixture declares 'format_version: 1' explicitly
        af = eos.AnalysisFile(_TESTD / 'valid-analysis-file.yaml')
        self.assertEqual(af.format_version, 1)

    def test_format_version_defaults_to_one(self):

        # a file without a 'format_version' key predates format versioning and is treated as v1
        af = eos.AnalysisFile(_TESTD / 'minimal-analysis-file.yaml')
        self.assertEqual(af.format_version, 1)


class TestAnalysisFileConstructionErrors(unittest.TestCase):
    """The mandatory-structure and uniqueness checks in ``__init__`` must raise."""

    def test_nonexistent_file(self):
        with self.assertRaises(RuntimeError):
            eos.AnalysisFile(_TESTD / 'this-file-does-not-exist.yaml')

    def test_not_a_file(self):
        # a directory is not a valid analysis file
        with self.assertRaises(RuntimeError):
            eos.AnalysisFile(_TESTD)

    def test_future_format_version(self):
        # a file declaring a newer format version than this EOS supports must be refused rather
        # than silently misinterpreted
        with self.assertRaises(RuntimeError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'future-format-version.yaml')

    def test_empty_file(self):
        # an empty file parses to None, which is not a top-level mapping; loading must raise a clear
        # RuntimeError rather than an opaque AttributeError from the subsequent .get()/** access
        with self.assertRaises(RuntimeError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'empty-file.yaml')

    def test_non_integer_format_version(self):
        # a 'format_version' that is not an integer (e.g. a quoted string) must be refused rather
        # than raising an opaque TypeError from the numeric comparison
        with self.assertRaises(RuntimeError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'string-format-version.yaml')

    def test_missing_priors(self):
        with self.assertRaises(RuntimeError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'no-priors.yaml')

    def test_missing_posteriors(self):
        with self.assertRaises(RuntimeError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'no-posteriors.yaml')

    def test_unknown_prior_reference(self):
        with self.assertRaises(RuntimeError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'unknown-prior-ref.yaml')

    def test_unknown_likelihood_reference(self):
        with self.assertRaises(RuntimeError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'unknown-likelihood-ref.yaml')

    def test_invalid_parameter_name(self):
        with self.assertRaises(ValueError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'bad-parameter.yaml')

    def test_duplicate_step_id(self):
        with self.assertRaises(RuntimeError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'duplicate-step-id.yaml')

    def test_duplicate_mask_name(self):
        with self.assertRaises(RuntimeError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'duplicate-mask-name.yaml')

    def test_unknown_mask_reference(self):
        with self.assertRaises(RuntimeError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'unknown-mask-ref.yaml')

    def test_invalid_file_local_component_name(self):
        with self.assertRaises(RuntimeError) as context:
            eos.AnalysisFile(_TESTD / 'invalid' / 'bad-component-name.yaml')

        self.assertIn(
            "likelihoods[1]/name: Invalid character '/' in likelihood name 'a/b'",
            str(context.exception),
        )

    def test_slash_in_observable_qualified_name_is_allowed(self):
        description = eos.analysis_file_description.AnalysisFileDescription.from_dict(
            priors=[{
                'name': 'prior',
                'descriptions': [{
                    'parameter': 'CKM::abs(V_ub)',
                    'min': 3.0e-3,
                    'max': 4.5e-3,
                    'type': 'uniform',
                }],
            }],
            likelihoods=[{
                'name': 'likelihood',
                'constraints': ['B^+->tau^+nu::BR@Belle:2014A;form-factors=BCL2008-4'],
            }],
            posteriors=[{
                'name': 'posterior',
                'prior': ['prior'],
                'likelihood': ['likelihood'],
            }],
            observables={
                'Test::ratio/a': {
                    'latex': r'R/a',
                    'unit': '1',
                    'expression': '1.0',
                },
            },
        )

        self.assertEqual([], list(description.validate_structure()))

    def test_multiple_structural_errors_are_reported_together(self):
        with self.assertRaises(RuntimeError) as context:
            eos.AnalysisFile(_TESTD / 'invalid' / 'multiple-structural-errors.yaml')

        message = str(context.exception)
        for fragment in (
            'priors/incomplete-prior/descriptions[0]/max',
            'likelihoods/empty-likelihood',
            'posteriors/broken-posterior/prior[0]',
            'posteriors/broken-posterior/likelihood[0]',
            'steps[0]/id',
            'steps[0]/tasks[0]/task',
            'masks/broken-mask/logical_combination',
            'masks/broken-mask/description[0]/mask_name',
        ):
            self.assertIn(fragment, message)


class TestAnalysisFileMethods(unittest.TestCase):
    """Exercises the accessors and factory methods on a fully-populated analysis file."""

    def setUp(self):
        # a fresh instance per test: analysis()/observables() read (and, for the latter, mutate in
        # place) the cached descriptions, so tests must not share a single instance
        self.af = eos.AnalysisFile(_TESTD / 'valid-analysis-file.yaml')

    def test_properties(self):
        self.assertIn('CKM', self.af.priors)
        self.assertIn('TH-pi', self.af.likelihoods)
        self.assertIn('CKM-all', self.af.posteriors)
        self.assertIn('leptonic-BR-CKM', self.af.predictions)
        self.assertIn('mask-A', self.af.masks)

    def test_repr_html(self):
        # the HTML representation renders every section present in the fixture
        html = self.af._repr_html_()
        for section in ('PRIORS', 'LIKELIHOODS', 'POSTERIORS', 'PREDICTIONS',
                        'OBSERVABLES', 'PARAMETERS', 'STEPS', 'MASKS'):
            self.assertIn(section, html)

    def test_dump(self):
        # dump() writes the input data back out as YAML
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            self.af.dump()
        self.assertIn('priors', buffer.getvalue())

    def test_analysis(self):
        analysis = self.af.analysis('CKM-all')
        self.assertIsNotNone(analysis)

    def test_analysis_unknown_posterior(self):
        with self.assertRaises(RuntimeError):
            self.af.analysis('does-not-exist')

    def test_observables(self):
        observables = self.af.observables('CKM-all', 'leptonic-BR-CKM', eos.Parameters())
        self.assertGreater(len(observables), 0)

    def test_observables_unknown_posterior(self):
        with self.assertRaises(RuntimeError):
            self.af.observables('does-not-exist', 'leptonic-BR-CKM', eos.Parameters())

    def test_observables_unknown_prediction(self):
        with self.assertRaises(RuntimeError):
            self.af.observables('CKM-all', 'does-not-exist', eos.Parameters())

    def test_observable(self):
        observable = self.af.observable('CKM-all', 'B_u->lnu::BR;l=e', eos.Parameters())
        self.assertIsNotNone(observable)

    def test_observable_applies_fixed_parameters(self):
        # 'WET-all' fixes 'CKM::abs(V_ub)', so observable() writes that value into the parameter set
        parameters = eos.Parameters()
        self.af.observable('WET-all', 'B_u->lnu::BR;l=e', parameters)
        self.assertAlmostEqual(float(parameters['CKM::abs(V_ub)']), 3.67e-3)

    def test_observable_detects_global_option_conflict(self):
        # the option-part (model=CKM) in the name conflicts with the posterior's global option
        # (model=WET); observable() must report this (regression test: the option-part keys are
        # qnpOptionKey objects that previously never matched the str-keyed option dict, so the
        # conflict went undetected). The observable is still created.
        with self.assertLogs('EOS', level='ERROR') as cm:
            observable = self.af.observable('WET-all', 'B_u->lnu::BR;l=e,model=CKM', eos.Parameters())
        self.assertIsNotNone(observable)
        self.assertTrue(any('overrides option part' in message for message in cm.output))

    def test_observables_detects_global_option_conflict(self):
        # the same conflict, but reported by observables() for a prediction observable whose name
        # (model=WET) contradicts the prediction's global option (model=CKM)
        with self.assertLogs('EOS', level='ERROR') as cm:
            observables = self.af.observables('CKM-all', 'option-conflict', eos.Parameters())
        self.assertGreater(len(observables), 0)
        self.assertTrue(any('overrides option part' in message for message in cm.output))

    def test_observable_unknown_posterior(self):
        with self.assertRaises(RuntimeError):
            self.af.observable('does-not-exist', 'B_u->lnu::BR;l=e', eos.Parameters())

    def test_observable_unknown_name(self):
        # an unknown observable name is reported with a clear message (regression test: the
        # 'if not observable' guard was dead because make() raises instead of returning None)
        with self.assertRaises(RuntimeError) as ctx:
            self.af.observable('CKM-all', 'B->pilnu::TOTALLY_UNKNOWN', eos.Parameters())
        self.assertIn('Unknown observable name', str(ctx.exception))


class TestAnalysisFileValidation(unittest.TestCase):
    """The reporting (non-raising) branches of ``validate``."""

    @staticmethod
    def _figure_description(observable, observables=None):
        return eos.analysis_file_description.AnalysisFileDescription.from_dict(
            figures=[{
                'name': 'semantic-figure',
                'type': 'single',
                'plot': {
                    'items': [{
                        'type': 'constraint',
                        'constraints': 'B^0->D^+e^-nu::BRs@Belle:2015A',
                        'variable': 'q2',
                        'observable': observable,
                    }],
                },
            }],
            observables=observables or {},
        )

    def test_figure_unknown_observable_is_located(self):
        description = self._figure_description('test::unknown-observable')
        diagnostics = list(description.validate_semantics(ValidationContext(description)))

        self.assertEqual(1, len(diagnostics))
        self.assertEqual(
            ('figures', 'semantic-figure', 'plot', 'items', 0, 'observable'),
            diagnostics[0].path,
        )
        self.assertIn("observable 'test::unknown-observable' is unknown", diagnostics[0].message)

    def test_figure_resolves_custom_observable_from_shadow_registry(self):
        name = 'test::figure-only'
        description = self._figure_description(name, {
            name: {'latex': 'x', 'unit': '1', 'expression': '1.0'},
        })
        context = ValidationContext(description)

        self.assertEqual([], list(description.figures[0].validate_semantics(context)))
        self.assertNotIn(name, context.unused('observable'))

    def test_parameter_alias_lookup_credits_declaration(self):
        name = 'test::aliased-parameter'
        alias = 'test::parameter-alias'
        description = eos.analysis_file_description.AnalysisFileDescription.from_dict(
            parameters={
                name: {
                    'latex': 'a',
                    'unit': '1',
                    'central': 0.0,
                    'min': -1.0,
                    'max': +1.0,
                    'alias_of': [alias],
                },
            },
        )
        context = ValidationContext(description)

        self.assertTrue(context.lookup('parameter', alias))
        self.assertNotIn(name, context.unused('parameter'))

    def test_unused_entities_are_warnings(self):
        af = eos.AnalysisFile(_TESTD / 'unused-entities.yaml')
        description = af._description

        structure = list(description.validate_structure())
        semantics = list(description.validate_semantics(ValidationContext(description)))
        self.assertFalse(any(
            diagnostic.severity is Severity.ERROR
            for diagnostic in structure + semantics
        ))

        warnings = [
            diagnostic
            for diagnostic in semantics
            if diagnostic.severity is Severity.WARNING
        ]
        self.assertEqual(
            [
                ('priors', 'unused-prior'),
                ('likelihoods', 'unused-likelihood'),
                ('masks', 'unused-mask'),
                ('observables', 'test::unused-observable'),
                ('parameters', 'test::unused-parameter'),
                ('masks', 'unused-mask', 'description', 'test::unused-mask-observable'),
                ('likelihoods', 'unused-likelihood', 'manual_constraints', 'test::unused-constraint'),
            ],
            [diagnostic.path for diagnostic in warnings],
        )

        warned_paths = {diagnostic.path for diagnostic in warnings}
        for path in (
            ('priors', 'used-prior'),
            ('likelihoods', 'used-likelihood'),
            ('masks', 'mask-used-by-mask'),
            ('masks', 'mask-composite'),
            ('masks', 'mask-used-by-task'),
            ('parameters', 'test::aliased-parameter'),
            ('parameters', 'test::expression-only-parameter'),
            ('observables', 'test::manual-constraint-observable'),
            ('observables', 'test::manual-constraint-observable-list'),
            ('observables', 'test::figure-only-observable'),
        ):
            self.assertNotIn(path, warned_paths)

    def test_observables_aggregates_unknown_names(self):
        # observables() reports every unknown observable in a prediction, not just the first
        # (regression test: make() raises on the first, so the aggregation must catch and continue),
        # as a comma-separated list in the order the observables appear in the analysis file
        af = eos.AnalysisFile(_TESTD / 'invalid' / 'semantic-errors.yaml')
        with self.assertRaises(RuntimeError) as ctx:
            af.observables('post-good', 'pred-bad-obs', eos.Parameters())
        message = str(ctx.exception)
        self.assertIn('B->pilnu::NOTREAL, B->pilnu::ALSO_NOTREAL', message)

    def test_semantic_errors_are_reported(self):
        af = eos.AnalysisFile(_TESTD / 'invalid' / 'semantic-errors.yaml')
        messages = af.validate()
        # validate() collects one message per problem rather than raising; the fixture contains
        # many independent problems, so it must return a non-empty list
        self.assertGreater(len(messages), 0)
        blob = '\n'.join(str(diagnostic) for diagnostic in messages)
        # a representative selection of the distinct problem classes in the fixture
        self.assertIn('Not::a-parameter', blob)                 # unknown prior parameter
        self.assertIn('Not::a-constraint@Nowhere:2000A', blob)  # unknown constraint
        self.assertIn('matches an already defined constraint', blob)  # colliding manual constraint
        self.assertIn('B->pilnu::NOTREAL', blob)                # unknown prediction observable
        self.assertIn('Not::a-fixed-parameter', blob)           # unknown posterior fixed parameter
        self.assertIn('Not::a-prediction-parameter', blob)      # unknown prediction fixed parameter
        self.assertIn('unknown steps', blob)                    # step dependency
        self.assertIn('unknown tasks', blob)                    # step default arguments
        self.assertIn("Posterior 'nonexistent-posterior'", blob)  # step task posterior
        self.assertIn('repeatedly', blob)                       # repeated mask expression name
        # every message is a Diagnostic, including those from the deep phase, so that callers can
        # filter by severity and location without special-casing the origin of a message
        self.assertTrue(all(isinstance(message, Diagnostic) for message in messages))
        self.assertTrue(any(
            diagnostic.severity is Severity.ERROR and diagnostic.path[0] in ('posteriors', 'predictions')
            for diagnostic in messages
        ))
        semantic_locations = {
            (diagnostic.path, diagnostic.severity)
            for diagnostic in messages
        }
        self.assertIn(
            (('priors', 'BAD-PARAM', 'descriptions', 0, 'parameter'), Severity.ERROR),
            semantic_locations,
        )
        self.assertIn(
            (
                ('predictions', 'pred-bad-obs', 'observables', 'B->pilnu::NOTREAL', 'name'),
                Severity.ERROR,
            ),
            semantic_locations,
        )
        self.assertIn(
            (('steps', 'bad-step', 'depends_on'), Severity.ERROR),
            semantic_locations,
        )

    def test_fixed_but_unused_parameter(self):
        # the fixture's posterior fixes one parameter its likelihood uses and one from an unrelated
        # decay; only the latter has no effect and must be reported
        af = eos.AnalysisFile(_TESTD / 'fixed-but-unused-parameter.yaml')

        reported = {
            diagnostic.path[-1]
            for diagnostic in af.validate(deep=True)
            if 'is fixed but used by neither' in diagnostic.message
        }
        self.assertEqual(reported, {'0->Kpi::M_(+,0)@KSvD2025'})
        self.assertNotIn('B->pi::b_+^1@BCL2008', reported)

        # the check needs the assembled likelihood, so the side-effect-free phases cannot see it
        self.assertFalse(any(
            'is fixed but used by neither' in diagnostic.message
            for diagnostic in af.validate(deep=False)
        ))


if __name__ == '__main__':
    unittest.main(verbosity=5)
