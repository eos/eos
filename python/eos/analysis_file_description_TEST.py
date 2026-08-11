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
from dataclasses import asdict

import eos

from eos.deserializable import InvalidComponent
from eos.diagnostic import Severity
from eos.validation_context import ValidationContext
from eos.analysis_file_description import (
    MetadataAuthorDescription,
    MetadataDescription,
    PriorDescription,
    PoissonPriorDescription,
    CurtailedGaussianDescription,
    GaussianPriorDescription,
    ScalePriorDescription,
    UniformPriorDescription,
    ConstraintPriorDescription,
    TransformPriorDescription,
    PriorComponent,
    ConstraintLikelihoodDescription,
    ManualConstraintDescription,
    PyHFConstraintDescription,
    LikelihoodComponent,
    PosteriorDescription,
    ObservableComponent,
    PredictionObservableComponent,
    PredictionDescription,
    ParameterComponent,
    TaskComponent,
    StepComponent,
    MaskDescription,
    MaskExpressionComponent,
    MaskObservableComponent,
    MaskNamedComponent,
    MaskComponent,
    AnalysisFileDescription,
)


class MetadataAuthorDescriptionTests(unittest.TestCase):

    def test_from_dict(self):
        author = MetadataAuthorDescription.from_dict(name='Jane Doe', affiliation='Some University', email='jane@example.org')
        self.assertEqual(author.name, 'Jane Doe')
        self.assertEqual(author.affiliation, 'Some University')
        self.assertEqual(author.email, 'jane@example.org')

    def test_defaults(self):
        author = MetadataAuthorDescription.from_dict(name='Jane Doe')
        self.assertEqual(author.affiliation, '')
        self.assertEqual(author.email, '')

    def test_fields_must_be_strings(self):
        for field_name in ('name', 'affiliation', 'email'):
            with self.subTest(field_name=field_name):
                values = {'name': 'Jane Doe', 'affiliation': '', 'email': ''}
                values[field_name] = 7
                author = MetadataAuthorDescription.from_dict(**values)
                diagnostics = list(author.validate_structure())
                self.assertEqual(len(diagnostics), 1)
                self.assertEqual(diagnostics[0].path, (field_name,))
                self.assertEqual(diagnostics[0].message, f'{field_name} must be a string')


class MetadataDescriptionTests(unittest.TestCase):

    def test_from_dict(self):
        md = MetadataDescription.from_dict(
            title='My Analysis',
            id='my-analysis',
            description='A detailed analysis description.',
            authors=[{'name': 'Jane Doe'}, {'name': 'John Roe', 'email': 'john@example.org'}],
        )
        self.assertEqual(md.title, 'My Analysis')
        self.assertEqual(md.id, 'my-analysis')
        self.assertEqual(md.description, 'A detailed analysis description.')
        # the authors are deserialized into MetadataAuthorDescription instances
        self.assertEqual(len(md.authors), 2)
        self.assertIsInstance(md.authors[0], MetadataAuthorDescription)
        self.assertEqual(md.authors[1].email, 'john@example.org')

    def test_defaults(self):
        md = MetadataDescription.from_dict()
        self.assertEqual(md.title, '')
        self.assertEqual(md.id, '')
        self.assertEqual(md.authors, [])
        self.assertEqual(md.description, '')


class PriorDescriptionTests(unittest.TestCase):

    def test_dispatch(self):
        "Check that the factory dispatches each prior type to the correct class."
        self.assertIsInstance(
            PriorDescription.from_dict(type='constraint', constraint='B->D::f_++f_0@HPQCD:2015A'),
            ConstraintPriorDescription)
        self.assertIsInstance(
            PriorDescription.from_dict(parameter='p', type='uniform', min=0.0, max=1.0),
            UniformPriorDescription)
        self.assertIsInstance(
            PriorDescription.from_dict(parameter='p', type='flat', min=0.0, max=1.0),
            UniformPriorDescription)
        self.assertIsInstance(
            PriorDescription.from_dict(parameter='p', type='gaussian', central=0.0, sigma=1.0),
            GaussianPriorDescription)
        # a gaussian with a 'min' key becomes a curtailed gaussian
        self.assertIsInstance(
            PriorDescription.from_dict(parameter='p', type='gaussian', central=0.0, sigma=1.0, min=-1.0, max=1.0),
            CurtailedGaussianDescription)
        self.assertIsInstance(
            PriorDescription.from_dict(parameter='p', type='poisson', k=3.0),
            PoissonPriorDescription)

    def test_unknown_type(self):
        unknown = PriorDescription.from_dict(parameter='p', type='not-a-type')
        missing = PriorDescription.from_dict()
        for desc in (unknown, missing):
            self.assertIsInstance(desc, InvalidComponent)
            diagnostics = list(desc.validate_structure())
            self.assertEqual(len(diagnostics), 1)
            self.assertEqual(diagnostics[0].path, ('type',))
            self.assertEqual(diagnostics[0].severity, Severity.ERROR)
            self.assertEqual(diagnostics[0].message, 'Unknown type of prior description')

    def test_gaussian_requires_both_bounds(self):
        "A gaussian prior with only one of 'min'/'max' is rejected with a clear error."
        missing_max = PriorDescription.from_dict(
            parameter='p', type='gaussian', central=0.0, sigma=1.0, min=-1.0)
        missing_min = PriorDescription.from_dict(
            parameter='p', type='gaussian', central=0.0, sigma=1.0, max=1.0)
        for desc, path in ((missing_max, ('max',)), (missing_min, ('min',))):
            self.assertIsInstance(desc, InvalidComponent)
            diagnostics = list(desc.validate_structure())
            self.assertEqual(len(diagnostics), 1)
            self.assertEqual(diagnostics[0].path, path)
            self.assertEqual(diagnostics[0].severity, Severity.ERROR)
            self.assertIn('must contain both', diagnostics[0].message)


class PoissonPriorDescriptionTests(unittest.TestCase):

    def test_from_dict(self):
        desc = PriorDescription.from_dict(parameter='mass::e', type='poisson', k=5.0)
        self.assertIsInstance(desc, PoissonPriorDescription)
        self.assertEqual(desc.parameter, 'mass::e')
        self.assertEqual(desc.k, 5.0)
        self.assertEqual(desc.type, 'poisson')


class CurtailedGaussianDescriptionTests(unittest.TestCase):

    def test_from_dict(self):
        desc = PriorDescription.from_dict(parameter='mass::e', type='gaussian', central=0.5, sigma=0.1, min=0.0, max=1.0)
        self.assertIsInstance(desc, CurtailedGaussianDescription)
        self.assertEqual(desc.parameter, 'mass::e')
        self.assertEqual(desc.central, 0.5)
        self.assertEqual(desc.sigma, 0.1)
        self.assertEqual(desc.min, 0.0)
        self.assertEqual(desc.max, 1.0)
        self.assertEqual(desc.type, 'gaussian')


class GaussianPriorDescriptionTests(unittest.TestCase):

    def test_from_dict(self):
        desc = PriorDescription.from_dict(parameter='decay-constant::B_u', type='gaussian', central=0.1894, sigma=0.0014)
        self.assertIsInstance(desc, GaussianPriorDescription)
        self.assertEqual(desc.parameter, 'decay-constant::B_u')
        self.assertEqual(desc.central, 0.1894)
        self.assertEqual(desc.sigma, 0.0014)
        self.assertEqual(desc.type, 'gaussian')


class ScalePriorDescriptionTests(unittest.TestCase):

    def test_from_dict(self):
        "Check that a scale prior description is deserialized correctly."

        desc = PriorDescription.from_dict(
            parameter='mass::b(MSbar)',
            type='scale',
            min=0.0,
            max=1.0,
            mu_0=1.0,
            lambda_scale=2.0,
        )

        # the factory dispatches the 'scale' type to ScalePriorDescription
        self.assertIsInstance(desc, ScalePriorDescription)
        self.assertEqual(desc.parameter, 'mass::b(MSbar)')
        self.assertEqual(desc.min, 0.0)
        self.assertEqual(desc.max, 1.0)
        self.assertEqual(desc.mu_0, 1.0)
        self.assertEqual(desc.lambda_scale, 2.0)
        self.assertEqual(desc.type, 'scale')

    def test_missing_field(self):
        "Check that an incomplete scale prior description is rejected."

        desc = PriorDescription.from_dict(parameter='mass::b(MSbar)', type='scale', min=0.0, max=1.0)
        self.assertIsInstance(desc, InvalidComponent)
        self.assertEqual(
            [(d.path, d.severity, d.message) for d in desc.validate_structure()],
            [(('lambda_scale',), Severity.ERROR, "Missing mandatory key 'lambda_scale'"),
             (('mu_0',), Severity.ERROR, "Missing mandatory key 'mu_0'")],
        )


class UniformPriorDescriptionTests(unittest.TestCase):

    def test_from_dict(self):
        desc = PriorDescription.from_dict(parameter='CKM::abs(V_ub)', type='uniform', min=3.0e-3, max=4.5e-3)
        self.assertIsInstance(desc, UniformPriorDescription)
        self.assertEqual(desc.parameter, 'CKM::abs(V_ub)')
        self.assertEqual(desc.min, 3.0e-3)
        self.assertEqual(desc.max, 4.5e-3)
        self.assertEqual(desc.type, 'uniform')


class ConstraintPriorDescriptionTests(unittest.TestCase):

    def test_from_dict(self):
        desc = PriorDescription.from_dict(type='constraint', constraint='B->D::f_++f_0@HPQCD:2015A')
        self.assertIsInstance(desc, ConstraintPriorDescription)
        self.assertEqual(desc.constraint, 'B->D::f_++f_0@HPQCD:2015A')
        self.assertEqual(desc.type, 'constraint')

    def test_from_dict_without_type_is_supported(self):
        "A constraint prior without an explicit 'type' key is still accepted (deprecated form)."
        desc = PriorDescription.from_dict(constraint='B->D::f_++f_0@HPQCD:2015A')
        self.assertIsInstance(desc, ConstraintPriorDescription)
        self.assertEqual(desc.constraint, 'B->D::f_++f_0@HPQCD:2015A')


class TransformPriorDescriptionTests(unittest.TestCase):

    def test_from_dict(self):
        "Check that a transform prior description is deserialized correctly."

        desc = PriorDescription.from_dict(
            parameters=['scnuee::Re{cVL}', 'scnuee::Re{cVR}'],
            shift=[0.0, 0.0],
            transform=[[0.707106, 0.707106], [-0.707106, 0.707106]],
            min=[-2.0, -2.0],
            max=[2.0, 2.0],
            type='transform',
        )

        # the factory dispatches the 'transform' type to TransformPriorDescription
        self.assertIsInstance(desc, TransformPriorDescription)
        self.assertEqual(desc.parameters, ['scnuee::Re{cVL}', 'scnuee::Re{cVR}'])
        self.assertEqual(desc.shift, [0.0, 0.0])
        self.assertEqual(desc.transform, [[0.707106, 0.707106], [-0.707106, 0.707106]])
        self.assertEqual(desc.min, [-2.0, -2.0])
        self.assertEqual(desc.max, [2.0, 2.0])
        self.assertEqual(desc.type, 'transform')


class PriorComponentTests(unittest.TestCase):

    def test_from_dict(self):
        comp = PriorComponent.from_dict(
            name='FF-pi',
            descriptions=[
                {'parameter': 'B->pi::f_+(0)@BCL2008', 'min': 0.21, 'max': 0.32, 'type': 'uniform'},
                {'parameter': 'B->pi::b_+^1@BCL2008', 'min': -2.96, 'max': -0.60, 'type': 'uniform'},
            ],
        )
        self.assertEqual(comp.name, 'FF-pi')
        self.assertEqual(len(comp.descriptions), 2)
        # each entry is deserialized through the PriorDescription factory
        self.assertIsInstance(comp.descriptions[0], UniformPriorDescription)


class ConstraintLikelihoodDescriptionTests(unittest.TestCase):

    def test_from_dict(self):
        desc = ConstraintLikelihoodDescription.from_dict(constraint='B->pi::f_++f_0@RBC+UKQCD:2015A')
        self.assertEqual(str(desc.constraint), 'B->pi::f_++f_0@RBC+UKQCD:2015A')


class ManualConstraintDescriptionTests(unittest.TestCase):

    def test_from_dict(self):
        info = {'type': 'Gaussian', 'observable': 'B->pilnu::BR'}
        desc = ManualConstraintDescription.from_dict(name='B->pilnu::BR@My:2026A', info=info)
        self.assertEqual(str(desc.name), 'B->pilnu::BR@My:2026A')
        self.assertEqual(desc.info, info)


class PyHFConstraintDescriptionTests(unittest.TestCase):

    def test_from_dict(self):
        desc = PyHFConstraintDescription.from_dict(file='workspace.json', parameter_map={'mu': 'EXP::mu'})
        self.assertEqual(desc.file, 'workspace.json')
        self.assertEqual(desc.parameter_map, {'mu': 'EXP::mu'})

    def test_default_parameter_map(self):
        desc = PyHFConstraintDescription.from_dict(file='workspace.json')
        self.assertEqual(desc.parameter_map, {})


class LikelihoodComponentTests(unittest.TestCase):

    def test_constraints(self):
        comp = LikelihoodComponent.from_dict(
            name='TH-pi',
            constraints=['B->pi::f_++f_0@RBC+UKQCD:2015A', 'B^0->pi^-l^+nu::BR@HFLAV:2019A'],
        )
        self.assertEqual(comp.name, 'TH-pi')
        self.assertEqual(len(comp.constraints), 2)
        self.assertIsInstance(comp.constraints[0], ConstraintLikelihoodDescription)

    def test_manual_constraints(self):
        comp = LikelihoodComponent.from_dict(
            name='manual',
            manual_constraints={'B->pilnu::BR@My:2026A': {'type': 'Gaussian'}},
        )
        self.assertEqual(len(comp.manual_constraints), 1)
        self.assertIsInstance(comp.manual_constraints[0], ManualConstraintDescription)

    def test_pyhf(self):
        comp = LikelihoodComponent.from_dict(name='pyhf', pyhf={'file': 'workspace.json'})
        self.assertIsInstance(comp.pyhf, PyHFConstraintDescription)

    def test_requires_at_least_one_source(self):
        comp = LikelihoodComponent.from_dict(name='empty')
        diagnostics = list(comp._diagnostics())
        self.assertEqual(len(diagnostics), 1)
        self.assertEqual(diagnostics[0].severity, Severity.ERROR)
        self.assertIn('must have at least one', diagnostics[0].message)


class PosteriorDescriptionTests(unittest.TestCase):

    def test_from_dict(self):
        desc = PosteriorDescription.from_dict(
            name='CKM-all',
            prior=['CKM', 'FF-pi'],
            likelihood=['TH-pi', 'EXP-pi'],
            global_options={'model': 'CKM'},
        )
        self.assertEqual(desc.name, 'CKM-all')
        self.assertEqual(desc.prior, ['CKM', 'FF-pi'])
        self.assertEqual(desc.likelihood, ['TH-pi', 'EXP-pi'])
        self.assertEqual(desc.global_options, {'model': 'CKM'})
        # fixed_parameters defaults to an empty dict
        self.assertEqual(desc.fixed_parameters, {})


class ObservableComponentTests(unittest.TestCase):

    def test_from_dict(self):
        desc = ObservableComponent.from_dict(
            name='B->pilnu::R_pi',
            latex='$R_{\\pi}$',
            unit='1',
            expression='<<B->pilnu::BR;l=tau>> / <<B->pilnu::BR;l=e>>',
        )
        self.assertEqual(desc.name, 'B->pilnu::R_pi')
        self.assertEqual(desc.unit, '1')
        self.assertTrue(desc.expression.startswith('<<'))
        # options defaults to an empty dict
        self.assertEqual(desc.options, {})


class PredictionObservableComponentTests(unittest.TestCase):

    def test_from_dict(self):
        desc = PredictionObservableComponent.from_dict(
            name='B->pilnu::dBR/dq2',
            kinematics={'q2': 1.0},
            options={'l': 'e'},
        )
        self.assertEqual(desc.name, 'B->pilnu::dBR/dq2')
        self.assertEqual(desc.kinematics, {'q2': 1.0})
        self.assertEqual(desc.options, {'l': 'e'})

    def test_defaults(self):
        desc = PredictionObservableComponent.from_dict(name='B_u->lnu::BR;l=e')
        self.assertEqual(desc.kinematics, {})
        self.assertEqual(desc.options, {})


class PredictionDescriptionTests(unittest.TestCase):

    def test_from_dict(self):
        desc = PredictionDescription.from_dict(
            name='leptonic-BR-CKM',
            global_options={'model': 'CKM'},
            observables=[
                {'name': 'B_u->lnu::BR;l=e'},
                {'name': 'B_u->lnu::BR;l=mu'},
            ],
        )
        self.assertEqual(desc.name, 'leptonic-BR-CKM')
        self.assertEqual(len(desc.observables), 2)
        # the observables are deserialized into PredictionObservableComponent instances
        self.assertIsInstance(desc.observables[0], PredictionObservableComponent)


class ParameterComponentTests(unittest.TestCase):

    def test_from_dict(self):
        desc = ParameterComponent.from_dict(
            name='ublnul::Re{cVL}',
            latex='$\\mathrm{Re}\\, \\mathcal{C}_{V_L}$',
            unit='1',
            central=1.0,
            min=-2.0,
            max=2.0,
            alias_of=['ubenue::Re{cVL}', 'ubmunumu::Re{cVL}'],
        )
        self.assertEqual(desc.name, 'ublnul::Re{cVL}')
        self.assertEqual(desc.central, 1.0)
        self.assertEqual(desc.min, -2.0)
        self.assertEqual(desc.max, 2.0)
        self.assertEqual(desc.alias_of, ['ubenue::Re{cVL}', 'ubmunumu::Re{cVL}'])

    def test_default_alias_of(self):
        desc = ParameterComponent.from_dict(
            name='my::parameter', latex='x', unit='1', central=0.0, min=-1.0, max=1.0)
        self.assertEqual(desc.alias_of, [])


class TaskComponentTests(unittest.TestCase):

    def test_from_dict(self):
        comp = TaskComponent.from_dict(task='corner-plot', arguments={'posterior': 'CKM-all', 'format': ['pdf']})
        self.assertEqual(comp.task, 'corner-plot')
        self.assertIn('posterior', comp.arguments)

    def test_argument_alias_mapping(self):
        "Check that CLI argument aliases are mapped to their internal names."
        comp = TaskComponent.from_dict(task='corner-plot', arguments={'posterior': 'CKM-all', 'F': ['pdf']})
        # the alias 'F' is mapped to 'format'
        self.assertIn('format', comp.arguments)
        self.assertNotIn('F', comp.arguments)

    def test_invalid_task(self):
        comp = TaskComponent.from_dict(task='not-a-real-task', arguments={})
        diagnostics = list(comp._diagnostics())
        self.assertEqual(len(diagnostics), 1)
        self.assertEqual(diagnostics[0].path, ('task',))
        self.assertEqual(diagnostics[0].severity, Severity.ERROR)

    def test_unknown_argument(self):
        comp = TaskComponent.from_dict(
            task='corner-plot',
            arguments={'posterior': 'CKM-all', 'not_an_argument': 1},
        )
        diagnostics = list(comp._diagnostics())
        self.assertTrue(any(
            d.path == ('arguments', 'not_an_argument') and d.severity == Severity.ERROR
            for d in diagnostics
        ))

    def test_missing_required_argument(self):
        comp = TaskComponent.from_dict(task='corner-plot', arguments={})
        diagnostics = list(comp._diagnostics())
        self.assertTrue(any(
            d.path == ('arguments', 'posterior') and d.severity == Severity.ERROR
            for d in diagnostics
        ))


class StepComponentTests(unittest.TestCase):

    def _tasks(self):
        return [{'task': 'corner-plot', 'arguments': {'posterior': 'CKM-all', 'format': ['pdf']}}]

    def test_from_dict(self):
        comp = StepComponent.from_dict(title='Corner plot', id='ckm.corner-plot', tasks=self._tasks())
        self.assertEqual(comp.title, 'Corner plot')
        self.assertEqual(comp.id, 'ckm.corner-plot')
        self.assertEqual(len(comp.tasks), 1)
        self.assertIsInstance(comp.tasks[0], TaskComponent)

    def test_invalid_id(self):
        slash = StepComponent.from_dict(title='t', id='has/slash', tasks=self._tasks())
        whitespace = StepComponent.from_dict(title='t', id='has space', tasks=self._tasks())
        self.assertTrue(any(d.path == ('id',) for d in slash._diagnostics()))
        self.assertTrue(any(d.path == ('id',) for d in whitespace._diagnostics()))

    def test_empty_tasks(self):
        comp = StepComponent.from_dict(title='t', id='no-tasks', tasks=[])
        diagnostics = list(comp._diagnostics())
        self.assertEqual(len(diagnostics), 1)
        self.assertEqual(diagnostics[0].path, ('tasks',))

    def test_task_and_default_arguments_are_checked_in_separate_phases(self):
        comp = StepComponent.from_dict(
            title='t',
            id='separate-argument-checks',
            tasks=[{
                'task': 'corner-plot',
                'arguments': {
                    'posterior': 'CKM-all',
                    'not_a_task_argument': 1,
                },
            }],
            default_arguments={
                'corner-plot': {
                    'not_a_default_argument': 2,
                },
            },
        )

        structural = list(comp.validate_structure())
        semantic = list(comp.validate_semantics(None))
        self.assertTrue(any('not_a_task_argument' in diagnostic.message for diagnostic in structural))
        self.assertFalse(any('not_a_task_argument' in diagnostic.message for diagnostic in semantic))
        self.assertTrue(any('not_a_default_argument' in diagnostic.message for diagnostic in semantic))


class MaskDescriptionTests(unittest.TestCase):

    def test_dispatch(self):
        "Check that the factory dispatches to the correct mask component."
        self.assertIsInstance(
            MaskDescription.from_dict(expression='1 > 0', name='B->pilnu::BR'),
            MaskExpressionComponent)
        self.assertIsInstance(
            MaskDescription.from_dict(mask_name='previous-mask'),
            MaskNamedComponent)
        self.assertIsInstance(
            MaskDescription.from_dict(name='B->pilnu::BR'),
            MaskObservableComponent)


class MaskExpressionComponentTests(unittest.TestCase):

    def test_from_dict(self):
        desc = MaskDescription.from_dict(expression='1 > 0', name='B->pilnu::BR')
        self.assertIsInstance(desc, MaskExpressionComponent)
        self.assertEqual(desc.expression, '1 > 0')
        self.assertEqual(desc.name, 'B->pilnu::BR')


class MaskObservableComponentTests(unittest.TestCase):

    def test_from_dict(self):
        desc = MaskDescription.from_dict(name='B->pilnu::BR')
        self.assertIsInstance(desc, MaskObservableComponent)
        self.assertEqual(desc.name, 'B->pilnu::BR')


class MaskNamedComponentTests(unittest.TestCase):

    def test_from_dict(self):
        desc = MaskDescription.from_dict(mask_name='previous-mask')
        self.assertIsInstance(desc, MaskNamedComponent)
        self.assertEqual(desc.mask_name, 'previous-mask')


class MaskComponentTests(unittest.TestCase):

    def test_from_dict(self):
        comp = MaskComponent.from_dict(
            name='my-mask',
            description=[{'name': 'B->pilnu::BR'}, {'mask_name': 'previous-mask'}],
        )
        self.assertEqual(comp.name, 'my-mask')
        self.assertEqual(comp.logical_combination, 'and')
        self.assertEqual(len(comp.description), 2)
        self.assertIsInstance(comp.description[0], MaskObservableComponent)
        self.assertIsInstance(comp.description[1], MaskNamedComponent)

    def test_invalid_logical_combination(self):
        comp = MaskComponent.from_dict(
            name='my-mask',
            description=[{'name': 'B->pilnu::BR'}],
            logical_combination='xor',
        )
        diagnostics = list(comp._diagnostics())
        self.assertEqual(len(diagnostics), 1)
        self.assertEqual(diagnostics[0].path, ('logical_combination',))
        self.assertEqual(diagnostics[0].severity, Severity.ERROR)


class AnalysisFileDescriptionTests(unittest.TestCase):

    def test_defaults(self):
        "An empty mapping yields the format-version baseline and empty sections."
        desc = AnalysisFileDescription.from_dict()
        self.assertEqual(desc.format_version, 1)
        self.assertIsInstance(desc.metadata, MetadataDescription)
        for section in (desc.priors, desc.likelihoods, desc.posteriors, desc.predictions,
                        desc.figures, desc.observables, desc.parameters, desc.steps, desc.masks):
            self.assertEqual(section, [])

    def test_from_dict(self):
        "Every recognized section is deserialized into its corresponding description objects."
        desc = AnalysisFileDescription.from_dict(
            format_version=1,
            metadata={'title': 'My Analysis', 'authors': [{'name': 'Jane Doe'}]},
            priors=[{'name': 'CKM', 'descriptions': [
                {'parameter': 'CKM::abs(V_ub)', 'min': 3.0e-3, 'max': 4.5e-3, 'type': 'uniform'}]}],
            likelihoods=[{'name': 'TH-pi', 'constraints': ['B->pi::f_++f_0@RBC+UKQCD:2015A']}],
            posteriors=[{'name': 'CKM-all', 'prior': ['CKM'], 'likelihood': ['TH-pi']}],
            predictions=[{'name': 'leptonic', 'observables': [{'name': 'B_u->lnu::BR;l=e'}]}],
            figures=[{'name': 'fig-A', 'type': 'single', 'plot': {'items': []}}],
            observables={'B->pilnu::R_pi': {'latex': '$R_\\pi$', 'unit': '1',
                                            'expression': '<<B->pilnu::BR;l=tau>> / <<B->pilnu::BR;l=e>>'}},
            parameters={'my::parameter': {'latex': 'x', 'unit': '1', 'central': 0.0, 'min': -1.0, 'max': 1.0}},
            steps=[{'title': 'Corner', 'id': 'ckm.corner', 'tasks': [
                {'task': 'corner-plot', 'arguments': {'posterior': 'CKM-all', 'format': ['pdf']}}]}],
            masks=[{'name': 'mask-A', 'description': [{'name': 'B->pilnu::BR'}]}],
        )
        self.assertEqual(desc.format_version, 1)
        self.assertEqual(desc.metadata.title, 'My Analysis')
        self.assertIsInstance(desc.metadata.authors[0], MetadataAuthorDescription)
        self.assertIsInstance(desc.priors[0], PriorComponent)
        self.assertIsInstance(desc.likelihoods[0], LikelihoodComponent)
        self.assertIsInstance(desc.posteriors[0], PosteriorDescription)
        self.assertIsInstance(desc.predictions[0], PredictionDescription)
        # figures are deserialized through the FigureFactory and carry their name
        self.assertEqual(desc.figures[0].name, 'fig-A')
        self.assertIsInstance(desc.steps[0], StepComponent)
        self.assertIsInstance(desc.masks[0], MaskComponent)

    def test_mapping_sections_inject_name(self):
        "The observables/parameters mappings become lists with the name injected into each entry."
        desc = AnalysisFileDescription.from_dict(
            observables={'B->pilnu::R_pi': {'latex': '$R_\\pi$', 'unit': '1', 'expression': '1'}},
            parameters={'my::parameter': {'latex': 'x', 'unit': '1', 'central': 0.0, 'min': -1.0, 'max': 1.0}},
        )
        self.assertEqual(len(desc.observables), 1)
        self.assertIsInstance(desc.observables[0], ObservableComponent)
        self.assertEqual(desc.observables[0].name, 'B->pilnu::R_pi')
        self.assertEqual(len(desc.parameters), 1)
        self.assertIsInstance(desc.parameters[0], ParameterComponent)
        self.assertEqual(desc.parameters[0].name, 'my::parameter')

    def test_unknown_top_level_key_is_rejected(self):
        "An unrecognized top-level key is rejected rather than silently ignored."
        desc = AnalysisFileDescription.from_dict(not_a_section=42)
        self.assertIsInstance(desc, InvalidComponent)
        diagnostics = list(desc.validate_structure())
        self.assertEqual(len(diagnostics), 1)
        self.assertEqual(diagnostics[0].path, ('not_a_section',))

    def test_validate_recurses_and_prefixes_raw_child_segments(self):
        desc = AnalysisFileDescription.from_dict(
            priors=[{
                'name': 'FF',
                'descriptions': [{'type': 'uniform', 'parameter': 'x', 'min': 0.0}],
            }],
            likelihoods=[{'name': 'empty'}],
            steps=[{
                'title': 'Broken step',
                'id': 'broken',
                'tasks': [{'arguments': {}}],
            }],
        )

        diagnostics = list(desc.validate_structure())
        self.assertIn(('priors', 'FF', 'descriptions', 0, 'max'), [d.path for d in diagnostics])
        self.assertIn(('likelihoods', 'empty'), [d.path for d in diagnostics])
        self.assertIn(('steps', 'broken', 'tasks', 0, 'task'), [d.path for d in diagnostics])

    def test_raw_child_segments_are_not_serialized(self):
        desc = AnalysisFileDescription.from_dict(
            priors=[{
                'name': 'P',
                'descriptions': [
                    {'type': 'uniform', 'parameter': 'x', 'min': 0.0, 'max': 1.0},
                ],
            }],
        )

        serialized = asdict(desc)
        self.assertFalse(any(key.startswith('_') for key in serialized))
        self.assertFalse(any(key.startswith('_') for key in serialized['priors'][0]))

    def test_validate_structure_checks_whole_file(self):
        desc = AnalysisFileDescription.from_dict(
            priors=[{
                'name': 'P',
                'descriptions': [
                    {'type': 'uniform', 'parameter': 'x', 'min': 0.0, 'max': 1.0},
                ],
            }],
            likelihoods=[{'name': 'L', 'constraints': ['x::constraint']}],
            posteriors=[{
                'name': 'posterior',
                'prior': ['missing-prior'],
                'likelihood': ['missing-likelihood'],
            }],
            steps=[
                {'title': 'First', 'id': 'duplicate', 'tasks': []},
                {'title': 'Second', 'id': 'duplicate', 'tasks': []},
            ],
            masks=[
                {'name': 'duplicate-mask', 'description': [{'mask_name': 'missing-mask'}]},
                {'name': 'duplicate-mask', 'description': []},
            ],
        )

        diagnostics = list(desc.validate_structure())
        paths = [diagnostic.path for diagnostic in diagnostics]
        self.assertIn(('posteriors', 'posterior', 'prior', 0), paths)
        self.assertIn(('posteriors', 'posterior', 'likelihood', 0), paths)
        self.assertIn(('steps', 'duplicate', 'id'), paths)
        self.assertIn(('masks', 'duplicate-mask', 'name'), paths)
        self.assertIn(
            ('masks', 'duplicate-mask', 'description', 0, 'mask_name'),
            paths,
        )

    def test_validate_structure_checks_section_presence(self):
        desc = AnalysisFileDescription.from_dict()
        diagnostics = list(desc.validate_structure())

        by_path = {diagnostic.path: diagnostic for diagnostic in diagnostics}
        self.assertEqual(by_path[('priors',)].severity, Severity.ERROR)
        self.assertEqual(by_path[('likelihoods',)].severity, Severity.WARNING)
        self.assertEqual(by_path[('posteriors',)].severity, Severity.ERROR)

    def test_fixed_parameter_varied_by_own_prior(self):
        "A posterior that fixes a parameter one of its own priors varies is contradictory."
        desc = AnalysisFileDescription.from_dict(
            priors=[
                {'name': 'varies-mB', 'descriptions': [
                    {'type': 'uniform', 'parameter': 'mass::B_d', 'min': 5.2, 'max': 5.4}]},
                {'name': 'varies-two', 'descriptions': [
                    {'type': 'transform', 'parameters': ['mass::B_u', 'mass::B_s'],
                     'shift': [0.0, 0.0], 'transform': [[1.0, 0.0], [0.0, 1.0]],
                     'min': [5.0, 5.0], 'max': [5.5, 5.5]}]},
                {'name': 'varies-elsewhere', 'descriptions': [
                    {'type': 'uniform', 'parameter': 'mass::e', 'min': 0.0, 'max': 1.0}]},
            ],
            likelihoods=[{'name': 'll', 'constraints': ['B->pi::f_++f_0@RBC+UKQCD:2015A']}],
            posteriors=[{
                'name': 'post',
                'prior': ['varies-mB', 'varies-two'],
                'likelihood': ['ll'],
                'fixed_parameters': {
                    'mass::B_d': 5.28,   # varied by 'varies-mB'          -> reported
                    'mass::B_s': 5.37,   # varied by 'varies-two'         -> reported
                    'mass::e': 0.0005,   # varied by a prior this posterior does not use -> silent
                    'mass::tau': 1.777,  # varied by no prior at all      -> silent
                },
            }],
        )

        reported = {
            diagnostic.path[-1]: diagnostic
            for diagnostic in desc.validate_semantics(ValidationContext(desc))
            if 'fixed by posterior' in diagnostic.message
        }
        self.assertEqual(set(reported), {'mass::B_d', 'mass::B_s'})
        for diagnostic in reported.values():
            # a warning, not an error: an error would abort the deep phase and hide every deep
            # finding elsewhere in the file
            self.assertEqual(diagnostic.severity, Severity.WARNING)
        self.assertIn("varied by its prior 'varies-mB'", reported['mass::B_d'].message)
        self.assertIn("varied by its prior 'varies-two'", reported['mass::B_s'].message)


if __name__ == '__main__':
    unittest.main(verbosity=5)
