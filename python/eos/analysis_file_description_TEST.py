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
import os

import eos

from eos.analysis_file_description import (
    AnalysisFileContext,
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
)


class AnalysisFileContextTests(unittest.TestCase):

    def test_data_path(self):
        "Check that data_path resolves a relative path against the base directory."
        ctx = AnalysisFileContext()
        path = ctx.data_path('some/relative/path')
        self.assertTrue(os.path.isabs(path))
        self.assertTrue(path.endswith(os.path.join('some', 'relative', 'path')))

    def test_invalid_base_directory(self):
        "Check that a non-existent or non-directory base directory is rejected."
        with self.assertRaises(ValueError):
            AnalysisFileContext(base_directory='/this/does/not/exist/at/all')
        # a path that exists but is a file, not a directory
        with self.assertRaises(ValueError):
            AnalysisFileContext(base_directory=__file__)


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


class MetadataDescriptionTests(unittest.TestCase):

    def test_from_dict(self):
        md = MetadataDescription.from_dict(
            title='My Analysis',
            id='my-analysis',
            authors=[{'name': 'Jane Doe'}, {'name': 'John Roe', 'email': 'john@example.org'}],
        )
        self.assertEqual(md.title, 'My Analysis')
        self.assertEqual(md.id, 'my-analysis')
        # the authors are deserialized into MetadataAuthorDescription instances
        self.assertEqual(len(md.authors), 2)
        self.assertIsInstance(md.authors[0], MetadataAuthorDescription)
        self.assertEqual(md.authors[1].email, 'john@example.org')

    def test_defaults(self):
        md = MetadataDescription.from_dict()
        self.assertEqual(md.title, '')
        self.assertEqual(md.id, '')
        self.assertEqual(md.authors, [])


class PriorDescriptionTests(unittest.TestCase):

    def test_dispatch(self):
        "Check that the factory dispatches each prior type to the correct class."
        self.assertIsInstance(
            PriorDescription.from_dict(constraint='B->D::f_++f_0@HPQCD:2015A'),
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
        with self.assertRaises(ValueError):
            PriorDescription.from_dict(parameter='p', type='not-a-type')
        with self.assertRaises(ValueError):
            PriorDescription.from_dict()


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

        with self.assertRaises(ValueError):
            PriorDescription.from_dict(parameter='mass::b(MSbar)', type='scale', min=0.0, max=1.0)


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
        with self.assertRaises(ValueError):
            LikelihoodComponent.from_dict(name='empty')


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
        with self.assertRaises(ValueError):
            TaskComponent.from_dict(task='not-a-real-task', arguments={})

    def test_unknown_argument(self):
        with self.assertRaises(ValueError):
            TaskComponent.from_dict(task='corner-plot', arguments={'posterior': 'CKM-all', 'not_an_argument': 1})


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
        with self.assertRaises(ValueError):
            StepComponent.from_dict(title='t', id='has/slash', tasks=self._tasks())
        with self.assertRaises(ValueError):
            StepComponent.from_dict(title='t', id='has space', tasks=self._tasks())

    def test_empty_tasks(self):
        with self.assertRaises(ValueError):
            StepComponent.from_dict(title='t', id='no-tasks', tasks=[])


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
        with self.assertRaises(ValueError):
            MaskComponent.from_dict(
                name='my-mask',
                description=[{'name': 'B->pilnu::BR'}],
                logical_combination='xor',
            )


if __name__ == '__main__':
    unittest.main(verbosity=5)
