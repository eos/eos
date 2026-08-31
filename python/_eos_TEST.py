#!/usr/bin/env python3
# vim: set sw=4 sts=4 et tw=120 :

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

import _eos
import atexit
import eos
import inspect
import os
import re
import tempfile
import unittest

# The module-level name under which free functions are registered below.
MODULE = '<module>'

# Members that Boost.Python adds to every exported class and that carry no binding of our own.
IGNORED_MEMBERS = { '__init__', '__doc__', '__module__', '__reduce__', '__instance_size__', '__slots__' }


# Shared fixtures. Each is created once and handed to the factories below.
_PARAMETERS  = eos.Parameters.Defaults()
_OPTIONS     = eos.Options({ 'form-factors': 'BSZ2015' })
_KINEMATICS  = eos.Kinematics(q2=1.0)
_OBSERVABLE  = eos.Observable.make('B->D::f_+(q2)', _PARAMETERS, _KINEMATICS, _OPTIONS)
_MODEL       = eos.Model.make('SM', _PARAMETERS, _OPTIONS)
_CONSTRAINT  = eos.Constraint.make('B->D::f_++f_0@HPQCD:2015A', _OPTIONS)
_SIGNAL_PDF  = _eos._SignalPDF.make('TestLegendre1D::P(z)', _PARAMETERS, eos.Kinematics(z=2.0, z_min=0.0, z_max=4.0), eos.Options())

# The Wilson-polynomial bindings need an observable that depends on the b -> s ll Wilson coefficients.
_WET_OPTIONS      = eos.Options({ 'model': 'WET', 'l': 'mu', 'tag': 'GvDV2020', 'form-factors': 'BSZ2015' })
_WET_OBSERVABLE   = eos.Observable.make('B->Kll::BR', _PARAMETERS, eos.Kinematics(q2_min=1.0, q2_max=6.0), _WET_OPTIONS)
_WET_COEFFICIENTS = [ 'b->smumu::Re{c9}', 'b->smumu::Re{c10}' ]

# Only a handful of constraint entries constrain parameters rather than observables, and only those
# can be turned into a prior.
_CONSTRAINT_ENTRY_NAME = 'B->K::FormFactors[parametric,LCSRLattice]@GKvD:2018A'


def _observable_cache():
    cache = _eos.ObservableCache(_PARAMETERS)
    cache.add(_OBSERVABLE)
    cache.update()
    return cache


def _log_likelihood():
    result = _eos.LogLikelihood(_PARAMETERS)
    result.add(_CONSTRAINT)
    result.evaluate() # updates the underlying cache, so that the blocks can be evaluated
    return result


def _log_posterior():
    result = _eos.LogPosterior(_log_likelihood())
    result.add(_eos.LogPrior.Flat(_PARAMETERS, 'mass::b(MSbar)', 4.0, 4.4), False)
    return result


_UNBINNED_GRID         = [ eos.Kinematics(z=z, z_min=0.0, z_max=4.0) for z in (0.5, 1.5, 2.5, 3.5) ]
_UNBINNED_RESOLUTION   = [ 0.0, 1.0, 0.0, 0.0 ]
_UNBINNED_OBSERVATIONS = [ eos.Kinematics(z=2.0, z_min=0.0, z_max=4.0) ]


def _wilson_coefficients():
    return _MODEL.wilson_coefficients_b_to_s(4.2, 'mu', False)


def _parameter_file():
    handle, path = tempfile.mkstemp(suffix='.yaml')
    with os.fdopen(handle, 'w') as f:
        f.write('mass::b(MSbar):\n    central: 4.18\n    min: 4.0\n    max: 4.4\n')

    atexit.register(lambda: os.path.exists(path) and os.remove(path))

    return path


# For each class exported by the _eos module, a callable returning an instance of that class.
BINDINGS_FACTORIES = {
    'BToSWilsonCoefficients':   _wilson_coefficients,
    'Constraint':               lambda: _CONSTRAINT,
    'ConstraintEntry':          lambda: eos.Constraints()[_CONSTRAINT_ENTRY_NAME],
    'ExpressionReferences':     lambda: _eos.analyze_expression('<<B->D::f_+(q2)>>'),
    'GoodnessOfFit':            lambda: _eos.GoodnessOfFit(_log_posterior()),
    'Kinematics':               lambda: _KINEMATICS,
    'KinematicVariable':        lambda: _KINEMATICS['q2'],
    'LogLikelihood':            _log_likelihood,
    'LogLikelihoodBlock':       lambda: list(_CONSTRAINT.blocks())[0],
    'LogPosterior':             _log_posterior,
    'LogPrior':                 lambda: _eos.LogPrior.Flat(_PARAMETERS, 'mass::b(MSbar)', 4.0, 4.4),
    'Model':                    lambda: _MODEL,
    'Observable':               lambda: _OBSERVABLE,
    'ObservableCache':          _observable_cache,
    'ObservableEntry':          lambda: eos.Observables()['B->D::f_+(q2)'],
    'ObservableGroup':          lambda: list(list(eos.Observables().sections())[0])[0],
    'ObservableId':             lambda: _eos.ObservableCache(_PARAMETERS).add(_OBSERVABLE),
    'ObservableSection':        lambda: list(eos.Observables().sections())[0],
    'Options':                  lambda: _OPTIONS,
    'OptionSpecification':      lambda: list(eos.Observables()['B->D::f_+(q2)'].options())[0],
    'Parameter':                lambda: _PARAMETERS['mass::b(MSbar)'],
    'ParameterGroup':           lambda: list(list(_PARAMETERS.sections())[0])[0],
    'ParameterSection':         lambda: list(_PARAMETERS.sections())[0],
    'ParameterUser':            lambda: _OBSERVABLE,
    'QualifiedName':            lambda: eos.QualifiedName('B->D::f_+(q2);form-factors=BSZ2015'),
    'Reference':                lambda: eos.References()['BSZ:2015A'],
    'ReferenceName':            lambda: eos.ReferenceName('BSZ:2015A'),
    'ReferenceUser':            lambda: _OBSERVABLE,
    'SignalPDFEntry':           lambda: eos.SignalPDFs()['B->Dlnu::P(q2)'],
    'SignalPDFGroup':           lambda: list(list(eos.SignalPDFs().sections())[0])[0],
    'SignalPDFSection':         lambda: list(eos.SignalPDFs().sections())[0],
    'Unit':                     lambda: eos.Unit.GeV(),
    '_Constraints':             lambda: eos.Constraints(),
    '_Observables':             lambda: eos.Observables(),
    '_Parameters':              lambda: _PARAMETERS,
    '_References':              lambda: eos.References(),
    '_SignalPDF':               lambda: _SIGNAL_PDF,
    '_SignalPDFs':              lambda: eos.SignalPDFs(),
    'qnpName':                  lambda: _eos.qnpName('f_+(q2)'),
    'qnpOptionKey':             lambda: _eos.qnpOptionKey('form-factors'),
    'qnpOptionValue':           lambda: _eos.qnpOptionValue('BSZ2015'),
    'qnpPrefix':                lambda: _eos.qnpPrefix('B->D'),
    'qnpSuffix':                lambda: _eos.qnpSuffix('FOO'),
    'rnpIndex':                 lambda: _eos.rnpIndex('A'),
    'rnpName':                  lambda: _eos.rnpName('BSZ'),
    'rnpYear':                  lambda: _eos.rnpYear('2015'),
    'test_statisticsChiSquare': lambda: list(_eos.GoodnessOfFit(_log_posterior()))[0][1],
}


# For each member of each exported class, either a tuple of arguments to call it with, or a
# callable that receives an instance of the class and exercises the member itself.
BINDINGS_TESTS = {
    ('BToSWilsonCoefficients', 'c1'):            (),
    ('BToSWilsonCoefficients', 'c2'):            (),

    ('Constraint', 'blocks'):                    (),
    ('Constraint', 'make'):                      ('B->D::f_++f_0@HPQCD:2015A', _OPTIONS),
    ('Constraint', 'name'):                      (),
    ('Constraint', 'observables'):               (),

    ('ConstraintEntry', 'deserialize'):          lambda e: _eos.ConstraintEntry.deserialize(e.name(), e.serialize()),
    ('ConstraintEntry', 'make'):                 (_CONSTRAINT_ENTRY_NAME, _OPTIONS),
    ('ConstraintEntry', 'make_prior'):           (_PARAMETERS, _OPTIONS),
    ('ConstraintEntry', 'name'):                 (),
    ('ConstraintEntry', 'observables'):          (),
    ('ConstraintEntry', 'references'):           (),
    ('ConstraintEntry', 'serialize'):            (),
    ('ConstraintEntry', 'type'):                 (),

    ('ExpressionReferences', 'observables'):     (),
    ('ExpressionReferences', 'parameters'):      (),

    ('GoodnessOfFit', '__iter__'):               lambda g: list(g),
    ('GoodnessOfFit', 'total_chi_square'):       (),
    ('GoodnessOfFit', 'total_degrees_of_freedom'): (),

    ('Kinematics', '__add__'):                   lambda k: k + eos.Kinematics(q2_min=0.0),
    ('Kinematics', '__contains__'):              lambda k: 'q2' in k,
    ('Kinematics', '__getitem__'):               ('q2',),
    ('Kinematics', '__iter__'):                  lambda k: list(k),
    ('Kinematics', '__str__'):                   lambda k: str(k),
    ('Kinematics', 'declare'):                   ('q2_test', 1.0),

    ('KinematicVariable', '__float__'):          lambda v: float(v),
    ('KinematicVariable', 'evaluate'):           (),
    ('KinematicVariable', 'name'):               (),
    ('KinematicVariable', 'set'):                (1.0,),

    ('LogLikelihood', '__iter__'):               lambda l: list(l),
    ('LogLikelihood', 'add'):                    (_CONSTRAINT,),
    ('LogLikelihood', 'evaluate'):               (),
    ('LogLikelihood', 'observable_cache'):       (),

    ('LogLikelihoodBlock', 'External'):          lambda b: _eos.LogLikelihoodBlock.External(_observable_cache(), _external_block_factory),
    ('LogLikelihoodBlock', 'Unbinned1D'):        lambda b: _eos.LogLikelihoodBlock.Unbinned1D(
                                                     _observable_cache(), 'TestLegendre1D::P(z)', _UNBINNED_GRID,
                                                     eos.Options(), _UNBINNED_RESOLUTION, _UNBINNED_OBSERVATIONS),
    ('LogLikelihoodBlock', '_Unbinned1D'):       lambda b: _eos.LogLikelihoodBlock._Unbinned1D(
                                                     _observable_cache(), 'TestLegendre1D::P(z)', _UNBINNED_GRID,
                                                     eos.Options(), [ 1.0, 0.0, 0.0, 0.0 ], _UNBINNED_OBSERVATIONS),
    ('LogLikelihoodBlock', '__str__'):           lambda b: str(b),
    ('LogLikelihoodBlock', 'evaluate'):          (),
    ('LogLikelihoodBlock', 'number_of_observations'): (),

    ('LogPosterior', 'add'):                     lambda p: p.add(_eos.LogPrior.Flat(_PARAMETERS, 'mass::c', 1.0, 1.5), False),
    ('LogPosterior', 'evaluate'):                (),
    ('LogPosterior', 'log_likelihood'):          (),
    ('LogPosterior', 'log_priors'):              lambda p: list(p.log_priors()),

    ('LogPrior', 'CurtailedGauss'):              (_PARAMETERS, 'mass::b(MSbar)', 4.0, 4.4, 4.15, 4.18, 4.21),
    ('LogPrior', 'External'):                    lambda p: _eos.LogPrior.External(_PARAMETERS, _external_prior_factory),
    ('LogPrior', 'Flat'):                        (_PARAMETERS, 'mass::b(MSbar)', 4.0, 4.4),
    ('LogPrior', 'Gaussian'):                    (_PARAMETERS, 'mass::b(MSbar)', 4.18, 0.03),
    ('LogPrior', 'Poisson'):                     (_PARAMETERS, 'mass::b(MSbar)', 3.0),
    ('LogPrior', 'Scale'):                       (_PARAMETERS, 'sbmumu::mu', 2.0, 8.0, 4.18, 2.0),
    ('LogPrior', 'Transform'):                   (_PARAMETERS, ['mass::b(MSbar)'], [0.0], [[1.0]], [4.0], [4.4]),
    ('LogPrior', 'Uniform'):                     (_PARAMETERS, 'mass::b(MSbar)', 4.0, 4.4),
    ('LogPrior', 'compute_cdf'):                 (),
    ('LogPrior', 'evaluate'):                    (),
    ('LogPrior', 'sample'):                      (),
    ('LogPrior', 'varied_parameters'):           lambda p: list(p.varied_parameters()),

    ('Model', 'alpha_s'):                        (4.2,),
    ('Model', 'ckm_cb'):                         (),
    ('Model', 'ckm_cd'):                         (),
    ('Model', 'ckm_cs'):                         (),
    ('Model', 'ckm_tb'):                         (),
    ('Model', 'ckm_td'):                         (),
    ('Model', 'ckm_ts'):                         (),
    ('Model', 'ckm_ub'):                         (),
    ('Model', 'ckm_ud'):                         (),
    ('Model', 'ckm_us'):                         (),
    ('Model', 'm_b_kin'):                        (1.0,),
    ('Model', 'm_b_msbar'):                      (4.2,),
    ('Model', 'm_b_pole'):                       (),
    ('Model', 'm_c_kin'):                        (1.0,),
    ('Model', 'm_c_msbar'):                      (1.5,),
    ('Model', 'm_c_pole'):                       (),
    ('Model', 'm_s_msbar'):                      (2.0,),
    ('Model', 'm_t_msbar'):                      (170.0,),
    ('Model', 'm_t_pole'):                       (),
    ('Model', 'm_ud_msbar'):                     (2.0,),
    ('Model', 'make'):                           ('SM', _PARAMETERS, _OPTIONS),
    ('Model', 'wilson_coefficients_b_to_s'):     (4.2, 'mu', False),

    ('Observable', 'evaluate'):                  (),
    ('Observable', 'kinematics'):                (),
    ('Observable', 'make'):                      ('B->D::f_+(q2)', _PARAMETERS, _KINEMATICS, _OPTIONS),
    ('Observable', 'name'):                      (),
    ('Observable', 'options'):                   (),
    ('Observable', 'parameters'):                (),

    ('ObservableCache', '__getitem__'):          lambda c: c[c.add(_OBSERVABLE)],
    ('ObservableCache', '__iter__'):             lambda c: list(c),
    ('ObservableCache', 'add'):                  (_OBSERVABLE,),
    ('ObservableCache', 'parameters'):           (),
    ('ObservableCache', 'update'):               (),

    ('ObservableEntry', 'kinematic_variables'):  (),
    ('ObservableEntry', 'latex'):                (),
    ('ObservableEntry', 'name'):                 (),
    ('ObservableEntry', 'options'):              (),
    ('ObservableEntry', 'unit'):                 (),

    ('ObservableGroup', '__iter__'):             lambda g: list(g),
    ('ObservableGroup', 'description'):          (),
    ('ObservableGroup', 'name'):                 (),

    ('ObservableId', '__eq__'):                  lambda i: i == i,
    ('ObservableId', '__hash__'):                lambda i: hash(i),
    ('ObservableId', '__repr__'):                lambda i: repr(i),
    ('ObservableId', 'value'):                   (),

    ('ObservableSection', '__iter__'):           lambda s: list(s),
    ('ObservableSection', 'description'):        (),
    ('ObservableSection', 'name'):               (),

    ('OptionSpecification', 'allowed_values'):   (),
    ('OptionSpecification', 'default_value'):    (),
    ('OptionSpecification', 'key'):              (),

    ('Options', '__contains__'):                 lambda o: 'form-factors' in o,
    ('Options', '__eq__'):                       lambda o: o == o,
    ('Options', '__getitem__'):                  ('form-factors',),
    ('Options', '__iter__'):                     lambda o: list(o),
    ('Options', '__str__'):                      lambda o: str(o),
    ('Options', 'declare'):                      ('l', 'mu'),

    ('Parameter', '__float__'):                  lambda p: float(p),
    ('Parameter', 'central'):                    (),
    ('Parameter', 'evaluate'):                   (),
    ('Parameter', 'evaluate_generator'):         (),
    ('Parameter', 'latex'):                      (),
    ('Parameter', 'max'):                        (),
    ('Parameter', 'min'):                        (),
    ('Parameter', 'name'):                       (),
    ('Parameter', 'set'):                        (4.18,),
    ('Parameter', 'set_generator'):              (0.5,),
    ('Parameter', 'set_max'):                    (4.4,),
    ('Parameter', 'set_min'):                    (4.0,),
    ('Parameter', 'unit'):                       (),

    ('ParameterGroup', '__iter__'):              lambda g: list(g),
    ('ParameterGroup', 'description'):           (),
    ('ParameterGroup', 'name'):                  (),

    ('ParameterSection', '__iter__'):            lambda s: list(s),
    ('ParameterSection', 'description'):         (),
    ('ParameterSection', 'name'):                (),

    ('ParameterUser', 'used_parameter_ids'):     lambda u: list(u.used_parameter_ids()),

    ('QualifiedName', '__eq__'):                 lambda q: q == q,
    ('QualifiedName', '__lt__'):                 lambda q: q < q,
    ('QualifiedName', '__ne__'):                 lambda q: q != q,
    ('QualifiedName', '__repr__'):               lambda q: repr(q),
    ('QualifiedName', '__str__'):                lambda q: str(q),
    ('QualifiedName', 'full'):                   (),
    ('QualifiedName', 'name_part'):              (),
    ('QualifiedName', 'options_part'):           (),
    ('QualifiedName', 'prefix_part'):            (),
    ('QualifiedName', 'suffix_part'):            (),

    ('Reference', 'authors'):                    (),
    ('Reference', 'eprint_archive'):             (),
    ('Reference', 'eprint_id'):                  (),
    ('Reference', 'inspire_id'):                 (),
    ('Reference', 'name'):                       (),
    ('Reference', 'title'):                      (),

    ('ReferenceName', '__eq__'):                 lambda r: r == r,
    ('ReferenceName', '__lt__'):                 lambda r: r < r,
    ('ReferenceName', '__ne__'):                 lambda r: r != r,
    ('ReferenceName', '__str__'):                lambda r: str(r),
    ('ReferenceName', 'index_part'):             (),
    ('ReferenceName', 'name_part'):              (),
    ('ReferenceName', 'year_part'):              (),

    ('ReferenceUser', 'references'):             lambda u: list(u.references()),

    ('SignalPDFEntry', 'description'):           (),
    ('SignalPDFEntry', 'name'):                  (),

    ('SignalPDFGroup', '__iter__'):              lambda g: list(g),
    ('SignalPDFGroup', 'description'):           (),
    ('SignalPDFGroup', 'name'):                  (),

    ('SignalPDFSection', '__iter__'):            lambda s: list(s),
    ('SignalPDFSection', 'description'):         (),
    ('SignalPDFSection', 'name'):                (),

    ('Unit', 'Femtometer2'):                     (),
    ('Unit', 'GeV'):                             (),
    ('Unit', 'GeV2'):                            (),
    ('Unit', 'GeV3'):                            (),
    ('Unit', 'GeVSecond'):                       (),
    ('Unit', 'InverseGeV'):                      (),
    ('Unit', 'InverseGeV2'):                     (),
    ('Unit', 'InverseGeV4'):                     (),
    ('Unit', 'InversePicoSecond'):               (),
    ('Unit', 'InverseSecond'):                   (),
    ('Unit', 'Second'):                          (),
    ('Unit', 'Undefined'):                       (),
    ('Unit', 'Unity'):                           (),
    ('Unit', '__eq__'):                          lambda u: u == u,
    ('Unit', '__str__'):                         lambda u: str(u),
    ('Unit', 'latex'):                           (),

    ('_Constraints', '__contains__'):            lambda c: 'B->D::f_++f_0@HPQCD:2015A' in c,
    ('_Constraints', '__getitem__'):             ('B->D::f_++f_0@HPQCD:2015A',),
    ('_Constraints', '__iter__'):                lambda c: list(c),
    ('_Constraints', 'insert'):                  lambda c: c.insert('B->D::f_+@Test:2026A',
                                                     eos.Constraints()[_CONSTRAINT_ENTRY_NAME].serialize()),

    ('_Observables', '__contains__'):            lambda o: 'B->D::f_+(q2)' in o,
    ('_Observables', '__getitem__'):             ('B->D::f_+(q2)',),
    ('_Observables', '__iter__'):                lambda o: list(o),
    ('_Observables', 'insert'):                  lambda o: o.insert('B->D::f_+_test(q2)', '', eos.Unit.Unity(),
                                                     eos.Options(), '<<B->D::f_+(q2)>>'),
    ('_Observables', 'sections'):                lambda o: list(o.sections()),

    ('_Parameters', 'Defaults'):                 (),
    ('_Parameters', '__contains__'):             lambda p: 'mass::b(MSbar)' in p,
    ('_Parameters', '__getitem__'):              ('mass::b(MSbar)',),
    ('_Parameters', '__iter__'):                 lambda p: list(p),
    ('_Parameters', 'by_id'):                    (0,),
    ('_Parameters', 'declare'):                  ('test::redirect-source', '', eos.Unit.Unity(), 1.0, 0.0, 2.0),
    ('_Parameters', 'declare_and_insert'):       ('test::declare_and_insert', '', eos.Unit.Unity(), 1.0, 0.0, 2.0),
    ('_Parameters', 'has'):                      ('mass::b(MSbar)',),
    ('_Parameters', 'override_from_file'):       lambda p: p.override_from_file(_parameter_file()),
    ('_Parameters', 'redirect'):                 lambda p: _eos._Parameters.redirect('test::redirect-source',
                                                     _eos._Parameters.declare('test::redirect-target', '', eos.Unit.Unity(), 2.0, 0.0, 4.0)),
    ('_Parameters', 'sections'):                 lambda p: list(p.sections()),
    ('_Parameters', 'set'):                      ('mass::b(MSbar)', 4.18),

    ('_References', '__contains__'):             lambda r: 'BSZ:2015A' in r,
    ('_References', '__getitem__'):              ('BSZ:2015A',),
    ('_References', '__iter__'):                 lambda r: list(r),

    ('_SignalPDF', 'evaluate'):                  (),
    ('_SignalPDF', 'kinematics'):                (),
    ('_SignalPDF', 'make'):                      ('TestLegendre1D::P(z)', _PARAMETERS,
                                                  eos.Kinematics(z=2.0, z_min=0.0, z_max=4.0), eos.Options()),
    ('_SignalPDF', 'name'):                      (),
    ('_SignalPDF', 'normalization'):             (),
    ('_SignalPDF', 'options'):                   (),
    ('_SignalPDF', 'parameters'):                (),

    ('_SignalPDFs', '__contains__'):             lambda s: 'B->Dlnu::P(q2)' in s,
    ('_SignalPDFs', '__getitem__'):              ('B->Dlnu::P(q2)',),
    ('_SignalPDFs', '__iter__'):                 lambda s: list(s),
    ('_SignalPDFs', 'insert'):                   lambda s: s.insert('B->Dlnu::P_test(q2)', '', eos.Options(),
                                                     'B->Dlnu::dBR/dq2', ['q2'], 'B->Dlnu::BR', ['q2_min', 'q2_max']),
    ('_SignalPDFs', 'sections'):                 lambda s: list(s.sections()),

    ('qnpName', '__lt__'):                       lambda n: n < n,
    ('qnpName', '__repr__'):                     lambda n: repr(n),
    ('qnpName', '__str__'):                      lambda n: str(n),

    ('qnpOptionKey', '__repr__'):                lambda k: repr(k),
    ('qnpOptionKey', '__str__'):                 lambda k: str(k),

    ('qnpOptionValue', '__eq__'):                lambda v: v == v,
    ('qnpOptionValue', '__lt__'):                lambda v: v < v,
    ('qnpOptionValue', '__repr__'):              lambda v: repr(v),
    ('qnpOptionValue', '__str__'):               lambda v: str(v),

    ('qnpPrefix', '__lt__'):                     lambda p: p < p,
    ('qnpPrefix', '__repr__'):                   lambda p: repr(p),
    ('qnpPrefix', '__str__'):                    lambda p: str(p),

    ('qnpSuffix', '__lt__'):                     lambda s: s < s,
    ('qnpSuffix', '__repr__'):                   lambda s: repr(s),
    ('qnpSuffix', '__str__'):                    lambda s: str(s),

    ('rnpIndex', '__lt__'):                      lambda i: i < i,
    ('rnpIndex', '__repr__'):                    lambda i: repr(i),
    ('rnpIndex', '__str__'):                     lambda i: str(i),

    ('rnpName', '__lt__'):                       lambda n: n < n,
    ('rnpName', '__repr__'):                     lambda n: repr(n),
    ('rnpName', '__str__'):                      lambda n: str(n),

    ('rnpYear', '__lt__'):                       lambda y: y < y,
    ('rnpYear', '__repr__'):                     lambda y: repr(y),
    ('rnpYear', '__str__'):                      lambda y: str(y),

    ('test_statisticsChiSquare', 'chi2'):        (),
    ('test_statisticsChiSquare', 'dof'):         (),
    ('test_statisticsChiSquare', 'signed_chi'):  (),

    (MODULE, '_emit_native_log'):                ('eos._eos_TEST', _eos._NativeLogLevel.DEBUG, 'test message'),
    (MODULE, '_register_log_callback'):          lambda _: _eos._register_log_callback(lambda *args: None),
    (MODULE, '_release_python_observables'):     (),
    (MODULE, '_set_native_log_level'):           (_eos._NativeLogLevel.SILENT,),
    (MODULE, 'analyze_expression'):              ('<<B->D::f_+(q2)>>',),
    (MODULE, 'compute_wilson_polynomial_coefficients'): (_WET_OBSERVABLE, _WET_COEFFICIENTS),
    (MODULE, 'delta_c7'):                        lambda _: _eos.delta_c7(complex(1.0, 0.0), 4.2, 0.22, 1.27, 4.18, _wilson_coefficients(), True),
    (MODULE, 'delta_c7_Qc'):                     lambda _: _eos.delta_c7_Qc(complex(1.0, 0.0), 4.2, 0.22, 1.27, 4.18, _wilson_coefficients(), True),
    (MODULE, 'delta_c9'):                        lambda _: _eos.delta_c9(complex(1.0, 0.0), 4.2, 0.22, 1.27, 4.18, _wilson_coefficients(), True),
    (MODULE, 'delta_c9_Qc'):                     lambda _: _eos.delta_c9_Qc(complex(1.0, 0.0), 4.2, 0.22, 1.27, 4.18, _wilson_coefficients(), True),
    (MODULE, 'make_wilson_polynomial_observable'): ('test::wilson-polynomial', _WET_OBSERVABLE, _WET_COEFFICIENTS),
    (MODULE, 'make_wilson_polynomial_ratio_observable'): ('test::wilson-polynomial-ratio', _WET_OBSERVABLE,
                                                     _WET_OBSERVABLE, _WET_COEFFICIENTS),
    (MODULE, 'register_python_observable'):      lambda _: _eos.register_python_observable('test::python-observable',
                                                     _ExternalObservableProvider(), '', eos.Unit.Unity()),
}


# Members that cannot presently be exercised, together with the reason why.
BINDINGS_XFAILS_UNTESTABLE = {
    ('Model', 'm_t_msbar'): 'SMComponent<components::QCD>::m_t_msbar is not implemented for any reachable scale',
}


_TODO         = 'not yet documented'


# Classes and members that carry no docstring, together with the reason why. A class is keyed by its
# name alone, a member by the pair of its class name and its own name. Dunders are exempt throughout,
# since Python's own conventions describe them.
BINDINGS_XFAILS_UNDOCUMENTED = {
    'BToSWilsonCoefficients':                     _TODO,

    ('BToSWilsonCoefficients', 'c1'):             _TODO,
    ('BToSWilsonCoefficients', 'c2'):             _TODO,
}


# Members whose docstring does not describe every argument that the binding declares, together with
# the reason why.
BINDINGS_XFAILS_UNDESCRIBED_ARGUMENTS = {
    ('LogPrior', 'Flat'):                  'an alias for LogPrior.Uniform, which describes the arguments',
    ('LogLikelihoodBlock', '_Unbinned1D'): 'the unbinned log-likelihood is not fully supported yet',
}


def _external_block_factory(cache):
    class Block:
        number_of_observations = 0

        def evaluate(self):
            return 0.0

    return Block()


class _ExternalObservableProvider:
    kinematic_variables = [ 'q2' ]

    def __call__(self, parameters, kinematics, options):
        class Observable:
            def evaluate(self):
                return 1.0

        return Observable()


def _external_prior_factory(parameters):
    class Prior:
        varied_parameters = [ 'mass::b(MSbar)' ]
        informative       = True

        def evaluate(self):
            return 0.0

        def sample(self):
            pass

        def compute_cdf(self):
            pass

    return Prior()


def is_dunder(member):
    "Returns whether the given member name has both leading and trailing double underscores."
    return member.startswith('__') and member.endswith('__')


def declared_arguments(subject, member):
    """Returns the arguments that the binding declares, taken from its automatic signature lines.

    Arguments only appear under their own name if the binding passes ``args(...)``; without it they
    appear as ``arg1``, ``arg2``, and so on, and carry nothing to compare a docstring against.
    """
    attribute = vars(subject)[member]
    text      = getattr(attribute.__func__ if isinstance(attribute, staticmethod) else attribute, '__doc__', None) or ''
    result    = []

    for line in text.split('\n'):
        if not re.match(rf'\s*{re.escape(member)}\s*\(', line) or '->' not in line:
            continue

        for _, name in re.findall(r'\(([^()]+)\)(\w+)', line):
            if name == 'self' or re.fullmatch(r'arg\d+', name) or name in result:
                continue

            result.append(name)

    return result


def described_arguments(subject, member):
    "Returns the arguments that the docstring of the given member describes."
    return [ match.group(1).strip() for match in re.finditer(r':param\s+([^:]+):', docstring(subject, member)) ]


def docstring(subject, member=None):
    "Returns the docstring of an exported class or member, with the automatic signature lines removed."
    if member is None:
        return (subject.__doc__ or '').strip()

    attribute = vars(subject)[member]
    text      = getattr(attribute.__func__ if isinstance(attribute, staticmethod) else attribute, '__doc__', None)
    lines     = [ line for line in (text or '').split('\n') if not re.match(rf'\s*{re.escape(member)}\s*\(', line) ]

    return '\n'.join(lines).strip()


def call_member(subject, member, test):
    "Exercises a single member, either through its own callable or by calling it with the given arguments."
    if callable(test):
        return test(subject)

    attribute = getattr(subject, member)

    return attribute(*test) if callable(attribute) else attribute


def exported_classes():
    "Returns the name and type of every class exported by the _eos module."
    return { name: value for name, value in vars(_eos).items() if inspect.isclass(value) }


def exported_members(cls):
    "Returns the name of every member bound on the given exported class."
    return { name for name in vars(cls) if name not in IGNORED_MEMBERS }


def exported_functions():
    "Returns the name of every free function exported by the _eos module."
    return { name for name, value in vars(_eos).items() if type(value).__name__ == 'function' and not inspect.isclass(value) }


def is_enumeration(cls):
    return any(isinstance(value, cls) for value in vars(cls).values())


class BindingsCoverageTests(unittest.TestCase):

    def test_every_class_has_a_factory(self):
        "Check that every class exported by the _eos module has an entry in BINDINGS_FACTORIES."
        missing = { name for name, cls in exported_classes().items() if not is_enumeration(cls) } - set(BINDINGS_FACTORIES)
        self.assertEqual(set(), missing,
                f'classes exported by _eos without an entry in BINDINGS_FACTORIES: {sorted(missing)}')

    def test_no_stale_factories(self):
        "Check that BINDINGS_FACTORIES does not describe classes that the _eos module no longer exports."
        stale = set(BINDINGS_FACTORIES) - set(exported_classes())
        self.assertEqual(set(), stale, f'entries in BINDINGS_FACTORIES without a class exported by _eos: {sorted(stale)}')

    def test_every_member_is_described(self):
        "Check that every member of every exported class has an entry in BINDINGS_TESTS or BINDINGS_XFAILS_UNTESTABLE."
        described = set(BINDINGS_TESTS) | set(BINDINGS_XFAILS_UNTESTABLE)
        missing   = set()

        for name, cls in exported_classes().items():
            if is_enumeration(cls):
                continue
            missing |= { (name, member) for member in exported_members(cls) if (name, member) not in described }

        missing |= { (MODULE, f) for f in exported_functions() if (MODULE, f) not in described }

        self.assertEqual(set(), missing,
                f'members exported by _eos without an entry in BINDINGS_TESTS or BINDINGS_XFAILS_UNTESTABLE: {sorted(missing)}')

    def test_no_stale_members(self):
        "Check that BINDINGS_TESTS and BINDINGS_XFAILS_UNTESTABLE do not describe members that the _eos module no longer exports."
        exported = { (MODULE, f) for f in exported_functions() }

        for name, cls in exported_classes().items():
            exported |= { (name, member) for member in exported_members(cls) }

        stale = (set(BINDINGS_TESTS) | set(BINDINGS_XFAILS_UNTESTABLE)) - exported
        self.assertEqual(set(), stale,
                f'entries in BINDINGS_TESTS or BINDINGS_XFAILS_UNTESTABLE without a member exported by _eos: {sorted(stale)}')

    def test_factories(self):
        "Check that every entry in BINDINGS_FACTORIES creates an instance of the class it describes."
        for name, factory in sorted(BINDINGS_FACTORIES.items()):
            with self.subTest(cls=name):
                try:
                    instance = factory()
                except Exception as e:
                    self.fail(f'cannot create an instance of {name}: {e}')

                self.assertIsInstance(instance, getattr(_eos, name), f'factory for {name} returns the wrong type')

    def test_members(self):
        "Check that every member described by BINDINGS_TESTS can be called."
        for (name, member), test in sorted(BINDINGS_TESTS.items(), key=lambda item: item[0]):
            if (name, member) in BINDINGS_XFAILS_UNTESTABLE:
                continue

            with self.subTest(cls=name, member=member):
                subject = _eos if name == MODULE else BINDINGS_FACTORIES[name]()

                try:
                    call_member(subject, member, test)
                except Exception as e:
                    self.fail(f'cannot call {name}.{member}: {type(e).__name__}: {e}')

    def test_docstrings(self):
        "Check that every exported class and member is documented, unless BINDINGS_XFAILS_UNDOCUMENTED excuses it."
        undocumented = set()

        for name, cls in exported_classes().items():
            if is_enumeration(cls):
                continue

            if not docstring(cls) and name not in BINDINGS_XFAILS_UNDOCUMENTED:
                undocumented.add(name)

            for member in exported_members(cls):
                if is_dunder(member): # Python's own conventions describe the 'dunders' (i.e., __*__ methods), so we do not require a docstring for them
                    continue

                if not docstring(cls, member) and (name, member) not in BINDINGS_XFAILS_UNDOCUMENTED:
                    undocumented.add((name, member))

        self.assertEqual(set(), undocumented,
                f'classes and members exported by _eos without a docstring: {sorted(undocumented, key=str)}')

    def test_arguments_are_described(self):
        "Check that a docstring describes every argument that its binding declares."
        undescribed = set()

        for name, cls in exported_classes().items():
            if is_enumeration(cls):
                continue

            for member in exported_members(cls):
                if (name, member) in BINDINGS_XFAILS_UNDESCRIBED_ARGUMENTS:
                    continue

                described = described_arguments(cls, member)

                if any(argument not in described for argument in declared_arguments(cls, member)):
                    undescribed.add((name, member))

        self.assertEqual(set(), undescribed,
                f'members whose docstring does not describe every declared argument: {sorted(undescribed)}')

    def test_arguments_are_described_under_their_own_name(self):
        "Check that a docstring does not describe an argument that its binding does not declare."
        misnamed = set()

        for name, cls in exported_classes().items():
            if is_enumeration(cls):
                continue

            for member in exported_members(cls):
                declared = declared_arguments(cls, member)

                if not declared:
                    continue

                if any(argument not in declared for argument in described_arguments(cls, member)):
                    misnamed.add((name, member))

        self.assertEqual(set(), misnamed,
                f'members whose docstring describes an argument that the binding does not declare: {sorted(misnamed)}')

    def test_no_stale_undescribed_arguments(self):
        "Check that BINDINGS_XFAILS_UNDESCRIBED_ARGUMENTS does not excuse a member that describes its arguments."
        classes = exported_classes()
        stale   = set()

        for (name, member), reason in BINDINGS_XFAILS_UNDESCRIBED_ARGUMENTS.items():
            cls = classes.get(name)

            if cls is None or member not in vars(cls):
                stale.add((name, member))
                continue

            described = described_arguments(cls, member)

            if all(argument in described for argument in declared_arguments(cls, member)):
                stale.add((name, member))

        self.assertEqual(set(), stale,
                f'entries in BINDINGS_XFAILS_UNDESCRIBED_ARGUMENTS that describe their arguments or are no longer exported: {sorted(stale)}')

    def test_no_stale_undocumented(self):
        "Check that BINDINGS_XFAILS_UNDOCUMENTED does not excuse anything that is documented or no longer exported."
        classes = exported_classes()
        stale   = set()

        for key, reason in BINDINGS_XFAILS_UNDOCUMENTED.items():
            name, member = (key, None) if isinstance(key, str) else key
            cls          = classes.get(name)

            if cls is None or (member is not None and member not in vars(cls)):
                stale.add(key)
            elif docstring(cls, member):
                stale.add(key)

        self.assertEqual(set(), stale,
                f'entries in BINDINGS_XFAILS_UNDOCUMENTED that are documented or no longer exported: {sorted(stale, key=str)}')

    def test_xfails(self):
        "Check that every member described by BINDINGS_XFAILS_UNTESTABLE still fails, so that stale entries do not survive a fix."
        for (name, member), reason in sorted(BINDINGS_XFAILS_UNTESTABLE.items()):
            with self.subTest(cls=name, member=member):
                self.assertIn((name, member), BINDINGS_TESTS,
                        f'{name}.{member} is listed in BINDINGS_XFAILS_UNTESTABLE but has no entry in BINDINGS_TESTS')

                subject = _eos if name == MODULE else BINDINGS_FACTORIES[name]()

                with self.assertRaises(Exception, msg=f'{name}.{member} no longer fails; remove it from BINDINGS_XFAILS_UNTESTABLE ({reason})'):
                    call_member(subject, member, BINDINGS_TESTS[(name, member)])


if __name__ == '__main__':
    unittest.main(verbosity=2)
