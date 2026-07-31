#!/usr/bin/python
# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2024-2026 Danny van Dyk
# Copyright (c) 2024-2025 Matthew Kirk
# Copyright (c) 2024      Carolina Bolognani
# Copyright (c) 2024      Méril Reboud
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

from .deserializable import Deserializable, InvalidComponent
from .diagnostic import Diagnostic, Severity
from eos.figure import FigureFactory
from dataclasses import dataclass, field
from collections import defaultdict
import copy as _copy
import eos
import inspect
from collections import Counter


class _AnalysisFileDeserializable(Deserializable):
    @classmethod
    def from_dict(cls, **kwargs):
        return Deserializable.make_with_diagnostics(cls, **kwargs)

    def _diagnostics(self):
        yield from ()

    def validate_structure(self):
        yield from self._diagnostics()


def _segments(raw_children, identifier='name'):
    return [
        child[identifier]
        if (
            isinstance(child, dict)
            and identifier in child
            and isinstance(child[identifier], str)
            and '/' not in child[identifier]
            and not any(character.isspace() for character in child[identifier])
        )
        else index
        for index, child in enumerate(raw_children)
    ]


def _validate_children(children, segments, prefix):
    assert len(children) == len(segments)
    for child, segment in zip(children, segments):
        validator = getattr(child, 'validate_structure', None)
        if validator is None:
            validator = getattr(child, 'validate', None)
        if validator is not None:
            yield from (diagnostic.prefixed(prefix, segment) for diagnostic in validator())


def _check_qualified(context, value, kind, path):
    try:
        qn = eos.QualifiedName(value)
    except RuntimeError as e:
        yield Diagnostic(path, Severity.ERROR, f"'{value}' is not a valid qualified name: {e}")
        return
    if not context.lookup(kind, qn):
        yield Diagnostic(path, Severity.ERROR, f"{kind} '{value}' is unknown to EOS")


def _check_file_local_name(value, kind, path):
    if '/' in value:
        yield Diagnostic(path, Severity.ERROR, f"Invalid character '/' in {kind} '{value}'")
    if any(character.isspace() for character in value):
        yield Diagnostic(path, Severity.ERROR, f"Invalid whitespace in {kind} '{value}'")


@dataclass
class MetadataAuthorDescription(_AnalysisFileDeserializable):
    """Describes a single author in an analysis file's metadata.

    :param name: The author's name.
    :type name: str
    :param affiliation: The author's affiliation. Optional.
    :type affiliation: str
    :param email: The author's email address. Optional.
    :type email: str
    """
    name:str
    affiliation:str=''
    email:str=''


@dataclass
class MetadataDescription(_AnalysisFileDeserializable):
    """Describes the metadata of an analysis file.

    :param title: A human-readable title for the analysis. Optional.
    :type title: str
    :param id: A unique identifier for the analysis. Optional.
    :type id: str
    :param authors: The list of authors of the analysis. Optional.
    :type authors: list[MetadataAuthorDescription]
    """
    title:str=''
    id:str=''
    authors:list[MetadataAuthorDescription]=field(default_factory=list)

    def validate_structure(self):
        yield from self._diagnostics()
        segments = getattr(self, '_author_segments', list(range(len(self.authors))))
        yield from _validate_children(self.authors, segments, 'authors')

    @staticmethod
    def from_dict(**kwargs):
        """Create a :class:`MetadataDescription` from its keyword description.

        Deserializes the nested ``authors`` list into :class:`MetadataAuthorDescription` instances.

        :returns: The instantiated metadata description.
        :rtype: MetadataDescription
        """
        _kwargs = _copy.deepcopy(kwargs)
        diagnostics = Deserializable.check_keys(MetadataDescription, _kwargs)
        if diagnostics:
            return InvalidComponent(diagnostics)
        if 'authors' in kwargs:
            author_segments = _segments(kwargs['authors'])
            _kwargs['authors'] = [MetadataAuthorDescription.from_dict(**a) for a in kwargs['authors']]
        else:
            author_segments = []
        result = MetadataDescription(**_kwargs)
        result._author_segments = author_segments
        return result


class PriorDescription:
    """Polymorphic description of a single prior on one or more parameters.

    This is a dispatcher rather than a concrete description: :meth:`from_dict` selects the concrete type
    from the ``type`` key. The recognized types are ``constraint`` (:class:`ConstraintPriorDescription`),
    ``uniform``/``flat`` (:class:`UniformPriorDescription`), ``scale`` (:class:`ScalePriorDescription`),
    ``gauss``/``gaussian`` (:class:`GaussianPriorDescription`, or :class:`CurtailedGaussianDescription`
    if ``min``/``max`` are given), ``poisson`` (:class:`PoissonPriorDescription`), and ``transform``
    (:class:`TransformPriorDescription`).

    For backwards compatibility, a description that omits the ``type`` key but provides a ``constraint``
    key is still interpreted as a :class:`ConstraintPriorDescription`; this form is deprecated.
    """

    @staticmethod
    def from_dict(**kwargs):
        """Create the concrete prior description matching the given keyword description.

        :returns: An instance of the concrete :class:`PriorDescription` subtype selected from ``type``.
        """
        if 'type' in kwargs:
            _kwargs = _copy.deepcopy(kwargs)
            _kwargs.pop('type')
            if kwargs['type'] in ("constraint",):
                return Deserializable.make_with_diagnostics(ConstraintPriorDescription, **_kwargs)
            elif kwargs['type'] in ("uniform", "flat"):
                return Deserializable.make_with_diagnostics(UniformPriorDescription, **_kwargs)
            elif kwargs['type'] in ("scale",):
                return Deserializable.make_with_diagnostics(ScalePriorDescription, **_kwargs)
            elif kwargs['type'] in ("gauss", "gaussian"):
                if ('min' in kwargs) != ('max' in kwargs):
                    missing_bound = 'max' if 'min' in kwargs else 'min'
                    return InvalidComponent([
                        Diagnostic(
                            (missing_bound,),
                            Severity.ERROR,
                            'A Gaussian prior description must contain both \'min\' and \'max\' or neither',
                        ),
                    ])
                if 'min' in kwargs:
                    return Deserializable.make_with_diagnostics(CurtailedGaussianDescription, **_kwargs)
                return Deserializable.make_with_diagnostics(GaussianPriorDescription, **_kwargs)
            elif kwargs['type'] in ("poisson",):
                return Deserializable.make_with_diagnostics(PoissonPriorDescription, **_kwargs)
            elif kwargs['type'] in ("transform",):
                return Deserializable.make_with_diagnostics(TransformPriorDescription, **_kwargs)
        elif 'constraint' in kwargs:
            eos.warn('A constraint prior description without a \'type\' key is deprecated and will be removed in a future version of EOS; add \'type: constraint\' instead')
            return Deserializable.make_with_diagnostics(ConstraintPriorDescription, **kwargs)

        return InvalidComponent([
            Diagnostic(('type',), Severity.ERROR, 'Unknown type of prior description'),
        ])

@dataclass
class PoissonPriorDescription(_AnalysisFileDeserializable):
    r"""Describes a Poisson prior on a single parameter.

    :param parameter: The qualified name of the parameter.
    :type parameter: str
    :param k: The number of observed counts that parametrizes the Poisson distribution.
    :type k: float
    """
    parameter:str
    k:float
    type:str=field(repr=False, init=False, default="poisson")

@dataclass
class CurtailedGaussianDescription(_AnalysisFileDeserializable):
    r"""Describes a Gaussian prior truncated to a finite interval (``type: gauss`` with ``min``/``max``).

    .. deprecated:: 1.0.21
       The curtailed Gaussian prior description is deprecated and will be removed in a future version
       of EOS; use ``type: gaussian`` without the ``min``/``max`` keys instead.

    :param parameter: The qualified name of the parameter.
    :type parameter: str
    :param central: The central value (mean) of the Gaussian.
    :type central: float
    :param sigma: The standard deviation of the Gaussian.
    :type sigma: float
    :param min: The lower bound to which the prior is truncated.
    :type min: float
    :param max: The upper bound to which the prior is truncated.
    :type max: float
    """
    parameter:str
    central:float
    sigma:float
    min:float
    max:float
    type:str=field(repr=False, init=False, default="gaussian")

    def __post_init__(self):
        pass

@dataclass
class GaussianPriorDescription(_AnalysisFileDeserializable):
    r"""Describes a Gaussian prior on a single parameter.

    :param parameter: The qualified name of the parameter.
    :type parameter: str
    :param central: The central value (mean) of the Gaussian.
    :type central: float
    :param sigma: The standard deviation of the Gaussian.
    :type sigma: float
    """
    parameter:str
    central:float
    sigma:float
    type:str=field(repr=False, init=False, default="gaussian")


@dataclass
class ScalePriorDescription(_AnalysisFileDeserializable):
    r"""Describes a prior on a renormalization scale parameter.

    :param parameter: The qualified name of the scale parameter.
    :type parameter: str
    :param min: The lower bound of the scale variation.
    :type min: float
    :param max: The upper bound of the scale variation.
    :type max: float
    :param mu_0: The default (central) scale.
    :type mu_0: float
    :param lambda_scale: The multiplicative factor that controls the width of the log-uniform variation.
    :type lambda_scale: float
    """
    parameter:str
    min:float
    max:float
    mu_0:float
    lambda_scale:float
    type:str=field(repr=False, init=False, default="scale")

@dataclass
class UniformPriorDescription(_AnalysisFileDeserializable):
    r"""Describes a uniform (flat) prior on a single parameter.

    :param parameter: The qualified name of the parameter.
    :type parameter: str
    :param min: The lower bound of the prior's support.
    :type min: float
    :param max: The upper bound of the prior's support.
    :type max: float
    """
    parameter:str
    min:float
    max:float
    type:str=field(repr=False, init=False, default="uniform")

@dataclass
class ConstraintPriorDescription(_AnalysisFileDeserializable):
    r"""Describes a (possibly correlated, multivariate) prior taken from a built-in EOS constraint.

    :param constraint: The qualified name of the EOS constraint used as the prior.
    :type constraint: str
    """
    constraint:str
    type:str=field(repr=False, init=False, default="constraint")

@dataclass
class TransformPriorDescription(_AnalysisFileDeserializable):
    r"""Describes a prior on a linear transformation of several parameters.

    :param parameters: The qualified names of the parameters that enter the transformation.
    :type parameters: list[str]
    :param shift: The constant shift applied to the parameters before the transformation.
    :type shift: list[float]
    :param transform: The transformation matrix applied to the (shifted) parameters.
    :type transform: list[list[float]]
    :param min: The lower bounds of the transformed parameters' support.
    :type min: list[float]
    :param max: The upper bounds of the transformed parameters' support.
    :type max: list[float]
    """
    parameters:list[str]
    shift:list[float]
    transform:list[list[float]]
    min:list[float]
    max:list[float]
    type:str=field(repr=False, init=False, default="transform")

# Maps the YAML selector of a prior description to the corresponding concrete description class.
# The canonical dispatch logic lives in PriorDescription.from_dict; this mapping mirrors it and is
# used to generate the reference documentation of the recognized prior types. It is assigned here,
# rather than in the class body, because the concrete classes are only defined above.
PriorDescription.registry = {
    'constraint':           ConstraintPriorDescription,
    'uniform':              UniformPriorDescription,
    'flat':                 UniformPriorDescription,
    'scale':                ScalePriorDescription,
    'gauss':                GaussianPriorDescription,
    'gaussian':             GaussianPriorDescription,
    'gauss (with min/max)': CurtailedGaussianDescription,
    'poisson':              PoissonPriorDescription,
    'transform':            TransformPriorDescription,
}

@dataclass
class PriorComponent(_AnalysisFileDeserializable):
    r"""Describes a single named prior, i.e. one entry of an analysis file's ``priors`` list.

    A named prior bundles one or more prior descriptions: a single description for a univariate prior, or
    several for an (uncorrelated) multivariate prior.

    :param name: The unique name of the prior, by which posteriors refer to it.
    :type name: str
    :param descriptions: The list of prior descriptions that make up this prior.
    :type descriptions: list[PriorDescription]
    """
    name:str
    descriptions:list

    def _diagnostics(self):
        yield from _check_file_local_name(self.name, 'prior name', ('name',))

    def validate_structure(self):
        yield from self._diagnostics()
        segments = getattr(self, '_description_segments', list(range(len(self.descriptions))))
        yield from _validate_children(self.descriptions, segments, 'descriptions')

    def validate_semantics(self, context):
        segments = getattr(self, '_description_segments', list(range(len(self.descriptions))))
        assert len(self.descriptions) == len(segments)
        for description, segment in zip(self.descriptions, segments):
            # The concrete prior descriptions reference EOS entities through one of three fields: a
            # single parameter (uniform, scale, gauss, poisson), a single constraint (constraint
            # priors), or a list of parameters (transform priors). An InvalidComponent has none of
            # them and is skipped here; its own diagnostics come from the structural phase.
            if hasattr(description, 'parameter'):
                yield from _check_qualified(
                    context,
                    description.parameter,
                    'parameter',
                    ('descriptions', segment, 'parameter'),
                )
            if hasattr(description, 'constraint'):
                yield from _check_qualified(
                    context,
                    description.constraint,
                    'constraint',
                    ('descriptions', segment, 'constraint'),
                )
            for index, parameter in enumerate(getattr(description, 'parameters', ())):
                yield from _check_qualified(
                    context,
                    parameter,
                    'parameter',
                    ('descriptions', segment, 'parameters', index),
                )

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`PriorComponent` from its keyword description.

        Deserializes each entry of ``descriptions`` via :meth:`PriorDescription.from_dict`. The deprecated
        ``parameters`` key is still accepted as a fall-back alias for ``descriptions``.

        :returns: The instantiated prior component.
        :rtype: PriorComponent
        """
        _kwargs = _copy.deepcopy(kwargs)
        if "descriptions" in kwargs:
            if "parameters" in kwargs:
                eos.error(f'Both \'descriptions\' and \'parameters\' are provided for prior component \'{kwargs["name"]}\', ignoring legacy support for \'parameters\'')
        if "parameters" in kwargs:
            eos.warn(f'\'parameters\' is in the description of prior component \'{kwargs["name"]}\', use \'descriptions\' instead')
            _kwargs['descriptions'] = _kwargs.pop("parameters")
        diagnostics = Deserializable.check_keys(cls, _kwargs)
        if diagnostics:
            return InvalidComponent(diagnostics)
        description_segments = _segments(_kwargs['descriptions'])
        _kwargs['descriptions'] = [PriorDescription.from_dict(**d) for d in _kwargs['descriptions']]
        result = cls(**_kwargs)
        result._description_segments = description_segments
        return result



@dataclass
class ConstraintLikelihoodDescription(_AnalysisFileDeserializable):
    r"""Describes a likelihood contribution taken from a built-in EOS constraint.

    :param constraint: The qualified name of the EOS constraint.
    :type constraint: eos.QualifiedName
    """
    constraint:eos.QualifiedName

@dataclass
class ManualConstraintDescription(_AnalysisFileDeserializable):
    r"""Describes a likelihood contribution from a constraint defined inline in the analysis file.

    :param name: The qualified name under which the manual constraint is registered.
    :type name: eos.QualifiedName
    :param info: The constraint definition, in the same format as a built-in EOS constraint entry.
    :type info: dict
    """
    name:eos.QualifiedName
    info:dict

@dataclass
class PyHFConstraintDescription(_AnalysisFileDeserializable):
    r"""Describes a likelihood contribution imported from a pyhf workspace.

    :param file: The path to the JSON file specifying the pyhf workspace.
    :type file: str
    :param parameter_map: An optional mapping from pyhf parameter names to EOS observables or parameters; see :class:`eos.PyhfLogLikelihood`.
    :type parameter_map: dict
    """
    file:str
    parameter_map:dict=field(default_factory=dict)

@dataclass
class LikelihoodComponent(_AnalysisFileDeserializable):
    r"""Describes a single named likelihood, i.e. one entry of an analysis file's ``likelihoods`` list.

    A named likelihood may combine any number of built-in constraints, inline (manual) constraints, and
    pyhf-based contributions; at least one of the three must be present.

    :param name: The unique name of the likelihood, by which posteriors refer to it.
    :type name: str
    :param constraints: The built-in EOS constraints contributing to the likelihood. Optional.
    :type constraints: list[ConstraintLikelihoodDescription]
    :param manual_constraints: The inline constraint definitions contributing to the likelihood. Optional.
    :type manual_constraints: list[ManualConstraintDescription]
    :param pyhf: The pyhf-based contribution to the likelihood, or None if there is none. Optional.
    :type pyhf: PyHFConstraintDescription | None
    """
    name:str
    constraints:list=field(default_factory=list)
    manual_constraints:list=field(default_factory=list)
    # unlike 'constraints' and 'manual_constraints', a likelihood carries at most one pyhf
    # contribution; from_dict deserializes it into a single PyHFConstraintDescription
    pyhf:PyHFConstraintDescription|None=None

    def _diagnostics(self):
        yield from _check_file_local_name(self.name, 'likelihood name', ('name',))
        if not (self.constraints or self.manual_constraints or self.pyhf):
            yield Diagnostic(
                (),
                Severity.ERROR,
                'LikelihoodComponent must have at least one of constraints, manual_constraints, or pyhf',
            )

    def validate_structure(self):
        yield from self._diagnostics()
        constraint_segments = getattr(self, '_constraint_segments', list(range(len(self.constraints))))
        yield from _validate_children(self.constraints, constraint_segments, 'constraints')
        manual_constraint_segments = getattr(
            self,
            '_manual_constraint_segments',
            list(range(len(self.manual_constraints))),
        )
        yield from _validate_children(
            self.manual_constraints,
            manual_constraint_segments,
            'manual_constraints',
        )
        if self.pyhf:
            yield from (diagnostic.prefixed('pyhf') for diagnostic in self.pyhf.validate_structure())

    def validate_semantics(self, context):
        constraint_segments = getattr(self, '_constraint_segments', list(range(len(self.constraints))))
        assert len(self.constraints) == len(constraint_segments)
        for description, segment in zip(self.constraints, constraint_segments):
            yield from _check_qualified(
                context,
                description.constraint,
                'constraint',
                ('constraints', segment),
            )

        known_constraints = eos.Constraints()
        manual_segments = getattr(
            self,
            '_manual_constraint_segments',
            list(range(len(self.manual_constraints))),
        )
        assert len(self.manual_constraints) == len(manual_segments)
        for manual, segment in zip(self.manual_constraints, manual_segments):
            yield from _check_qualified(
                context,
                manual.name,
                'constraint',
                ('manual_constraints', segment, 'name'),
            )
            try:
                known_constraints[manual.name]
                yield Diagnostic(
                    ('manual_constraints', segment),
                    Severity.ERROR,
                    f"Error in likelihood {self.name}: The manual constraint named '{manual.name}' matches an already defined constraint name",
                )
            except RuntimeError:
                pass

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`LikelihoodComponent` from its keyword description.

        Deserializes the ``constraints``, ``manual_constraints``, and ``pyhf`` entries into their
        respective description types.

        :returns: The instantiated likelihood component.
        :rtype: LikelihoodComponent
        """
        _kwargs = _copy.deepcopy(kwargs)
        diagnostics = Deserializable.check_keys(cls, _kwargs)
        if diagnostics:
            return InvalidComponent(diagnostics)
        if 'constraints' in kwargs:
            constraint_segments = _segments(kwargs['constraints'])
            _kwargs['constraints'] = [ConstraintLikelihoodDescription.from_dict(constraint=c) for c in kwargs['constraints']]
        else:
            constraint_segments = []
        if 'manual_constraints' in kwargs:
            raw_manual_constraints = [
                {'name': name, 'info': info}
                for name, info in kwargs['manual_constraints'].items()
            ]
            manual_constraint_segments = _segments(raw_manual_constraints)
            _kwargs['manual_constraints'] = [
                ManualConstraintDescription.from_dict(**description)
                for description in raw_manual_constraints
            ]
        else:
            manual_constraint_segments = []
        if 'pyhf' in kwargs:
            _kwargs['pyhf'] = PyHFConstraintDescription.from_dict(**kwargs['pyhf'])
        result = cls(**_kwargs)
        result._constraint_segments = constraint_segments
        result._manual_constraint_segments = manual_constraint_segments
        return result



@dataclass
class PosteriorDescription(_AnalysisFileDeserializable):
    r"""Describes a single named posterior, i.e. one entry of an analysis file's ``posteriors`` list.

    A posterior combines one or more named priors with one or more named likelihoods, optionally fixing
    parameters and setting global options for all of its theory predictions.

    :param name: The unique name of the posterior.
    :type name: str
    :param prior: The names of the priors that make up the posterior.
    :type prior: list[str]
    :param likelihood: The names of the likelihoods that make up the posterior.
    :type likelihood: list[str]
    :param global_options: Options forwarded to all theory predictions of the posterior. Optional.
    :type global_options: dict
    :param fixed_parameters: Parameters to fix, given as a mapping from qualified name to value. Optional.
    :type fixed_parameters: dict
    """
    name:str
    prior:list
    likelihood:list
    global_options:dict=field(default_factory=dict)
    fixed_parameters:dict=field(default_factory=dict)

    def _diagnostics(self):
        yield from _check_file_local_name(self.name, 'posterior name', ('name',))

    def validate_semantics(self, context):
        for parameter in self.fixed_parameters:
            yield from _check_qualified(
                context,
                parameter,
                'parameter',
                ('fixed_parameters', parameter),
            )



@dataclass
class ObservableComponent(_AnalysisFileDeserializable):
    r"""Describes a custom observable defined in an analysis file's ``observables`` section.

    :param name: The qualified name under which the new observable is registered.
    :type name: str
    :param latex: The LaTeX representation of the observable, used in figures and tables.
    :type latex: str
    :param unit: The unit of the observable.
    :type unit: str
    :param expression: The EOS expression that defines the observable.
    :type expression: str
    :param options: Default options applied when evaluating the observable. Optional.
    :type options: dict
    """
    name:str
    latex:str
    unit:str
    expression:str
    options:dict=field(default_factory=dict)

    def validate_semantics(self, context):
        try:
            references = eos.analyze_expression(self.expression)
        except RuntimeError as error:
            message = str(error)
            if message.startswith('Parsing Error:'):
                yield Diagnostic(
                    ('expression',),
                    Severity.ERROR,
                    f'Expression does not parse: {message}',
                )
                return
            elif 'is not a valid qualified name:' in message:
                yield Diagnostic(
                    ('expression',),
                    Severity.ERROR,
                    f'Expression contains a malformed qualified name: {message}',
                )
                return
            else:
                raise

        for observable in references.observables:
            yield from _check_qualified(context, observable.full(), 'observable', ('expression',))
        for parameter in references.parameters:
            yield from _check_qualified(context, parameter.full(), 'parameter', ('expression',))



@dataclass
class PredictionObservableComponent(_AnalysisFileDeserializable):
    r"""Describes a single observable to be predicted, i.e. one entry of a prediction's ``observables`` list.

    :param name: The qualified name of the observable to predict.
    :type name: str
    :param kinematics: The kinematic configuration(s) at which to evaluate the observable. Optional.
    :type kinematics: dict
    :param options: Options applied when evaluating the observable. Optional.
    :type options: dict
    """
    name:str
    kinematics:dict=field(default_factory=dict) # TODO: once we bump the minimum python version to 3.10, use dict | list[dict] instead
    options:dict=field(default_factory=dict)

@dataclass
class PredictionDescription(_AnalysisFileDeserializable):
    r"""Describes a single named prediction, i.e. one entry of an analysis file's ``predictions`` list.

    A prediction lists the observables to be evaluated on previously obtained importance samples,
    optionally fixing parameters and setting global options.

    :param name: The unique name of the prediction.
    :type name: str
    :param observables: The observables to be predicted.
    :type observables: list[PredictionObservableComponent]
    :param global_options: Options forwarded to all observables of the prediction. Optional.
    :type global_options: dict
    :param fixed_parameters: Parameters to fix, given as a mapping from qualified name to value. Optional.
    :type fixed_parameters: dict
    """
    name:str
    observables:list
    global_options:dict=field(default_factory=dict)
    fixed_parameters:dict=field(default_factory=dict)

    def _diagnostics(self):
        yield from _check_file_local_name(self.name, 'prediction name', ('name',))

    def validate_structure(self):
        yield from self._diagnostics()
        segments = getattr(self, '_observable_segments', list(range(len(self.observables))))
        yield from _validate_children(self.observables, segments, 'observables')

    def validate_semantics(self, context):
        segments = getattr(self, '_observable_segments', list(range(len(self.observables))))
        assert len(self.observables) == len(segments)
        for observable, segment in zip(self.observables, segments):
            yield from _check_qualified(
                context,
                observable.name,
                'observable',
                ('observables', segment, 'name'),
            )

        for parameter in self.fixed_parameters:
            yield from _check_qualified(
                context,
                parameter,
                'parameter',
                ('fixed_parameters', parameter),
            )

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`PredictionDescription` from its keyword description.

        Deserializes each entry of ``observables`` into a :class:`PredictionObservableComponent`.

        :returns: The instantiated prediction description.
        :rtype: PredictionDescription
        """
        _kwargs = _copy.deepcopy(kwargs)
        diagnostics = Deserializable.check_keys(cls, _kwargs)
        if diagnostics:
            return InvalidComponent(diagnostics)
        observable_segments = _segments(kwargs["observables"])
        _kwargs["observables"] = [PredictionObservableComponent.from_dict(**o) for o in kwargs["observables"]]
        result = cls(**_kwargs)
        result._observable_segments = observable_segments
        return result



@dataclass
class ParameterComponent(_AnalysisFileDeserializable):
    r"""Describes a custom parameter defined in an analysis file's ``parameters`` section.

    :param name: The new EOS qualified name of the parameter.
    :type name: str
    :param latex: The LaTeX representation of the parameter, used in figures and tables.
    :type latex: str
    :param unit: The unit of the parameter.
    :type unit: str
    :param central: The default (central) value of the parameter.
    :type central: float
    :param min: The lower end of the parameter's allowed range.
    :type min: float
    :param max: The upper end of the parameter's allowed range.
    :type max: float
    :param alias_of: The qualified names of existing parameters for which this parameter is an alias. Optional.
    :type alias_of: list[str]
    """
    name:str
    latex:str
    unit:str
    central:float
    min:float
    max:float
    alias_of:list=field(default_factory=list)

    def validate_semantics(self, context):
        for index, alias in enumerate(self.alias_of):
            yield from _check_qualified(
                context,
                alias,
                'parameter',
                ('alias_of', index),
            )



# Mapping between CLI argument names and internal names for individual tasks
_task_argument_map = {
    # list-likelihoods
    ('list-likelihoods', 'd'): 'display', ('list-likelihoods', 'display-details'): 'display', ('list-likelihoods', 'DISPLAY'): 'display',
    # sample-mcmc
    ('sample-mcmc', 'POSTERIOR'): 'posterior',
    ('sample-mcmc', 'CHAIN-IDX'): 'chain',
    ('sample-mcmc', 'number-of-samples'): 'N',
    ('sample-mcmc', 'S'): 'stride', ('sample-mcmc', 'STRIDE'): 'stride',
    ('sample-mcmc', 'p'): 'preruns', ('sample-mcmc', 'number-of-preruns'): 'preruns', ('sample-mcmc', 'PRERUNS'): 'preruns',
    ('sample-mcmc', 'n'): 'pre_N', ('sample-mcmc', 'number-of-prerun-samples'): 'pre_N', ('sample-mcmc', 'PRE_N'): 'pre_N',
    ('sample-mcmc', 's'): 'start_point', ('sample-mcmc', 'start-point'): 'start_point', ('sample-mcmc', 'START_POINT'): 'start_point',
    ('sample-mcmc', 'c'): 'cov_scale', ('sample-mcmc', 'cov-scale'): 'cov_scale', ('sample-mcmc', 'COV_SCALE'): 'cov_scale',
    # sample-pmc
    ('sample-pmc', 'POSTERIOR'): 'posterior',
    ('sample-pmc', 'n'): 'step_N', ('sample-pmc', 'number-of-adaptation-samples'): 'step_N', ('sample-pmc', 'STEP_N'): 'step_N',
    ('sample-pmc', 's'): 'steps', ('sample-pmc', 'number-of-adaptation-steps'): 'steps', ('sample-pmc', 'STEPS'): 'steps',
    ('sample-pmc', 't'): 'perplexity_threshold', ('sample-pmc', 'perplexity-threshold'): 'perplexity_threshold', ('sample-pmc', 'PERPLEXITY_THRESHOLD'): 'perplexity_threshold',
    ('sample-pmc', 'w'): 'weight_thresold', ('sample-pmc', 'weight-threshold'): 'weight_thresold', ('sample-pmc', 'WEIGHT_THRESHOLD'): 'weight_thresold',
    ('sample-pmc', 'N'): 'final_N', ('sample-pmc', 'numerb-of-final-samples'): 'final_N', ('sample-pmc', 'FINAL_N'): 'final_N',
    ('sample-pmc', 'PMC_ITERATIONS'): 'pmc_iterations',
    ('sample-pmc', 'pmc-rel-tol'): 'pmc_rel_tol', ('sample-pmc', 'PMC_REL_TOL'): 'pmc_rel_tol',
    ('sample-pmc', 'pmc-abs-tol'): 'pmc_abs_tol', ('sample-pmc', 'PMC_ABS_TOL'): 'pmc_abs_tol',
    ('sample-pmc', 'l'): 'pmc_lookback', ('sample-pmc', 'pmc-lookback'): 'pmc_lookback', ('sample-pmc', 'PMC_LOOKBACK'): 'pmc_lookback',
    ('sample-pmc', 'p'): 'initial_proposal', ('sample-pmc', 'initial-proposal'): 'initial_proposal', ('sample-pmc', 'INITIAL_PROPOSAL'): 'initial_proposal',
    ('sample-pmc', 'S'): 'sigma_stat_test', ('sample-pmc', 'sigma-stat-test'): 'sigma_stat_test', ('sample-pmc', 'SIGMA_STAT_TEST'): 'sigma_stat_test',
    # sample-nested
    ('sample-nested', 'B'): 'bound', ('sample-nested', 'target-bound'): 'bound', ('sample-nested', 'BOUND'): 'bound',
    ('sample-nested', 'n'): 'nlive', ('sample-nested', 'number-of-live-points'): 'nlive', ('sample-nested', 'NLIVE'): 'nlive',
    ('sample-nested', 'd'): 'dlogz', ('sample-nested', 'evidence-tolerance'): 'dlogz', ('sample-nested', 'DLOGZ'): 'dlogz',
    ('sample-nested', 'm'): 'maxiter', ('sample-nested', 'max-number-iterations'): 'maxiter', ('sample-nested', 'MAXITER'): 'maxiter',
    ('sample-nested', 'min-number-iterations'): 'miniter', ('sample-nested', 'MINITER'): 'miniter',
    ('sample-nested', 's'): 'seed', ('sample-nested', 'use-random-seed'): 'seed', ('sample-nested', 'SEED'): 'seed',
    ('sample-nested', 'M'): 'sample', ('sample-nested', 'sampling-method'): 'sample', ('sample-nested', 'SAMPLE'): 'sample',
    # plot-samples
    ('plot-samples', 'POSTERIOR'): 'posterior',
    ('plot-samples', 'B'): 'bins', ('plot-samples', 'BINS'): 'bins',
    # find-mode
    ('find-mode', 'POSTERIOR'): 'posterior',
    ('find-mode', 'o'): 'optimizations', ('find-mode', 'OPTIMIZATIONS'): 'optimizations',
    ('find-mode', 'c'): 'chain', ('find-mode', 'from-mcmc'): 'chain', ('find-mode', 'CHAIN'): 'chain',
    ('find-mode', 'S'): 'importance_samples', ('find-mode', 'from-samples'): 'importance_samples',
    ('find-mode', 'p'): 'start_point', ('find-mode', 'from-point'): 'start_point', ('find-mode', 'START_POINT'): 'start_point',
    ('find-mode', 's'): 'seed', ('find-mode', 'use-random-seed'): 'seed', ('find-mode', 'SEED'): 'seed',
    ('find-mode', 'L'): 'label', ('find-mode', 'LABEL'): 'label',
    # mixture-product
    ('mixture-product', 'POSTERIOR'): 'posterior',
    ('mixture-product', 'POSTERIORS'): 'posteriors',
    # find-clusters
    ('find-clusters', 'POSTERIOR'): 'posterior',
    ('find-clusters', 't'): 'threshold', ('find-clusters', 'THRESHOLD'): 'threshold',
    ('find-clusters', 'c'): 'K_g', ('find-clusters', 'clusters-per-group'): 'cluster', ('find-clusters', 'CLUSTER'): 'cluster', ('find-clusters', 'K_G'): 'K_g',
    # predict-observables
    ('predict-observables', 'POSTERIOR'): 'posterior',
    ('predict-obserables', 'PREDICTION'): 'prediction',
    ('predict-observables', 'B'): 'begin', ('predict-observables', 'begin-index'): 'begin', ('predict-observables', 'BEGIN'): 'begin',
    ('predict-observables', 'E'): 'end', ('predict-observables', 'end-index'): 'end', ('predict-observables', 'END'): 'end',
    # corner-plot
    ('corner-plot', 'B'): 'begin', ('corner-plot', 'begin-parameter'): 'begin', ('corner-plot', 'BEGIN'): 'begin',
    ('corner-plot', 'E'): 'end', ('corner-plot', 'end-parameter'): 'end', ('corner-plot', 'END'): 'end',
    ('corner-plot', 'F'): 'format', ('corner-plot', 'FORMAT'): 'format',
    # list-step-dependencies
    ('list-step-dependencies', 'ID'): 'id',
    # report
    ('report', 'TEMPLATE'): 'template',
    ('report', 'o'): 'output_file', ('report', 'output-file'): 'output_file', ('report', 'OUTPUT_FILE'): 'output_file',
    # draw-figure
    ('draw-figure', 'FIGURE'): 'figure_name', ('draw-figure', 'figure'): 'figure_name',
    ('draw-figure', 'F'): 'format', ('draw-figure', 'FORMAT'): 'format',
}

@dataclass
class TaskComponent(_AnalysisFileDeserializable):
    r"""Describes a single task invocation within a step, i.e. one entry of a step's ``tasks`` list.

    :param task: The name of the task to run; must be a registered EOS task.
    :type task: str
    :param arguments: The arguments passed to the task, given as a mapping from argument name to value. Optional.
    :type arguments: dict
    """
    task:str
    arguments:dict=field(default_factory=dict)

    def __post_init__(self):
        self.arguments = {
            (_task_argument_map[(self.task, k)] if (self.task, k) in _task_argument_map else k): v
            for k,v in self.arguments.items()
        }

    def _diagnostics(self):
        if self.task not in eos.tasks._tasks:
            yield Diagnostic(('task',), Severity.ERROR, f'Task \'{self.task}\' is not a valid task')
            return

        task = eos.tasks._tasks[self.task]
        provided_arguments = set(self.arguments)
        # 'analysis_file' and 'base_directory' are provided implicitly when the task is executed.
        provided_arguments.update(('analysis_file', 'base_directory'))

        known_arguments = set(inspect.signature(task).parameters)
        default_arguments = {
            key for key, value in inspect.signature(task).parameters.items()
            if value.default is not inspect.Parameter.empty
        }
        required_arguments = known_arguments - default_arguments

        # This checks TaskComponent.arguments only. StepComponent.default_arguments is checked in
        # StepComponent.validate_semantics().
        for argument in sorted(provided_arguments - known_arguments):
            yield Diagnostic(
                ('arguments', argument),
                Severity.ERROR,
                f'Task \'{self.task}\' does not recognize argument \'{argument}\'',
            )

        for argument in sorted(required_arguments - provided_arguments):
            yield Diagnostic(
                ('arguments', argument),
                Severity.ERROR,
                f'Task \'{self.task}\' requires provision of argument \'{argument}\'',
            )


@dataclass
class StepComponent(_AnalysisFileDeserializable):
    r"""Describes a single step of a reproducible analysis, i.e. one entry of an analysis file's ``steps`` list.

    A step bundles one or more task invocations, may declare dependencies on other steps, and may supply
    default arguments shared by its tasks.

    :param title: A human-readable title for the step.
    :type title: str
    :param id: The unique identifier of the step. Must not contain ``/`` or whitespace.
    :type id: str
    :param tasks: The tasks executed by the step.
    :type tasks: list[TaskComponent]
    :param depends_on: The ids of steps that must complete before this step runs. Optional.
    :type depends_on: list[str]
    :param default_arguments: Default arguments per task name, applied to this step's task invocations. Optional.
    :type default_arguments: dict
    """
    title:str
    id:str
    tasks:list
    depends_on:list=field(default_factory=list)
    default_arguments:defaultdict=field(default_factory=lambda: defaultdict(dict))

    def __post_init__(self):
        self.default_arguments = defaultdict(dict, self.default_arguments)

    def _diagnostics(self):
        yield from _check_file_local_name(self.id, 'step id', ('id',))
        if not self.tasks:
            yield Diagnostic(('tasks',), Severity.ERROR, f'Step \'{self.id}\' has no tasks')

    def validate_structure(self):
        yield from self._diagnostics()
        segments = getattr(self, '_task_segments', list(range(len(self.tasks))))
        yield from _validate_children(self.tasks, segments, 'tasks')

    def validate_semantics(self, context):
        unknown_tasks = self.default_arguments.keys() - eos.tasks._tasks.keys()
        if unknown_tasks:
            yield Diagnostic(
                ('default_arguments',),
                Severity.ERROR,
                f"Step '{self.id}' has default arguments for unknown tasks: {unknown_tasks}",
            )

        for task_component in self.tasks:
            if task_component.task not in eos.tasks._tasks:
                continue
            task = eos.tasks._tasks[task_component.task]
            # TaskComponent.arguments is checked by TaskComponent._diagnostics() in phase 1.
            provided_arguments = self.default_arguments[task_component.task].keys()
            known_arguments = set(inspect.signature(task).parameters)
            for argument in provided_arguments - known_arguments:
                yield Diagnostic(
                    ('default_arguments', task_component.task, argument),
                    Severity.ERROR,
                    f"Task '{task_component.task}' does not recognize argument '{argument}'",
                )

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`StepComponent` from its keyword description.

        Deserializes each entry of ``tasks`` into a :class:`TaskComponent` and normalizes the
        ``default_arguments`` keys against the per-task command-line argument map.

        :returns: The instantiated step component.
        :rtype: StepComponent
        """
        _kwargs = _copy.deepcopy(kwargs)
        if "default_arguments" in kwargs:
            for task, args in kwargs["default_arguments"].items():
                _kwargs["default_arguments"][task] = {
                    (_task_argument_map[(task, k)] if (task, k) in _task_argument_map else k): v
                    for k, v in args.items()
                }
        diagnostics = Deserializable.check_keys(cls, _kwargs)
        if diagnostics:
            return InvalidComponent(diagnostics)
        task_segments = _segments(kwargs["tasks"])
        _kwargs["tasks"] = [TaskComponent.from_dict(**t) for t in kwargs["tasks"]]
        result = cls(**_kwargs)
        result._task_segments = task_segments
        return result


class MaskDescription(_AnalysisFileDeserializable):
    """Polymorphic description of a single entry in a mask's ``description`` list.

    This is a dispatcher rather than a concrete description: :meth:`from_dict` inspects the keys of an
    entry and returns the matching concrete type. An entry containing an ``expression`` key yields a
    :class:`MaskExpressionComponent`, one containing a ``mask_name`` key yields a
    :class:`MaskNamedComponent`, and otherwise the entry yields a :class:`MaskObservableComponent`.
    """

    @staticmethod
    def from_dict(**kwargs):
        """Create the concrete mask-description entry matching the given keyword description.

        :returns: An instance of the concrete :class:`MaskDescription` subtype selected from the keys.
        :rtype: MaskExpressionComponent | MaskNamedComponent | MaskObservableComponent
        """
        if 'expression' in kwargs:
            return Deserializable.make_with_diagnostics(MaskExpressionComponent, **kwargs)
        if 'mask_name' in kwargs:
            return Deserializable.make_with_diagnostics(MaskNamedComponent, **kwargs)
        else:
            return Deserializable.make_with_diagnostics(MaskObservableComponent, **kwargs)

@dataclass
class MaskExpressionComponent(_AnalysisFileDeserializable):
    r"""Describes a mask entry given by a new observable defined through an expression.

    :param expression: The EOS expression that defines the (pseudo-)observable used for masking.
    :type expression: str
    :param name: The qualified name under which the new observable is registered.
    :type name: str
    """
    expression:str
    name:str

    def validate_semantics(self, context):
        try:
            references = eos.analyze_expression(self.expression)
        except RuntimeError as error:
            message = str(error)
            if message.startswith('Parsing Error:'):
                yield Diagnostic(
                    ('expression',),
                    Severity.ERROR,
                    f'Expression does not parse: {message}',
                )
                return
            elif 'is not a valid qualified name:' in message:
                yield Diagnostic(
                    ('expression',),
                    Severity.ERROR,
                    f'Expression contains a malformed qualified name: {message}',
                )
                return
            else:
                raise

        for observable in references.observables:
            yield from _check_qualified(context, observable.full(), 'observable', ('expression',))
        for parameter in references.parameters:
            yield from _check_qualified(context, parameter.full(), 'parameter', ('expression',))

@dataclass
class MaskObservableComponent(_AnalysisFileDeserializable):
    r"""Describes a mask entry given by the name of an existing EOS observable.

    :param name: The qualified name of the observable used for masking.
    :type name: str
    """
    name:str

    def validate_semantics(self, context):
        yield from _check_qualified(context, self.name, 'observable', ('name',))

@dataclass
class MaskNamedComponent(_AnalysisFileDeserializable):
    r"""Describes a mask entry that refers to another, previously defined mask.

    :param mask_name: The name of the previously defined mask to include.
    :type mask_name: str
    """
    mask_name:str

@dataclass
class MaskComponent(_AnalysisFileDeserializable):
    r"""Describes a single named mask, i.e. one entry of an analysis file's ``masks`` list.

    A mask selects a subset of posterior samples by combining one or more (pseudo-)observables and/or
    other masks, keeping samples for which all (``and``) or any (``or``) of them evaluate positive.

    :param name: The unique name of the mask.
    :type name: str
    :param description: The entries that make up the mask.
    :type description: list[MaskDescription]
    :param logical_combination: How the entries are combined, either ``'and'`` or ``'or'``. Defaults to ``'and'``.
    :type logical_combination: str
    """
    name:str
    description:list
    logical_combination: str = "and"

    def __post_init__(self):
        pass

    def _diagnostics(self):
        yield from _check_file_local_name(self.name, 'mask name', ('name',))
        if self.logical_combination not in ('and', 'or'):
            yield Diagnostic(
                ('logical_combination',),
                Severity.ERROR,
                f'Invalid logical combination \'{self.logical_combination}\' for mask \'{self.name}\'',
            )

    def validate_structure(self):
        yield from self._diagnostics()
        segments = getattr(self, '_description_segments', list(range(len(self.description))))
        yield from _validate_children(self.description, segments, 'description')

    def validate_semantics(self, context):
        segments = getattr(self, '_description_segments', list(range(len(self.description))))
        assert len(self.description) == len(segments)
        for description, segment in zip(self.description, segments):
            if hasattr(description, 'validate_semantics'):
                yield from (
                    diagnostic.prefixed('description', segment)
                    for diagnostic in description.validate_semantics(context)
                )

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`MaskComponent` from its keyword description.

        Deserializes each entry of ``description`` via :meth:`MaskDescription.from_dict`.

        :returns: The instantiated mask component.
        :rtype: MaskComponent
        """
        _kwargs = _copy.deepcopy(kwargs)
        diagnostics = Deserializable.check_keys(cls, _kwargs)
        if diagnostics:
            return InvalidComponent(diagnostics)
        description_segments = _segments(kwargs['description'])
        _kwargs['description'] = [MaskDescription.from_dict(**d) for d in kwargs['description']]
        result = cls(**_kwargs)
        result._description_segments = description_segments
        return result

# Maps the discriminating key of a mask-description entry to the corresponding concrete description
# class. The canonical dispatch logic lives in MaskDescription.from_dict; this mapping mirrors it and
# is used to generate the reference documentation of the recognized mask-description types.
MaskDescription.registry = {
    'observable': MaskObservableComponent,
    'expression': MaskExpressionComponent,
    'mask_name':  MaskNamedComponent,
}


@dataclass
class AnalysisFileDescription(_AnalysisFileDeserializable):
    r"""Structured, deserialized representation of a complete analysis file.

    This is the in-memory counterpart of an analysis file's top-level YAML mapping: each field mirrors
    one top-level section, and :meth:`from_dict` deserializes every section into the corresponding
    description objects. It performs only the pure structural deserialization; the mandatory-section
    and cross-reference checks (e.g. that a posterior refers only to defined priors and likelihoods),
    the enforcement of unique step ids and mask names, and the EOS-runtime side effects (declaring
    custom parameters, inserting custom observables) remain the responsibility of :class:`eos.AnalysisFile`.

    Unknown top-level keys are rejected: adding a new section requires adding a field here (and, for a
    backwards-incompatible change, bumping the analysis file format version).

    :param format_version: The schema version of the analysis file. Defaults to 1, the version predating format versioning.
    :type format_version: int
    :param metadata: The analysis file's metadata.
    :type metadata: MetadataDescription
    :param priors: The named priors.
    :type priors: list[PriorComponent]
    :param likelihoods: The named likelihoods.
    :type likelihoods: list[LikelihoodComponent]
    :param posteriors: The named posteriors.
    :type posteriors: list[PosteriorDescription]
    :param predictions: The named predictions.
    :type predictions: list[PredictionDescription]
    :param figures: The figures.
    :type figures: list
    :param observables: The custom observables (from the YAML mapping keyed by observable name).
    :type observables: list[ObservableComponent]
    :param parameters: The custom parameters (from the YAML mapping keyed by parameter name).
    :type parameters: list[ParameterComponent]
    :param steps: The reproducible-analysis steps.
    :type steps: list[StepComponent]
    :param masks: The named sample masks.
    :type masks: list[MaskComponent]
    """
    format_version:int           = 1
    metadata:MetadataDescription = field(default_factory=MetadataDescription)
    priors:list                  = field(default_factory=list)
    likelihoods:list             = field(default_factory=list)
    posteriors:list              = field(default_factory=list)
    predictions:list             = field(default_factory=list)
    figures:list                 = field(default_factory=list)
    observables:list             = field(default_factory=list)
    parameters:list              = field(default_factory=list)
    steps:list                   = field(default_factory=list)
    masks:list                   = field(default_factory=list)

    def validate_structure_children(self):
        yield from self._diagnostics()
        if self.metadata:
            yield from (diagnostic.prefixed('metadata') for diagnostic in self.metadata.validate_structure())

        child_sections = (
            ('priors', self.priors, '_prior_segments'),
            ('likelihoods', self.likelihoods, '_likelihood_segments'),
            ('posteriors', self.posteriors, '_posterior_segments'),
            ('predictions', self.predictions, '_prediction_segments'),
            ('figures', self.figures, '_figure_segments'),
            ('observables', self.observables, '_observable_segments'),
            ('parameters', self.parameters, '_parameter_segments'),
            ('steps', self.steps, '_step_segments'),
            ('masks', self.masks, '_mask_segments'),
        )
        for section, children, attribute in child_sections:
            segments = getattr(self, attribute, list(range(len(children))))
            yield from _validate_children(children, segments, section)

    def validate_structure(self):
        yield from self.validate_structure_children()

        present_sections = getattr(self, '_present_sections', set())
        if 'priors' not in present_sections:
            yield Diagnostic(
                ('priors',),
                Severity.ERROR,
                'Cannot load analysis file: need at least one prior component',
            )
        if 'likelihoods' not in present_sections:
            yield Diagnostic(
                ('likelihoods',),
                Severity.WARNING,
                'No likelihood components found in analysis file',
            )
        if 'posteriors' not in present_sections:
            yield Diagnostic(
                ('posteriors',),
                Severity.ERROR,
                'Cannot load analysis file: need at least one posterior',
            )

        prior_names = {
            prior.name for prior in self.priors
            if not isinstance(prior, InvalidComponent)
        }
        likelihood_names = {
            likelihood.name for likelihood in self.likelihoods
            if not isinstance(likelihood, InvalidComponent)
        }
        posterior_segments = getattr(self, '_posterior_segments', list(range(len(self.posteriors))))
        assert len(self.posteriors) == len(posterior_segments)
        for posterior, segment in zip(self.posteriors, posterior_segments):
            if isinstance(posterior, InvalidComponent):
                continue
            for index, prior in enumerate(posterior.prior):
                if prior not in prior_names:
                    yield Diagnostic(
                        ('posteriors', segment, 'prior', index),
                        Severity.ERROR,
                        f'Posterior \'{posterior.name}\' references prior \'{prior}\' which is not defined',
                    )
            for index, likelihood in enumerate(posterior.likelihood):
                if likelihood not in likelihood_names:
                    yield Diagnostic(
                        ('posteriors', segment, 'likelihood', index),
                        Severity.ERROR,
                        f'Posterior \'{posterior.name}\' references likelihood \'{likelihood}\' which is not defined',
                    )

        step_ids = [
            step.id for step in self.steps
            if not isinstance(step, InvalidComponent)
        ]
        duplicate_step_ids = {
            step_id for step_id, count in Counter(step_ids).items()
            if count > 1
        }
        for step_id in sorted(duplicate_step_ids):
            yield Diagnostic(
                ('steps', step_id, 'id'),
                Severity.ERROR,
                f'Duplicate step id \'{step_id}\'; all steps must have a unique id',
            )

        mask_names = [
            mask.name for mask in self.masks
            if not isinstance(mask, InvalidComponent)
        ]
        duplicate_mask_names = {
            mask_name for mask_name, count in Counter(mask_names).items()
            if count > 1
        }
        for mask_name in sorted(duplicate_mask_names):
            yield Diagnostic(
                ('masks', mask_name, 'name'),
                Severity.ERROR,
                f'Duplicate mask name \'{mask_name}\'; all masks must have a unique name',
            )

        known_masks = set(mask_names)
        mask_segments = getattr(self, '_mask_segments', list(range(len(self.masks))))
        assert len(self.masks) == len(mask_segments)
        for mask, mask_segment in zip(self.masks, mask_segments):
            if isinstance(mask, InvalidComponent):
                continue
            description_segments = getattr(
                mask,
                '_description_segments',
                list(range(len(mask.description))),
            )
            assert len(mask.description) == len(description_segments)
            for description, description_segment in zip(mask.description, description_segments):
                if isinstance(description, MaskNamedComponent) and description.mask_name not in known_masks:
                    yield Diagnostic(
                        ('masks', mask_segment, 'description', description_segment, 'mask_name'),
                        Severity.ERROR,
                        f'Mask {mask.name} references unknown mask {description.mask_name}',
                    )

    def validate_semantics(self, context):
        child_sections = (
            ('priors', self.priors, '_prior_segments'),
            ('likelihoods', self.likelihoods, '_likelihood_segments'),
            ('posteriors', self.posteriors, '_posterior_segments'),
            ('predictions', self.predictions, '_prediction_segments'),
            ('observables', self.observables, '_observable_segments'),
            ('parameters', self.parameters, '_parameter_segments'),
            ('steps', self.steps, '_step_segments'),
            ('masks', self.masks, '_mask_segments'),
        )
        for section, children, attribute in child_sections:
            segments = getattr(self, attribute, list(range(len(children))))
            assert len(children) == len(segments)
            for child, segment in zip(children, segments):
                if hasattr(child, 'validate_semantics'):
                    yield from (
                        diagnostic.prefixed(section, segment)
                        for diagnostic in child.validate_semantics(context)
                    )

        step_ids = {step.id for step in self.steps}
        step_segments = getattr(self, '_step_segments', list(range(len(self.steps))))
        assert len(self.steps) == len(step_segments)
        for step, segment in zip(self.steps, step_segments):
            unknown_dependencies = set(step.depends_on) - step_ids
            if unknown_dependencies:
                yield Diagnostic(
                    ('steps', segment, 'depends_on'),
                    Severity.ERROR,
                    f"Step '{step.id}' depends on unknown steps: {unknown_dependencies}",
                )

            task_segments = getattr(step, '_task_segments', list(range(len(step.tasks))))
            assert len(step.tasks) == len(task_segments)
            for task, task_segment in zip(step.tasks, task_segments):
                if 'posterior' in task.arguments:
                    posterior = task.arguments['posterior']
                    if not context.lookup('posterior', posterior):
                        yield Diagnostic(
                            ('steps', segment, 'tasks', task_segment, 'arguments', 'posterior'),
                            Severity.ERROR,
                            f"Error in step {step.id}: Posterior '{posterior}' not known to EOS",
                        )

        expression_names = [
            description.name
            for mask in self.masks
            for description in mask.description
            if isinstance(description, MaskExpressionComponent)
        ]
        repeated_names = {
            name for name, count in Counter(expression_names).items()
            if count > 1
        }
        for name in repeated_names:
            yield Diagnostic(
                ('masks',),
                Severity.ERROR,
                f"Error in masks: Name '{name}' is used repeatedly",
            )

    @classmethod
    def from_dict(cls, **kwargs):
        r"""Create an :class:`AnalysisFileDescription` from the top-level mapping of an analysis file.

        Each recognized section is deserialized into its corresponding description objects via that
        type's own ``from_dict``. The ``observables`` and ``parameters`` sections are YAML mappings
        keyed by name; the name is injected into each entry so the resulting objects match those built
        elsewhere. Unknown top-level keys are rejected by the dataclass constructor.

        :returns: The deserialized analysis file description.
        :rtype: AnalysisFileDescription
        """
        _kwargs = dict(kwargs) # shallow copy: only the top-level keys are reassigned below
        diagnostics = Deserializable.check_keys(cls, _kwargs)
        if diagnostics:
            return InvalidComponent(diagnostics)
        if 'metadata' in kwargs:
            _kwargs['metadata'] = MetadataDescription.from_dict(**kwargs['metadata'])
        if 'priors' in kwargs:
            prior_segments = _segments(kwargs['priors'])
            _kwargs['priors'] = [PriorComponent.from_dict(**p) for p in kwargs['priors']]
        else:
            prior_segments = []
        if 'likelihoods' in kwargs:
            likelihood_segments = _segments(kwargs['likelihoods'])
            _kwargs['likelihoods'] = [LikelihoodComponent.from_dict(**ll) for ll in kwargs['likelihoods']]
        else:
            likelihood_segments = []
        if 'posteriors' in kwargs:
            posterior_segments = _segments(kwargs['posteriors'])
            _kwargs['posteriors'] = [PosteriorDescription.from_dict(**p) for p in kwargs['posteriors']]
        else:
            posterior_segments = []
        if 'predictions' in kwargs:
            prediction_segments = _segments(kwargs['predictions'])
            _kwargs['predictions'] = [PredictionDescription.from_dict(**p) for p in kwargs['predictions']]
        else:
            prediction_segments = []
        if 'figures' in kwargs:
            figure_segments = _segments(kwargs['figures'])
            _kwargs['figures'] = [FigureFactory.from_dict(**f) for f in kwargs['figures']]
        else:
            figure_segments = []
        if 'observables' in kwargs:
            raw_observables = [
                {'name': name, **description}
                for name, description in kwargs['observables'].items()
            ]
            observable_segments = _segments(raw_observables)
            _kwargs['observables'] = [ObservableComponent.from_dict(**o) for o in raw_observables]
        else:
            observable_segments = []
        if 'parameters' in kwargs:
            raw_parameters = [
                {'name': name, **description}
                for name, description in kwargs['parameters'].items()
            ]
            parameter_segments = _segments(raw_parameters)
            _kwargs['parameters'] = [ParameterComponent.from_dict(**p) for p in raw_parameters]
        else:
            parameter_segments = []
        if 'steps' in kwargs:
            step_segments = _segments(kwargs['steps'], identifier='id')
            _kwargs['steps'] = [StepComponent.from_dict(**s) for s in kwargs['steps']]
        else:
            step_segments = []
        if 'masks' in kwargs:
            mask_segments = _segments(kwargs['masks'])
            _kwargs['masks'] = [MaskComponent.from_dict(**m) for m in kwargs['masks']]
        else:
            mask_segments = []
        result = cls(**_kwargs)
        result._prior_segments = prior_segments
        result._likelihood_segments = likelihood_segments
        result._posterior_segments = posterior_segments
        result._prediction_segments = prediction_segments
        result._figure_segments = figure_segments
        result._observable_segments = observable_segments
        result._parameter_segments = parameter_segments
        result._step_segments = step_segments
        result._mask_segments = mask_segments
        result._present_sections = set(kwargs)
        return result

# AnalysisFile schema

# dict with keys:
#   format_version (optional): int, the schema version of the file; defaults to 1 (the version
#       predating format versioning). EOS refuses to load a file declaring a newer version than
#       it supports.
#   metadata (optional)
#   priors (mandatory)
#   likelihoods (mandatory)
#   posteriors (mandatory)
#   figures (optional)
#   observables (optional)
#   predictions (optional)
#   parameters (optional)
#   steps (optional)
#   masks (optional)


# metadata schema:
# dict with keys:
#   title (optional): string
#   id (optional): string
#   authors (optional): list of dicts, each with keys:
#       name (mandatory): string
#       affiliation (optional): string
#       email (optional): string


# priors schema:
# list of dicts, each with keys:
#   name (mandatory): string
#   descriptions (mandatory) (or parameters in depreciated style)
#       list of dicts, each with keys:
#           either:
#             parameter (mandatory): string
#             type (mandatory): string (uniform, flat, scale, gauss, gaussian, poisson)
#             .... other keys depend on type
#           or:
#             constraint (mandatory): string


# likelihoods schema:
# list of dicts, each with keys:
#   name (mandatory): string
#   constraints (optional): list of strings
#   manual_constraints (optional): dict where each key is the name and the corresponding value is a dict
#   pyhf (optional): dict with keys:
        #  file (mandatory): string
        #  parameter_map (optional): dict
        #  ???????


# posteriors schema:
# list of dicts, each with keys:
#  name (mandatory): string
#  prior (mandatory): list of strings (corresponding to prior names)
#  likelihood (mandatory): list of strings (corresponding to likelihood names)
#  global_options (optional): dict
#  fixed_parameters (optional): dict


# figures schema:
# list of dicts, each dict describing a valid type of figure;
# see eos/figure/figure.py for the schema of each figure type


# observables schema:
# dict, where each key is the observable name and the corresponding value is a dict with keys:
#    latex (mandatory), str
#    unit (mandatory), str
#    options (optional), dict
#    expression (mandatory), str


# predictions schema:
# list of dicts, each with keys:
#  name (mandatory): string
#  observables (mandatory): list of dicts, each with keys:
#       name (mandatory): string
#       kinematics (optional): dict or list of dicts
#       options (optional): dict
#  global_options (optional): dict
#  fixed_parameters (optional): dict


# parameters schema:
# dict, where each key is a new EOS qualified name and corresponding value is a dict with keys:
#   latex (mandatory), str
#   unit (mandatory), str
#   central (mandatory), float
#   min (mandatory), float
#   max (mandatory), float
#   alias_of (optional), list of strings (each string is a valid EOS qualified name)


# steps schema:
# list of dicts, each with keys:
#  title (mandatory): string
#  id (mandatory): unique string
#  depends_on (optional): list of strings (each string is a step id)
#  default_arguments (optional): dict, whose keys are task names, and values are dicts with argument names and values for that task
#  tasks (mandatory): list of dicts, each with keys:
    #  task (mandatory) : string
    #  arguments (optional): dict, whose keys are arguments for the task


# masks schema:
# list of dicts, each with keys:
#  name (mandatory): unique string
#  logical_combination (optional): str, either 'and' or 'or', defaults to 'and'
#  description (mandatory): list of dicts, with keys:
#    either:
#      name: str, valid EOS observable name
#    or:
#      name: str, valid EOS observable name
#      expression: str, valid EOS observable expression
#    or:
#      mask_name: a previously defined mask name
