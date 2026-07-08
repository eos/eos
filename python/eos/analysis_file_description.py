from .deserializable import Deserializable
from dataclasses import dataclass, field
from collections import defaultdict
import copy as _copy
import eos
import inspect
import os

@dataclass(kw_only=True)
class AnalysisFileContext:
    """
    Auxiliary class used to resolve paths relative to an analysis file's base directory.

    Instances are passed to consumers of an analysis file (e.g. figures and tasks) so that relative data
    paths can be resolved against a common base directory.

    :param base_directory: The base directory against which relative paths are resolved. Defaults to the current working directory.
    :type base_directory: str
    """
    base_directory:str='./'

    def __post_init__(self):
        if not os.path.exists(self.base_directory):
            raise ValueError(f'Base directory \'{self.base_directory}\' does not exist')

        if not os.path.isdir(self.base_directory):
            raise ValueError(f'Base directory \'{self.base_directory}\' is not a directory')

    def data_path(self, relative_path:str):
        """Resolve a path relative to :attr:`base_directory` into an absolute path.

        :param relative_path: A path relative to the base directory.
        :type relative_path: str
        :returns: The corresponding absolute path.
        :rtype: str
        """
        return os.path.abspath(os.path.join(self.base_directory, relative_path))


@dataclass
class MetadataAuthorDescription(Deserializable):
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
class MetadataDescription(Deserializable):
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

    @staticmethod
    def from_dict(**kwargs):
        """Create a :class:`MetadataDescription` from its keyword description.

        Deserializes the nested ``authors`` list into :class:`MetadataAuthorDescription` instances.

        :returns: The instantiated metadata description.
        :rtype: MetadataDescription
        """
        _kwargs = _copy.deepcopy(kwargs)
        if 'authors' in kwargs:
            _kwargs['authors'] = [MetadataAuthorDescription.from_dict(**a) for a in kwargs['authors']]
        return Deserializable.make(MetadataDescription, **_kwargs)


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
        :raises ValueError: If the keys do not identify a known type of prior description.
        """
        if 'type' in kwargs:
            _kwargs = _copy.deepcopy(kwargs)
            _kwargs.pop('type')
            if kwargs['type'] in ("constraint",):
                return Deserializable.make(ConstraintPriorDescription, **_kwargs)
            elif kwargs['type'] in ("uniform", "flat"):
                return Deserializable.make(UniformPriorDescription, **_kwargs)
            elif kwargs['type'] in ("scale",):
                return Deserializable.make(ScalePriorDescription, **_kwargs)
            elif kwargs['type'] in ("gauss", "gaussian"):
                if ('min' in kwargs) != ('max' in kwargs):
                    raise ValueError('A Gaussian prior description must contain both \'min\' and \'max\' or neither')
                if 'min' in kwargs:
                    return Deserializable.make(CurtailedGaussianDescription, **_kwargs)
                return Deserializable.make(GaussianPriorDescription, **_kwargs)
            elif kwargs['type'] in ("poisson",):
                return Deserializable.make(PoissonPriorDescription, **_kwargs)
            elif kwargs['type'] in ("transform",):
                return Deserializable.make(TransformPriorDescription, **_kwargs)
        elif 'constraint' in kwargs:
            eos.warn('A constraint prior description without a \'type\' key is deprecated and will be removed in a future version of EOS; add \'type: constraint\' instead')
            return Deserializable.make(ConstraintPriorDescription, **kwargs)

        raise ValueError('Unknown type of prior description')

@dataclass
class PoissonPriorDescription(Deserializable):
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
class CurtailedGaussianDescription(Deserializable):
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
        eos.warn('The curtailed Gaussian prior description (a \'gauss\'/\'gaussian\' prior with \'min\' and \'max\') is deprecated and will be removed in a future version of EOS; use \'type: gaussian\' without the \'min\'/\'max\' keys instead')

@dataclass
class GaussianPriorDescription(Deserializable):
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
class ScalePriorDescription(Deserializable):
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
class UniformPriorDescription(Deserializable):
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
class ConstraintPriorDescription(Deserializable):
    r"""Describes a (possibly correlated, multivariate) prior taken from a built-in EOS constraint.

    :param constraint: The qualified name of the EOS constraint used as the prior.
    :type constraint: str
    """
    constraint:str
    type:str=field(repr=False, init=False, default="constraint")

@dataclass
class TransformPriorDescription(Deserializable):
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
class PriorComponent(Deserializable):
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
            _kwargs['descriptions'] = [PriorDescription.from_dict(**d) for d in kwargs['descriptions']]
        if "parameters" in kwargs:
            eos.warn(f'\'parameters\' is in the description of prior component \'{kwargs["name"]}\', use \'descriptions\' instead')
            _kwargs.pop("parameters")
            _kwargs['descriptions'] = [PriorDescription.from_dict(**d) for d in kwargs['parameters']]
        return Deserializable.make(cls, **_kwargs)



@dataclass
class ConstraintLikelihoodDescription(Deserializable):
    r"""Describes a likelihood contribution taken from a built-in EOS constraint.

    :param constraint: The qualified name of the EOS constraint.
    :type constraint: eos.QualifiedName
    """
    constraint:eos.QualifiedName

@dataclass
class ManualConstraintDescription(Deserializable):
    r"""Describes a likelihood contribution from a constraint defined inline in the analysis file.

    :param name: The qualified name under which the manual constraint is registered.
    :type name: eos.QualifiedName
    :param info: The constraint definition, in the same format as a built-in EOS constraint entry.
    :type info: dict
    """
    name:eos.QualifiedName
    info:dict

@dataclass
class PyHFConstraintDescription(Deserializable):
    r"""Describes a likelihood contribution imported from a pyhf workspace.

    :param file: The path to the JSON file specifying the pyhf workspace.
    :type file: str
    :param parameter_map: An optional mapping from pyhf parameter names to EOS observables or parameters; see :class:`eos.PyhfLogLikelihood`.
    :type parameter_map: dict
    """
    file:str
    parameter_map:dict=field(default_factory=dict)

@dataclass
class LikelihoodComponent(Deserializable):
    r"""Describes a single named likelihood, i.e. one entry of an analysis file's ``likelihoods`` list.

    A named likelihood may combine any number of built-in constraints, inline (manual) constraints, and
    pyhf-based contributions; at least one of the three must be present.

    :param name: The unique name of the likelihood, by which posteriors refer to it.
    :type name: str
    :param constraints: The built-in EOS constraints contributing to the likelihood. Optional.
    :type constraints: list[ConstraintLikelihoodDescription]
    :param manual_constraints: The inline constraint definitions contributing to the likelihood. Optional.
    :type manual_constraints: list[ManualConstraintDescription]
    :param pyhf: The pyhf-based contribution to the likelihood. Optional.
    :type pyhf: PyHFConstraintDescription
    """
    name:str
    constraints:list=field(default_factory=list)
    manual_constraints:list=field(default_factory=list)
    pyhf:list=field(default_factory=list)

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`LikelihoodComponent` from its keyword description.

        Deserializes the ``constraints``, ``manual_constraints``, and ``pyhf`` entries into their
        respective description types.

        :returns: The instantiated likelihood component.
        :rtype: LikelihoodComponent
        :raises ValueError: If none of ``constraints``, ``manual_constraints``, or ``pyhf`` is provided.
        """
        if not ('constraints' in kwargs or 'manual_constraints' in kwargs or 'pyhf' in kwargs):
            raise ValueError('LikelihoodComponent must have at least one of constraints, manual_constraints, or pyhf')
        _kwargs = _copy.deepcopy(kwargs)
        if 'constraints' in kwargs:
            _kwargs['constraints'] = [ConstraintLikelihoodDescription.from_dict(constraint=c) for c in kwargs['constraints']]
        if 'manual_constraints' in kwargs:
            _kwargs['manual_constraints'] = [ManualConstraintDescription.from_dict(name=n, info=d) for n, d in kwargs['manual_constraints'].items()]
        if 'pyhf' in kwargs:
            _kwargs['pyhf'] = PyHFConstraintDescription.from_dict(**kwargs['pyhf'])
        return Deserializable.make(cls, **_kwargs)



@dataclass
class PosteriorDescription(Deserializable):
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



@dataclass
class ObservableComponent(Deserializable):
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



@dataclass
class PredictionObservableComponent(Deserializable):
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
class PredictionDescription(Deserializable):
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

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`PredictionDescription` from its keyword description.

        Deserializes each entry of ``observables`` into a :class:`PredictionObservableComponent`.

        :returns: The instantiated prediction description.
        :rtype: PredictionDescription
        """
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs["observables"] = [PredictionObservableComponent.from_dict(**o) for o in kwargs["observables"]]
        return Deserializable.make(cls, **_kwargs)



@dataclass
class ParameterComponent(Deserializable):
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
class TaskComponent(Deserializable):
    r"""Describes a single task invocation within a step, i.e. one entry of a step's ``tasks`` list.

    :param task: The name of the task to run; must be a registered EOS task.
    :type task: str
    :param arguments: The arguments passed to the task, given as a mapping from argument name to value. Optional.
    :type arguments: dict
    """
    task:str
    arguments:dict=field(default_factory=dict)

    def __post_init__(self):
        if self.task not in eos.tasks._tasks.keys():
            raise ValueError(f'Task \'{self.task}\' is not a valid task')
        task = eos.tasks._tasks[self.task]

        self.arguments = {
            (_task_argument_map[(self.task, k)] if (self.task, k) in _task_argument_map else k): v
            for k,v in self.arguments.items()
        }

        provided_arguments = set(self.arguments.keys())
        # inject 'analysis_file' and 'base_directory' into known arguments
        # these will be provided implicitly when the task is executed
        provided_arguments.add('analysis_file')
        provided_arguments.add('base_directory')

        known_arguments = set(inspect.signature(task).parameters.keys())
        default_arguments = { k for k, v in inspect.signature(task).parameters.items() if v.default is not inspect.Parameter.empty }
        required_arguments = known_arguments - default_arguments

        for arg in provided_arguments - known_arguments:
            raise ValueError(f'Task \'{self.task}\' does not recognize argument \'{arg}\'')

        for arg in required_arguments - provided_arguments:
            raise ValueError(f'Task \'{self.task}\' requires provision of argument \'{arg}\'')

        for arg in default_arguments - provided_arguments:
            eos.warn(f'Task \'{self.task}\' has a default value for argument \'{arg}\', which can change across versions; consider providing a value explicitly')


@dataclass
class StepComponent(Deserializable):
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
        if '/' in self.id:
            raise ValueError(f'Invalid character \'/\' in step id \'{self.id}\'')
        if ' ' in self.id:
            raise ValueError(f'Invalid character \' \' in step id \'{self.id}\'')
        if len(self.tasks) == 0:
            raise ValueError(f'Step \'{self.id}\' has no tasks')
        self.default_arguments = defaultdict(dict, self.default_arguments)

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
        _kwargs["tasks"] = [TaskComponent.from_dict(**t) for t in kwargs["tasks"]]
        return Deserializable.make(cls, **_kwargs)


class MaskDescription(Deserializable):
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
            return Deserializable.make(MaskExpressionComponent, **kwargs)
        if 'mask_name' in kwargs:
            return Deserializable.make(MaskNamedComponent, **kwargs)
        else:
            return Deserializable.make(MaskObservableComponent, **kwargs)

@dataclass
class MaskExpressionComponent(Deserializable):
    r"""Describes a mask entry given by a new observable defined through an expression.

    :param expression: The EOS expression that defines the (pseudo-)observable used for masking.
    :type expression: str
    :param name: The qualified name under which the new observable is registered.
    :type name: str
    """
    expression:str
    name:str

@dataclass
class MaskObservableComponent(Deserializable):
    r"""Describes a mask entry given by the name of an existing EOS observable.

    :param name: The qualified name of the observable used for masking.
    :type name: str
    """
    name:str

@dataclass
class MaskNamedComponent(Deserializable):
    r"""Describes a mask entry that refers to another, previously defined mask.

    :param mask_name: The name of the previously defined mask to include.
    :type mask_name: str
    """
    mask_name:str

@dataclass
class MaskComponent(Deserializable):
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
        if self.logical_combination not in ['and', 'or']:
            raise ValueError(f'Invalid logical combination \'{self.logical_combination}\' for mask \'{self.name}\'')

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`MaskComponent` from its keyword description.

        Deserializes each entry of ``description`` via :meth:`MaskDescription.from_dict`.

        :returns: The instantiated mask component.
        :rtype: MaskComponent
        """
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs['description'] = [MaskDescription.from_dict(**d) for d in kwargs['description']]
        return Deserializable.make(cls, **_kwargs)

# Maps the discriminating key of a mask-description entry to the corresponding concrete description
# class. The canonical dispatch logic lives in MaskDescription.from_dict; this mapping mirrors it and
# is used to generate the reference documentation of the recognized mask-description types.
MaskDescription.registry = {
    'observable': MaskObservableComponent,
    'expression': MaskExpressionComponent,
    'mask_name':  MaskNamedComponent,
}

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
