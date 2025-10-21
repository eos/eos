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
    Auxiliary class to transfer information to consumers of the analysis file.
    """
    base_directory:str='./'

    def __post_init__(self):
        if not os.path.exists(self.base_directory):
            raise ValueError(f'Base directory \'{self.base_directory}\' does not exist')

        if not os.path.isdir(self.base_directory):
            raise ValueError(f'Base directory \'{self.base_directory}\' is not a directory')

    def data_path(self, relative_path:str):
        return os.path.abspath(os.path.join(self.base_directory, relative_path))


@dataclass
class MetadataAuthorDescription(Deserializable):
    name:str
    affiliation:str=''
    email:str=''


@dataclass
class MetadataDescription(Deserializable):
    title:str=''
    id:str=''
    authors:list[MetadataAuthorDescription]=field(default_factory=list)

    @staticmethod
    def from_dict(**kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        if 'authors' in kwargs:
            _kwargs['authors'] = [MetadataAuthorDescription.from_dict(**a) for a in kwargs['authors']]
        return Deserializable.make(MetadataDescription, **_kwargs)


class PriorDescription:
    @staticmethod
    def from_dict(**kwargs):
        if 'constraint' in kwargs:
            return Deserializable.make(ConstraintPriorDescription, **kwargs)
        elif 'parameter' in kwargs and 'type' in kwargs:
            _kwargs = _copy.deepcopy(kwargs)
            _kwargs.pop('type')
            if kwargs['type'] in ("uniform", "flat"):
                return Deserializable.make(UniformPriorDescription, **_kwargs)
            elif kwargs['type'] in ("scale",):
                return Deserializable.make(ScalePriorDescription, **_kwargs)
            elif kwargs['type'] in ("gauss", "gaussian"):
                if "min" in kwargs:
                    return Deserializable.make(CurtailedGaussianDescription, **_kwargs)
                return Deserializable.make(GaussianPriorDescription, **_kwargs)
            elif kwargs['type'] in ("poisson",):
                return Deserializable.make(PoissonPriorDescription, **_kwargs)
        elif 'parameters' in kwargs:
            _kwargs = _copy.deepcopy(kwargs)
            _kwargs.pop('type')
            if kwargs['type'] in ("transform"):
                return Deserializable.make(TransformPriorDescription, **_kwargs)

        raise ValueError('Unknown type of prior description')

@dataclass
class PoissonPriorDescription(Deserializable):
    parameter:str
    k:float
    type:str=field(repr=False, init=False, default="poisson")

@dataclass
class CurtailedGaussianDescription(Deserializable):
    parameter:str
    central:float
    sigma:float
    min:float
    max:float
    type:str=field(repr=False, init=False, default="gaussian")

@dataclass
class GaussianPriorDescription(Deserializable):
    parameter:str
    central:float
    sigma:float
    type:str=field(repr=False, init=False, default="gaussian")


@dataclass
class ScalePriorDescription(Deserializable):
    parameter:str
    min:float
    max:float
    mu_0:float
    lambda_scale:float
    type:str=field(repr=False, init=False, default="scale")

@dataclass
class UniformPriorDescription(Deserializable):
    parameter:str
    min:float
    max:float
    type:str=field(repr=False, init=False, default="uniform")

@dataclass
class ConstraintPriorDescription(Deserializable):
     constraint:str

@dataclass
class TransformPriorDescription(Deserializable):
    parameters:list[str]
    shift:list[float]
    transform:list[list[float]]
    min:list[float]
    max:list[float]
    type:str=field(repr=False, init=False, default="transform")

@dataclass
class PriorComponent(Deserializable):
    name:str
    descriptions:list

    @classmethod
    def from_dict(cls, **kwargs):
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
    constraint:eos.QualifiedName

@dataclass
class ManualConstraintDescription(Deserializable):
    name:eos.QualifiedName
    info:dict

@dataclass
class PyHFConstraintDescription(Deserializable):
    file:str
    parameter_map:dict=field(default_factory=dict)

@dataclass
class LikelihoodComponent(Deserializable):
    name:str
    constraints:list=field(default_factory=list)
    manual_constraints:list=field(default_factory=list)
    pyhf:list=field(default_factory=list)

    @classmethod
    def from_dict(cls, **kwargs):
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
    name:str
    prior:list
    likelihood:list
    global_options:dict=field(default_factory=dict)
    fixed_parameters:dict=field(default_factory=dict)



@dataclass
class ObservableComponent(Deserializable):
    name:str
    latex:str
    unit:str
    expression:str
    options:dict=field(default_factory=dict)



@dataclass
class PredictionObservableComponent(Deserializable):
    name:str
    kinematics:dict=field(default_factory=dict) # TODO: once we bump the minimum python version to 3.10, use dict | list[dict] instead
    options:dict=field(default_factory=dict)

@dataclass
class PredictionDescription(Deserializable):
    name:str
    observables:list
    global_options:dict=field(default_factory=dict)
    fixed_parameters:dict=field(default_factory=dict)

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs["observables"] = [PredictionObservableComponent.from_dict(**o) for o in kwargs["observables"]]
        return Deserializable.make(cls, **_kwargs)



@dataclass
class ParameterComponent(Deserializable):
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
    @staticmethod
    def from_dict(**kwargs):
        if 'expression' in kwargs:
            return Deserializable.make(MaskExpressionComponent, **kwargs)
        if 'mask_name' in kwargs:
            return Deserializable.make(MaskNamedComponent, **kwargs)
        else:
            return Deserializable.make(MaskObservableComponent, **kwargs)

@dataclass
class MaskExpressionComponent(Deserializable):
    expression:str
    name:str

@dataclass
class MaskObservableComponent(Deserializable):
    name:str

@dataclass
class MaskNamedComponent(Deserializable):
    mask_name:str

@dataclass
class MaskComponent(Deserializable):
    name:str
    description:dict
    logical_combination: str = "and"

    def __post_init__(self):
        if self.logical_combination not in ['and', 'or']:
            raise ValueError(f'Invalid logical combination \'{self.logical_combination}\' for mask \'{self.name}\'')

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs['description'] = [MaskDescription.from_dict(**d) for d in kwargs['description']]
        return Deserializable.make(cls, **_kwargs)

# AnalysisFile schema

# dict with keys:
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
