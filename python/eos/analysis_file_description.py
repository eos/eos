from .deserializable import Deserializable
from dataclasses import dataclass, field
import copy as _copy
import eos

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



@dataclass
class TaskComponent(Deserializable):
    task:str
    arguments:dict=field(default_factory=dict)

@dataclass
class StepComponent(Deserializable):
    name:str
    tasks:list
    iterations:list=field(default_factory=list)
    depends_on:list=field(default_factory=list)

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs["tasks"] = [TaskComponent.from_dict(**t) for t in kwargs["tasks"]]
        if "depends_on" in kwargs:
            _kwargs["depends_on"] = kwargs["depends-on"]
        return Deserializable.make(cls, **_kwargs)


@dataclass
class SampleNestedMatchDescription(Deserializable):
    posterior:str
    bound:str = None
    nlive:int = None
    dlogz:float = None
    maxiter:int = None
    seed:int = None
    base_directory:str = None

@dataclass
class SampleNestedSetDescription(Deserializable):
    bound:str
    nlive:int
    dlogz:float
    maxiter:int
    seed:int
    base_directory:str

@dataclass
class SampleMCMCMatchDescription(Deserializable):
    posterior:str
    chain:int
    N:int = None
    stride:int = None
    preruns:int = None
    pre_N:int = None
    start_point:tuple[float] = None
    cov_scale:float = None
    base_directory:str = None

@dataclass
class SampleMCMCSetDescription(Deserializable):
    N:int
    stride:int
    preruns:int
    pre_N:int
    start_point:tuple[float]
    cov_scale:float
    base_directory:str

@dataclass
class MixtureProductMatchDescription(Deserializable):
    posterior:str
    posteriors:tuple[str]
    base_directory:str = None

@dataclass
class MixtureProductSetDescription(Deserializable):
    base_directory:str

@dataclass
class FindClustersMatchDescription(Deserializable):
    posterior:str
    threshold:float = None
    K_g:int = None
    base_directory:str = None

@dataclass
class FindClustersSetDescription(Deserializable):
    threshold:float
    K_g:int
    base_directory:str

@dataclass
class PredictObservablesMatchDescription(Deserializable):
    posterior:str
    prediction:str
    begin:int = None
    end:int = None
    base_directory:str = None

@dataclass
class PredictObservablesSetDescription(Deserializable):
    begin:int
    end:int
    base_directory:str

dataclass
class CornerPlotMatchDescription(Deserializable):
    posterior:str
    begin:int = None
    end:int = None
    base_directory:str = None

@dataclass
class CornerPlotSetDescription(Deserializable):
    begin:int
    end:int
    base_directory:str

@dataclass
class SamplePMCMatchDescription(Deserializable):
    posterior:str
    step_N:int = None
    steps:int = None
    perplexity_threshold:float = None
    weight_threshold:float = None
    final_N:int = None
    pmc_iterations:int = None
    pmc_rel_tol:float = None
    pmc_abs_tol:float = None
    pmc_lookback:int = None
    initial_proposal:str = None
    sigma_test_stat:tuple[float] = None
    base_directory:str = None

@dataclass
class SamplePMCSetDescription(Deserializable):
    step_N:int
    steps:int
    perplexity_threshold:float
    weight_threshold:float
    final_N:int
    pmc_iterations:int
    pmc_rel_tol:float
    pmc_abs_tol:float
    pmc_lookback:int
    initial_proposal:str
    sigma_test_stat:tuple[float]
    base_directory:str

@dataclass
class PlotSamplesMatchDescription(Deserializable):
    posterior:str
    bins:int = None
    base_directory:str = None

@dataclass
class PlotSamplesSetDescription(Deserializable):
    bins:int
    base_directory:str

@dataclass
class FindModeMatchDescription(Deserializable):
    posterior:str
    optimizations:int = None
    importance_samples:bool = None
    label:str = None
    base_directory:str = None

@dataclass
class FindModeSetDescription(Deserializable):
    optimizations:int
    importance_samples:bool
    label:str
    base_directory:str

@dataclass
class ConfigurationComponent(Deserializable):
    task:str
    match:dict
    set:dict

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)

        if 'task' not in kwargs or 'match' not in kwargs or 'set' not in kwargs:
            raise ValueError('ConfigurationComponent must have all of task, match and set')

        task = kwargs["task"]

        # Check: nothing in match is in set
        if (overlap := kwargs["match"].keys() & kwargs["set"].keys()): # Check for any intersection
            raise ValueError(f'Arguments {overlap} are in both match and set of the configuration of task {task}')

        if task == "sample-nested":
            _kwargs["match"] = SampleNestedMatchDescription.from_dict(**kwargs["match"])
            # Remove the mandatory argument, then add anything else to "set" so it is complete
            kwargs["match"].pop("posterior")
            kwargs["set"].update(kwargs["match"])
            _kwargs["set"] = SampleNestedSetDescription.from_dict(**kwargs["set"])
        elif task == "sample-pmc":
            if "sigma_test_stat" in kwargs["match"]:
                kwargs["match"]["sigma_test_stat"] = tuple(kwargs["match"]["sigma_test_stat"])
            _kwargs["match"] = SamplePMCMatchDescription.from_dict(**kwargs["match"])
            # Remove the mandatory argument, then add anything else to "set" so it is complete
            kwargs["match"].pop("posterior")
            kwargs["set"].update(kwargs["match"])
            _kwargs["set"] = SamplePMCSetDescription.from_dict(**kwargs["set"])
        elif task == "plot-samples":
            _kwargs["match"] = PlotSamplesMatchDescription.from_dict(**kwargs["match"])
            # Remove the mandatory argument, then add anything else to "set" so it is complete
            kwargs["match"].pop("posterior")
            kwargs["set"].update(kwargs["match"])
            _kwargs["set"] = PlotSamplesSetDescription.from_dict(**kwargs["set"])
        elif task == "find-mode":
            _kwargs["match"] = FindModeMatchDescription.from_dict(**kwargs["match"])
            # Remove the mandatory argument, then add anything else to "set" so it is complete
            kwargs["match"].pop("posterior")
            kwargs["set"].update(kwargs["match"])
            _kwargs["set"] = FindModeSetDescription.from_dict(**kwargs["set"])
        elif task == "corner-plot":
            _kwargs["match"] = CornerPlotMatchDescription.from_dict(**kwargs["match"])
            # Remove the mandatory argument, then add anything else to "set" so it is complete
            kwargs["match"].pop("posterior")
            kwargs["set"].update(kwargs["match"])
            _kwargs["set"] = CornerPlotSetDescription.from_dict(**kwargs["set"])
        elif task == "sample-mcmc":
            if "start_point" in kwargs["match"]:
                kwargs["match"]["start_point"] = tuple(kwargs["match"]["start_point"])
            _kwargs["match"] = SampleMCMCMatchDescription.from_dict(**kwargs["match"])
            # Remove the mandatory argument, then add anything else to "set" so it is complete
            kwargs["match"].pop("posterior")
            kwargs["match"].pop("chain")
            kwargs["set"].update(kwargs["match"])
            _kwargs["set"] = SampleMCMCSetDescription.from_dict(**kwargs["set"])
        elif task == "mixture-product":
            kwargs["match"]["posteriors"] = tuple(kwargs["match"]["posteriors"])
            _kwargs["match"] = MixtureProductMatchDescription.from_dict(**kwargs["match"])
            # Remove the mandatory argument, then add anything else to "set" so it is complete
            kwargs["match"].pop("posterior")
            kwargs["match"].pop("posteriors")
            kwargs["set"].update(kwargs["match"])
            _kwargs["set"] = MixtureProductSetDescription.from_dict(**kwargs["set"])
        elif task == "find-clusters":
            _kwargs["match"] = FindClustersMatchDescription.from_dict(**kwargs["match"])
            # Remove the mandatory argument, then add anything else to "set" so it is complete
            kwargs["match"].pop("posterior")
            kwargs["set"].update(kwargs["match"])
            _kwargs["set"] = FindClustersSetDescription.from_dict(**kwargs["set"])
        elif task == "predict-observables":
            _kwargs["match"] = PredictObservablesMatchDescription.from_dict(**kwargs["match"])
            # Remove the mandatory argument, then add anything else to "set" so it is complete
            kwargs["match"].pop("posterior")
            kwargs["match"].pop("prediction")
            kwargs["set"].update(kwargs["match"])
            _kwargs["set"] = PredictObservablesSetDescription.from_dict(**kwargs["set"])
        else:
            raise ValueError(f'Configuration for unknown task : {task}')

        return Deserializable.make(cls, **_kwargs)


# AnalysisFile schema

# dict with keys:
#   priors (mandatory)
#   likelihoods (mandatory)
#   posteriors (mandatory)
#   observables (optional)
#   predictions (optional)
#   parameters (optional)
#   steps (optional)
#   configuration (optional)


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
#  name (mandatory): string
#  tasks (mandatory): list of dicts, each with keys:
    #  task (mandatory) : string
    #  arguments (optional): dict
#  iterations (optional): list of dicts
#  depends-on (optional): list of strings ???????
#  ???????

# configuration schema:
# dict, where each key is a EOS task, and the corresponding value is a list of dicts each with keys:
#   match (mandatory): dict, whose keys are arguments of the task and values are regex-ish strings
#   set (mandatory): dict, whose keys are options of the task and values are the corresponding values
