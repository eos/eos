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


# AnalysisFile schema

# dict with keys:
#   priors (mandatory)
#   likelihoods (mandatory)
#   posteriors (mandatory)
#   observables (optional)
#   predictions (optional)
#   parameters (optional)
#   steps (optional)


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
