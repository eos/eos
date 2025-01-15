from .deserializable import Deserializable
from dataclasses import dataclass, field
import copy as _copy
import eos


@dataclass
class PublicLikelihoodDescription(Deserializable):
    """
    Represents a single EOS likelihood.
    """

    filename:str
    filetype:str

    def __post_init__(self):
        VALID_FILETYPES = ['MixtureDensity', 'NabuLikelihood']
        if not self.filetype in VALID_FILETYPES:
            raise ValueError(f'Invalid likelihood file type: {self.filetype}')


@dataclass
class DataSetDescription(Deserializable):
    """
    Represents a single EOS data set.
    """

    authors:list[str]
    title:str
    keywords:list[str]
    eos_version:str
    likelihoods:dict[str, PublicLikelihoodDescription]
    doi:str=None

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs['likelihoods'] = { k: PublicLikelihoodDescription.from_dict(**v) for k, v in kwargs['likelihoods'].items() }
        return Deserializable.make(cls, **_kwargs)
