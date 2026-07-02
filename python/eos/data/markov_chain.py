# Copyright (c) 2020-2026 Danny van Dyk
# Copyright (c) 2023 Stephan Kürten
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

from dataclasses import dataclass, field, asdict
from eos.data._common import ParameterDescription, load_array
from eos.deserializable import Deserializable
from eos.serializable import Serializable

import copy as _copy
import eos
import os
import numpy as _np


@dataclass(kw_only=True)
class MarkovChainDescription(Serializable, Deserializable):
    r"""Schema of the ``description.yaml`` written for a :class:`MarkovChain`.

    This is the single source of truth for the metadata stored alongside a Markov chain, and it backs
    both the read path (:meth:`from_yaml_file`) and the write path (:meth:`to_yaml_file`). The
    ``type`` discriminator is validated on deserialization and is not part of the constructor.

    :param version: The version of EOS that wrote the data object.
    :type version: str
    :param parameters: The descriptions of the varied parameters.
    :type parameters: list[ParameterDescription]
    :param has_weights: Whether the chain stores importance weights (on-disk key ``has-weights``).
    :type has_weights: bool
    """
    version:str
    parameters:list[ParameterDescription]
    has_weights:bool = field(default=False)
    type:str = field(init=False, default='MarkovChain')

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`MarkovChainDescription` from its on-disk keyword description.

        Validates the ``type`` discriminator, normalizes the on-disk ``has-weights`` key to the
        ``has_weights`` field, and deserializes each entry of ``parameters`` into a
        :class:`ParameterDescription`.

        :raises ValueError: If the description does not identify a Markov chain.
        """
        _kwargs = _copy.deepcopy(kwargs)

        _type = _kwargs.pop('type', None)
        if _type != 'MarkovChain':
            raise ValueError(f'Expected a description of type \'MarkovChain\', got \'{_type}\'')

        if 'has-weights' in _kwargs:
            _kwargs['has_weights'] = _kwargs.pop('has-weights')

        if 'parameters' in _kwargs:
            _kwargs['parameters'] = [ParameterDescription.from_dict(**p) for p in _kwargs['parameters']]

        return Deserializable.make(cls, **_kwargs)

    def to_dict(self):
        """Serialize this description into the on-disk mapping written to ``description.yaml``.

        Emits the on-disk ``has-weights`` key and the ``type`` discriminator, inverting
        :meth:`from_dict`.
        """
        return {
            'version':     self.version,
            'type':        self.type,
            'parameters':  [asdict(p) for p in self.parameters],
            'has-weights': self.has_weights,
        }


class MarkovChain:
    r"""Represents the parameter samples of a single Markov chain stored on disk.

    A Markov chain bundles samples drawn in parameter space together with their counterparts in the
    unit-hypercube ``u`` space and, optionally, importance weights. Instances are created either by
    reading an existing chain from disk (passing its ``path`` to the constructor) or by writing a new
    chain with :meth:`create`.

    :ivar type: The type identifier of the data object, always ``'MarkovChain'``.
    :ivar varied_parameters: The descriptions (name, min, max) of the varied parameters.
    :ivar lookup_table: A mapping from each parameter name to its column index in :attr:`samples`.
    :ivar samples: The samples in parameter space as a 2D array of shape (N, P).
    :ivar usamples: The samples in unit-hypercube ``u`` space as a 2D array of shape (N, P).
    :ivar weights: The importance weights on a linear scale as a 1D array of shape (N, ), or ``None`` if the chain is unweighted.
    """

    def __init__(self, path):
        """ Read a MarkovChain object from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError(f'Path {path} does not exist or is not a directory')

        description = MarkovChainDescription.from_yaml_file(os.path.join(path, 'description.yaml'))

        self.type = description.type
        self.varied_parameters = [asdict(p) for p in description.parameters]
        self.lookup_table = { p.name: idx for idx, p in enumerate(description.parameters) }

        ncols = len(description.parameters)
        self.samples  = load_array(path, 'samples.npy',  ncols=ncols)
        self.usamples = load_array(path, 'usamples.npy', ncols=ncols)

        if description.has_weights:
            self.weights = load_array(path, 'weights.npy', nrows=self.samples.shape[0])
        else:
            self.weights = None


    @staticmethod
    def create(path, parameters, samples, usamples, weights=None):
        """ Write a new MarkovChain object to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param parameters: Parameter descriptions as a 1D array of shape (N, ).
        :type parameters: list or iterable of eos.Parameter
        :param samples: Samples in parameter space as a 2D array of shape (N, P).
        :type samples: 2D numpy array
        :param usamples: Samples in u space as a 2D array of shape (N, P).
        :type usamples: 2D numpy array
        :param weights: Weights on a linear scale as a 1D array of shape (N, ).
        :type weights: 1D numpy array, optional
        """
        if not samples.shape[1] == len(parameters):
            raise RuntimeError(f'Shape of samples {samples.shape} incompatible with number of parameters {len(parameters)}')

        if not usamples.shape[1] == len(parameters):
            raise RuntimeError(f'Shape of usamples {usamples.shape} incompatible with number of parameters {len(parameters)}')

        if not weights is None and not samples.shape[0] == weights.shape[0]:
            raise RuntimeError(f'Shape of weights {weights.shape} incompatible with shape of samples {samples.shape}')

        description = MarkovChainDescription(
            version     = eos.__version__,
            parameters  = [ParameterDescription(name=p.name(), min=p.min(), max=p.max()) for p in parameters],
            has_weights = weights is not None,
        )

        os.makedirs(path, exist_ok=True)
        description.to_yaml_file(os.path.join(path, 'description.yaml'))
        _np.save(os.path.join(path, 'samples.npy'), samples)
        _np.save(os.path.join(path, 'usamples.npy'), usamples)

        if not weights is None:
            _np.save(os.path.join(path, 'weights.npy'), weights)
