# Copyright (c) 2022-2026 Danny van Dyk
# Copyright (c) 2021-2024 Méril Reboud
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
class ImportanceSamplesDescription(Serializable, Deserializable):
    r"""Schema of the ``description.yaml`` written for an :class:`ImportanceSamples` object.

    This is the single source of truth for the metadata stored alongside a set of importance samples,
    and it backs both the read path (:meth:`from_yaml_file`) and the write path (:meth:`to_yaml_file`).
    The ``type`` discriminator is validated on deserialization and is not part of the constructor.

    :param version: The version of EOS that wrote the data object.
    :type version: str
    :param parameters: The descriptions of the varied parameters.
    :type parameters: list[ParameterDescription]
    """
    version:str
    parameters:list[ParameterDescription]
    type:str = field(init=False, default='ImportanceSamples')

    @classmethod
    def from_dict(cls, **kwargs):
        """Create an :class:`ImportanceSamplesDescription` from its on-disk keyword description.

        Validates the ``type`` discriminator and deserializes each entry of ``parameters`` into a
        :class:`ParameterDescription`.

        :raises ValueError: If the description does not identify an importance-samples object.
        """
        _kwargs = _copy.deepcopy(kwargs)

        _type = _kwargs.pop('type', None)
        if _type != 'ImportanceSamples':
            raise ValueError(f'Expected a description of type \'ImportanceSamples\', got \'{_type}\'')

        if 'parameters' in _kwargs:
            _kwargs['parameters'] = [ParameterDescription.from_dict(**p) for p in _kwargs['parameters']]

        return Deserializable.make(cls, **_kwargs)

    def to_dict(self):
        """Serialize this description into the on-disk mapping written to ``description.yaml``.

        Emits the ``type`` discriminator, inverting :meth:`from_dict`.
        """
        return {
            'version':    self.version,
            'type':       self.type,
            'parameters': [asdict(p) for p in self.parameters],
        }


class ImportanceSamples:
    r"""Represents a set of weighted importance samples stored on disk.

    Bundles samples in parameter space with their importance weights and, optionally, the values of the
    posterior density at each sample. Instances are created either by reading an existing data set from
    disk (passing its ``path`` to the constructor) or by writing a new data set with :meth:`create`.

    :ivar type: The type identifier of the data object, always ``'ImportanceSamples'``.
    :ivar varied_parameters: The descriptions (name, min, max) of the varied parameters.
    :ivar lookup_table: A mapping from each parameter name to its column index in :attr:`samples`.
    :ivar samples: The samples in parameter space as a 2D array of shape (N, P).
    :ivar weights: The importance weights on a linear scale as a 1D array of shape (N, ).
    :ivar posterior_values: The posterior density values at each sample as a 1D array of shape (N, ), or ``None`` if not stored.
    """

    def __init__(self, path):
        """ Read an ImportanceSamples object from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError(f'Path {path} does not exist or is not a directory')

        description = ImportanceSamplesDescription.from_yaml_file(os.path.join(path, 'description.yaml'))

        self.type = description.type
        self.varied_parameters = [asdict(p) for p in description.parameters]
        self.lookup_table = { p.name: idx for idx, p in enumerate(description.parameters) }

        ncols = len(description.parameters)
        self.samples = load_array(path, 'samples.npy', ncols=ncols)
        self.weights = load_array(path, 'weights.npy', nrows=self.samples.shape[0])

        # posterior_values are optional and detected by the presence of their file
        if os.path.isfile(os.path.join(path, 'posterior_values.npy')):
            self.posterior_values = load_array(path, 'posterior_values.npy', nrows=self.samples.shape[0])
        else:
            self.posterior_values = None


    @staticmethod
    def create(path, parameters, samples, weights, posterior_values=None):
        """ Write a new ImportanceSamples object to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param parameters: Parameter descriptions as a 1D array of shape (P, ).
        :type parameters: list or iterable of eos.Parameter
        :param samples: Samples as a 2D array of shape (N, P).
        :type samples: 2D numpy array
        :param weights: Weights on a linear scale as a 1D array of shape (N, ).
        :type weights: 1D numpy array, optional
        """
        if not samples.shape[1] == len(parameters):
            raise RuntimeError(f'Shape of samples {samples.shape} incompatible with number of parameters {len(parameters)}')

        if not weights is None and not samples.shape[0] == weights.shape[0]:
            raise RuntimeError(f'Shape of weights {weights.shape} incompatible with shape of samples {samples.shape}')

        if not posterior_values is None and not samples.shape[0] == posterior_values.shape[0]:
            raise RuntimeError(f'Shape of posterior values {posterior_values.shape} incompatible with shape of samples {samples.shape}')

        description = ImportanceSamplesDescription(
            version    = eos.__version__,
            parameters = [ParameterDescription(
                name = p.name(),
                min  = p.min() if 'min' in dir(p) else -_np.inf,
                max  = p.max() if 'max' in dir(p) else +_np.inf,
            ) for p in parameters],
        )

        os.makedirs(path, exist_ok=True)
        description.to_yaml_file(os.path.join(path, 'description.yaml'))
        _np.save(os.path.join(path, 'samples.npy'), samples)
        _np.save(os.path.join(path, 'weights.npy'), weights)
        if not posterior_values is None:
            _np.save(os.path.join(path, 'posterior_values.npy'), posterior_values)
