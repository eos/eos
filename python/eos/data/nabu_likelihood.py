# Copyright (c) 2024-2026 Danny van Dyk
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
from eos.data._common import ParameterDescription
from eos.deserializable import Deserializable
from eos.serializable import Serializable

import copy as _copy
import eos
import os


@dataclass(kw_only=True)
class NabuLikelihoodDescription(Serializable, Deserializable):
    r"""Schema of the ``description.yaml`` written for a :class:`NabuLikelihood` object.

    This is the single source of truth for the metadata stored alongside a nabu-serialized likelihood,
    and it backs both the read path (:meth:`from_yaml_file`) and the write path (:meth:`to_yaml_file`).
    The ``type`` discriminator is validated on deserialization and is not part of the constructor.

    :param version: The version of EOS that wrote the data object.
    :type version: str
    :param parameters: The descriptions of the varied parameters.
    :type parameters: list[ParameterDescription]
    """
    version:str
    parameters:list[ParameterDescription]
    type:str = field(init=False, default='NabuLikelihood')

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`NabuLikelihoodDescription` from its on-disk keyword description.

        Validates the ``type`` discriminator and deserializes each entry of ``parameters`` into a
        :class:`ParameterDescription`.

        :raises ValueError: If the description does not identify a nabu-likelihood object.
        """
        _kwargs = _copy.deepcopy(kwargs)

        _type = _kwargs.pop('type', None)
        if _type != 'NabuLikelihood':
            raise ValueError(f'Expected a description of type \'NabuLikelihood\', got \'{_type}\'')

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


class NabuLikelihood:
    r"""Represents a likelihood serialized with the nabu package, stored on disk.

    Wraps a :class:`nabu.Likelihood` object together with the descriptions of the varied parameters.
    Instances are created either by reading an existing likelihood from disk (passing its ``path`` to
    the constructor) or by writing a new likelihood with :meth:`create`. Reading requires the optional
    nabu module.

    :ivar type: The type identifier of the data object, always ``'NabuLikelihood'``.
    :ivar varied_parameters: The descriptions (name, min, max) of the varied parameters.
    :ivar lookup_table: A mapping from each parameter name to its index in :attr:`varied_parameters`.
    :ivar likelihood: The underlying :class:`nabu.Likelihood` object.
    """

    def __init__(self, path):
        """ Read a nabu serialized likelihood from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError(f'Path {path} does not exist or is not a directory')

        description = NabuLikelihoodDescription.from_yaml_file(os.path.join(path, 'description.yaml'))

        self.type = description.type
        self.varied_parameters = [asdict(p) for p in description.parameters]
        self.lookup_table = { p.name: idx for idx, p in enumerate(description.parameters) }

        f = os.path.join(path, 'likelihood.nabu')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Nabu likelihood file {f} does not exist or is not a file')

        import nabu
        self.likelihood = nabu.Likelihood.load(f)


    @staticmethod
    def create(path, parameters, likelihood):
        """ Write a new NabuLikelihood object to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param parameters: Parameter descriptions as a 1D array of shape (N, ).
        :type parameters: list or iterable of eos.Parameter
        :param likelihood: The likelihood object created by nabu.
        :type likelihood: nabu.Likelihood or descendant
        """
        description = NabuLikelihoodDescription(
            version    = eos.__version__,
            parameters = [ParameterDescription(name=p.name(), min=p.min(), max=p.max()) for p in parameters],
        )

        os.makedirs(path, exist_ok=True)
        description.to_yaml_file(os.path.join(path, 'description.yaml'))

        likelihood.save(os.path.join(path, 'likelihood.nabu'))
