# Copyright (c) 2025 Matthew Kirk
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

from dataclasses import dataclass, field
from eos.data._common import load_array
from eos.deserializable import Deserializable
from eos.serializable import Serializable

import copy as _copy
import eos
import os
import numpy as _np


@dataclass(kw_only=True)
class SampleMaskDescription(Serializable, Deserializable):
    r"""Schema of the ``description.yaml`` written for a :class:`SampleMask`.

    This is the single source of truth for the metadata stored alongside a sample mask, and it backs
    both the read path (:meth:`from_yaml_file`) and the write path (:meth:`to_yaml_file`).

    The ``type`` discriminator is always written as ``'SampleMask'``. For backward compatibility the
    legacy discriminator ``'Mask'`` is also accepted on deserialization (with a warning), since older
    files were written with it.

    :param version: The version of EOS that wrote the data object.
    :type version: str
    :param observables: The names of the observables to which the mask applies.
    :type observables: list
    """
    version:str
    observables:list
    type:str = field(init=False, default='SampleMask')

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`SampleMaskDescription` from its on-disk keyword description.

        Accepts either the current ``'SampleMask'`` discriminator or the legacy ``'Mask'`` one; the
        latter is deprecated and triggers a warning.

        :raises ValueError: If the description does not identify a sample mask.
        """
        _kwargs = _copy.deepcopy(kwargs)

        _type = _kwargs.pop('type', None)
        if _type == 'Mask':
            eos.warn('The \'Mask\' data type is deprecated; it will be written as \'SampleMask\' when this data set is re-saved')
        elif _type != 'SampleMask':
            raise ValueError(f'Expected a description of type \'SampleMask\' (or the legacy \'Mask\'), got \'{_type}\'')

        return Deserializable.make(cls, **_kwargs)

    def to_dict(self):
        """Serialize this description into the on-disk mapping written to ``description.yaml``.

        Always emits the current ``'SampleMask'`` discriminator, inverting :meth:`from_dict`.
        """
        return {
            'version':     self.version,
            'type':        self.type,
            'observables': self.observables,
        }


class SampleMask:
    r"""Represents a boolean mask over a set of observables, stored on disk.

    A mask selects a subset of samples (e.g. to remove unphysical points) and records the observables it
    applies to. Instances are created either by reading an existing mask from disk (passing its ``path``
    to the constructor) or by writing a new mask with :meth:`create`.

    :ivar type: The type identifier of the data object, always ``'SampleMask'``.
    :ivar observables: The observables to which the mask applies.
    :ivar mask: The boolean mask as a 1D array.
    """

    def __init__(self, path):
        """ Read a mask from a file.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError(f'Path {path} does not exist or is not a directory')

        description = SampleMaskDescription.from_yaml_file(os.path.join(path, 'description.yaml'))

        self.type = description.type
        self.observables = description.observables
        self.mask = load_array(path, 'mask.npy')

    @staticmethod
    def create(path, mask, observables):
        """ Write a new SampleMask object to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param mask: Mask as a 1D boolean array
        :type mask: list or iterable of bool
        :param observables: Observables as a 1D array
        """
        description = SampleMaskDescription(version=eos.__version__, observables=observables)

        os.makedirs(path, exist_ok=True)
        description.to_yaml_file(os.path.join(path, 'description.yaml'))
        _np.save(os.path.join(path, 'mask.npy'), mask)
