# Copyright (c) 2022-2024 Méril Reboud
# Copyright (c) 2023-2026 Danny van Dyk
# Copyright (c) 2023 Stephan Kürten
# Copyright (c) 2024 Lorenz Gärtner
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
class ModeDescription(Serializable, Deserializable):
    r"""Schema of the ``description.yaml`` written for a :class:`Mode`.

    This is the single source of truth for the metadata stored for a (local) mode, and it backs both
    the read path (:meth:`from_yaml_file`) and the write path (:meth:`to_yaml_file`). The ``type``
    discriminator is validated on deserialization and is not part of the constructor. The
    ``global_chi2`` and ``dof`` keys are optional for backward compatibility with older files.

    :param version: The version of EOS that wrote the data object.
    :type version: str
    :param parameters: The descriptions of the varied parameters.
    :type parameters: list[ParameterDescription]
    :param mode: The location of the mode in parameter space.
    :type mode: list[float]
    :param pvalue: The global p-value at the mode, or ``None``.
    :type pvalue: float | None
    :param local_pvalues: The local p-values at the mode.
    :type local_pvalues: dict
    :param global_chi2: The global :math:`\chi^2` value at the mode, or ``None`` if not stored.
    :type global_chi2: float | None
    :param dof: The number of degrees of freedom, or ``None`` if not stored.
    :type dof: float | None
    """
    version:str
    parameters:list[ParameterDescription]
    mode:list[float]
    pvalue:float|None
    local_pvalues:dict
    global_chi2:float|None = field(default=None)
    dof:float|None = field(default=None)
    type:str = field(init=False, default='Mode')

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`ModeDescription` from its on-disk keyword description.

        Validates the ``type`` discriminator and deserializes each entry of ``parameters`` into a
        :class:`ParameterDescription`.

        :raises ValueError: If the description does not identify a mode.
        """
        _kwargs = _copy.deepcopy(kwargs)

        _type = _kwargs.pop('type', None)
        if _type != 'Mode':
            raise ValueError(f'Expected a description of type \'Mode\', got \'{_type}\'')

        if 'parameters' in _kwargs:
            _kwargs['parameters'] = [ParameterDescription.from_dict(**p) for p in _kwargs['parameters']]

        return Deserializable.make(cls, **_kwargs)

    def to_dict(self):
        """Serialize this description into the on-disk mapping written to ``description.yaml``.

        Emits the ``type`` discriminator, inverting :meth:`from_dict`.
        """
        return {
            'version':       self.version,
            'type':          self.type,
            'parameters':    [asdict(p) for p in self.parameters],
            'mode':          self.mode,
            'pvalue':        self.pvalue,
            'local_pvalues': self.local_pvalues,
            'global_chi2':   self.global_chi2,
            'dof':           self.dof,
        }


class Mode:
    r"""Represents a (local) mode of a posterior density stored on disk.

    Stores the location of the mode in parameter space together with goodness-of-fit information.
    Instances are created either by reading an existing mode from disk (passing its ``path`` to the
    constructor) or by writing a new mode with :meth:`create`.

    :ivar type: The type identifier of the data object, always ``'Mode'``.
    :ivar varied_parameters: The descriptions (name, min, max) of the varied parameters.
    :ivar mode: The location of the mode in parameter space.
    :ivar pvalue: The global p-value at the mode.
    :ivar local_pvalues: The local p-values at the mode.
    :ivar global_chi2: The global :math:`\chi^2` value at the mode, or ``None`` if not stored.
    :ivar dof: The number of degrees of freedom, or ``None`` if not stored.
    """

    def __init__(self, path):
        """ Read a posterior's (local) mode from a file.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError(f'Path {path} does not exist or is not a directory')

        description = ModeDescription.from_yaml_file(os.path.join(path, 'description.yaml'))

        self.type = description.type
        self.varied_parameters = [asdict(p) for p in description.parameters]
        self.mode = description.mode
        self.pvalue = description.pvalue
        self.local_pvalues = description.local_pvalues
        self.global_chi2 = description.global_chi2
        self.dof = description.dof

    @staticmethod
    def create(path, parameters, mode, pvalue, local_pvalues, global_chi2, dof):
        """ Write a new Mode object to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param parameters: Parameter descriptions as a 1D array of shape (N, ).
        :type parameters: list or iterable of eos.Parameter
        :param mode: The mode to be stored.
        :type mode: numpy.ndarray
        :param pvalue: The p-value of the mode.
        :type pvalue: float
        :param local_pvalues: The local p-values of the mode.
        :type local_pvalues: dict
        :param global_chi2: The global chi2 value of the mode.
        :type global_chi2: float
        :param dof: The degrees of freedom of the mode.
        :type dof: float
        """
        description = ModeDescription(
            version       = eos.__version__,
            parameters    = [ParameterDescription(name=p.name(), min=p.min(), max=p.max()) for p in parameters],
            mode          = mode.tolist(),
            pvalue        = float(pvalue) if pvalue is not None else None,
            local_pvalues = local_pvalues,
            global_chi2   = global_chi2,
            dof           = dof,
        )

        os.makedirs(path, exist_ok=True)
        description.to_yaml_file(os.path.join(path, 'description.yaml'))
