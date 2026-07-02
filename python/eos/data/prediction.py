# Copyright (c) 2020-2026 Danny van Dyk
# Copyright (c) 2024-2025 Méril Reboud
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
from eos.data._common import load_array
from eos.deserializable import Deserializable
from eos.serializable import Serializable

import copy as _copy
import eos
import os
import numpy as _np


# The current on-disk format version written by :meth:`Prediction.create`. Files without an explicit
# 'format' key predate this versioning and are treated as format 1 (no recorded column 'kind').
_PREDICTION_FORMAT = 2


def _classify_column(name, options, kinematics, parameters, observables):
    r"""Classify a prediction column as an ``'observable'`` or a ``'parameter'``.

    Mirrors the precedence of :cpp:func:`eos.Observable.make`: a column is a parameter (clothed as an
    observable) only if it carries no options and no kinematics, its qualified name has an empty option
    part, its name does **not** match a registered observable (a registered observable always shadows a
    parameter of the same name), and its name **does** match a parameter.

    :param name: The qualified name of the column.
    :param options: The column's options as a mapping (empty for a parameter).
    :param kinematics: The column's kinematics as a mapping (empty for a parameter).
    :param parameters: The parameter set to test membership against.
    :type parameters: eos.Parameters
    :param observables: The observable registry to test the shadowing rule against.
    :type observables: eos.Observables
    :rtype: str
    """
    if options or kinematics:
        return 'observable'
    if list(eos.QualifiedName(name).options_part()):
        return 'observable'
    # a registered observable of the same name shadows any parameter
    try:
        observables[name]
        return 'observable'
    except RuntimeError:
        pass
    try:
        parameters[name]
        return 'parameter'
    except RuntimeError:
        return 'observable'


@dataclass(kw_only=True)
class ObservableDescription(Serializable, Deserializable):
    r"""Describes a single column of a :class:`Prediction`.

    A column is either a genuine observable (with its options and kinematics) or a parameter clothed as
    an observable (with empty options and kinematics), as recorded by :attr:`kind`.

    :param name: The qualified name of the observable (or parameter).
    :type name: str
    :param kind: Whether the column is an ``'observable'`` or a ``'parameter'``.
    :type kind: str
    :param kinematics: The kinematic variables and their values. Empty for a parameter.
    :type kinematics: dict
    :param options: The options and their values. Empty for a parameter.
    :type options: dict
    """
    name:str
    kind:str = field(default='observable')
    kinematics:dict = field(default_factory=dict)
    options:dict = field(default_factory=dict)

    def __post_init__(self):
        if self.kind not in ('observable', 'parameter'):
            raise ValueError(f'Invalid column kind \'{self.kind}\'; must be \'observable\' or \'parameter\'')


@dataclass(kw_only=True)
class PredictionDescription(Serializable, Deserializable):
    r"""Schema of the ``description.yaml`` written for a :class:`Prediction`.

    This is the single source of truth for the metadata stored alongside a set of posterior-predictive
    samples, and it backs both the read path (:meth:`from_yaml_file`) and the write path
    (:meth:`to_yaml_file`). The ``type`` discriminator is validated on deserialization and is not part
    of the constructor.

    The ``format`` key versions the on-disk layout. Files without it predate the versioning and are
    read as format 1, which does not record whether each column is an observable or a parameter; that
    classification is then re-derived on read. New files are always written as format
    :data:`_PREDICTION_FORMAT`, recording each column's :attr:`~ObservableDescription.kind`.

    :param version: The version of EOS that wrote the data object.
    :type version: str
    :param observables: The descriptions of the predicted columns.
    :type observables: list[ObservableDescription]
    """
    version:str
    observables:list[ObservableDescription]
    format:int = field(default=_PREDICTION_FORMAT)
    type:str = field(init=False, default='Prediction')

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`PredictionDescription` from its on-disk keyword description.

        Validates the ``type`` discriminator and deserializes each entry of ``observables`` into an
        :class:`ObservableDescription`.

        :raises ValueError: If the description does not identify a prediction.
        """
        _kwargs = _copy.deepcopy(kwargs)

        _type = _kwargs.pop('type', None)
        if _type != 'Prediction':
            raise ValueError(f'Expected a description of type \'Prediction\', got \'{_type}\'')

        _format = _kwargs.pop('format', 1)

        if 'observables' in _kwargs:
            _kwargs['observables'] = [ObservableDescription.from_dict(**o) for o in _kwargs['observables']]

        # Record the format the description was read from (provenance); freshly created
        # descriptions default to the current format.
        _kwargs['format'] = _format

        return Deserializable.make(cls, **_kwargs)

    def to_dict(self):
        """Serialize this description into the on-disk mapping written to ``description.yaml``.

        Always emits the current on-disk format (:data:`_PREDICTION_FORMAT`), independent of
        :attr:`format` (which records the provenance of a read description): the in-memory structure
        is always the current one, so any file we write is a current-format file.
        """
        return {
            'version':     self.version,
            'type':        self.type,
            'format':      _PREDICTION_FORMAT,
            'observables': [o.to_dict() for o in self.observables],
        }


class Prediction:
    r"""Represents weighted samples of theory predictions for one or more observables, stored on disk.

    Bundles the predictive samples of a set of observables with their importance weights. Instances are
    created either by reading an existing data set from disk (passing its ``path`` to the constructor) or
    by writing a new data set with :meth:`create`.

    Each column is either a genuine observable or a parameter clothed as an observable; the two are
    distinguished by the ``kind`` key of the corresponding entry in :attr:`varied_parameters`.

    :ivar type: The type identifier of the data object, always ``'Prediction'``.
    :ivar format: The on-disk format version the prediction was read from (1 for legacy files).
    :ivar varied_parameters: The descriptions (name, kind, kinematics, options) of the predicted columns.
    :ivar lookup_table: A mapping from each column's qualified name (including options and kinematics) to its column index in :attr:`samples`.
    :ivar samples: The predictive samples as a 2D array of shape (N, O).
    :ivar weights: The importance weights on a linear scale as a 1D array of shape (N, ).
    """

    def __init__(self, path):
        """ Read a Prediction object from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError(f'Path {path} does not exist or is not a directory')

        description = PredictionDescription.from_yaml_file(os.path.join(path, 'description.yaml'))

        self.type = description.type
        self.format = description.format
        observables = description.observables

        # For legacy files (no recorded kind), re-derive each column's kind, mirroring
        # eos.Observable.make (best-effort, against the default registries).
        if description.format < 2:
            registry   = eos.Observables()
            parameters = eos.Parameters()
            for o in observables:
                o.kind = _classify_column(o.name, o.options, o.kinematics, parameters, registry)

        self.varied_parameters = [asdict(o) for o in observables]

        self.lookup_table = {}
        for idx, o in enumerate(observables):
            id = o.name
            id += ';' + str(eos.Options(o.options)).replace(" ", "")
            id += '[' + str(eos.Kinematics(o.kinematics)).replace(" ", "") + ']'
            self.lookup_table[id] = idx

        ncols = len(observables)
        self.samples = load_array(path, 'samples.npy', ncols=ncols)
        self.weights = load_array(path, 'weights.npy', nrows=self.samples.shape[0])


    @staticmethod
    def create(path, observables, samples, weights):
        """ Write a new Prediction object to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param observables: Observables as a 1D array of shape (O, ).
        :type observables: list or iterable of eos.Observable
        :param samples: Samples as a 2D array of shape (N, O).
        :type samples: 2D numpy array
        :param weights: Weights on a linear scale as a 1D array of shape (N, ).
        :type weights: 1D numpy array
        """
        if not samples.shape[1] == len(observables):
            raise RuntimeError(f'Shape of samples {samples.shape} incompatible with number of observables {len(observables)}')

        if not samples.shape[0] == weights.shape[0]:
            raise RuntimeError(f'Shape of weights {weights.shape} incompatible with shape of samples {samples.shape}')

        registry = eos.Observables()
        observable_descriptions = []
        for o in observables:
            name       = o.name().full()
            kinematics = { k.name(): float(k) for k in o.kinematics() }
            options    = { str(k): str(v) for k, v in o.options() }
            kind       = _classify_column(name, options, kinematics, o.parameters(), registry)
            observable_descriptions.append(ObservableDescription(
                name=name, kind=kind, kinematics=kinematics, options=options))

        description = PredictionDescription(version=eos.__version__, observables=observable_descriptions)

        os.makedirs(path, exist_ok=True)
        description.to_yaml_file(os.path.join(path, 'description.yaml'))
        _np.save(os.path.join(path, 'samples.npy'), samples)
        _np.save(os.path.join(path, 'weights.npy'), weights)
