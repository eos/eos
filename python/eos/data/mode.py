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
import numpy as _np


# The current on-disk format version written by :meth:`Mode.create`. Files without an explicit
# 'format_version' key predate this versioning and are treated as format v1.
_MODE_FORMAT_VERSION = 2


@dataclass(kw_only=True)
class TestStatisticDescription(Serializable, Deserializable):
    r"""Describes the test statistic of a single element of the likelihood evaluated at a mode.

    :param name: The name of the likelihood element (e.g. the constraint name).
    :type name: str
    :param local_pvalue: The local p-value of the element at the mode.
    :type local_pvalue: float
    :param type: The type of the primary test statistic, typically ``'chi^2'``.
    :type type: str
    :param value: The value of the primary test statistic, or ``None`` if not recorded (e.g. read
        from a format-v1 file, which stored only the local p-value).
    :type value: float | None
    """
    name:str
    local_pvalue:float
    type:str = field(default='chi^2')
    value:float|None = field(default=None)


@dataclass(kw_only=True)
class ModeDescription(Serializable, Deserializable):
    r"""Schema of the ``description.yaml`` written for a :class:`Mode`.

    This is the single source of truth for the metadata stored for a (local) mode, and it backs both
    the read path (:meth:`from_yaml_file`) and the write path (:meth:`to_yaml_file`). The ``type``
    discriminator is validated on deserialization and is not part of the constructor.

    The ``format_version`` key versions the on-disk layout. Files without it predate the versioning and are
    read as format v1: their flat ``local_pvalues`` map is upgraded into :attr:`test_statistics`
    entries (with ``type='chi^2'`` and no recorded statistic value). New files are always written as
    format :data:`_MODE_FORMAT_VERSION`. The ``global_chi2`` and ``dof`` keys remain optional.

    :param version: The version of EOS that wrote the data object.
    :type version: str
    :param parameters: The descriptions of the varied parameters.
    :type parameters: list[ParameterDescription]
    :param mode: The location of the mode in parameter space.
    :type mode: list[float]
    :param pvalue: The global p-value at the mode, or ``None``.
    :type pvalue: float | None
    :param test_statistics: The per-element test statistics at the mode.
    :type test_statistics: list[TestStatisticDescription]
    :param global_chi2: The global :math:`\chi^2` value at the mode, or ``None`` if not stored.
    :type global_chi2: float | None
    :param dof: The number of degrees of freedom, or ``None`` if not stored.
    :type dof: float | None
    """
    version:str
    parameters:list[ParameterDescription]
    mode:list[float]
    pvalue:float|None
    test_statistics:list[TestStatisticDescription]
    global_chi2:float|None = field(default=None)
    dof:float|None = field(default=None)
    format_version:int = field(default=_MODE_FORMAT_VERSION)
    type:str = field(init=False, default='Mode')

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`ModeDescription` from its on-disk keyword description.

        Validates the ``type`` discriminator, deserializes each entry of ``parameters``, and reads the
        per-element test statistics. Format-v1 files (those without a ``format_version`` key) are upgraded: the
        flat ``local_pvalues`` map is converted into :class:`TestStatisticDescription` entries.

        :raises ValueError: If the description does not identify a mode, or a format-v1 description is
            missing its ``local_pvalues`` map.
        """
        _kwargs = _copy.deepcopy(kwargs)

        _type = _kwargs.pop('type', None)
        if _type != 'Mode':
            raise ValueError(f'Expected a description of type \'Mode\', got \'{_type}\'')

        _format = _kwargs.pop('format_version', 1)

        if 'parameters' in _kwargs:
            _kwargs['parameters'] = [ParameterDescription.from_dict(**p) for p in _kwargs['parameters']]

        if _format >= 2:
            if 'test_statistics' in _kwargs:
                _kwargs['test_statistics'] = [TestStatisticDescription.from_dict(**t) for t in _kwargs['test_statistics']]
        else:
            # Upgrade a format-v1 description: synthesize test statistics from the flat local_pvalues map.
            local_pvalues = _kwargs.pop('local_pvalues', None)
            if local_pvalues is None:
                raise ValueError('A format-v1 Mode description must contain \'local_pvalues\'')
            _kwargs['test_statistics'] = [
                TestStatisticDescription.from_dict(name=name, local_pvalue=local_pvalue)
                for name, local_pvalue in local_pvalues.items()
            ]

        # Record the format the description was read from (provenance); freshly created
        # descriptions default to the current format.
        _kwargs['format_version'] = _format

        return Deserializable.make(cls, **_kwargs)

    def to_dict(self):
        """Serialize this description into the on-disk mapping written to ``description.yaml``.

        Always emits the current on-disk format (:data:`_MODE_FORMAT_VERSION`), independent of
        :attr:`format_version` (which records the provenance of a read description): the in-memory structure
        is always the current one, so any file we write is a current-format file.
        """
        return {
            'version':         self.version,
            'type':            self.type,
            'format_version':  _MODE_FORMAT_VERSION,
            'parameters':      [asdict(p) for p in self.parameters],
            'mode':            self.mode,
            'pvalue':          self.pvalue,
            'global_chi2':     self.global_chi2,
            'dof':             self.dof,
            'test_statistics': [asdict(t) for t in self.test_statistics],
        }


class Mode:
    r"""Represents a (local) mode of a posterior density stored on disk.

    Stores the location of the mode in parameter space together with goodness-of-fit information.
    Instances are created either by reading an existing mode from disk (passing its ``path`` to the
    constructor) or by writing a new mode with :meth:`create` (low-level) or :meth:`from_bfp_and_gof`
    (from an optimization result).

    :ivar type: The type identifier of the data object, always ``'Mode'``.
    :ivar format_version: The on-disk format version the mode was read from (1 for legacy files).
    :ivar varied_parameters: The descriptions (name, min, max) of the varied parameters.
    :ivar mode: The location of the mode in parameter space.
    :ivar pvalue: The global p-value at the mode.
    :ivar test_statistics: The per-element test statistics at the mode (name, type, value, local p-value).
    :ivar local_pvalues: A mapping from each element's name to its local p-value, derived from :attr:`test_statistics`.
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
        self.format_version = description.format_version
        self.varied_parameters = [asdict(p) for p in description.parameters]
        self.mode = description.mode
        self.pvalue = description.pvalue
        self.global_chi2 = description.global_chi2
        self.dof = description.dof
        self.test_statistics = [asdict(t) for t in description.test_statistics]
        self.local_pvalues = { t['name']: t['local_pvalue'] for t in self.test_statistics }

    @staticmethod
    def create(path, parameters, mode, pvalue, test_statistics, global_chi2, dof):
        """ Write a new Mode object to disk.

        This is the low-level serializer. To record a mode obtained from an optimization, prefer
        :meth:`from_bfp_and_gof`, which derives the arguments below from the best-fit point and the
        goodness-of-fit object.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param parameters: Parameter descriptions as a 1D array of shape (N, ).
        :type parameters: list or iterable of eos.Parameter
        :param mode: The mode to be stored.
        :type mode: numpy.ndarray
        :param pvalue: The global p-value of the mode.
        :type pvalue: float
        :param test_statistics: The per-element test statistics, each a mapping with the keys
            ``name``, ``local_pvalue`` and, optionally, ``type`` (default ``'chi^2'``) and ``value``.
        :type test_statistics: list of dict
        :param global_chi2: The global chi2 value of the mode.
        :type global_chi2: float
        :param dof: The degrees of freedom of the mode.
        :type dof: float
        """
        description = ModeDescription(
            version         = eos.__version__,
            parameters      = [ParameterDescription(name=p.name(), min=p.min(), max=p.max()) for p in parameters],
            mode            = _np.asarray(mode).tolist(),
            pvalue          = float(pvalue) if pvalue is not None else None,
            test_statistics = [TestStatisticDescription.from_dict(**t) for t in test_statistics],
            global_chi2     = global_chi2,
            dof             = dof,
        )

        os.makedirs(path, exist_ok=True)
        description.to_yaml_file(os.path.join(path, 'description.yaml'))

    @staticmethod
    def from_bfp_and_gof(path, bfp, gof):
        """ Write a new Mode object to disk from an optimization result.

        Derives the varied parameters and the mode from the best-fit point ``bfp`` and the global and
        per-element test statistics (chi^2 values, degrees of freedom, and p-values) from the
        goodness-of-fit object ``gof``, then delegates to :meth:`create`.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param bfp: The best-fit point, as returned by :meth:`eos.Analysis.optimize`.
        :type bfp: eos.BestFitPoint
        :param gof: The goodness-of-fit at the best-fit point.
        :type gof: eos.GoodnessOfFit
        """
        import scipy.stats

        total_chi2 = gof.total_chi_square()
        dof        = gof.total_degrees_of_freedom()
        pvalue     = float(1.0 - scipy.stats.chi2(dof).cdf(total_chi2))

        test_statistics = []
        for name, entry in gof:
            local_pvalue = float(1.0 - scipy.stats.chi2(entry.dof).cdf(entry.chi2))
            test_statistics.append({
                'name':         str(name),
                'type':         'chi^2',
                'value':        float(entry.chi2),
                'local_pvalue': local_pvalue,
            })

        Mode.create(path, bfp.analysis.varied_parameters, bfp.point, pvalue, test_statistics, total_chi2, dof)
