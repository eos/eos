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

from dataclasses import dataclass, field, asdict
from eos.data._common import load_array
from eos.deserializable import Deserializable
from eos.serializable import Serializable

import copy as _copy
import eos
import hashlib as _hashlib
import os
import numpy as _np


@dataclass(kw_only=True)
class UnbinnedAxisDescription(Serializable, Deserializable):
    r"""Describes one native sampling axis of an :class:`UnbinnedLikelihood` data object.

    :param variable: The kinematic variable sampled along this axis.
    :type variable: str
    :param min: The lower bound of the native grid.
    :type min: float
    :param max: The upper bound of the native grid.
    :type max: float
    :param points: The number of native grid points (even, >= 2).
    :type points: int
    """
    variable:str
    min:float
    max:float
    points:int


@dataclass(kw_only=True)
class UnbinnedOffsetAxisDescription(Serializable, Deserializable):
    r"""Describes one axis of a sampled resolution kernel's offset grid.

    :param variable: The offset variable sampled along this axis.
    :type variable: str
    :param min: The lower bound of the declared offset grid (the most negative offset sampled).
    :type min: float
    :param spacing: The node spacing of the declared offset grid.
    :type spacing: float
    :param points: The number of declared offset-grid nodes.
    :type points: int
    """
    variable:str
    min:float
    spacing:float
    points:int


@dataclass(kw_only=True)
class UnbinnedSampledResolutionDescription(Serializable, Deserializable):
    r"""Describes a resolution kernel given as numerical samples (``kind: samples``).

    The samples are stored in ``resolution.npy``, in centred order relative to the declared
    :attr:`offsets` grid.

    :param offsets: The axis descriptions of the declared offset grid, in axis order.
    :type offsets: list[UnbinnedOffsetAxisDescription]
    """
    offsets:list
    kind:str = field(init=False, default='samples')

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs.pop('kind', None)
        if 'offsets' in _kwargs:
            _kwargs['offsets'] = [UnbinnedOffsetAxisDescription.from_dict(**o) for o in _kwargs['offsets']]
        return Deserializable.make(cls, **_kwargs)

    def to_dict(self):
        return {
            'kind':    self.kind,
            'offsets': [asdict(o) for o in self.offsets],
        }


@dataclass(kw_only=True)
class UnbinnedExpressionResolutionDescription(Serializable, Deserializable):
    r"""Describes a resolution kernel given as an EOS expression (``kind: expression``).

    :param expression: The EOS expression that defines the (unnormalized) resolution kernel over the
        per-axis offset variables.
    :type expression: str
    """
    expression:str
    kind:str = field(init=False, default='expression')

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs.pop('kind', None)
        return Deserializable.make(cls, **_kwargs)

    def to_dict(self):
        return {
            'kind':       self.kind,
            'expression': self.expression,
        }


@dataclass(kw_only=True)
class UnbinnedPDFResolutionDescription(Serializable, Deserializable):
    r"""Describes a resolution kernel given as a registered SignalPDF (``kind: pdf``).

    :param pdf: The qualified name of the resolution SignalPDF, a density over the per-axis offset
        variables.
    :type pdf: str
    """
    pdf:str
    kind:str = field(init=False, default='pdf')

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs.pop('kind', None)
        return Deserializable.make(cls, **_kwargs)

    def to_dict(self):
        return {
            'kind': self.kind,
            'pdf':  self.pdf,
        }


class UnbinnedResolutionDescription:
    """Polymorphic description of an :class:`UnbinnedLikelihood`'s resolution.

    This is a dispatcher rather than a concrete description: :meth:`from_dict` selects the concrete
    type from the ``kind`` key: ``samples`` (:class:`UnbinnedSampledResolutionDescription`),
    ``expression`` (:class:`UnbinnedExpressionResolutionDescription`), or ``pdf``
    (:class:`UnbinnedPDFResolutionDescription`).
    """

    @staticmethod
    def from_dict(**kwargs):
        """Create the concrete resolution description matching the given keyword description.

        :returns: An instance of the concrete :class:`UnbinnedResolutionDescription` subtype.
        :raises ValueError: If ``kind`` is missing or unrecognized.
        """
        kind = kwargs.get('kind')
        if kind == 'samples':
            return UnbinnedSampledResolutionDescription.from_dict(**kwargs)
        elif kind == 'expression':
            return UnbinnedExpressionResolutionDescription.from_dict(**kwargs)
        elif kind == 'pdf':
            return UnbinnedPDFResolutionDescription.from_dict(**kwargs)

        raise ValueError(f"Unknown resolution kind '{kind}'; expected one of 'samples', 'expression', 'pdf'")


@dataclass(kw_only=True)
class UnbinnedLikelihoodDescription(Serializable, Deserializable):
    r"""Schema of the ``description.yaml`` written for an :class:`UnbinnedLikelihood` object.

    This is the single source of truth for the metadata stored alongside an unbinned-likelihood data
    object, and it backs both the read path (:meth:`from_yaml_file`) and the write path
    (:meth:`to_yaml_file`). The ``type`` discriminator is validated on deserialization and is not part
    of the constructor. The rank of the data object is ``len(axes)``, detected rather than declared.

    :param version: The version of EOS that wrote the data object.
    :type version: str
    :param axes: The native sampling axes, in axis order.
    :type axes: list[UnbinnedAxisDescription]
    :param resolution: The resolution description.
    :type resolution: UnbinnedSampledResolutionDescription | UnbinnedExpressionResolutionDescription | UnbinnedPDFResolutionDescription
    """
    version:str
    axes:list
    resolution:object
    type:str = field(init=False, default='UnbinnedLikelihood')

    @classmethod
    def from_dict(cls, **kwargs):
        """Create an :class:`UnbinnedLikelihoodDescription` from its on-disk keyword description.

        Validates the ``type`` discriminator and deserializes ``axes`` and ``resolution`` into their
        respective description types.

        :raises ValueError: If the description does not identify an unbinned-likelihood object.
        """
        _kwargs = _copy.deepcopy(kwargs)

        _type = _kwargs.pop('type', None)
        if _type != 'UnbinnedLikelihood':
            raise ValueError(f'Expected a description of type \'UnbinnedLikelihood\', got \'{_type}\'')

        if 'axes' in _kwargs:
            _kwargs['axes'] = [UnbinnedAxisDescription.from_dict(**a) for a in _kwargs['axes']]
        if 'resolution' in _kwargs:
            _kwargs['resolution'] = UnbinnedResolutionDescription.from_dict(**_kwargs['resolution'])

        return Deserializable.make(cls, **_kwargs)

    def to_dict(self):
        """Serialize this description into the on-disk mapping written to ``description.yaml``.

        Emits the ``type`` discriminator, inverting :meth:`from_dict`.
        """
        return {
            'version':    self.version,
            'type':       self.type,
            'axes':       [asdict(a) for a in self.axes],
            'resolution': self.resolution.to_dict(),
        }


class UnbinnedLikelihood:
    r"""Represents observed unbinned events and their detector resolution, stored on disk.

    Publishes the native grid granularity that the resolution supports (the axes' ``points``), and
    the resolution itself, in whatever form it was recorded (numerical samples, an EOS expression, or
    the qualified name of a registered SignalPDF). The signal PDF fitted to the events, and its
    options, are not part of this object: it describes the measured events, not the model fitted to
    them. Instances are created either by reading an existing data object from disk (passing its
    ``path`` to the constructor) or by writing a new one with :meth:`create`.

    :ivar type: The type identifier of the data object, always ``'UnbinnedLikelihood'``.
    :ivar axes: The native sampling axes (:class:`UnbinnedAxisDescription`), in axis order.
    :ivar resolution: The resolution description.
    :ivar observations: The observed events, shape ``(n_events, rank)``, columns in axis order.
    :ivar resolution_values: The sampled resolution kernel, shape equal to the native grid, or
        ``None`` unless :attr:`resolution` is a sampled (``kind: samples``) resolution.
    """

    def __init__(self, path):
        """Read an unbinned-likelihood data object from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError(f'Path {path} does not exist or is not a directory')

        description = UnbinnedLikelihoodDescription.from_yaml_file(os.path.join(path, 'description.yaml'))

        self.type = description.type
        self.axes = description.axes
        self.resolution = description.resolution

        rank = len(self.axes)
        self.observations = load_array(path, 'observations.npy', ncols=rank)

        if isinstance(self.resolution, UnbinnedSampledResolutionDescription):
            expected_shape = tuple(axis.points for axis in self.resolution.offsets)
            values = load_array(path, 'resolution.npy')
            if tuple(values.shape) != expected_shape:
                raise RuntimeError(
                    f'resolution.npy has shape {values.shape}, expected {expected_shape} from the declared offset grid'
                )
            self.resolution_values = values
        else:
            self.resolution_values = None

    @property
    def rank(self):
        """The dimensionality of the data object, detected as ``len(axes)``."""
        return len(self.axes)

    def cropped_axes(self, crop=None):
        r"""Compute the per-axis grid used by a :class:`DetectorLevelPDF`, applying an analysis-file crop.

        Each axis is cropped by snapping the requested bounds onto the native grid's nodes; the node
        spacing itself never changes, only the range and the point count. An axis absent from ``crop``
        keeps its full native range.

        :param crop: A mapping from variable name to a ``(min, max)`` pair, or ``None`` to use the
            native range for every axis.
        :type crop: dict | None
        :returns: One dictionary per native axis, in axis order, with the keys ``variable``, ``min``,
            ``max``, ``points``.
        :rtype: list[dict]
        :raises ValueError: If a cropped axis spans fewer than two native nodes, or its point count is
            odd (the real DFT that the convolution relies on requires an even count per axis).
        """
        crop = crop or {}
        result = []
        for axis in self.axes:
            spacing = (axis.max - axis.min) / (axis.points - 1)
            user_min, user_max = crop.get(axis.variable, (axis.min, axis.max))

            index_min = round((user_min - axis.min) / spacing)
            index_max = round((user_max - axis.min) / spacing)
            index_min = max(0, min(axis.points - 1, index_min))
            index_max = max(0, min(axis.points - 1, index_max))

            if index_max <= index_min:
                raise ValueError(
                    f"Axis '{axis.variable}': cropped range [{user_min}, {user_max}] contains fewer than "
                    'two native grid nodes'
                )

            points = index_max - index_min + 1
            if points % 2 != 0:
                raise ValueError(
                    f"Axis '{axis.variable}': the cropped range yields {points} points, which is odd; the "
                    'convolution requires an even point count per axis'
                )

            result.append({
                'variable': axis.variable,
                'min':      axis.min + index_min * spacing,
                'max':      axis.min + index_max * spacing,
                'points':   points,
            })

        return result

    def resolution_kernel(self, cropped_axes, *, tolerance=1.0e-9):
        r"""Slice the stored sampled-resolution kernel onto the grid required by ``cropped_axes``.

        Compares the declared offset spacing against the spacing implied by ``cropped_axes`` within
        ``tolerance``, checks that every required offset lands on a declared node within ``tolerance``,
        and slices the stored values onto the required offsets. The required offsets are
        :math:`(i - N/2) \cdot \mathrm{spacing}` for :math:`i = 0 \ldots N-1`, i.e. numpy's
        ``fftshift`` convention for even :math:`N`. A declared grid wider than required (e.g. a
        symmetric, odd-count sample set) simply leaves its outermost samples unused.

        :param cropped_axes: The per-axis grid, as returned by :meth:`cropped_axes`.
        :type cropped_axes: list[dict]
        :param tolerance: The absolute tolerance used for the spacing and offset-alignment checks.
        :type tolerance: float
        :returns: The resolution kernel, in centred order, with shape matching ``cropped_axes``.
        :rtype: numpy.ndarray
        :raises RuntimeError: If :attr:`resolution` is not a sampled resolution, if the declared
            spacing does not match the required spacing, or if a required offset is not among the
            declared samples.
        """
        if not isinstance(self.resolution, UnbinnedSampledResolutionDescription):
            raise RuntimeError('resolution_kernel() requires a sampled (kind: samples) resolution')

        values = self.resolution_values
        for axis_index, (offset_axis, axis) in enumerate(zip(self.resolution.offsets, cropped_axes)):
            points = axis['points']
            spacing = (axis['max'] - axis['min']) / (points - 1)

            if abs(offset_axis.spacing - spacing) > tolerance:
                raise RuntimeError(
                    f"Axis '{offset_axis.variable}': the declared resolution spacing {offset_axis.spacing} "
                    f'does not match the required spacing {spacing} within tolerance {tolerance}'
                )

            required_offset_0 = -(points // 2) * spacing
            index_0 = round((required_offset_0 - offset_axis.min) / offset_axis.spacing)
            declared_node = offset_axis.min + index_0 * offset_axis.spacing

            if (
                abs(declared_node - required_offset_0) > tolerance
                or index_0 < 0
                or index_0 + points > offset_axis.points
            ):
                raise RuntimeError(
                    f"Axis '{offset_axis.variable}': the required offset {required_offset_0} is not among "
                    f'the declared resolution samples within tolerance {tolerance}'
                )

            values = _np.take(values, range(index_0, index_0 + points), axis=axis_index)

        return values

    def sampled_resolution_kernel(self, cropped_axes, *, parameters=None, options=None):
        r"""Sample an evaluatable resolution (``kind: expression`` or ``kind: pdf``) onto ``cropped_axes``.

        Evaluates the resolution once at every node of the offset grid required by ``cropped_axes``,
        in the same centred order that :meth:`resolution_kernel` slices a sampled resolution into.

        :param cropped_axes: The per-axis grid, as returned by :meth:`cropped_axes`.
        :type cropped_axes: list[dict]
        :param parameters: The parameters the resolution is evaluated with. Defaults to a fresh
            :class:`eos.Parameters` object.
        :type parameters: eos.Parameters | None
        :param options: The options the resolution is evaluated with. Defaults to no options.
        :type options: eos.Options | None
        :returns: The resolution kernel, in centred order, with shape matching ``cropped_axes``.
        :rtype: numpy.ndarray
        :raises RuntimeError: If :attr:`resolution` is a sampled resolution.
        """
        parameters = parameters if parameters is not None else eos.Parameters()
        options = options if options is not None else eos.Options()

        if isinstance(self.resolution, UnbinnedPDFResolutionDescription):
            pdf_name = self.resolution.pdf

            def _evaluate(kinematics):
                return eos.SignalPDF.make(pdf_name, parameters, kinematics, options).evaluate_linear()
        elif isinstance(self.resolution, UnbinnedExpressionResolutionDescription):
            digest = _hashlib.sha1(self.resolution.expression.encode()).hexdigest()[:16]
            name = eos.QualifiedName(f'Unbinned::resolution_{digest}@kernel')
            if name not in eos.Observables():
                eos.Observables().insert(name, '', eos.Unit.Unity(), options, self.resolution.expression)

            def _evaluate(kinematics):
                return eos.Observable.make(name, parameters, kinematics, options).evaluate()
        else:
            raise RuntimeError('sampled_resolution_kernel() requires an evaluatable (kind: expression or kind: pdf) resolution')

        shape = tuple(axis['points'] for axis in cropped_axes)
        values = _np.empty(shape, dtype=float)
        for flat_index in _np.ndindex(shape):
            kinematics = eos.Kinematics()
            for axis, index in zip(cropped_axes, flat_index):
                points = axis['points']
                spacing = (axis['max'] - axis['min']) / (points - 1)
                kinematics.declare(axis['variable'], (index - points // 2) * spacing)
            values[flat_index] = _evaluate(kinematics)

        return values

    @staticmethod
    def create(path, axes, resolution, observations, resolution_values=None):
        """Write a new UnbinnedLikelihood object to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param axes: The native sampling axes, in axis order; each an :class:`UnbinnedAxisDescription`
            or an equivalent keyword mapping.
        :type axes: list
        :param resolution: The resolution description; one of :class:`UnbinnedSampledResolutionDescription`,
            :class:`UnbinnedExpressionResolutionDescription`, or :class:`UnbinnedPDFResolutionDescription`.
        :type resolution: UnbinnedSampledResolutionDescription | UnbinnedExpressionResolutionDescription | UnbinnedPDFResolutionDescription
        :param observations: The observed events, shape ``(n_events, rank)``, columns in axis order.
        :type observations: numpy.ndarray or nested list of float
        :param resolution_values: The resolution kernel, in centred order relative to ``resolution.offsets``,
            with shape equal to the offsets' native grid. Mandatory for a sampled resolution, ignored
            otherwise.
        :type resolution_values: numpy.ndarray or nested list of float, optional
        """
        axes = [a if isinstance(a, UnbinnedAxisDescription) else UnbinnedAxisDescription(**a) for a in axes]
        rank = len(axes)

        observations = _np.asarray(observations, dtype=float)
        if observations.ndim != 2 or observations.shape[1] != rank:
            raise RuntimeError(f'observations has shape {observations.shape}, expected (n_events, {rank})')

        if isinstance(resolution, UnbinnedSampledResolutionDescription):
            if resolution_values is None:
                raise RuntimeError('resolution_values must be given for a sampled resolution')

            expected_shape = tuple(o.points for o in resolution.offsets)
            resolution_values = _np.asarray(resolution_values, dtype=float)
            if resolution_values.shape != expected_shape:
                raise RuntimeError(f'resolution_values has shape {resolution_values.shape}, expected {expected_shape} from the declared offsets')

        description = UnbinnedLikelihoodDescription(version=eos.__version__, axes=axes, resolution=resolution)

        os.makedirs(path, exist_ok=True)
        description.to_yaml_file(os.path.join(path, 'description.yaml'))
        _np.save(os.path.join(path, 'observations.npy'), observations)

        if isinstance(resolution, UnbinnedSampledResolutionDescription):
            _np.save(os.path.join(path, 'resolution.npy'), resolution_values)
