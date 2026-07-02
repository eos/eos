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
from eos.deserializable import Deserializable
from eos.serializable import Serializable

import numpy as _np
import os as _os


@dataclass(kw_only=True)
class ParameterDescription(Serializable, Deserializable):
    r"""Describes a single varied parameter recorded in a data object's ``description.yaml``.

    :param name: The qualified name of the parameter.
    :type name: str
    :param min: The lower bound of the parameter's prior range. Defaults to :math:`-\infty`.
    :type min: float
    :param max: The upper bound of the parameter's prior range. Defaults to :math:`+\infty`.
    :type max: float
    """
    name:str
    min:float = field(default=-_np.inf)
    max:float = field(default=+_np.inf)


def load_array(directory:str, filename:str, *, ncols:int=None, nrows:int=None):
    r"""Load a NumPy array stored beside a data object's description and validate its shape.

    Centralizes the presence and shape checks that the individual :class:`eos.data` objects would
    otherwise repeat when loading their ``.npy`` payloads. The optional ``ncols``/``nrows`` arguments
    let the caller assert that the payload is consistent with the deserialized description (e.g. that
    the number of sample columns matches the number of varied parameters).

    :param directory: The storage directory that contains the array file.
    :type directory: str
    :param filename: The name of the ``.npy`` file within ``directory``.
    :type filename: str
    :param ncols: If given, require a 2D array with exactly this many columns.
    :type ncols: int | None
    :param nrows: If given, require exactly this many entries along the first axis.
    :type nrows: int | None
    :returns: The loaded array.
    :rtype: numpy.ndarray
    :raises RuntimeError: If the file is missing, is not a file, or its shape does not match.
    """
    path = _os.path.join(directory, filename)
    if not _os.path.exists(path) or not _os.path.isfile(path):
        raise RuntimeError(f'Data file {path} does not exist or is not a file')

    array = _np.load(path)

    if ncols is not None and (array.ndim != 2 or array.shape[1] != ncols):
        raise RuntimeError(f'Data file {path} has shape {array.shape}, expected a 2D array with {ncols} columns')

    if nrows is not None and (array.ndim < 1 or array.shape[0] != nrows):
        raise RuntimeError(f'Data file {path} has shape {array.shape}, expected {nrows} entries along the first axis')

    return array
