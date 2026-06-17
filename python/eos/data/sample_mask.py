# Copyright (c) 2025 Matthew Kirk
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

import eos
import os
import numpy as _np
import yaml


class SampleMask:
    r"""Represents a boolean mask over a set of observables, stored on disk.

    A mask selects a subset of samples (e.g. to remove unphysical points) and records the observables it
    applies to. Instances are created either by reading an existing mask from disk (passing its ``path``
    to the constructor) or by writing a new mask with :meth:`create`.

    :ivar type: The type identifier of the data object, always ``'Mask'``.
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

        f = os.path.join(path, 'description.yaml')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Description file {f} does not exist or is not a file')

        with open(f) as df:
            description = yaml.load(df, Loader=yaml.SafeLoader)

        if not description['type'] == 'Mask':
            raise RuntimeError(f'Path {path} not pointing to a Mask file')

        self.type = 'Mask'
        self.observables = description['observables']

        f = os.path.join(path, 'mask.npy')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Mask file {f} does not exist or is not a file')
        self.mask = _np.load(f)

    @staticmethod
    def create(path, mask, observables):
        """ Write a new Mask object to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param mask: Mask as a 1D boolean array
        :type mask: list or iterable of bool
        :param observables: Observables as a 1D array
        """
        description = {}
        description['version'] = eos.__version__
        description['type'] = 'Mask'
        description['observables'] = observables

        os.makedirs(path, exist_ok=True)
        with open(os.path.join(path, 'description.yaml'), 'w') as description_file:
            yaml.dump(description, description_file, default_flow_style=False)
        _np.save(os.path.join(path, 'mask.npy'), mask)
