# Copyright (c) 2022-2025 Danny van Dyk
# Copyright (c) 2023 Méril Reboud
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


class ImportanceSamples:
    def __init__(self, path):
        """ Read an ImportanceSamples object from disk.

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

        if not description['type'] == 'ImportanceSamples':
            raise RuntimeError(f'Path {path} not pointing to an ImportanceSamples object')

        self.type = 'ImportanceSamples'
        self.varied_parameters = description['parameters']
        self.lookup_table = { item['name']: idx for idx, item in enumerate(self.varied_parameters) }

        f = os.path.join(path, 'samples.npy')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Samples file {f} does not exist or is not a file')
        self.samples = _np.load(f)

        f = os.path.join(path, 'weights.npy')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Weights file {f} does not exist or is not a file')
        self.weights = _np.load(f)

        f = os.path.join(path, 'posterior_values.npy')
        if not os.path.exists(f) and not os.path.isfile(f):
            self.posterior_values = None
        else:
            self.posterior_values = _np.load(f)


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
        description = {}
        description['version'] = eos.__version__
        description['type'] = 'ImportanceSamples'
        description['parameters'] = [{
            'name': p.name(),
            'min': p.min() if 'min' in dir(p) else -_np.inf,
            'max': p.max() if 'max' in dir(p) else +_np.inf
        } for p in parameters]

        if not samples.shape[1] == len(parameters):
            raise RuntimeError(f'Shape of samples {samples.shape} incompatible with number of parameters {len(parameters)}')

        if not weights is None and not samples.shape[0] == weights.shape[0]:
            raise RuntimeError(f'Shape of weights {weights.shape} incompatible with shape of samples {samples.shape}')

        if not posterior_values is None and not samples.shape[0] == posterior_values.shape[0]:
            raise RuntimeError(f'Shape of posterior values {posterior_values.shape} incompatible with shape of samples {samples.shape}')

        os.makedirs(path, exist_ok=True)
        with open(os.path.join(path, 'description.yaml'), 'w') as description_file:
            yaml.dump(description, description_file, default_flow_style=False)
        _np.save(os.path.join(path, 'samples.npy'), samples)
        _np.save(os.path.join(path, 'weights.npy'), weights)
        if not posterior_values is None:
            _np.save(os.path.join(path, 'posterior_values.npy'), posterior_values)
