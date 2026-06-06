# Copyright (c) 2020-2023 Danny van Dyk
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

import eos
import os
import numpy as _np
import yaml


class Prediction:
    def __init__(self, path):
        """ Read a Prediction object from disk.

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

        if not description['type'] == 'Prediction':
            raise RuntimeError(f'Path {path} not pointing to a Prediction')

        self.type = 'Prediction'
        self.varied_parameters = description['observables']
        self.lookup_table = {}
        for idx, item in enumerate(self.varied_parameters):
            id = item['name']
            if 'options' in item:
                id += ';' + str(eos.Options(item['options'])).replace(" ", "")
            if 'kinematics' in item:
                id += '[' + str(eos.Kinematics(item['kinematics'])).replace(" ", "") + ']'
            self.lookup_table[id] = idx

        f = os.path.join(path, 'samples.npy')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Samples file {f} does not exist or is not a file')
        self.samples = _np.load(f)

        f = os.path.join(path, 'weights.npy')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Weights file {f} does not exist or is not a file')
        self.weights = _np.load(f)


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
        description = {}
        description['version'] = eos.__version__
        description['type'] = 'Prediction'
        description['observables'] = [{
            'name': o.name().full(),
            'kinematics': { k.name(): float(k) for k in o.kinematics() },
            'options': { str(k): str(v) for k, v in o.options() }
        } for o in observables]

        if not samples.shape[1] == len(observables):
            raise RuntimeError(f'Shape of samples {samples.shape} incompatible with number of observables {len(observables)}')

        if not samples.shape[0] == weights.shape[0]:
            raise RuntimeError(f'Shape of weights {weights.shape} incompatible with shape of samples {samples.shape}')

        os.makedirs(path, exist_ok=True)
        with open(os.path.join(path, 'description.yaml'), 'w') as description_file:
            yaml.dump(description, description_file, default_flow_style=False)
        _np.save(os.path.join(path, 'samples.npy'), samples)
        _np.save(os.path.join(path, 'weights.npy'), weights)
