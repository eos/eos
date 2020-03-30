# Copyright (c) 2019 Danny van Dyk
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
import pypmc
import yaml

class MarkovChain:
    def __init__(self, path):
        """ Read a MarkovChain object from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError('Path {} does not exist or is not a directory'.formt(path))

        f = os.path.join(path, 'description.yaml')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError('Description file {} does not exist or is not a file'.format(f))

        with open(f, 'r') as df:
            description = yaml.load(df, Loader=yaml.SafeLoader)

        if not description['type'] == 'MarkovChain':
            raise RuntimeError('Path {} not pointing to a MarkovChain'.format(path))

        self.type = 'MarkovChain'
        self.varied_parameters = description['parameters']

        f = os.path.join(path, 'samples.npy')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError('Samples file {} does not exist or is not a file'.format(f))
        self.samples = _np.load(f)

        if description['has-weights']:
            f = os.path.join(path, 'weights.npy')
            if not os.path.exists(f) or not os.path.isfile(f):
                raise RuntimeError('Weights file {} does not exist or is not a file'.format(f))
            self.weights = _np.load(f)
        else:
            self.weights = None


    @staticmethod 
    def create(path, parameters, samples, weights=None):
        """ Write a new MarkovChain object to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param parameters: Parameter descriptions as a 1D array of shape (N, ).
        :type parameters: list or iterable of eos.Parameter
        :param samples: Samples as a 2D array of shape (N, P).
        :type samples: 2D numpy array
        :param weights: Weights on a linear scale as a 2D array of shape (N, 1).
        :type weights: 2D numpy array, optional
        """
        description = {}
        description['version'] = eos.version()
        description['type'] = 'MarkovChain'
        description['parameters'] = [{
            'name': p.name(),
            'min': p.min(),
            'max': p.max()
        } for p in parameters]
        description['has-weights'] = (not weights is None)

        if not samples.shape[1] == len(parameters):
            raise RuntimeError('Shape of samples {} incompatible with number of parameters {}'.format(samples.shape, len(parameters)))

        if not weights is None and not samples.shape[0] == weights.shape[0]:
            raise RuntimeError('Shape of weights {} incompatible with shape of samples {}'.format(weights.shape, samples.shape))

        os.makedirs(path, exist_ok=True)
        with open(os.path.join(path, 'description.yaml'), 'w') as description_file:
            yaml.dump(description, description_file, default_flow_style=False)
        _np.save(os.path.join(path, 'samples.npy'), samples)

        if not weights is None:
            _np.save(os.path.join(path, 'weights.npy'), weights)


class MixtureDensity:
    def __init__(self, path):
        """ Read a MixtureDensity object from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError('Path {} does not exist or is not a directory'.formt(path))

        f = os.path.join(path, 'description.yaml')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError('Description file {} does not exist or is not a file'.format(f))

        with open(f, 'r') as df:
            description = yaml.load(df, Loader=yaml.SafeLoader)

        if not description['type'] == 'MixtureDensity':
            raise RuntimeError('Path {} not pointing to a MixtureDensity'.format(path))

        self.type = 'MixtureDensity'
        self.components = description['components']
        self.weights    = description['weights']


    def density(self):
        """ Return a pypmc.density.MixtureDensity. """
        components = [pypmc.density.gauss.Gauss(c['mu'], c['sigma']) for c in self.components]
        return pypmc.density.mixture.MixtureDensity(components, self.weights)


    @staticmethod 
    def create(path, density):
        """ Write a new MixtureDensity object to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param density: Mixture density.
        :type density: pypmc.density.MixtureDensity
        """
        description = {}
        description['version'] = eos.version()
        description['type'] = 'MixtureDensity'
        description['components'] = []
        for c in density.components:
            if type(c) is pypmc.density.gauss.Gauss:
                description['components'].append({
                    'type': 'gauss',
                    'mu': c.mu.tolist(),
                    'sigma': c.sigma.tolist()
                })
            else:
                raise RuntimeError('Unsupported type of MixtureDensity component: {}'.format(type(c)))
        description['weights'] = density.weights.tolist()

        os.makedirs(path, exist_ok=True)
        with open(os.path.join(path, 'description.yaml'), 'w') as description_file:
            yaml.dump(description, description_file, default_flow_style=False)


class PMCSampler:
    def __init__(self, path):
        """ Read a PMCSampler object from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError('Path {} does not exist or is not a directory'.formt(path))

        f = os.path.join(path, 'description.yaml')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError('Description file {} does not exist or is not a file'.format(f))

        with open(f, 'r') as df:
            description = yaml.load(df, Loader=yaml.SafeLoader)

        if not description['type'] == 'PMCSampler':
            raise RuntimeError('Path {} not pointing to a PMCSampler'.format(path))

        self.type = 'PMCSampler'
        self.varied_parameters = description['parameters']
        self.components = _np.array([pypmc.density.gauss.Gauss(_np.array(c['mu']), _np.array(c['sigma'])) for c in description['proposal']['components']])
        self.weights    = _np.array(description['proposal']['weights'])

        f = os.path.join(path, 'samples.npy')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError('Samples file {} does not exist or is not a file'.format(f))
        self.samples = _np.load(f)

        f = os.path.join(path, 'weights.npy')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError('Weights file {} does not exist or is not a file'.format(f))
        self.weights = _np.load(f)


    @staticmethod 
    def create(path, parameters, samples, weights, proposal):
        """ Write a new PMCSampler object to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param parameters: Parameter descriptions as a 1D array of shape (N, ).
        :type parameters: list or iterable of eos.Parameter
        :param samples: Samples as a 2D array of shape (N, P).
        :type samples: 2D numpy array
        :param weights: Weights on a linear scale as a 2D array of shape (N, 1).
        :type weights: 1D numpy array, optional
        """
        description = {}
        description['version'] = eos.version()
        description['type'] = 'PMCSampler'
        description['parameters'] = [{
            'name': p.name(),
            'min': p.min(),
            'max': p.max()
        } for p in parameters]
        description['proposal'] = {
            'components': [
                { 'mu': c.mu.tolist(), 'sigma': c.sigma.tolist() } for c in proposal.components
            ],
            'weights': proposal.weights.tolist()
        }

        if not samples.shape[1] == len(parameters):
            raise RuntimeError('Shape of samples {} incompatible with number of parameters {}'.format(samples.shape, len(parameters)))

        if not weights is None and not samples.shape[0] == weights.shape[0]:
            raise RuntimeError('Shape of weights {} incompatible with shape of samples {}'.format(weights.shape, samples.shape))

        os.makedirs(path, exist_ok=True)
        with open(os.path.join(path, 'description.yaml'), 'w') as description_file:
            yaml.dump(description, description_file, default_flow_style=False)
        _np.save(os.path.join(path, 'samples.npy'), samples)
        _np.save(os.path.join(path, 'weights.npy'), weights)


class Prediction:
    def __init__(self, path):
        """ Read a Prediction object from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError('Path {} does not exist or is not a directory'.formt(path))

        f = os.path.join(path, 'description.yaml')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError('Description file {} does not exist or is not a file'.format(f))

        with open(f, 'r') as df:
            description = yaml.load(df, Loader=yaml.SafeLoader)

        if not description['type'] == 'Prediction':
            raise RuntimeError('Path {} not pointing to a Prediction'.format(path))

        self.type = 'Prediction'
        self.varied_parameters = description['observables']

        f = os.path.join(path, 'samples.npy')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError('Samples file {} does not exist or is not a file'.format(f))
        self.samples = _np.load(f)

        f = os.path.join(path, 'weights.npy')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError('Weights file {} does not exist or is not a file'.format(f))
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
        description['version'] = eos.version()
        description['type'] = 'Prediction'
        description['observables'] = [{
            'name': str(o.name()),
            'kinematics': { k.name(): float(k) for k in o.kinematics() }
        } for o in observables]

        if not samples.shape[1] == len(observables):
            raise RuntimeError('Shape of samples {} incompatible with number of observables {}'.format(samples.shape, len(observables)))

        if not samples.shape[0] == weights.shape[0]:
            raise RuntimeError('Shape of weights {} incompatible with shape of samples {}'.format(weights.shape, samples.shape))

        os.makedirs(path, exist_ok=True)
        with open(os.path.join(path, 'description.yaml'), 'w') as description_file:
            yaml.dump(description, description_file, default_flow_style=False)
        _np.save(os.path.join(path, 'samples.npy'), samples)
        _np.save(os.path.join(path, 'weights.npy'), weights)
