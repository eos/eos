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
from scipy.special import erf
from scipy.linalg import block_diag

class MarkovChain:
    def __init__(self, path):
        """ Read a MarkovChain object from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError('Path {} does not exist or is not a directory'.format(path))

        f = os.path.join(path, 'description.yaml')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError('Description file {} does not exist or is not a file'.format(f))

        with open(f, 'r') as df:
            description = yaml.load(df, Loader=yaml.SafeLoader)

        if not description['type'] == 'MarkovChain':
            raise RuntimeError('Path {} not pointing to a MarkovChain'.format(path))

        self.type = 'MarkovChain'
        self.varied_parameters = description['parameters']
        self.lookup_table = { item['name']: idx for idx, item in enumerate(self.varied_parameters) }

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
        description['version'] = eos.__version__
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
            raise RuntimeError('Path {} does not exist or is not a directory'.format(path))

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
        description['version'] = eos.__version__
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

    @staticmethod
    def _cartesian_product(densityA, densityB):
        """ Construct the cartesian product of two mixture densities. The means and covariances
        of the components are copied and concatenated and the weights are multiplied. Note that
        this product is not commutative, the order of densities affects the order of parameters
        in the output density. Only Gaussian mixtures are supported.

        :param densityA: First mixture density.
        :type densityA: pypmc.density.MixtureDensity
        :param densityB: Second mixture density.
        :type densityB: pypmc.density.MixtureDensity
        """
        components = []
        weights = []
        for cA, wA in zip(densityA.components, densityA.weights):
            if type(cA) is pypmc.density.gauss.Gauss:
                for cB, wB in zip(densityB.components, densityB.weights):
                    if type(cB) is pypmc.density.gauss.Gauss:
                        cAB_mu = _np.concatenate((cA.mu, cB.mu))
                        cAB_sigma = block_diag(cA.sigma, cB.sigma)
                        components.append(pypmc.density.gauss.Gauss(cAB_mu, cAB_sigma))
                        weights.append(wA * wB)
                    else:
                        raise RuntimeError('Unsupported type of MixtureDensity component: {}'.format(type(cB)))
            else:
                raise RuntimeError('Unsupported type of MixtureDensity component: {}'.format(type(cA)))
        return pypmc.density.mixture.MixtureDensity(components, weights / _np.sum(weights))

    @staticmethod
    def cartesian_product(densities):
        """ Construct the cartesian product of a list of mixture densities. The means and covariances
        of the components are copied and concatenated and the weights are multiplied. Note that
        this product is not commutative, the order of densities affects the order of parameters
        in the output density. Only Gaussian mixtures are supported.

        :param densities: List of mixture densities.
        :type densities: iterable of pypmc.density.MixtureDensity
        """

        product_density = densities[0]
        for density in densities[1:]:
            product_density = MixtureDensity._cartesian_product(product_density, density)

        return product_density


class PMCSampler:
    def __init__(self, path):
        """ Read a PMCSampler object from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError('Path {} does not exist or is not a directory'.format(path))

        f = os.path.join(path, 'description.yaml')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError('Description file {} does not exist or is not a file'.format(f))

        with open(f, 'r') as df:
            description = yaml.load(df, Loader=yaml.SafeLoader)

        if not description['type'] == 'PMCSampler':
            raise RuntimeError('Path {} not pointing to a PMCSampler'.format(path))

        self.type = 'PMCSampler'
        self.varied_parameters = description['parameters']
        self.lookup_table = { item['name']: idx for idx, item in enumerate(self.varied_parameters) }
        self.components        = _np.array([pypmc.density.gauss.Gauss(_np.array(c['mu']), _np.array(c['sigma'])) for c in description['proposal']['components']])
        self.component_weights = _np.array(description['proposal']['weights'])


    def density(self):
        """ Return a pypmc.density.MixtureDensity. """
        return pypmc.density.mixture.MixtureDensity(self.components, self.component_weights)


    @staticmethod
    def _evaluate_mixture_pdf(density, x):
        """ Internal function that evaluates the PDF of the mixture density at point x."""
        def evaluate_component(comp, x):
            # To limit numerical instabilities, the covariance matrix is devided by the median of its diagonal entries
            median_variance = _np.median(_np.diag(comp.sigma))
            chi2 = _np.dot(
                _np.transpose(x - comp.mu) / median_variance,
                _np.dot(_np.linalg.inv(comp.sigma / median_variance), x - comp.mu)
            )
            return (2 * _np.pi * median_variance)**(-0.5 * len(x)) * _np.exp(-0.5 * chi2) \
                / _np.sqrt(_np.linalg.det(comp.sigma / median_variance))

        return _np.sum(
            _np.array(density.weights)
            * _np.array([evaluate_component(comp, _np.array(x)) for comp in density.components])
        )


    @staticmethod
    def create(path, parameters, proposal, sigma_test_stat=None, samples=None, weights=None):
        """ Write a new PMCSampler object to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param parameters: Parameter descriptions as a 1D array of shape (P, ).
        :type parameters: list or iterable of eos.Parameter
        :param sigma_test_stat: (optional) If provided, the inverse CDF of -2*log(PDF) will be evaluated, using the provided values as the respective significance.
        :type sigma_test_stat: list or iterable
        :param samples: Samples as a 2D array of shape (N, P). Needed to generate the test statistic.
        :type samples: 2D numpy array, optional
        :param weights: Weights on a linear scale as a 2D array of shape (N, 1). Needed to generate the test statistic.
        :type weights: 1D numpy array, optional
        """
        description = {}
        description['version'] = eos.__version__
        description['type'] = 'PMCSampler'
        description['parameters'] = [{
            'name': p.name(),
            'min': p.min(),
            'max': p.max()
        } for p in parameters]

        # Don't write components that have a 0 weight
        purged_components = []
        purged_weights = []
        for comp, w in zip(proposal.components, proposal.weights.tolist()):
            if w != 0:
                purged_components.append({ 'mu': comp.mu.tolist(), 'sigma': comp.sigma.tolist() })
                purged_weights.append(w)
        description['proposal'] = { 'components': purged_components, 'weights': purged_weights }

        # The test statistics defaults to two empty lists
        description['test statistics'] = { "sigma": [], "densities": [] }
        if sigma_test_stat and samples and weights:
            sigma_test_stat = _np.array(sigma_test_stat)
            samplesPDF = [x for x in map(lambda x: -2.0 * _np.log(PMCSampler._evaluate_mixture_pdf(proposal, x)), samples)]
            ind = _np.argsort(samplesPDF)
            sorted_samplesPDF = _np.array(samplesPDF)[ind]
            sorted_weights = weights[ind]
            cumulant = 1.*sorted_weights.cumsum()/sorted_weights.sum()
            percents = erf(sigma_test_stat/_np.sqrt(2))
            weighted_percentile = _np.interp(percents, cumulant, sorted_samplesPDF)

            description['test statistics'] = {
                "sigma": sigma_test_stat.tolist(),
                "densities": weighted_percentile.tolist()
            }

        os.makedirs(path, exist_ok=True)
        with open(os.path.join(path, 'description.yaml'), 'w') as description_file:
            yaml.dump(description, description_file, default_flow_style=False)


class ImportanceSamples:
    def __init__(self, path):
        """ Read an ImportanceSamples object from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError('Path {} does not exist or is not a directory'.format(path))

        f = os.path.join(path, 'description.yaml')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError('Description file {} does not exist or is not a file'.format(f))

        with open(f, 'r') as df:
            description = yaml.load(df, Loader=yaml.SafeLoader)

        if not description['type'] == 'ImportanceSamples':
            raise RuntimeError('Path {} not pointing to an ImportanceSamples object'.format(path))

        self.type = 'ImportanceSamples'
        self.varied_parameters = description['parameters']
        self.lookup_table = { item['name']: idx for idx, item in enumerate(self.varied_parameters) }

        f = os.path.join(path, 'samples.npy')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError('Samples file {} does not exist or is not a file'.format(f))
        self.samples = _np.load(f)

        f = os.path.join(path, 'weights.npy')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError('Weights file {} does not exist or is not a file'.format(f))
        self.weights = _np.load(f)


    @staticmethod
    def create(path, parameters, samples, weights):
        """ Write a new ImportanceSamples object to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param parameters: Parameter descriptions as a 1D array of shape (P, ).
        :type parameters: list or iterable of eos.Parameter
        :param samples: Samples as a 2D array of shape (N, P).
        :type samples: 2D numpy array
        :param weights: Weights on a linear scale as a 2D array of shape (N, 1).
        :type weights: 1D numpy array, optional
        """
        description = {}
        description['version'] = eos.__version__
        description['type'] = 'ImportanceSamples'
        description['parameters'] = [{
            'name': p.name(),
            'min': p.min(),
            'max': p.max()
        } for p in parameters]

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
            raise RuntimeError('Path {} does not exist or is not a directory'.format(path))

        f = os.path.join(path, 'description.yaml')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError('Description file {} does not exist or is not a file'.format(f))

        with open(f, 'r') as df:
            description = yaml.load(df, Loader=yaml.SafeLoader)

        if not description['type'] == 'Prediction':
            raise RuntimeError('Path {} not pointing to a Prediction'.format(path))

        self.type = 'Prediction'
        self.varied_parameters = description['observables']
        self.lookup_table = { item['name']: idx for idx, item in enumerate(self.varied_parameters) }

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
        description['version'] = eos.__version__
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
