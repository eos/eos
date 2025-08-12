# Copyright (c) 2019-2024 Danny van Dyk
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
import dynesty
from scipy.special import erf
from scipy.linalg import block_diag

class Mode:
    def __init__(self, path):
        """ Read a posterior's (local) mode from a file.

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

        if not description['type'] == 'Mode':
            raise RuntimeError(f'Path {path} not pointing to a Mode file')

        self.type = 'Mode'
        self.varied_parameters = description['parameters']
        self.mode = description['mode']
        self.pvalue = description['pvalue']
        self.local_pvalues = description['local_pvalues']
        self.global_chi2 = None
        if 'global_chi2' in description:
            self.global_chi2 = description['global_chi2']
        self.dof = None
        if 'dof' in description:
            self.dof = description['dof']

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
        description = {}
        description['version'] = eos.__version__
        description['type'] = 'Mode'
        description['parameters'] = [{
            'name': p.name(),
            'min': p.min(),
            'max': p.max()
        } for p in parameters]
        description['mode'] = mode.tolist()
        description['pvalue'] = float(pvalue) if pvalue is not None else None
        description['local_pvalues'] = local_pvalues
        description['global_chi2'] = global_chi2
        description['dof'] = dof

        os.makedirs(path, exist_ok=True)
        with open(os.path.join(path, 'description.yaml'), 'w') as description_file:
            yaml.dump(description, description_file, default_flow_style=False)


class MarkovChain:
    def __init__(self, path):
        """ Read a MarkovChain object from disk.

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

        if not description['type'] == 'MarkovChain':
            raise RuntimeError(f'Path {path} not pointing to a MarkovChain')

        self.type = 'MarkovChain'
        self.varied_parameters = description['parameters']
        self.lookup_table = { item['name']: idx for idx, item in enumerate(self.varied_parameters) }

        f = os.path.join(path, 'samples.npy')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Samples file {f} does not exist or is not a file')
        self.samples = _np.load(f)

        f = os.path.join(path, 'usamples.npy')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'U-space samples file {f} does not exist or is not a file')
        self.usamples = _np.load(f)

        if description['has-weights']:
            f = os.path.join(path, 'weights.npy')
            if not os.path.exists(f) or not os.path.isfile(f):
                raise RuntimeError(f'Weights file {f} does not exist or is not a file')
            self.weights = _np.load(f)
        else:
            self.weights = None


    @staticmethod
    def create(path, parameters, samples, usamples, weights=None):
        """ Write a new MarkovChain object to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param parameters: Parameter descriptions as a 1D array of shape (N, ).
        :type parameters: list or iterable of eos.Parameter
        :param samples: Samples in parameter space as a 2D array of shape (N, P).
        :type samples: 2D numpy array
        :param usamples: Samples in u space as a 2D array of shape (N, P).
        :type usamples: 2D numpy array
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
            raise RuntimeError(f'Shape of samples {samples.shape} incompatible with number of parameters {len(parameters)}')

        if not usamples.shape[1] == len(parameters):
            raise RuntimeError(f'Shape of usamples {usamples.shape} incompatible with number of parameters {len(parameters)}')

        if not weights is None and not samples.shape[0] == weights.shape[0]:
            raise RuntimeError(f'Shape of weights {weights.shape} incompatible with shape of samples {samples.shape}')

        os.makedirs(path, exist_ok=True)
        with open(os.path.join(path, 'description.yaml'), 'w') as description_file:
            yaml.dump(description, description_file, default_flow_style=False)
        _np.save(os.path.join(path, 'samples.npy'), samples)
        _np.save(os.path.join(path, 'usamples.npy'), usamples)

        if not weights is None:
            _np.save(os.path.join(path, 'weights.npy'), weights)


class MixtureDensity:
    def __init__(self, path):
        """ Read a MixtureDensity object from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path):
            raise RuntimeError(f'Path {path} does not exist')

        if os.path.isfile(path):
            f = path
        elif os.path.isdir(path):
            f = os.path.join(path, 'description.yaml')
        else:
            raise RuntimeError(f'Path {path} is neither a file nor a directory')

        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Description file {f} does not exist or is not a file')

        with open(f) as df:
            description = yaml.load(df, Loader=yaml.SafeLoader)

        if not description['type'] == 'MixtureDensity':
            raise RuntimeError(f'Path {path} not pointing to a MixtureDensity')

        self.type = 'MixtureDensity'
        self.components = description['components']
        self.weights    = description['weights']
        self.varied_parameters = description['varied_parameters'] if 'varied_parameters' in description else None
        self.test_statistics = description['test_statistics'] if 'test_statistics' in description else None

    def density(self):
        """ Return a pypmc.density.MixtureDensity. """
        components = [pypmc.density.gauss.Gauss(c['mu'], c['sigma']) for c in self.components]
        return pypmc.density.mixture.MixtureDensity(components, self.weights)


    @staticmethod
    def create(path, density, varied_parameters=None, sigma_test_stat=None, samples=None, weights=None):
        """ Write a new MixtureDensity object to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param density: Mixture density.
        :type density: pypmc.density.MixtureDensity
        :param varied_parameters: List of the qualified names of varied parameters.
        :type varied_parameters: ``numpy.array`` of str, optional
        :param sigma_test_stat: (optional) If provided, the inverse CDF of -2*log(PDF) will be evaluated, using the provided values as the respective significance.
        :type sigma_test_stat: list or iterable
        :param samples: Samples as a 2D array of shape (N, P). Needed to generate the test statistic.
        :type samples: 2D numpy array, optional
        :param weights: Weights on a linear scale as a 2D array of shape (N, 1). Needed to generate the test statistic.
        :type weights: 1D numpy array, optional

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
                raise RuntimeError(f'Unsupported type of MixtureDensity component: {type(c)}')
        description['weights'] = density.weights.tolist()
        if varied_parameters is not None:
            description['varied_parameters'] = list(varied_parameters)

        # The test statistics defaults to two empty lists
        description['test_statistics'] = { "sigma": [], "densities": [] }
        if sigma_test_stat is not None and samples is not None and weights is not None:
            sigma_test_stat = _np.array(sigma_test_stat)
            samples_pdf = list(map(lambda x: -density.evaluate(x), samples))
            ind = _np.argsort(samples_pdf)
            samples_pdf_sorted = _np.array(samples_pdf)[ind]
            sorted_weights = weights[ind]
            cumulant = 1.0 * sorted_weights.cumsum() / sorted_weights.sum()
            percents = erf(sigma_test_stat/_np.sqrt(2))
            weighted_percentile = _np.interp(percents, cumulant, samples_pdf_sorted)

            description['test_statistics'] = {
                "sigma": sigma_test_stat.tolist(),
                "densities": weighted_percentile.tolist()
            }

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
                        raise RuntimeError(f'Unsupported type of MixtureDensity component: {type(cB)}')
            else:
                raise RuntimeError(f'Unsupported type of MixtureDensity component: {type(cA)}')
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
            raise RuntimeError(f'Path {path} does not exist or is not a directory')

        f = os.path.join(path, 'description.yaml')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Description file {f} does not exist or is not a file')

        with open(f) as df:
            description = yaml.load(df, Loader=yaml.SafeLoader)

        if not description['type'] == 'PMCSampler':
            raise RuntimeError(f'Path {path} not pointing to a PMCSampler')

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
        if sigma_test_stat is not None and samples is not None and weights is not None:
            sigma_test_stat = _np.array(sigma_test_stat)
            samplesPDF = list(map(lambda x: -2.0 * _np.log(PMCSampler._evaluate_mixture_pdf(proposal, x)), samples))
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
        :param weights: Weights on a linear scale as a 2D array of shape (N, 1).
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


class DynestyResults:
    def __init__(self, path):
        """ Read Results object (in the dynesty.results module) from disk.

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

        if not description['type'] == 'DynestyResults':
            raise RuntimeError(f'Path {path} not pointing to a DynestyResults file')

        self.type = 'DynestyResults'
        self.varied_parameters = description['parameters']
        self.lookup_table = { item['name']: idx for idx, item in enumerate(self.varied_parameters) }

        f = os.path.join(path, 'dynesty_results.npy')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Dynesty results file {f} does not exist or is not a file')

        res_dict = _np.load(f, allow_pickle=True).item()
        if "blob" in res_dict:
            res_dict.pop('blob')
        self.results = dynesty.results.Results(res_dict)
        self.samples = self.results.samples
        self.weights = _np.exp(self.results.logwt - self.results.logz[-1])


    @staticmethod
    def create(path, parameters, results):
        """ Write a new Results object (in the dynesty.results module) to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param parameters: Parameter descriptions as a 1D array of shape (N, ).
        :type parameters: list or iterable of eos.Parameter
        :param results: The results of a nested sampling run.
        :type results: dynesty.results.Results
        """
        description = {}
        description['version'] = eos.__version__
        description['type'] = 'DynestyResults'
        description['parameters'] = [{
            'name': p.name(),
            'min': p.min(),
            'max': p.max()
        } for p in parameters]

        os.makedirs(path, exist_ok=True)
        with open(os.path.join(path, 'description.yaml'), 'w') as description_file:
            yaml.dump(description, description_file, default_flow_style=False)

        res_dict = results.asdict()
        _np.save(os.path.join(path, 'dynesty_results.npy'), res_dict)


class NabuLikelihood:
    def __init__(self, path):
        """ Read a nabu serialized likelihood from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        import nabu

        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError(f'Path {path} does not exist or is not a directory')

        f = os.path.join(path, 'description.yaml')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Description file {f} does not exist or is not a file')

        with open(f) as df:
            description = yaml.load(df, Loader=yaml.SafeLoader)

        if not description['type'] == 'NabuLikelihood':
            raise RuntimeError(f'Path {path} not pointing to a NabuLikelihood file')

        self.type = 'NabuLikelihood'
        self.varied_parameters = description['parameters']
        self.lookup_table = { item['name']: idx for idx, item in enumerate(self.varied_parameters) }

        f = os.path.join(path, 'likelihood.nabu')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Nabu likelihood file {f} does not exist or is not a file')

        self.likelihood = nabu.Likelihood.load(f)


    @staticmethod
    def create(path, parameters, likelihood):
        """ Write a new Results object (in the dynesty.results module) to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param parameters: Parameter descriptions as a 1D array of shape (N, ).
        :type parameters: list or iterable of eos.Parameter
        :param likelihood: The likelihood object created by nabu.
        :type results: nabu.Likelihood or descendant
        """
        description = {}
        description['version'] = eos.__version__
        description['type'] = 'NabuLikelihood'
        description['parameters'] = [{
            'name': p.name(),
            'min': p.min(),
            'max': p.max()
        } for p in parameters]

        os.makedirs(path, exist_ok=True)
        with open(os.path.join(path, 'description.yaml'), 'w') as description_file:
            yaml.dump(description, description_file, default_flow_style=False)

        likelihood.save(os.path.join(path, 'likelihood.nabu'))


class SampleMask:
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
            raise RuntimeError(f'Path {path} not pointing to a Mode file')

        self.type = 'Mask'
        self.observables = description['observables']

        f = os.path.join(path, 'mask.npy')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Weights file {f} does not exist or is not a file')
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
