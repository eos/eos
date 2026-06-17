# Copyright (c) 2020 Danny van Dyk
# Copyright (c) 2024 Méril Reboud
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
from scipy.special import erf
from scipy.linalg import block_diag


class MixtureDensity:
    r"""Represents a mixture density (e.g. a PMC proposal or a fit to a posterior) stored on disk.

    Stores the components and weights of a mixture of (Gaussian) densities, optionally together with the
    qualified names of the varied parameters and precomputed test statistics. Instances are created
    either by reading an existing density from disk (passing its ``path`` to the constructor) or by
    writing a new density with :meth:`create`. Use :meth:`density` to obtain the corresponding
    :class:`pypmc.density.mixture.MixtureDensity`.

    :ivar type: The type identifier of the data object, always ``'MixtureDensity'``.
    :ivar components: The descriptions of the mixture components.
    :ivar weights: The component weights of the mixture.
    :ivar varied_parameters: The qualified names of the varied parameters, or ``None`` if not stored.
    :ivar test_statistics: The precomputed test statistics, or ``None`` if not stored.
    """

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
        """Construct the corresponding PyPMC mixture density.

        :returns: A Gaussian mixture density built from the stored components and weights.
        :rtype: pypmc.density.mixture.MixtureDensity
        :raises ImportError: If the optional PyPMC module is not installed.
        """
        try:
            import pypmc
        except ImportError as e:
            raise ImportError('eos.data.MixtureDensity.density requires the PyPMC python module, which can be installed from PyPI.') from e
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
        try:
            import pypmc
        except ImportError as e:
            raise ImportError('eos.data.MixtureDensity.create requires the PyPMC python module, which can be installed from PyPI.') from e
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
        try:
            import pypmc
        except ImportError as e:
            raise ImportError('eos.data.MixtureDensity._cartesian_product requires the PyPMC python module, which can be installed from PyPI.') from e
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
