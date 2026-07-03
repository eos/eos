# Copyright (c) 2020-2026 Danny van Dyk
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

from dataclasses import dataclass, field
from eos.data._common import GaussianComponentDescription
from eos.deserializable import Deserializable
from eos.serializable import Serializable

import copy as _copy
import eos
import os
import numpy as _np
from scipy.special import erf
from scipy.linalg import block_diag


@dataclass(kw_only=True)
class MixtureDensityDescription(Serializable, Deserializable):
    r"""Schema of the ``description.yaml`` written for a :class:`MixtureDensity`.

    This is the single source of truth for the metadata stored for a mixture density, and it backs both
    the read path (:meth:`from_yaml_file`) and the write path (:meth:`to_yaml_file`). The ``type``
    discriminator is validated on deserialization and is not part of the constructor.

    :param version: The version of EOS that wrote the data object.
    :type version: str
    :param components: The Gaussian components of the mixture.
    :type components: list[GaussianComponentDescription]
    :param weights: The component weights of the mixture.
    :type weights: list
    :param varied_parameters: The qualified names of the varied parameters, or ``None`` if not stored.
    :type varied_parameters: list | None
    :param test_statistics: The precomputed test statistics, or ``None`` if not stored.
    :type test_statistics: dict | None
    """
    version:str
    components:list[GaussianComponentDescription]
    weights:list
    varied_parameters:list|None = field(default=None)
    test_statistics:dict|None = field(default=None)
    type:str = field(init=False, default='MixtureDensity')

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`MixtureDensityDescription` from its on-disk keyword description.

        Validates the ``type`` discriminator and deserializes each Gaussian component (which carries a
        ``type: 'gauss'`` marker on disk).

        :raises ValueError: If the description does not identify a mixture density, or a component has
            an unsupported type.
        """
        _kwargs = _copy.deepcopy(kwargs)

        _type = _kwargs.pop('type', None)
        if _type != 'MixtureDensity':
            raise ValueError(f'Expected a description of type \'MixtureDensity\', got \'{_type}\'')

        if 'components' in _kwargs:
            components = []
            for c in _kwargs['components']:
                _c = dict(c)
                ctype = _c.pop('type', 'gauss')
                if ctype != 'gauss':
                    raise ValueError(f'Unsupported MixtureDensity component type \'{ctype}\'; only \'gauss\' is supported')
                components.append(GaussianComponentDescription.from_dict(**_c))
            _kwargs['components'] = components

        return Deserializable.make(cls, **_kwargs)

    def to_dict(self):
        """Serialize this description into the on-disk mapping written to ``description.yaml``.

        Re-emits the ``type: 'gauss'`` marker on each component, inverting :meth:`from_dict`.
        """
        result = {
            'version':    self.version,
            'type':       self.type,
            'components': [{'type': 'gauss', 'mu': c.mu, 'sigma': c.sigma} for c in self.components],
            'weights':    self.weights,
        }
        if self.varied_parameters is not None:
            result['varied_parameters'] = self.varied_parameters
        result['test_statistics'] = self.test_statistics if self.test_statistics is not None else {'sigma': [], 'densities': []}
        return result


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

        description = MixtureDensityDescription.from_yaml_file(f)

        self.type = description.type
        self.components = [{'type': 'gauss', 'mu': c.mu, 'sigma': c.sigma} for c in description.components]
        self.weights    = description.weights
        self.varied_parameters = description.varied_parameters
        self.test_statistics = description.test_statistics

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

        components = []
        for c in density.components:
            if type(c) is pypmc.density.gauss.Gauss:
                components.append(GaussianComponentDescription(mu=c.mu.tolist(), sigma=c.sigma.tolist()))
            else:
                raise RuntimeError(f'Unsupported type of MixtureDensity component: {type(c)}')

        varied_parameters = list(varied_parameters) if varied_parameters is not None else None

        # The test statistics defaults to two empty lists
        test_statistics = { "sigma": [], "densities": [] }
        if sigma_test_stat is not None and samples is not None and weights is not None:
            sigma_test_stat = _np.array(sigma_test_stat)
            samples_pdf = list(map(lambda x: -density.evaluate(x), samples))
            ind = _np.argsort(samples_pdf)
            samples_pdf_sorted = _np.array(samples_pdf)[ind]
            sorted_weights = weights[ind]
            cumulant = 1.0 * sorted_weights.cumsum() / sorted_weights.sum()
            percents = erf(sigma_test_stat/_np.sqrt(2))
            weighted_percentile = _np.interp(percents, cumulant, samples_pdf_sorted)

            test_statistics = {
                "sigma": sigma_test_stat.tolist(),
                "densities": weighted_percentile.tolist()
            }

        description = MixtureDensityDescription(
            version           = eos.__version__,
            components        = components,
            weights           = density.weights.tolist(),
            varied_parameters = varied_parameters,
            test_statistics   = test_statistics,
        )

        os.makedirs(path, exist_ok=True)
        description.to_yaml_file(os.path.join(path, 'description.yaml'))

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
