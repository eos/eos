# Copyright (c) 2021-2022 Danny van Dyk
# Copyright (c) 2021-2024 Méril Reboud
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

        try:
            import pypmc
        except ImportError as e:
            raise ImportError('eos.data.PMCSampler requires the PyPMC python module, which can be installed from PyPI.') from e
        self.type = 'PMCSampler'
        self.varied_parameters = description['parameters']
        self.lookup_table = { item['name']: idx for idx, item in enumerate(self.varied_parameters) }
        self.components        = _np.array([pypmc.density.gauss.Gauss(_np.array(c['mu']), _np.array(c['sigma'])) for c in description['proposal']['components']])
        self.component_weights = _np.array(description['proposal']['weights'])


    def density(self):
        """ Return a pypmc.density.MixtureDensity. """
        try:
            import pypmc
        except ImportError as e:
            raise ImportError('eos.data.PMCSampler.density requires the PyPMC python module, which can be installed from PyPI.') from e
        return pypmc.density.mixture.MixtureDensity(self.components, self.component_weights)


    @staticmethod
    def _evaluate_mixture_pdf(density, x):
        """ Internal function that evaluates the PDF of the mixture density at point x."""
        def evaluate_component(comp, x):
            # To limit numerical instabilities, the covariance matrix is divided by the median of its diagonal entries
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
