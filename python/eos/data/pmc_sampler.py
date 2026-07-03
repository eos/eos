# Copyright (c) 2021-2026 Danny van Dyk
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

from dataclasses import dataclass, field, asdict
from eos.data._common import GaussianComponentDescription, ParameterDescription
from eos.deserializable import Deserializable
from eos.serializable import Serializable

import copy as _copy
import eos
import os
import numpy as _np
from scipy.special import erf


@dataclass(kw_only=True)
class ProposalDescription(Serializable, Deserializable):
    r"""Describes the Gaussian-mixture proposal density of a :class:`PMCSampler`.

    :param components: The Gaussian components of the proposal mixture.
    :type components: list[GaussianComponentDescription]
    :param weights: The component weights of the proposal mixture.
    :type weights: list
    """
    components:list[GaussianComponentDescription]
    weights:list

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`ProposalDescription`, deserializing its Gaussian components."""
        _kwargs = _copy.deepcopy(kwargs)
        if 'components' in _kwargs:
            _kwargs['components'] = [GaussianComponentDescription.from_dict(**c) for c in _kwargs['components']]
        return Deserializable.make(cls, **_kwargs)


@dataclass(kw_only=True)
class PMCSamplerDescription(Serializable, Deserializable):
    r"""Schema of the ``description.yaml`` written for a :class:`PMCSampler`.

    This is the single source of truth for the metadata stored alongside a PMC proposal, and it backs
    both the read path (:meth:`from_yaml_file`) and the write path (:meth:`to_yaml_file`). The ``type``
    discriminator is validated on deserialization and is not part of the constructor.

    The test statistics are written under the ``test_statistics`` key. For backward compatibility the
    legacy misspelled key ``'test statistics'`` (with a space) is also accepted on deserialization
    (with a warning), since older files were written with it.

    :param version: The version of EOS that wrote the data object.
    :type version: str
    :param parameters: The descriptions of the varied parameters.
    :type parameters: list[ParameterDescription]
    :param proposal: The Gaussian-mixture proposal density.
    :type proposal: ProposalDescription
    :param test_statistics: The precomputed test statistics.
    :type test_statistics: dict
    """
    version:str
    parameters:list[ParameterDescription]
    proposal:ProposalDescription
    test_statistics:dict = field(default_factory=lambda: {'sigma': [], 'densities': []})
    type:str = field(init=False, default='PMCSampler')

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`PMCSamplerDescription` from its on-disk keyword description.

        Validates the ``type`` discriminator, deserializes the parameters and the proposal, and
        normalizes the legacy ``'test statistics'`` key to ``test_statistics``.

        :raises ValueError: If the description does not identify a PMC sampler.
        """
        _kwargs = _copy.deepcopy(kwargs)

        _type = _kwargs.pop('type', None)
        if _type != 'PMCSampler':
            raise ValueError(f'Expected a description of type \'PMCSampler\', got \'{_type}\'')

        if 'test statistics' in _kwargs:
            eos.warn('The \'test statistics\' key is deprecated; it will be written as \'test_statistics\' when this data set is re-saved')
            legacy = _kwargs.pop('test statistics')
            _kwargs.setdefault('test_statistics', legacy)

        if 'parameters' in _kwargs:
            _kwargs['parameters'] = [ParameterDescription.from_dict(**p) for p in _kwargs['parameters']]

        if 'proposal' in _kwargs:
            _kwargs['proposal'] = ProposalDescription.from_dict(**_kwargs['proposal'])

        return Deserializable.make(cls, **_kwargs)

    def to_dict(self):
        """Serialize this description into the on-disk mapping written to ``description.yaml``.

        Emits the ``type`` discriminator and the (normalized) ``test_statistics`` key, inverting
        :meth:`from_dict`.
        """
        return {
            'version':         self.version,
            'type':            self.type,
            'parameters':      [asdict(p) for p in self.parameters],
            'proposal':        self.proposal.to_dict(),
            'test_statistics': self.test_statistics,
        }


class PMCSampler:
    r"""Represents the proposal density of a population Monte Carlo (PMC) sampler stored on disk.

    Stores the Gaussian components and weights of the PMC proposal mixture together with the descriptions
    of the varied parameters. Instances are created either by reading an existing sampler from disk
    (passing its ``path`` to the constructor) or by writing a new sampler with :meth:`create`. Use
    :meth:`density` to obtain the corresponding :class:`pypmc.density.mixture.MixtureDensity`. Reading
    requires the optional PyPMC module.

    :ivar type: The type identifier of the data object, always ``'PMCSampler'``.
    :ivar varied_parameters: The descriptions (name, min, max) of the varied parameters.
    :ivar lookup_table: A mapping from each parameter name to its index in :attr:`varied_parameters`.
    :ivar components: The Gaussian components of the proposal mixture.
    :ivar component_weights: The weights of the proposal mixture components.
    :ivar test_statistics: The precomputed test statistics.
    """

    def __init__(self, path):
        """ Read a PMCSampler object from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError(f'Path {path} does not exist or is not a directory')

        description = PMCSamplerDescription.from_yaml_file(os.path.join(path, 'description.yaml'))

        try:
            import pypmc
        except ImportError as e:
            raise ImportError('eos.data.PMCSampler requires the PyPMC python module, which can be installed from PyPI.') from e

        self.type = description.type
        self.varied_parameters = [asdict(p) for p in description.parameters]
        self.lookup_table = { p.name: idx for idx, p in enumerate(description.parameters) }
        self.components        = _np.array([pypmc.density.gauss.Gauss(_np.array(c.mu), _np.array(c.sigma)) for c in description.proposal.components])
        self.component_weights = _np.array(description.proposal.weights)
        self.test_statistics   = description.test_statistics


    def density(self):
        """Construct the corresponding PyPMC mixture density.

        :returns: A Gaussian mixture density built from the stored proposal components and weights.
        :rtype: pypmc.density.mixture.MixtureDensity
        :raises ImportError: If the optional PyPMC module is not installed.
        """
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
        :param weights: Weights on a linear scale as a 1D array of shape (N, ). Needed to generate the test statistic.
        :type weights: 1D numpy array, optional
        """
        # Don't write components that have a 0 weight
        purged_components = []
        purged_weights = []
        for comp, w in zip(proposal.components, proposal.weights.tolist()):
            if w != 0:
                purged_components.append(GaussianComponentDescription(mu=comp.mu.tolist(), sigma=comp.sigma.tolist()))
                purged_weights.append(w)
        proposal_description = ProposalDescription(components=purged_components, weights=purged_weights)

        # The test statistics default to two empty lists
        test_statistics = { "sigma": [], "densities": [] }
        if sigma_test_stat is not None and samples is not None and weights is not None:
            sigma_test_stat = _np.array(sigma_test_stat)
            samplesPDF = list(map(lambda x: -2.0 * _np.log(PMCSampler._evaluate_mixture_pdf(proposal, x)), samples))
            ind = _np.argsort(samplesPDF)
            sorted_samplesPDF = _np.array(samplesPDF)[ind]
            sorted_weights = weights.flatten()[ind]
            cumulant = 1.*sorted_weights.cumsum()/sorted_weights.sum()
            percents = erf(sigma_test_stat/_np.sqrt(2))
            weighted_percentile = _np.interp(percents, cumulant, sorted_samplesPDF)

            test_statistics = {
                "sigma": sigma_test_stat.tolist(),
                "densities": weighted_percentile.tolist()
            }

        description = PMCSamplerDescription(
            version         = eos.__version__,
            parameters      = [ParameterDescription(name=p.name(), min=p.min(), max=p.max()) for p in parameters],
            proposal        = proposal_description,
            test_statistics = test_statistics,
        )

        os.makedirs(path, exist_ok=True)
        description.to_yaml_file(os.path.join(path, 'description.yaml'))
