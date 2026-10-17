# Copyright (c) 2024-2026 Danny van Dyk
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
from eos.data._common import ParameterDescription
from eos.deserializable import Deserializable
from eos.serializable import Serializable

import copy as _copy
import eos
import numpy as _np
import os


_jax_initialized = False

def _init_jax_single_threaded():
    # XLA sizes its CPU thread pool from the affinity at backend start-up; a large pool slows
    # down the single-point evaluations used by NabuLogPrior about eightfold on many-core hosts.
    global _jax_initialized
    if _jax_initialized:
        return
    import jax
    _jax_initialized = True
    if not hasattr(os, 'sched_setaffinity'):
        return
    affinity = os.sched_getaffinity(0)
    threads = set(os.listdir('/proc/self/task'))
    os.sched_setaffinity(0, {min(affinity)})
    try:
        jax.devices()
    finally:
        # the threads started by XLA inherit the restricted affinity, but keep their small pool once released
        for tid in ['0'] + sorted(set(os.listdir('/proc/self/task')) - threads):
            try:
                os.sched_setaffinity(int(tid), affinity)
            except ProcessLookupError:
                pass


@dataclass(kw_only=True)
class NabuLikelihoodDescription(Serializable, Deserializable):
    r"""Schema of the ``description.yaml`` written for a :class:`NabuLikelihood` object.

    This is the single source of truth for the metadata stored alongside a nabu-serialized likelihood,
    and it backs both the read path (:meth:`from_yaml_file`) and the write path (:meth:`to_yaml_file`).
    The ``type`` discriminator is validated on deserialization and is not part of the constructor.

    :param version: The version of EOS that wrote the data object.
    :type version: str
    :param parameters: The descriptions of the varied parameters.
    :type parameters: list[ParameterDescription]
    """
    version:str
    parameters:list[ParameterDescription]
    type:str = field(init=False, default='NabuLikelihood')

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`NabuLikelihoodDescription` from its on-disk keyword description.

        Validates the ``type`` discriminator and deserializes each entry of ``parameters`` into a
        :class:`ParameterDescription`.

        :raises ValueError: If the description does not identify a nabu-likelihood object.
        """
        _kwargs = _copy.deepcopy(kwargs)

        _type = _kwargs.pop('type', None)
        if _type != 'NabuLikelihood':
            raise ValueError(f'Expected a description of type \'NabuLikelihood\', got \'{_type}\'')

        if 'parameters' in _kwargs:
            _kwargs['parameters'] = [ParameterDescription.from_dict(**p) for p in _kwargs['parameters']]

        return Deserializable.make(cls, **_kwargs)

    def to_dict(self):
        """Serialize this description into the on-disk mapping written to ``description.yaml``.

        Emits the ``type`` discriminator, inverting :meth:`from_dict`.
        """
        return {
            'version':    self.version,
            'type':       self.type,
            'parameters': [asdict(p) for p in self.parameters],
        }


class NabuLikelihood:
    r"""Represents a likelihood serialized with the nabu package, stored on disk.

    Wraps a :class:`nabu.Likelihood` object together with the descriptions of the varied parameters.
    Instances are created either by reading an existing likelihood from disk (passing its ``path`` to
    the constructor) or by writing a new likelihood with :meth:`create`. Reading requires the optional
    nabu module. The likelihood can be used as a prior through :meth:`log_prior`.

    Reading a likelihood initializes the CPU backend of JAX, unless already initialized, with a
    thread pool sized for a single core; no thread of the process remains pinned to that core.
    This speeds up the evaluation of single points about eightfold on many-core hosts, but makes
    JAX computations on large batches in the same process (e.g. fitting a flow) slower.

    :ivar type: The type identifier of the data object, always ``'NabuLikelihood'``.
    :ivar varied_parameters: The descriptions (name, min, max) of the varied parameters.
    :ivar lookup_table: A mapping from each parameter name to its index in :attr:`varied_parameters`.
    :ivar likelihood: The underlying :class:`nabu.Likelihood` object.
    """

    def __init__(self, path):
        """ Read a nabu serialized likelihood from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError(f'Path {path} does not exist or is not a directory')

        description = NabuLikelihoodDescription.from_yaml_file(os.path.join(path, 'description.yaml'))

        self.type = description.type
        self.varied_parameters = [asdict(p) for p in description.parameters]
        self.lookup_table = { p.name: idx for idx, p in enumerate(description.parameters) }

        f = os.path.join(path, 'likelihood.nabu')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Nabu likelihood file {f} does not exist or is not a file')

        _init_jax_single_threaded()
        import nabu
        self.likelihood = nabu.Likelihood.load(f)


    def log_prior(self, parameters):
        """Return this likelihood as an external prior on the parameters listed in ``description.yaml``.

        :param parameters: The parameters to which the prior is bound.
        :type parameters: eos.Parameters
        :returns: The prior, varying the parameters in the order of :attr:`varied_parameters`.
        :rtype: eos.LogPrior
        """
        return eos.LogPrior.External(parameters, lambda p: NabuLogPrior(p, self))


    @staticmethod
    def make_log_prior(parameters, path):
        """Read the likelihood at ``path`` and return it as a prior; see :meth:`log_prior`."""
        return NabuLikelihood(path).log_prior(parameters)


    @staticmethod
    def create(path, parameters, likelihood):
        """ Write a new NabuLikelihood object to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param parameters: Parameter descriptions as a 1D array of shape (N, ).
        :type parameters: list or collections.abc.Iterable of eos.Parameter
        :param likelihood: The likelihood object created by nabu.
        :type likelihood: nabu.Likelihood or descendant
        """
        description = NabuLikelihoodDescription(
            version    = eos.__version__,
            parameters = [ParameterDescription(name=p.name(), min=p.min(), max=p.max()) for p in parameters],
        )

        os.makedirs(path, exist_ok=True)
        description.to_yaml_file(os.path.join(path, 'description.yaml'))

        likelihood.save(os.path.join(path, 'likelihood.nabu'))


class NabuLogPrior:
    r"""Provider of an external prior backed by a :class:`NabuLikelihood`, for use with
    :meth:`eos.LogPrior.External`.

    The prior density is the nabu likelihood's density in the physical basis. Inverse transform
    sampling uses the likelihood's map from and to the unit hypercube. Outside the likelihood's
    support, the prior evaluates to :math:`-\infty` and :meth:`compute_cdf` yields NaN.

    :param parameters: The parameters to which the prior is bound.
    :type parameters: eos.Parameters
    :param likelihood: The likelihood that provides the prior density.
    :type likelihood: NabuLikelihood
    """
    informative = True

    def __init__(self, parameters, likelihood):
        self._likelihood = likelihood.likelihood
        self.varied_parameters = [p['name'] for p in likelihood.varied_parameters]
        self._parameters = [parameters[n] for n in self.varied_parameters]

    def _values(self):
        return _np.array([p.evaluate() for p in self._parameters])

    def evaluate(self):
        return float(self._likelihood.log_prob(self._values()))

    def sample(self):
        u = _np.array([p.evaluate_generator() for p in self._parameters])
        for p, v in zip(self._parameters, self._likelihood.from_unit_cube(u)):
            p.set(float(v))

    def compute_cdf(self):
        for p, u in zip(self._parameters, self._likelihood.to_unit_cube(self._values())):
            p.set_generator(float(u))
