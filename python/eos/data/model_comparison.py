# Copyright (c) 2026 Danny van Dyk
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
from eos.deserializable import Deserializable
from eos.serializable import Serializable

import copy as _copy
import eos
import os


@dataclass(kw_only=True)
class ModelComparisonEntryDescription(Serializable, Deserializable):
    r"""Describes the evidence of a single posterior within a model comparison.

    :param posterior: The name of the posterior.
    :type posterior: str
    :param log_evidence: The natural logarithm of the evidence, as estimated by nested sampling.
    :type log_evidence: float
    :param log_evidence_uncertainty: The uncertainty of :attr:`log_evidence`.
    :type log_evidence_uncertainty: float
    :param log_prior_volume_adjustment: The correction for differing prior ranges of shared parameters.
    :type log_prior_volume_adjustment: float
    :param adjusted_log_evidence: The sum of :attr:`log_evidence` and :attr:`log_prior_volume_adjustment`.
    :type adjusted_log_evidence: float
    :param log_bayes_factor: The adjusted log Bayes factor of this posterior with respect to the reference.
    :type log_bayes_factor: float
    :param log_bayes_factor_uncertainty: The uncertainty of :attr:`log_bayes_factor`.
    :type log_bayes_factor_uncertainty: float
    :param strength: The strength of the evidence against this posterior, relative to the reference.
    :type strength: str
    """
    posterior:str
    log_evidence:float
    log_evidence_uncertainty:float
    log_prior_volume_adjustment:float
    adjusted_log_evidence:float
    log_bayes_factor:float
    log_bayes_factor_uncertainty:float
    strength:str


@dataclass(kw_only=True)
class PairwiseComparisonDescription(Serializable, Deserializable):
    r"""Describes the comparison of two posteriors, ordered such that ``first`` has the larger adjusted evidence.

    :param first: The name of the favoured posterior.
    :type first: str
    :param second: The name of the disfavoured posterior.
    :type second: str
    :param log_bayes_factor: The adjusted log Bayes factor of ``first`` over ``second``.
    :type log_bayes_factor: float
    :param log_bayes_factor_uncertainty: The uncertainty of :attr:`log_bayes_factor`.
    :type log_bayes_factor_uncertainty: float
    :param unadjusted_log_bayes_factor: The log Bayes factor without the prior-volume adjustment.
    :type unadjusted_log_bayes_factor: float
    :param strength: The strength of the evidence in favour of ``first``.
    :type strength: str
    :param failed_checks: The names of the stability checks that this pair fails.
    :type failed_checks: list[str]
    """
    first:str
    second:str
    log_bayes_factor:float
    log_bayes_factor_uncertainty:float
    unadjusted_log_bayes_factor:float
    strength:str
    failed_checks:list[str] = field(default_factory=list)


@dataclass(kw_only=True)
class ModelComparisonCheckDescription(Serializable, Deserializable):
    r"""Describes the outcome of a named stability check applied to all pairs of posteriors.

    :param name: The name of the check.
    :type name: str
    :param status: Either ``'passed'`` or ``'failed'``.
    :type status: str
    :param description: A human-readable description of what the check tests.
    :type description: str
    :param failed_pairs: The pairs of posteriors, as ``[first, second]``, for which the check fails.
    :type failed_pairs: list[list[str]]
    """
    name:str
    status:str
    description:str
    failed_pairs:list[list[str]] = field(default_factory=list)


@dataclass(kw_only=True)
class ModelComparisonDescription(Serializable, Deserializable):
    r"""Schema of the ``description.yaml`` written for a :class:`ModelComparison`.

    :param version: The version of EOS that wrote the data object.
    :type version: str
    :param group: The name of the group of posteriors that is compared.
    :type group: str
    :param reference: The name of the posterior with the largest adjusted evidence.
    :type reference: str
    :param posteriors: The per-posterior evidences, in the order in which the posteriors were given.
    :type posteriors: list[ModelComparisonEntryDescription]
    :param comparisons: The comparisons of all pairs of posteriors.
    :type comparisons: list[PairwiseComparisonDescription]
    :param checks: The outcomes of the stability checks.
    :type checks: list[ModelComparisonCheckDescription]
    """
    version:str
    group:str
    reference:str
    posteriors:list[ModelComparisonEntryDescription]
    comparisons:list[PairwiseComparisonDescription]
    checks:list[ModelComparisonCheckDescription]
    type:str = field(init=False, default='ModelComparison')

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`ModelComparisonDescription` from its on-disk keyword description.

        :raises ValueError: If the description does not identify a model comparison.
        """
        _kwargs = _copy.deepcopy(kwargs)

        _type = _kwargs.pop('type', None)
        if _type != 'ModelComparison':
            raise ValueError(f'Expected a description of type \'ModelComparison\', got \'{_type}\'')

        for key, schema in (('posteriors', ModelComparisonEntryDescription), ('comparisons', PairwiseComparisonDescription), ('checks', ModelComparisonCheckDescription)):
            if key in _kwargs:
                _kwargs[key] = [schema.from_dict(**e) for e in _kwargs[key]]

        return Deserializable.make(cls, **_kwargs)

    def to_dict(self):
        """Serialize this description into the on-disk mapping written to ``description.yaml``."""
        return {
            'version':     self.version,
            'type':        self.type,
            'group':       self.group,
            'reference':   self.reference,
            'posteriors':  [asdict(e) for e in self.posteriors],
            'comparisons': [asdict(c) for c in self.comparisons],
            'checks':      [asdict(c) for c in self.checks],
        }


class ModelComparison:
    r"""Represents the Bayesian comparison of a group of posteriors, stored on disk.

    :ivar type: The type identifier of the data object, always ``'ModelComparison'``.
    :ivar group: The name of the group of posteriors that is compared.
    :ivar reference: The name of the posterior with the largest adjusted evidence.
    :ivar posteriors: The per-posterior evidences and log Bayes factors relative to ``reference``.
    :ivar comparisons: The comparisons of all pairs of posteriors.
    :ivar checks: The named stability checks and their statuses.
    """

    def __init__(self, path):
        """ Read a model comparison from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError(f'Path {path} does not exist or is not a directory')

        description = ModelComparisonDescription.from_yaml_file(os.path.join(path, 'description.yaml'))

        self.type        = description.type
        self.group       = description.group
        self.reference   = description.reference
        self.posteriors  = [asdict(e) for e in description.posteriors]
        self.comparisons = [asdict(c) for c in description.comparisons]
        self.checks      = [asdict(c) for c in description.checks]

    @property
    def stable(self):
        """Whether all stability checks passed.

        :rtype: bool
        """
        return all(c['status'] == 'passed' for c in self.checks)

    @staticmethod
    def create(path, group, reference, posteriors, comparisons, checks):
        """ Write a new ModelComparison object to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param group: The name of the group of posteriors that is compared.
        :type group: str
        :param reference: The name of the posterior with the largest adjusted evidence.
        :type reference: str
        :param posteriors: The per-posterior entries, as mappings of the fields of ``ModelComparisonEntryDescription``.
        :type posteriors: list of dict
        :param comparisons: The pairwise entries, as mappings of the fields of ``PairwiseComparisonDescription``.
        :type comparisons: list of dict
        :param checks: The check outcomes, as mappings of the fields of ``ModelComparisonCheckDescription``.
        :type checks: list of dict
        """
        description = ModelComparisonDescription(
            version     = eos.__version__,
            group       = group,
            reference   = reference,
            posteriors  = [ModelComparisonEntryDescription.from_dict(**e) for e in posteriors],
            comparisons = [PairwiseComparisonDescription.from_dict(**c) for c in comparisons],
            checks      = [ModelComparisonCheckDescription.from_dict(**c) for c in checks],
        )

        os.makedirs(path, exist_ok=True)
        description.to_yaml_file(os.path.join(path, 'description.yaml'))

    def _repr_html_(self):
        rows = '\n'.join(
            f'<tr><td><tt>{e["posterior"]}</tt></td><td>{e["adjusted_log_evidence"]:.2f} &pm; {e["log_evidence_uncertainty"]:.2f}</td>'
            f'<td>{e["log_bayes_factor"]:.2f} &pm; {e["log_bayes_factor_uncertainty"]:.2f}</td><td>{e["strength"]}</td></tr>'
            for e in self.posteriors
        )
        checks = '\n'.join(f'<tr><td>{c["name"]}</td><td colspan="3">{c["status"]}</td></tr>' for c in self.checks)
        return rf'''
        <table>
            <thead>
                <tr><th colspan="4">model comparison <tt>{self.group}</tt></th></tr>
                <tr><th>posterior</th><th>adjusted $\ln Z$</th><th>$\ln B$ vs. <tt>{self.reference}</tt></th><th>strength</th></tr>
            </thead>
            <tbody>
{rows}
{checks}
            </tbody>
        </table>'''
