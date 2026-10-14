# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2021-2026 Danny van Dyk
# Copyright (c) 2026      Mark E Smith
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
import collections
import itertools
import numpy as _np
import os

from dataclasses import asdict
from ._framework import task, _check_varied_parameters_match, _staged


_model_comparison_checks = {}

def _model_comparison_check(name, description):
    """Register a named stability check, which receives an :class:`eos.data.model_comparison.PairwiseComparisonDescription`
    and returns whether the conclusion for that pair of posteriors is stable."""
    def _register(func):
        _model_comparison_checks[name] = (description, func)
        return func
    return _register


# Jeffreys' scale (1961): thresholds at half-integer powers of 10 in B, given on |ln B|
_BAYES_FACTOR_STRENGTHS = [
    (2.0 * _np.log(10.0), 'decisive'),
    (1.5 * _np.log(10.0), 'very strong'),
    (1.0 * _np.log(10.0), 'strong'),
    (0.5 * _np.log(10.0), 'substantial'),
]
_BAYES_FACTOR_NEGLIGIBLE = 'barely worth mentioning'

def _bayes_factor_strength(log_bayes_factor):
    for threshold, strength in _BAYES_FACTOR_STRENGTHS:
        if abs(log_bayes_factor) >= threshold:
            return strength
    return _BAYES_FACTOR_NEGLIGIBLE


def _bayes_factor_conclusion(log_bayes_factor):
    strength = _bayes_factor_strength(log_bayes_factor)
    if strength == _BAYES_FACTOR_NEGLIGIBLE:
        return (0, strength)
    return (int(_np.sign(log_bayes_factor)), strength)


@_model_comparison_check('uncertainty', 'The conclusion does not change within one standard deviation of the log Bayes factor.')
def _check_uncertainty(pair):
    # the conclusion is monotonic in ln B, so comparing the interval's endpoints suffices
    lower = pair.log_bayes_factor - pair.log_bayes_factor_uncertainty
    upper = pair.log_bayes_factor + pair.log_bayes_factor_uncertainty
    return _bayes_factor_conclusion(lower) == _bayes_factor_conclusion(upper)


@_model_comparison_check('prior-volume-adjustment', 'The conclusion does not change when the prior-volume adjustment is omitted.')
def _check_prior_volume_adjustment(pair):
    return _bayes_factor_conclusion(pair.log_bayes_factor) == _bayes_factor_conclusion(pair.unadjusted_log_bayes_factor)


def _uniform_log_volume_and_box(description):
    from ..analysis_file_description import TransformPriorDescription

    if isinstance(description, TransformPriorDescription):
        transform = _np.asarray(description.transform, dtype=float)
        shift     = _np.asarray(description.shift, dtype=float)
        lower     = _np.asarray(description.min, dtype=float)
        upper     = _np.asarray(description.max, dtype=float)
        # the support is the image of the box [min, max] under x = transform . y + shift
        corners = _np.array([transform @ y + shift for y in itertools.product(*zip(lower, upper))])
        log_volume = float(_np.sum(_np.log(upper - lower)) + _np.linalg.slogdet(transform)[1])
        return log_volume, dict(zip(description.parameters, zip(corners.min(axis=0), corners.max(axis=0))))

    lower, upper = float(description.min), float(description.max)
    return float(_np.log(upper - lower)), { description.parameter: (lower, upper) }


def _log_prior_volume_adjustments(analysis_file, posteriors, results):
    from ..analysis_file_description import TransformPriorDescription, UniformPriorDescription

    def _parameters(description):
        if hasattr(description, 'parameter'):
            return [description.parameter]
        return list(getattr(description, 'parameters', []))

    descriptions = {
        p: [d for prior in analysis_file.posteriors[p].prior for d in analysis_file.priors[prior].descriptions]
        for p in posteriors
    }
    shared = set.intersection(*(set(results[p].lookup_table) for p in posteriors))

    # group the parameters into blocks that are linked by a common prior description in any posterior
    parent = {}
    def _find(name):
        parent.setdefault(name, name)
        while parent[name] != name:
            name = parent[name]
        return name
    for p in posteriors:
        for d in descriptions[p]:
            names = _parameters(d)
            for name in names:
                parent[_find(name)] = _find(names[0])
    blocks = {}
    for name in list(parent) + sorted(shared):
        blocks.setdefault(_find(name), set()).add(name)

    adjustments = { p: 0.0 for p in posteriors }
    log_prior_volumes = []
    for block in sorted(blocks.values(), key=min):
        if not block & shared:
            continue
        involving = { p: [d for d in descriptions[p] if block & set(_parameters(d))] for p in posteriors }
        if all(involving[p] == involving[posteriors[0]] for p in posteriors):
            continue

        if not block <= shared:
            raise ValueError(f'The priors of shared parameters {sorted(block & shared)} differ and also involve the model-specific parameters {sorted(block - shared)}; cannot adjust for the prior volume')
        if not all(isinstance(d, (UniformPriorDescription, TransformPriorDescription)) for p in posteriors for d in involving[p]):
            raise ValueError(f'Shared parameters {sorted(block)} have differing priors that are not all uniform or transform priors; cannot adjust for the prior volume')

        log_volumes, boxes = {}, {}
        for p in posteriors:
            log_volumes[p], boxes[p] = 0.0, {}
            for d in involving[p]:
                log_volume, box = _uniform_log_volume_and_box(d)
                log_volumes[p] += log_volume
                boxes[p].update(box)
            if sorted(boxes[p]) != sorted(block) or len(boxes[p]) != sum(len(_parameters(d)) for d in involving[p]):
                raise ValueError(f'The priors of posterior \'{p}\' do not cover the shared parameters {sorted(block)} exactly once; cannot adjust for the prior volume')

        log_reference_volume = float(sum(
            _np.log(max(boxes[p][name][1] for p in posteriors) - min(boxes[p][name][0] for p in posteriors))
            for name in block
        ))
        for p in posteriors:
            adjustments[p] += log_volumes[p] - log_reference_volume
        log_prior_volumes.append({ 'parameters': sorted(block), 'log_reference_volume': log_reference_volume, 'log_volumes': log_volumes })

    return adjustments, log_prior_volumes


@task('model-comparison', 'model-comparison/{group}')
def model_comparison(analysis_file:str, posteriors:list, base_directory:str='./', group:str='default'):
    r"""
    Compares a group of named posteriors that share the same likelihood by means of their Bayesian evidences.

    Bayesian model comparison ranks models by their evidences :math:`Z`, i.e., the likelihood averaged over
    the prior. The posteriors to be compared must share the same likelihood, i.e., the same list of named
    likelihoods; they may differ in their priors, global options, and fixed parameters, which define the models.
    The evidences are read from the nested-sampling results in EOS_BASE_DIRECTORY/data/POSTERIOR/nested,
    as produced by the ``sample-nested`` task.

    Each evidence includes an Occam factor for the prior volume of every varied parameter. For parameters that
    distinguish the models this is intended. For parameters varied by all posteriors whose priors differ,
    e.g. because a range was widened for sampling, it is an artefact. These shared parameters are grouped into
    blocks linked by common prior descriptions. For each block whose priors differ, all priors must be uniform or
    transform priors, and each evidence is corrected to a common prior on the smallest axis-aligned box that
    contains all of their supports. The log Bayes factors do not depend on the choice of that common box. This
    correction assumes that no posterior is truncated by its own prior on a shared parameter.

    The log Bayes factor :math:`\ln B` of each pair of posteriors is graded on Jeffreys' scale:

    =================================  =======================
    lower bound on :math:`|\ln B|`     strength
    =================================  =======================
    :math:`0`                          barely worth mentioning
    :math:`\ln 10^{1/2} \approx 1.15`  substantial
    :math:`\ln 10 \approx 2.30`        strong
    :math:`\ln 10^{3/2} \approx 3.45`  very strong
    :math:`\ln 10^{2} \approx 4.61`    decisive
    =================================  =======================

    The conclusion for each pair, i.e., the favoured posterior and the strength, is subjected to the stability checks

    - ``uncertainty``: the conclusion does not change within one standard deviation of :math:`\ln B`;
    - ``prior-volume-adjustment``: the conclusion does not change when the prior-volume correction is omitted.

    Failed checks are reported as warnings and recorded in the output.
    The output will be stored in EOS_BASE_DIRECTORY/model-comparison/GROUP; see :class:`eos.data.ModelComparison` for its content.

    :param analysis_file: The name of the analysis file that describes the named posteriors, or an object of class `eos.AnalysisFile`.
    :type analysis_file: str or `eos.AnalysisFile`
    :param posteriors: The names of the posteriors to compare.
    :type posteriors: list[str]
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    :param group: The name of the group of posteriors, used to label the output. Defaults to 'default'.
    :type group: str, optional
    :returns: The model comparison.
    :rtype: eos.data.ModelComparison
    """
    if len(posteriors) < 2 or len(set(posteriors)) != len(posteriors):
        raise ValueError('A model comparison requires at least two distinct posteriors')
    for p in posteriors:
        if p not in analysis_file.posteriors:
            raise ValueError(f'Unknown posterior \'{p}\'')

    first = analysis_file.posteriors[posteriors[0]]
    for p in posteriors[1:]:
        other = analysis_file.posteriors[p]
        a, b = collections.Counter(first.likelihood), collections.Counter(other.likelihood)
        if a != b:
            raise ValueError(f'Posteriors \'{first.name}\' and \'{p}\' do not share the same likelihood; differing likelihoods: {sorted(((a - b) + (b - a)).elements())}')
        if other.global_options != first.global_options or other.fixed_parameters != first.fixed_parameters:
            eos.info(f'Posteriors \'{first.name}\' and \'{p}\' differ in their global options or fixed parameters')

    results = {}
    for p in posteriors:
        path = os.path.join(base_directory, 'data', p, 'nested')
        if not os.path.isdir(path):
            raise RuntimeError(f'No nested-sampling results found for posterior \'{p}\' in \'{path}\'; run the \'sample-nested\' task first')
        results[p] = eos.data.DynestyResults(path)
        _check_varied_parameters_match(analysis_file.analysis(p, base_directory=base_directory), results[p])

    adjustments, log_prior_volumes = _log_prior_volume_adjustments(analysis_file, posteriors, results)
    log_z    = { p: float(results[p].results.logz[-1])    for p in posteriors }
    sigma    = { p: float(results[p].results.logzerr[-1]) for p in posteriors }
    adjusted = { p: log_z[p] + adjustments[p]             for p in posteriors }
    reference = max(posteriors, key=adjusted.get)

    entries = [{
        'posterior':                    p,
        'log_evidence':                 log_z[p],
        'log_evidence_uncertainty':     sigma[p],
        'log_prior_volume_adjustment':  adjustments[p],
        'adjusted_log_evidence':        adjusted[p],
        'log_bayes_factor':             adjusted[p] - adjusted[reference],
        'log_bayes_factor_uncertainty': float(_np.hypot(sigma[p], sigma[reference])) if p != reference else 0.0,
        'strength':                     _bayes_factor_strength(adjusted[p] - adjusted[reference]) if p != reference else 'reference',
    } for p in posteriors]

    from eos.data.model_comparison import PairwiseComparisonDescription

    pairs = []
    for a, b in itertools.combinations(posteriors, 2):
        a, b = (a, b) if adjusted[a] >= adjusted[b] else (b, a)
        pairs.append(PairwiseComparisonDescription(
            first                        = a,
            second                       = b,
            log_bayes_factor             = adjusted[a] - adjusted[b],
            log_bayes_factor_uncertainty = float(_np.hypot(sigma[a], sigma[b])),
            unadjusted_log_bayes_factor  = log_z[a] - log_z[b],
            strength                     = _bayes_factor_strength(adjusted[a] - adjusted[b]),
        ))

    checks = []
    for name, (description, check) in _model_comparison_checks.items():
        failed_pairs = []
        for pair in pairs:
            if check(pair):
                continue
            pair.failed_checks.append(name)
            failed_pairs.append([pair.first, pair.second])
            eos.warn(f'The conclusion for \'{pair.first}\' over \'{pair.second}\' fails the \'{name}\' check')
        checks.append({ 'name': name, 'status': 'failed' if failed_pairs else 'passed', 'description': description, 'failed_pairs': failed_pairs })

    comparisons = [asdict(pair) for pair in pairs]

    for e in entries:
        eos.info(f'{e["posterior"]}: adjusted ln Z = {e["adjusted_log_evidence"]:.2f} +/- {e["log_evidence_uncertainty"]:.2f}, '
                 f'ln B vs. \'{reference}\' = {e["log_bayes_factor"]:.2f} +/- {e["log_bayes_factor_uncertainty"]:.2f} ({e["strength"]})')

    path = _staged(os.path.join(base_directory, 'model-comparison', group))
    eos.data.ModelComparison.create(path, group, reference, entries, comparisons, checks, log_prior_volumes)
    return eos.data.ModelComparison(path)
