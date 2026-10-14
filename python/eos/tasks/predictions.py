# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2020-2026 Danny van Dyk
# Copyright (c) 2023      Méril Reboud
# Copyright (c) 2024-2025 Matthew Kirk
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
import numpy as _np
import os
import scipy
import yaml

from dataclasses import dataclass
from ._framework import task, _check_path_component, _check_varied_parameters_match, _mask_output_fields, _staged


# Predict observables
@task('predict-observables', 'data/{posterior}/pred-{prediction}{mask}', output_fields=_mask_output_fields)
def predict_observables(analysis_file:str, posterior:str, prediction:str, base_directory:str='./', begin:int=0, end:int=None, mask_name:str=None):
    '''
    Predicts a set of observables based on previously obtained importance samples.

    The input files are expected in EOS_BASE_DIRECTORY/data/POSTERIOR/samples.
    The output files will be stored in EOS_BASE_DIRECTORY/data/POSTERIOR/pred-PREDICTION.

    :param analysis_file: The name of the analysis file that describes the named posterior, or an object of class `eos.AnalysisFile`.
    :type analysis_file: str or `eos.AnalysisFile`
    :param posterior: The name of the posterior.
    :type posterior: str
    :param prediction: The name of the set of observables to predict.
    :type prediction: str
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    :param begin: The index of the first sample to use for the predictions. Defaults to 0.
    :type begin: int
    :param end: The index beyond the last sample to use for the predictions. Defaults to None.
    :type begin: int
    :param mask_name: The label of the mask to apply to the observables. Defaults to None.
    :type mask_name: str, optional
    '''
    _analysis      = analysis_file.analysis(posterior, base_directory=base_directory)
    _parameters    = _analysis.parameters
    cache          = eos.ObservableCache(_parameters)
    observables    = analysis_file.observables(posterior, prediction, _parameters)
    observable_ids = [cache.add(o) for o in observables]

    if mask_name is not None and (begin != 0 or end is not None):
        raise ValueError('The arguments mask-name and begin or end are mutually exclusive')
    if mask_name is not None:
        _check_path_component('mask_name', mask_name)
        mask = eos.data.SampleMask(os.path.join(base_directory, 'data', posterior, f'mask-{mask_name}')).mask

    data = eos.data.ImportanceSamples(os.path.join(base_directory, 'data', posterior, 'samples'))
    _check_varied_parameters_match(_analysis, data)

    try:
        from tqdm.auto import tqdm
        progressbar = tqdm
    except ImportError:
        progressbar = lambda x: x

    parameters = [_parameters[p['name']] for p in data.varied_parameters]
    observable_samples = []
    nsamples = len(data.samples[begin:end])
    eos.inprogress(f'Predicting observables from set \'{prediction}\' for {nsamples} samples')
    for i, sample in enumerate(progressbar(data.samples[begin:end])):
        for p, v in zip(parameters, sample):
            p.set(v)
        try:
            cache.update()
            observable_samples.append([cache[id] for id in observable_ids])
        except RuntimeError as e:
            eos.error(f'skipping prediction for sample {i} due to runtime error ({e}): {sample}')
            observable_samples.append([_np.nan for _ in observable_ids])
    observable_samples = _np.array(observable_samples)
    if mask_name is not None:
        eos.info(f'Applying mask {mask_name} to the observables')
        observable_samples = observable_samples[mask]
    eos.completed('... done')

    if mask_name is not None:
        weights = data.weights[mask]
    else:
        weights = data.weights[begin:end]
    filename = f'pred-{prediction}{_mask_output_fields(mask_name)["mask"]}'
    output_path = _staged(os.path.join(base_directory, 'data', posterior, filename))
    eos.data.Prediction.create(output_path, observables, observable_samples, weights)


def _calculate_mask(observable: eos.Observable, analysis_parameters: eos.Parameters, data: eos.data.ImportanceSamples):
    try:
        from tqdm.auto import tqdm
        progressbar = tqdm
    except ImportError:
        progressbar = lambda x: x

    parameters = [analysis_parameters[p['name']] for p in data.varied_parameters]
    nsamples = len(data.samples)
    observable_samples = _np.zeros(nsamples)
    eos.inprogress(f'Predicting observable \'{observable.name()}\' for {nsamples} samples')
    for i, sample in enumerate(progressbar(data.samples)):
        for p, v in zip(parameters, sample):
            p.set(v)
        try:
            observable_samples[i] = observable.evaluate()
        except RuntimeError as e:
            eos.error(f'skipping prediction for sample {i} due to runtime error ({e}): {sample}')
            observable_samples[i] = _np.nan
    return observable_samples > 0

# Create mask
@task('create-mask', 'data/{posterior}/mask-{mask_name}')
def create_mask(analysis_file:str, posterior:str, mask_name:str, base_directory:str='./'):
    """
    Creates a sample mask for a named posterior from a mask description in the analysis file.

    The named mask's description is a logical combination (``and`` or ``or``) of observables and/or other
    masks. Each observable contributes the subset of samples for which it evaluates positive, and nested
    masks are created recursively. The combined mask is stored as an :class:`eos.data.SampleMask` in
    EOS_BASE_DIRECTORY/data/POSTERIOR/mask-MASK_NAME, and is also returned as a boolean array.

    :param analysis_file: The name of the analysis file that describes the named posterior and mask, or an object of class `eos.AnalysisFile`.
    :type analysis_file: str or `eos.AnalysisFile`
    :param posterior: The name of the posterior whose samples are masked.
    :type posterior: str
    :param mask_name: The name of the mask to create, as defined in the analysis file.
    :type mask_name: str
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    :returns: The combined boolean mask, with one entry per posterior sample.
    :rtype: numpy.ndarray
    """
    _analysis = analysis_file.analysis(posterior, base_directory=base_directory)
    _parameters = _analysis.parameters
    data = eos.data.ImportanceSamples(os.path.join(base_directory, 'data', posterior, 'samples'))
    _check_varied_parameters_match(_analysis, data)

    mask_component = analysis_file._masks[mask_name]
    mask_logical_combination = mask_component.logical_combination
    if mask_logical_combination == "and":
        mask_combination_function = _np.all
    elif mask_logical_combination == "or":
        mask_combination_function = _np.any

    masks = []
    observables = []
    for d in mask_component.description:
        try:
            masks.append(create_mask(analysis_file, posterior, d.mask_name, base_directory))
            observables.append(d.mask_name)
        except AttributeError: # Not a mask but an observable
            masks.append(_calculate_mask(analysis_file.observable(posterior, d.name, _parameters), _parameters, data))
            observables.append(d.name)
    mask = mask_combination_function(_np.stack(masks), axis=0)
    eos.data.SampleMask.create(_staged(os.path.join(base_directory, 'data', posterior, f'mask-{mask_name}')), mask, observables)
    return mask


@dataclass
class GoodnessOfFitResult:
    r"""The result of a named goodness-of-fit test, as run by :func:`create_constraint` or
    called directly to diagnose why a fit failed one.

    :ivar name: The name of the test that produced this result.
    :ivar pvalue: The p-value; a fit is judged to fail the test if this is below the caller's threshold.
    :ivar statistic: The test's own statistic (e.g. the Kolmogorov-Smirnov distance for 'KS-test').
    :ivar squared_distances: The per-sample squared Mahalanobis distances against the fit, sorted ascending.
    """
    name:str
    pvalue:float
    statistic:float
    squared_distances:_np.ndarray


goodness_of_fit_tests = {}

def _goodness_of_fit_test(name):
    """Registers a goodness-of-fit test function under `name` in `goodness_of_fit_tests`.

    Private: this is registration machinery, not something a user calls; `goodness_of_fit_tests`
    itself (the registry it populates) and every test function in it are the public surface.
    """
    def _register(func):
        goodness_of_fit_tests[name] = func
        return func
    return _register

def _weighted_ecdf(values, weights):
    order = _np.argsort(values)
    sorted_values = values[order]
    cumulative_weights = _np.cumsum(weights[order])
    cumulative_weights = cumulative_weights / cumulative_weights[-1]
    return sorted_values, cumulative_weights

@_goodness_of_fit_test('KS-test')
def ks_test(samples, weights, mean, covariance):
    r"""
    Kolmogorov-Smirnov test of a multivariate Gaussian fit's squared Mahalanobis distances.

    For each sample x_i, the squared Mahalanobis distance to the fit,
    d_i^2 = (x_i - mean)^T covariance^-1 (x_i - mean), is expected to follow a chi2(nu)
    distribution (nu = the number of fitted dimensions) if the fit is a good description of
    the samples. This test compares the weighted empirical distribution of the d_i^2 against
    that expectation via a one-sample Kolmogorov-Smirnov statistic, converting the importance
    weights into an effective sample size (Kish's approximation) to obtain a p-value.

    Call this directly -- with the same samples/weights `create_constraint` used, and the
    mean/covariance it fitted -- to diagnose why a fit failed this test; the returned
    `squared_distances` can be compared against `scipy.stats.chi2(nu)` (nu = samples.shape[-1])
    to reproduce the diagnostic this test is based on, e.g. as a QQ-plot.

    :param samples: The samples to test, as a 2D array of shape (N, P).
    :type samples: numpy.ndarray
    :param weights: The importance weights of the samples, as a 1D array of shape (N, ).
    :type weights: numpy.ndarray
    :param mean: The fitted mean, as a 1D array of shape (P, ).
    :type mean: numpy.ndarray
    :param covariance: The fitted covariance, as a 2D array of shape (P, P).
    :type covariance: numpy.ndarray
    :returns: The test's result, including its p-value and the detail needed to diagnose a failure.
    :rtype: GoodnessOfFitResult
    """
    nu = samples.shape[-1]
    diffs = samples - mean
    inverse_covariance = _np.linalg.inv(covariance)
    squared_distances = _np.einsum('ij,jk,ik->i', diffs, inverse_covariance, diffs)

    sorted_distances, empirical_cdf = _weighted_ecdf(squared_distances, weights)
    theoretical_cdf = scipy.stats.chi2(nu).cdf(sorted_distances)
    ks_statistic = _np.max(_np.abs(empirical_cdf - theoretical_cdf))

    effective_sample_size = _np.sum(weights)**2 / _np.sum(weights**2)
    pvalue = scipy.stats.kstwo.sf(ks_statistic, max(1, round(effective_sample_size)))
    return GoodnessOfFitResult(name='KS-test', pvalue=pvalue, statistic=ks_statistic, squared_distances=sorted_distances)


@task('create-constraint', 'constraints/{constraint_name}')
def create_constraint(analysis_file:str, posterior:str, constraint_name:str, base_directory:str='./',
                       source:str='samples', parameters:list=None, tests:list=None,
                       goodness_of_fit_threshold:float=0.03, strict:bool=True):
    """
    Creates a MultivariateGaussian(Covariance) constraint from previously computed samples of a named posterior.

    The fit is the closed-form weighted maximum-likelihood estimate: the weighted sample mean and
    weighted sample covariance of the selected columns. One or more named goodness-of-fit tests
    are then run on the fit; see `tests` and `strict` below for how a failing test is handled.

    The input files are expected in EOS_BASE_DIRECTORY/data/POSTERIOR/SOURCE.
    The output file will be stored in EOS_BASE_DIRECTORY/constraints/CONSTRAINT_NAME/constraint.yaml.

    :param analysis_file: The name of the analysis file that describes the named posterior, or an object of class `eos.AnalysisFile`.
    :type analysis_file: str or `eos.AnalysisFile`
    :param posterior: The name of the posterior PDF from which to draw the samples.
    :type posterior: str
    :param constraint_name: The qualified name under which the constraint is registered; also used verbatim as the output subdirectory name.
    :type constraint_name: str
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    :param source: The name of the data file inside the posterior's data directory to fit. Defaults to 'samples'
        (the posterior's own importance samples). A value starting with 'pred-' selects the posterior-predictive
        importance samples of that name, previously written by the predict-observables task.
    :type source: str, optional
    :param parameters: The list of column names to constrain, matched against the loaded data's lookup table
        (bare parameter names for source='samples'; composite observable identifiers for a 'pred-' source).
        If None, every column of the loaded data is used.
    :type parameters: list of str, optional
    :param tests: The names of the goodness-of-fit tests that must be run before the constraint is accepted. Defaults to ['KS-test'].
        Each name must be a key of `goodness_of_fit_tests`; call one of those functions directly (e.g. `eos.ks_test`)
        on your own reloaded samples to diagnose why a test failed.
    :type tests: list of str, optional
    :param goodness_of_fit_threshold: The p-value threshold below which a test counts as failed. Either a single
        float applied to every test in `tests`, or a dict mapping each test name in `tests` to its own threshold.
    :type goodness_of_fit_threshold: float or dict[str, float], optional
    :param strict: If True (the default), a failed test logs an error and aborts the task -- no constraint file is written.
        If False, a failed test only logs a warning, and the constraint is written regardless.
    :type strict: bool, optional
    """
    # Validate tests/thresholds before touching any data, so a typo'd test name or a
    # threshold dict missing an entry fails immediately -- before any file access, fit, or
    # (for strict=True) before any earlier-in-the-list test has already run.
    _tests = tests if tests is not None else ['KS-test']
    unknown_tests = [t for t in _tests if t not in goodness_of_fit_tests]
    if unknown_tests:
        raise ValueError(
            f"Unknown goodness-of-fit test(s) {unknown_tests}; supported: {sorted(goodness_of_fit_tests)}"
        )
    if isinstance(goodness_of_fit_threshold, dict):
        missing_thresholds = [t for t in _tests if t not in goodness_of_fit_threshold]
        if missing_thresholds:
            raise ValueError(f"No goodness_of_fit_threshold given for test(s) {missing_thresholds}")
        thresholds = goodness_of_fit_threshold
    else:
        thresholds = {t: goodness_of_fit_threshold for t in _tests}

    analysis = analysis_file.analysis(posterior, base_directory=base_directory)

    path = os.path.join(base_directory, 'data', posterior, source)
    if source == 'samples':
        data = eos.data.ImportanceSamples(path)
        _check_varied_parameters_match(analysis, data)
    elif source.startswith('pred-'):
        data = eos.data.Prediction(path)
    else:
        raise ValueError(f"Unsupported source '{source}'; expected 'samples' or a 'pred-<name>' prediction directory")

    selected = parameters if parameters is not None else list(data.lookup_table.keys())
    columns = [data.lookup_table[name] for name in selected]
    samples = data.samples[:, columns]
    weights = data.weights

    eos.inprogress(f"Fitting a multivariate Gaussian to {len(selected)} column(s) of '{source}' from posterior '{posterior}'")
    mean = _np.average(samples, axis=0, weights=weights)
    covariance = _np.atleast_2d(_np.cov(samples.T, aweights=weights))

    for test_name in _tests:
        result = goodness_of_fit_tests[test_name](samples, weights, mean, covariance)
        threshold = thresholds[test_name]
        eos.info(
            f"Goodness-of-fit test '{test_name}' for constraint '{constraint_name}': "
            f"p-value = {result.pvalue:.3g}, statistic = {result.statistic:.3g}"
        )
        if result.pvalue < threshold:
            message = (
                f"Goodness-of-fit test '{test_name}' for constraint '{constraint_name}' failed "
                f"(p-value {result.pvalue:.3g} below threshold {threshold}); call "
                f"eos.goodness_of_fit_tests['{test_name}'](samples, weights, mean, covariance) "
                f"directly on these inputs to diagnose further"
            )
            if strict:
                eos.error(message)
                raise RuntimeError(message)
            eos.warn(message)

    if source == 'samples':
        names = list(selected)
        kinematics_list = [{} for _ in selected]
        options_list = [{} for _ in selected]
    else:
        names = [data.varied_parameters[i]['name'] for i in columns]
        kinematics_list = [data.varied_parameters[i]['kinematics'] for i in columns]
        options_list = [data.varied_parameters[i]['options'] for i in columns]

    body = {
        'type': 'MultivariateGaussian(Covariance)',
        'observables': names,
        'kinematics': kinematics_list,
        'options': options_list,
        'means': mean.tolist(),
        'covariance': covariance.tolist(),
    }
    output_path = os.path.join(base_directory, 'constraints', constraint_name, 'constraint.yaml')
    with open(os.path.join(_staged(os.path.dirname(output_path)), 'constraint.yaml'), 'w') as f:
        yaml.safe_dump({constraint_name: body}, f, sort_keys=False, default_flow_style=False)

    eos.success(f"Wrote constraint '{constraint_name}' to '{output_path}'")
    return body
