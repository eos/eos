# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2021-2026 Danny van Dyk
# Copyright (c) 2022-2023 Filip Novak
# Copyright (c) 2022-2025 Méril Reboud
# Copyright (c) 2023      Stephan Kürten
# Copyright (c) 2024      Ery McPartland
# Copyright (c) 2024-2025 Matthew Kirk
# Copyright (c) 2024      Stefan Meiser
# Copyright (c) 2025      Carolina Bolognani
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
import dynesty as _dynesty

from ._framework import task, _check_path_component, _check_varied_parameters_match, _mask_output_fields, _staged


@task('find-mode', 'data/{posterior}/mode-{label}{mask}', output_fields=_mask_output_fields)
def find_mode(analysis_file:str, posterior:str, base_directory:str='./', optimizations:int=3, start_point:list=None, chain:int=None,
              importance_samples:bool=None, seed:int=None, label:str='default', mask_name:str=None):
    '''
    Finds the mode of the named posterior using scipy SLSQP method by default.

    The optimization process can be initialized either with a random point,
    a provided parameter point, or by extracting the point with the largest posterior
    from among previously obtained MCMC or importance samples. The latter possibility usually gives the better performence.
    The optimization can also be iterated to increase the accuracy of the result.

    The output file will be stored in EOS_BASE_DIRECTORY/data/POSTERIOR/mode-LABEL.

    :param analysis_file: The name of the analysis file that describes the named posterior, or an object of class `eos.AnalysisFile`.
    :type analysis_file: str or `eos.AnalysisFile`
    :param posterior: The name of the posterior PDF from which to draw the samples.
    :type posterior: str
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    :param optimizations: The number of calls to the optimization algorithm.
    :type optimizations: int, optional
    :param start_point: If provided, the optimization will start at this point.
    :type start_point: `numpy.ndarray`, optional
    :param chain: If provided, the optimization will start at the point with the largest posterior of previously computed MCMC samples.
        The samples are expected to be stored in the `mcmc-XXXX` subdirectories of the `base_directory`, where 'XXXX' is the provided int.
    :type chain: int, optional
    :param importance_samples: If set to True, the optimization will start at the point with the largest posterior of previously computed importance samples.
    :type importance_samples: bool, optional
    :param seed: If provided, the optimization will start from a random point and the random number generator will be seeded with this value.
    :type seed: int, optional
    '''
    if optimizations < 1:
        raise ValueError('The number of optimizations should be larger than zero.')

    if not seed is None and not start_point is None:
        raise ValueError('The arguments seed and start_point are mutually exclusive')

    if not seed is None and not chain is None:
        raise ValueError('The arguments seed and chain are mutually exclusive')

    if not seed is None and not importance_samples is None:
        raise ValueError('The arguments seed and importance_samples are mutually exclusive')

    if not start_point is None and not chain is None:
        raise ValueError('The arguments chain and start_point are mutually exclusive')

    if not start_point is None and not importance_samples is None:
        raise ValueError('The arguments importance_samples and start_point are mutually exclusive')

    if not chain is None and not importance_samples is None:
        raise ValueError('The arguments importance_samples and chain are mutually exclusive')

    if mask_name is not None and importance_samples is None:
        raise ValueError('The argument mask-name can only be used with importance_samples')

    if mask_name is not None:
        _check_path_component('mask_name', mask_name)
        mask = eos.data.SampleMask(os.path.join(base_directory, 'data', posterior, f'mask-{mask_name}')).mask
    else:
        mask = slice(None) # Equivalent to mask = : but allowed

    analysis = analysis_file.analysis(posterior, base_directory=base_directory)

    eos.inprogress(f'Beginning minimization in {optimizations} points')
    if not start_point is None:
        _start_point = _np.array(start_point)
        eos.info('Using a user-provided starting point')

        _bfp = analysis.optimize(start_point=_start_point)
        _gof = eos.GoodnessOfFit(analysis._log_posterior)
        _chi2 = _gof.total_chi_square()
    elif not chain is None:
        eos.info(f'Using a starting point based on the existing MCMC data file (chain {chain:04})')
        _chain = eos.data.MarkovChain(os.path.join(base_directory, 'data', posterior, f'mcmc-{chain:04}'))
        idx_mode = _np.argmax(_chain.weights)
        _check_varied_parameters_match(analysis, _chain)
        analysis_varied_params = [p.name() for p in analysis.varied_parameters]
        samples_varied_params = [p["name"] for p in _chain.varied_parameters]
        if analysis_varied_params != samples_varied_params:
            raise ValueError("Parameters varied in the analysis file don't match those from the loaded sample")
        mode = '[ ' + ', '.join([f'{v:.4g}' for v in _chain.samples[idx_mode]]) + ' ]'
        eos.info(f'Using starting point {mode}')
        for p, v in zip(analysis.varied_parameters, _chain.samples[idx_mode]):
            p.set(v)

        _bfp = analysis.optimize()
        _gof = eos.GoodnessOfFit(analysis._log_posterior)
        _chi2 = _gof.total_chi_square()
    elif importance_samples:
        eos.info('Using a starting point based on the existing importance samples')
        _file = eos.data.ImportanceSamples(os.path.join(base_directory, 'data', posterior, 'samples'))
        if _file.posterior_values is None:
            FileNotFoundError("The argument importance_samples requires a valid 'posterior_values.npy' file.")
        idx_mode = _np.argmax(_file.posterior_values[mask])
        _check_varied_parameters_match(analysis, _file)
        mode = '[ ' + ', '.join([f'{v:.4g}' for v in _file.samples[mask][idx_mode]]) + ' ]'
        eos.info(f'Using starting point {mode}')
        for p, v in zip(analysis.varied_parameters, _file.samples[mask][idx_mode]):
            p.set(v)

        _bfp = analysis.optimize()
        _gof = eos.GoodnessOfFit(analysis._log_posterior)
        _chi2 = _gof.total_chi_square()
    else:
        eos.info('Using a random point')
        if seed is None:
            seed = 17
        _bfp = analysis.optimize(start_point='random', rng=_np.random.mtrand.RandomState(seed))
        _gof = eos.GoodnessOfFit(analysis._log_posterior)
        _chi2 = _gof.total_chi_square()

    eos.info('First optimization finished')
    # the first optimization provides the baseline, even if its chi^2 is not finite
    bfp      = _bfp
    gof      = _gof
    min_chi2 = _chi2

    for i in range(optimizations - 1):
        starting_point = [float(p) for p in analysis.varied_parameters]
        _bfp = analysis.optimize(start_point = starting_point)
        _gof = eos.GoodnessOfFit(analysis._log_posterior)
        _chi2 = _gof.total_chi_square()
        # a non-finite chi^2 never replaces the baseline, and never compares as an improvement
        if _np.isfinite(_chi2) and (not _np.isfinite(min_chi2) or _chi2 < min_chi2):
            gof = _gof
            bfp = _bfp
            min_chi2 = _chi2

    eos.completed('... optimization finished')
    eos.info('The best-fit point is:')
    for p, v in zip(analysis.varied_parameters, bfp.point):
        eos.info(f'  - {p.name()} -> {v}')
    eos.info(f'total chi^2 = {min_chi2:.2f}')
    eos.info(f'total dof   = {gof.total_degrees_of_freedom()}')
    pvalue = (1.0 - scipy.stats.chi2(gof.total_degrees_of_freedom()).cdf(gof.total_chi_square()))
    eos.info(f'p value     = {100 * pvalue:.2f}%')
    eos.info('individual test statistics:')
    for n, e in gof:
        local_pvalue = (1.0 - scipy.stats.chi2(e.dof).cdf(e.chi2))
        eos.info(f'  - {n}: chi^2 / dof = {e.chi2:.2f} / {e.dof}, local_pvalue = {100 * local_pvalue:.2f}%')

    eos.data.Mode.from_bfp_and_gof(
        _staged(os.path.join(base_directory, 'data', posterior, f'mode-{label}{_mask_output_fields(mask_name)["mask"]}')),
        bfp,
        gof
        )
    if (pvalue < 0.03):
        eos.warn(f'Final p value is {100 * pvalue:.2f}%, which is below the a-priori threshold 3%')
    else:
        eos.success(f'Final p value is {100 * pvalue:.2f}%, indicating a successful fit')

    return (bfp, gof)


@task('sample-prior', 'data/{posterior}/samples')
def sample_prior(analysis_file:str, posterior:str, base_directory:str='./', N:int=1000, seed:int=1701):
    """
    Samples from a named posterior PDF w/o likelihood information, an effective prior PDF.

    The output file will be stored in EOS_BASE_DIRECTORY/data/POSTERIOR/samples.

    :param analysis_file: The name of the analysis file that describes the named prior, or an object of class `eos.AnalysisFile`.
    :type analysis_file: str or `eos.AnalysisFile`
    :param posterior: The name of the posterior PDF (w/o any likelihood information) from which to draw the samples.
    :type posterior: str
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    :param N: The number of samples to be stored in the output file. Defaults to 1000.
    :type N: int, optional
    :param seed: The seed used to initialize the Mersenne Twister pseudo-random number generator.
    :type seed: int, optional
    """

    analysis = analysis_file.analysis(posterior, base_directory=base_directory)
    rng = _np.random.mtrand.RandomState(seed)
    eos.inprogress('Beginning prior sampling...')
    samples = analysis.sample_prior(N=N, rng=rng)
    weights =  _np.ones(N) / N
    eos.data.ImportanceSamples.create(_staged(os.path.join(base_directory, 'data', posterior, 'samples')), analysis.varied_parameters, samples, weights)
    eos.completed('...finished!')
    eos.info(f'Generated {N} samples of prior PDF for posterior {posterior}.')


class DynestyResultLogger:
    """Throttled progress logger for dynesty nested-sampling runs.

    An instance is passed as the ``print_function`` to a dynesty sampler. It forwards a one-line progress
    summary to the EOS logger, but at most once every 15 seconds, to avoid flooding the log.
    """

    def __init__(self):
        self.timer = _dynesty.utils.DelayTimer(15)

    def print_function(self, results, niter, ncall, add_live_it=None, dlogz=None, stop_val=None, nbatch=None, logl_min=-_np.inf, logl_max=_np.inf):
        """Log a one-line summary of the current nested-sampling state.

        Intended to be supplied as the ``print_function`` of a dynesty sampler. The summary is emitted to
        the EOS logger at most once per throttling interval; calls in between return without logging. The
        arguments mirror those passed by dynesty to its print function.

        :param results: The current sampling results as provided by dynesty.
        :param niter: The current iteration number.
        :type niter: int
        :param ncall: The cumulative number of likelihood evaluations.
        :type ncall: int
        :param add_live_it: The iteration at which the final live points were added, if applicable.
        :param dlogz: The current remaining-evidence estimate, if applicable.
        :param stop_val: The current value of the stopping criterion, if applicable.
        :param nbatch: The current batch number for dynamic nested sampling, if applicable.
        :param logl_min: The lower log-likelihood bound of the current batch.
        :param logl_max: The upper log-likelihood bound of the current batch.
        """
        if not self.timer.is_time():
            return

        fn_args = _dynesty.utils.get_print_fn_args(results, niter, ncall, add_live_it=add_live_it, dlogz=dlogz, stop_val=stop_val, nbatch=nbatch, logl_min=logl_min, logl_max=logl_max)
        eos.info(f'iteration {fn_args.niter} | {" | ".join(fn_args.long_str)}')


# Nested sampling
@task('sample-nested', 'data/{posterior}/nested')
def sample_nested(analysis_file:str, posterior:str, base_directory:str='./', bound:str='multi', nlive:int=250, dlogz:float=1.0, maxiter:int=None, miniter:int=0, min_ess:int=0, seed:int=10, sample:str='auto'):
    """
    Samples from a likelihood associated with a named posterior using dynamic nested sampling.

    The output will be stored in EOS_BASE_DIRECTORY/data/POSTERIOR/nested.
    In addition, an ImportanceSamples object is exported to EOS_BASE_DIRECTORY/data/POSTERIOR/samples.

    :param analysis_file: The name of the analysis file that describes the named posterior, or an object of class `eos.AnalysisFile`.
    :type analysis_file: str or `eos.AnalysisFile`
    :param posterior: The name of the posterior.
    :type posterior: str
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    :param bound: The option for bounding the target distribution. For valid values, see the dynesty documentation. Defaults to 'multi'.
    :type bound: str, optional
    :param nlive: The number of live points.
    :type nlive: int, optional
    :param dlogz: The relative tolerance for the remaining evidence. Defaults to 1.0.
    :type dlogz: float, optional
    :param maxiter: The maximum number of iterations. Iterations may stop earlier if the termination condition is reached.
    :type maxiter: int, optional
    :param miniter: The minimum number of iterations. If not provided, the sampler will run until the termination condition is reached. Defaults to 0.
    :type miniter: int, optional
    :param min_ess: The minimum effective sample size (ESS). Samples will be added until the estimated ESS reaches min_ess, even if the termination condition is reached earlier. If 0 (default), dynesty's default ESS target is used.
    :type min_ess: int, optional
    :param seed: The seed used to initialize the Mersenne Twister pseudo-random number generator.
    :type seed: int, optional
    :param sample: The method used for sampling within the likelihood constraints. For valid values, see dynesty documentation. Defaults to 'auto'.
    :type sample: str, optional
    """
    eos.inprogress('Beginning sampling...')
    analysis = analysis_file.analysis(posterior, base_directory=base_directory)
    logger = DynestyResultLogger()
    results = analysis.sample_nested(bound=bound, nlive=nlive, dlogz=dlogz, maxiter=maxiter, miniter=miniter, min_ess=min_ess, print_function=logger.print_function, seed=seed, sample=sample)
    samples = results.samples
    posterior_values = results.logwt - results.logz[-1]
    weights = _np.exp(posterior_values)
    ess = _dynesty.utils.get_neff_from_logwt(results.logwt)
    eos.completed('...finished!')
    eos.info(f'Finished sampling with {len(samples)} samples (effective sample size {ess:.0f}) and evidence estimate {results.logz[-1]:.2f} +/- {results.logzerr[-1]:.2f}')
    eos.data.DynestyResults.create(_staged(os.path.join(base_directory, 'data', posterior, 'nested')), analysis.varied_parameters, results)
    eos.data.ImportanceSamples.create(_staged(os.path.join(base_directory, 'data', posterior, 'samples')), analysis.varied_parameters,
                                      samples, weights, posterior_values=posterior_values)
