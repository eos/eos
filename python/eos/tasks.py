# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2020-2025 Danny van Dyk
# Copyright (c) 2023      Philip Lüghausen
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
import contextlib
import functools
import glob
import inspect
import logging
import eos.analysis_file_description
import numpy as _np
import os
import pypmc
import scipy
import sys
import copy as _copy
import warnings
import dynesty as _dynesty

from dataclasses import asdict
from .ipython import __ipython__

class LogfileHandler:
    def __init__(self, path, name='log', mode='w'):
        self.path = path
        self.name = name
        self.mode = mode
        self.formatter = logging.Formatter('%(asctime)-15s %(levelname)-10s %(message)s')
        self.handler = logging.FileHandler(os.path.join(path, name), mode=mode)
        self.handler.setLevel(logging.DEBUG)
        self.handler.setFormatter(self.formatter)

    def __enter__(self):
        eos.logger.addHandler(self.handler)

    def __exit__(self, type, value, traceback):
        eos.logger.removeHandler(self.handler)


_tasks = {}

def task(name, output, mode=lambda **kwargs: 'w', modules=[], logfile=True):
    def _task(func):
        @functools.wraps(func)
        def task_wrapper(*args, **kwargs):
            # import task-specific optional modules
            try:
                import importlib as _importlib
                for module in modules:
                    _importlib.import_module(module)
            except ModuleNotFoundError as e:
                eos.error(f'failed to import missing module \'{module}\': {e}')
                raise e
            except ImportError as e:
                eos.error(f'failed to import module \'{module}\': {e}')
                raise e

            # extract default arguments
            _args = {
                k: v.default
                for k, v in inspect.signature(func).parameters.items()
                if v.default is not inspect.Parameter.empty
            }
            _args.update(zip(func.__code__.co_varnames, args))
            _args.update(kwargs)
            if 'analysis_file' in _args and type(_args['analysis_file']) is str:
                _args.update({ 'analysis_file': eos.AnalysisFile(_args['analysis_file'])})
            # create output directory if needed directly or for logging
            if output or logfile:
                outputpath = ('{base_directory}/' + output).format(**_args)
                os.makedirs(outputpath, exist_ok=True)
            if logfile:
                # create invocation-specific log file handler
                handler = LogfileHandler(path=outputpath, mode=mode(**_args))
            else:
                handler = contextlib.nullcontext()
            # use invocation-specific log file handler
            with handler:
                iaccordion = None
                ioutput = contextlib.suppress() # can be replaced with contextlib.nullcontext once Python >=3.7 is ensured
                if __ipython__:
                    import ipywidgets as _ipywidgets
                    ioutput = _ipywidgets.Output(layout={'height': '200px', 'overflow': 'auto'})
                    iaccordion = _ipywidgets.Accordion(children=[ioutput])
                    iaccordion.set_title(0, output.format(**_args))
                    display(iaccordion)
                # use invocation-specific ipython output widget (if available)
                with ioutput:
                    result = func(**_args)
                    if iaccordion:
                        iaccordion.selected_index = None
                    return result
        _tasks[name] = task_wrapper
        return task_wrapper
    return _task


def _check_varied_parameters_match(analysis: eos.Analysis, data):
    # Check the parameters varied in the analysis match those of the loaded samples
    analysis_varied_params = [p.name() for p in analysis.varied_parameters]
    samples_varied_params = [p["name"] for p in data.varied_parameters]
    if analysis_varied_params != samples_varied_params:
        raise ValueError(f"Parameters varied in the analysis file don't match those from the loaded sample")


@task('find-mode', 'data/{posterior}/mode-{label}')
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
        mask = eos.data.SampleMask(os.path.join(base_directory, 'data', posterior, f'mask-{mask_name}')).mask
        label += f'_mask-{mask_name}'
    else:
        mask = slice(None) # Equivalent to mask = : but allowed

    analysis = analysis_file.analysis(posterior)
    min_chi2 = sys.float_info.max
    gof = None
    bfp = None

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
            raise ValueError(f"Parameters varied in the analysis file don't match those from the loaded sample")
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

    eos.info(f'First optimization finished')
    for i in range(optimizations - 1):
        starting_point = [float(p) for p in analysis.varied_parameters]
        _bfp = analysis.optimize(start_point = starting_point)
        _gof = eos.GoodnessOfFit(analysis._log_posterior)
        _chi2 = _gof.total_chi_square()
        if _chi2 < min_chi2:
            gof = _gof
            bfp = _bfp
            min_chi2 = _chi2

    eos.completed(f'... optimization finished')
    eos.info('The best-fit point is:')
    for p, v in zip(analysis.varied_parameters, _bfp.point):
        eos.info(f'  - {p.name()} -> {v}')
    eos.info(f'total chi^2 = {min_chi2:.2f}')
    eos.info(f'total dof   = {gof.total_degrees_of_freedom()}')
    pvalue = (1.0 - scipy.stats.chi2(gof.total_degrees_of_freedom()).cdf(gof.total_chi_square()))
    eos.info(f'p value     = {100 * pvalue:.2f}%')
    eos.info('individual test statistics:')
    local_pvalues = {}
    for n, e in gof:
        local_pvalue = (1.0 - scipy.stats.chi2(e.dof).cdf(e.chi2))
        local_pvalues[f'{n}'] = float(local_pvalue)
        eos.info(f'  - {n}: chi^2 / dof = {e.chi2:.2f} / {e.dof}, local_pvalue = {100 * local_pvalue:.2f}%')

    eos.data.Mode.create(
        os.path.join(base_directory, 'data', posterior, f'mode-{label}'),
        analysis.varied_parameters,
        bfp.point,
        pvalue,
        local_pvalues,
        min_chi2,
        gof.total_degrees_of_freedom()
        )
    if (pvalue < 0.03):
        eos.warn(f'Final p value is {100 * pvalue:.2f}%, which is below the a-priori threshold 3%')
    else:
        eos.success(f'Final p value is {100 * pvalue:.2f}%, indicating a successful fit')

    return (bfp, gof)


@task('sample-mcmc', 'data/{posterior}/mcmc-{chain:04}')
def sample_mcmc(analysis_file:str, posterior:str, chain:int, base_directory:str='./', pre_N:int=150, preruns:int=3, N:int=1000, stride:int=5, cov_scale:float=0.1, start_point:list=None):
    """
    Samples from a named posterior PDF using Markov Chain Monte Carlo (MCMC) methods.

    The output file will be stored in EOS_BASE_DIRECTORY/data/POSTERIOR/mcmc-CHAIN.

    :param analysis_file: The name of the analysis file that describes the named posterior, or an object of class `eos.AnalysisFile`.
    :type analysis_file: str or `eos.AnalysisFile`
    :param posterior: The name of the posterior PDF from which to draw the samples.
    :type posterior: str
    :param chain: The index assigned to the Markov chain. This value is used to seed the RNG for a reproducible analysis.
    :type chain: int >= 0
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    :param pre_N: The number of samples to be used for an adaptation in each prerun steps. These samples will be discarded.
    :type pre_N: int, optional
    :param preruns: The number of prerun steps, which are used to adapt the MCMC proposal to the posterior.
    :type preruns: int, optional
    :param N: The number of samples to be stored in the output file. Defaults to 1000.
    :type N: int, optional
    :param stride: The ratio of samples drawn over samples stored. For every S samples, S - 1 will be discarded. Defaults to 5.
    :type stride: int, optional
    :param cov_scale: Scale factor for the initial guess of the covariance matrix.
    :type cov_scale: float, optional
    :param start_point: Optional starting point for the chain
    :type start_point: list-like, optional
    """

    eos.inprogress(f'Beginning sampling...')

    analysis = analysis_file.analysis(posterior)
    rng = _np.random.mtrand.RandomState(int(chain) + 1701)
    try:
        samples, usamples, weights = analysis.sample(N=N, stride=stride, pre_N=pre_N, preruns=preruns, rng=rng, cov_scale=cov_scale, start_point=start_point, return_uspace=True)
        eos.data.MarkovChain.create(os.path.join(base_directory, 'data', posterior, f'mcmc-{chain:04}'), analysis.varied_parameters, samples, usamples, weights)
    except RuntimeError as e:
        eos.error(f'encountered run time error ({e}) in parameter point:')
        for p in analysis.varied_parameters:
            eos.error(f' - {p.name()}: {p.evaluate()}')
    eos.completed(f'...finished!')
    eos.info(f'Generated {N} samples from posterior {posterior}.')

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

    analysis = analysis_file.analysis(posterior)
    rng = _np.random.mtrand.RandomState(seed)
    eos.inprogress(f'Beginning prior sampling...')
    samples = analysis.sample_prior(N=N, rng=rng)
    weights =  _np.ones(N) / N
    eos.data.ImportanceSamples.create(os.path.join(base_directory, 'data', posterior, 'samples'), analysis.varied_parameters, samples, weights)
    eos.completed('...finished!')
    eos.info(f'Generated {N} samples of prior PDF for posterior {posterior}.')

@task('find-clusters', 'data/{posterior}/clusters')
def find_clusters(posterior:str, base_directory:str='./', threshold:float=2.0, K_g:int=1, analysis_file:str=None):
    r"""
    Finds clusters among posterior MCMC samples, grouped by Gelman-Rubin R value, and creates a Gaussian mixture density.

    Finding clusters and creating a Gaussian mixture density is a necessary intermediate step before using the sample-pmc subcommand.
    The input files are expected in EOS_BASE_DIRECTORY/data/POSTERIOR/mcmc-\*. All MCMC input files present will be used in the clustering.
    The output files will be stored in EOS_BASE_DIRECTORY/data/POSTERIOR/clusters.

    :param posterior: The name of the posterior.
    :type posterior: str
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    :param threshold: The R value threshold. If two sample subsets have an R value larger than this threshold, they will be treated as two distinct clusters. Defaults to 2.0.
    :type threshold: float > 1.0, optional
    :param K_g: The number of mixture components per cluster. Default to 1.
    :type K_g: int >= 1, optional
    """

    import pathlib
    input_paths = [str(p) for p in pathlib.Path(os.path.join(base_directory, 'data', posterior)).glob('mcmc-*')]
    chains    = [eos.data.MarkovChain(path).usamples for path in input_paths]
    n = len(chains[0])
    for chain in chains:
        assert len(chain) == n, 'Every chains must contain the same number of samples'

    groups = pypmc.mix_adapt.r_value.r_group([_np.mean(chain.T, axis=1) for chain in chains],
                           [_np.var (chain.T, axis=1, ddof=1) for chain in chains],
                           n, threshold)
    eos.info(f'Found {len(groups)} groups using an R value threshold of {threshold}')
    density   = pypmc.mix_adapt.r_value.make_r_gaussmix(chains, K_g=K_g, critical_r=threshold)
    eos.success(f'Created mixture density with {len(density.components)} components')
    eos.data.MixtureDensity.create(os.path.join(base_directory, 'data', posterior, 'clusters'), density)


@task('mixture-product', 'data/{posterior}/product')
def mixture_product(posterior:str, posteriors:list, base_directory:str='./', analysis_file:str=None):
    """
    Compute the cartesian product of the densities listed in posteriors. Note that this product is not commutative.

    The input densities are read from EOS_BASE_DIRECTORY/data/POSTERIOR_i/pmc, where POSTERIOR_i is listed in posteriors.
    The output density will be stored in EOS_BASE_DIRECTORY/data/POSTERIOR/product.

    :param posterior: The name of the posterior.
    :type posterior: str
    :param posteriors: The list of names of the posteriors whose mixture densities will be concatenated.
    :type posteriors: iterable of str
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    """

    eos.inprogress('Beginning computation.')
    densities = [eos.data.PMCSampler(os.path.join(base_directory, 'data', p, 'pmc')).density() for p in posteriors]
    output_path = os.path.join(base_directory, 'data', posterior, 'product')
    eos.data.MixtureDensity.create(output_path, eos.data.MixtureDensity.cartesian_product(densities))
    eos.completed('...finished!')

# Sample PMC
@task('sample-pmc', 'data/{posterior}/pmc', mode=lambda initial_proposal, **kwargs: 'a' if initial_proposal != 'clusters' else 'a')
def sample_pmc(analysis_file:str, posterior:str, base_directory:str='./', step_N:int=500, steps:int=10, final_N:int=5000,
               perplexity_threshold:float=1.0, weight_threshold:float=1e-10, sigma_test_stat:list=None, initial_proposal:str='clusters',
               pmc_iterations:int=1, pmc_rel_tol:float=1e-10, pmc_abs_tol:float=1e-05, pmc_lookback:int=1):
    """
    Samples from a named posterior using the Population Monte Carlo (PMC) methods.

    The results of the find-cluster command are expected in EOS_BASE_DIRECTORY/data/POSTERIOR/clusters.
    The output will be stored in EOS_BASE_DIRECTORY/data/POSTERIOR/pmc.
    In addition, an ImportanceSamples object is exported to EOS_BASE_DIRECTORY/data/POSTERIOR/samples.

    :param analysis_file: The name of the analysis file that describes the named posterior, or an object of class `eos.AnalysisFile`.
    :type analysis_file: str or `eos.AnalysisFile`
    :param posterior: The name of the posterior.
    :type posterior: str
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    :param step_N: The number of samples to be used in each adaptation step. These samples will be discarded. Defaults to 500.
    :type step_N: int > 0, optional
    :param steps: The number of adaptation steps, which are used to adapt the PMC proposal to the posterior. Defaults to 10.
    :type steps: int > 0, optional
    :param final_N: The number of samples to be stored in the output file. Defaults to 5000,
    :type final_N: int > 0, optional
    :param perplexity_threshold: The threshold for the perplexity in the last step after which further adaptation steps are to be skipped. Defaults to 1.0.
    :type perplexity_threshold: 0.0 < float <= 1.0, optional
    :param weight_threshold: Mixture components with a weight smaller than this threshold are pruned.
    :type weight_threshold: 0.0 < float <= 1.0, optional.
    :param sigma_test_stat: If provided, the inverse CDF of -2*log(PDF) will be evaluated, using the provided values as the respective significance.
    :type sigma_test_stat: list or iterable
    :param initial_proposal: Specify where the initial proposal should be taken from:

     * ``clusters``: use the proposal obtained using `find-clusters` (default)

     * ``product``: use the proposal obtained from `mixture_product`

     * ``pmc``: continue sampling from the previous `sample-pmc` results.

    :type initial_proposal: str, optional
    :param pmc_iterations: Maximum number of update of the PMC, changing this value may make the update unstable.
    :type pmc_iterations: int > 0, optional, advanced
    :param pmc_rel_tol: Relative tolerance of the PMC. If two consecutive values of the current density log-likelihood are relatively smaller than this value, the convergence is declared.
    :type pmc_rel_tol: float > 0.0, optional, advanced
    :param pmc_abs_tol: Absolute tolerance of the PMC. If two consecutive values of the current density log-likelihood are smaller than this value, the convergence is declared.
    :type pmc_abs_tol: float > 0.0, optional, advanced
    :param pmc_lookback: Use reweighted samples from the previous update steps when adjusting the mixture density. The parameter determines the number of update steps to "look back". The default value of 1 disables this feature, a value of 0 means that all previous steps are used.
    :type pmc_lookback: int >= 0, optional
    """

    analysis = analysis_file.analysis(posterior)
    rng = _np.random.mtrand.RandomState(1701)
    eos.inprogress('Beginning sampling...')
    if initial_proposal == 'clusters':
        initial_density = eos.data.MixtureDensity(os.path.join(base_directory, 'data', posterior, 'clusters')).density()
    elif initial_proposal == 'pmc':
        previous_sampler = eos.data.PMCSampler(os.path.join(base_directory, 'data', posterior, 'pmc'))
        initial_density = previous_sampler.density()
    elif initial_proposal == 'product':
        initial_density = eos.data.MixtureDensity(os.path.join(base_directory, 'data', posterior, 'product')).density()
    else:
        eos.error(f"Could not initialize proposal in sample_pmc: argument {initial_proposal} is not supported.")

    samples, weights, posterior_values, proposal = analysis.sample_pmc(initial_density, step_N=step_N, steps=steps, final_N=final_N,
                                                     rng=rng, final_perplexity_threshold=perplexity_threshold,
                                                     weight_threshold=weight_threshold, pmc_iterations=pmc_iterations,
                                                     pmc_rel_tol=pmc_rel_tol, pmc_abs_tol=pmc_abs_tol, pmc_lookback=pmc_lookback)

    if initial_proposal == 'pmc':
        samples = _np.concatenate((previous_sampler.samples, samples), axis=0)
        weights = _np.concatenate((previous_sampler.weights, weights), axis=0)

    eos.data.PMCSampler.create(os.path.join(base_directory, 'data', posterior, 'pmc'), analysis.varied_parameters, proposal,
                               sigma_test_stat=sigma_test_stat, samples=samples, weights=weights)
    eos.data.ImportanceSamples.create(os.path.join(base_directory, 'data', posterior, 'samples'), analysis.varied_parameters,
                                      samples, weights, posterior_values=posterior_values)
    eos.completed('...finished!')
    eos.info(f'Finished sampling with {len(samples)} samples.')

# Predict observables
@task('predict-observables', 'data/{posterior}/pred-{prediction}')
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
    _analysis      = analysis_file.analysis(posterior)
    _parameters    = _analysis.parameters
    cache          = eos.ObservableCache(_parameters)
    observables    = analysis_file.observables(posterior, prediction, _parameters)
    observable_ids = [cache.add(o) for o in observables]

    if mask_name is not None and (begin != 0 or end is not None):
        raise ValueError('The arguments mask-name and begin or end are mutually exclusive')
    if mask_name is not None:
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
    eos.completed(f'... done')

    filename = f'pred-{prediction}'
    if mask_name is not None:
        filename += f'_mask-{mask_name}'
        weights = data.weights[mask]
    else:
        weights = data.weights[begin:end]
    output_path = os.path.join(base_directory, 'data', posterior, filename)
    eos.data.Prediction.create(output_path, observables, observable_samples, weights)


# Run one analysis step
@task('run', '')
def run(analysis_file:str, id:str, base_directory:str='./', dry_run:bool=False):
    """
    Runs one out of a list of predefined steps recorded in the analysis file.

    Each step corresponds to a call to one or more of the common tasks, e.g.,
     - sample-mcmc
     - find-cluster
     - sample-pmc
     - predict-observables

    :param analysis_file: The name of the analysis file that describes the named posterior, or an object of class `eos.AnalysisFile`.
    :type analysis_file: str or `eos.AnalysisFile`
    :param id: The id of the step to run.
    :type id: str
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    :param dry_run: The flag that disables execution and instead prints the full information on the tasks that would be run to standard output. Defaults to `False`.
    :type dry_run: bool, optional
    """
    try:
        step = analysis_file._steps[id]
    except KeyError:
        raise ValueError(f'Step with id \'{id}\' not found in analysis file')

    for task in step.tasks:
        if dry_run:
            arguments = step.default_arguments[task.task] | task.arguments | { 'analysis_file': analysis_file.analysis_file, 'base_directory': base_directory }
            print(f"{task.task}({','.join([f'{k}={v}' for k, v in arguments.items()])})")
        else:
            arguments = step.default_arguments[task.task] | task.arguments | { 'analysis_file': analysis_file, 'base_directory': base_directory }
            task_function = _tasks[task.task]

            task_function(**arguments)


class DynestyResultLogger:
    def __init__(self):
        self.timer = _dynesty.utils.DelayTimer(15)

    def print_function(self, results, niter, ncall, add_live_it=None, dlogz=None, stop_val=None, nbatch=None, logl_min=-_np.inf, logl_max=_np.inf):
        if not self.timer.is_time():
            return

        fn_args = _dynesty.utils.get_print_fn_args(results, niter, ncall, add_live_it=add_live_it, dlogz=dlogz, stop_val=stop_val, nbatch=nbatch, logl_min=logl_min, logl_max=logl_max)
        eos.info(f'iteration {fn_args.niter} | {" | ".join(fn_args.long_str)}')


# Nested sampling
@task('sample-nested', 'data/{posterior}/nested')
def sample_nested(analysis_file:str, posterior:str, base_directory:str='./', bound:str='multi', nlive:int=250, dlogz:float=1.0, maxiter:int=None, miniter:int=0, seed:int=10, sample:str='auto'):
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
    :param dlogz: Relative tolerance for the remaining evidence. Defaults to 5%.
    :type dlogz: float, optional
    :param maxiter: The maximum number of iterations. Iterations may stop earlier if the termination condition is reached.
    :type maxiter: int, optional
    :param miniter: The minimum number of iterations. If not provided, the sampler will run until the termination condition is reached. Defaults to 0.
    :type miniter: int, optional
    :param seed: The seed used to initialize the Mersenne Twister pseudo-random number generator.
    :type seed: int, optional
    :param sample: The method used for sampling within the likelihood constraints. For valid values, see dynesty documentation. Defaults to 'auto'.
    :type sample: str, optional
    """
    eos.inprogress('Beginning sampling...')
    analysis = analysis_file.analysis(posterior)
    logger = DynestyResultLogger()
    results = analysis.sample_nested(bound=bound, nlive=nlive, dlogz=dlogz, maxiter=maxiter, miniter=miniter, print_function=logger.print_function, seed=seed, sample=sample)
    samples = results.samples
    posterior_values = results.logwt - results.logz[-1]
    weights = _np.exp(posterior_values)
    eos.completed('...finished!')
    eos.info(f'Finished sampling with {len(samples)} samples and evidence estimate {results.logz[-1]:.2f} +/- {results.logzerr[-1]:.2f}')
    eos.data.DynestyResults.create(os.path.join(base_directory, 'data', posterior, 'nested'), analysis.varied_parameters, results)
    eos.data.ImportanceSamples.create(os.path.join(base_directory, 'data', posterior, 'samples'), analysis.varied_parameters,
                                      samples, weights, posterior_values=posterior_values)


def _get_modes(posterior:str, base_directory:str='./'):
    result = []
    search_path = os.path.join(base_directory, 'data', posterior, 'mode-*')

    for mode_dir in glob.glob(search_path):
        name = mode_dir.split('mode-')[1]
        mode = eos.data.Mode(mode_dir)
        result.append((name, mode))

    return result

def _get_references(analysis_file):
    all_references = eos.References()
    all_constraints = eos.Constraints()
    result = []
    for likelihood in analysis_file.likelihoods.values():
        constraints         = likelihood.constraints
        for constraint in constraints:
            constraint_name = constraint.constraint.split(';')[0]
            for reference in all_constraints[constraint_name].references():

                try:
                    result.append(all_references[reference].inspire_id())
                except:
                    warnings.warn(f'No reference found for {reference}')

    result = _np.unique(result)
    return result

# Create a report
@task('report', 'reports')
def report(analysis_file:str, template_file:str, base_directory:str='./', generate_pdf:bool=True):
    """
    Generates a report from an analysis file and a Jinja2 template file.

    The processed file will be stored in the same directory as the template file, with the the extension '.jinja' or '.jinja2' stripped from the file name.
    If `generate_pdf` is set to True, the intermediate file will be converted to PDF using Pandoc and XeLaTeX. This requires a working installation of both programs.

    :param analysis_file: The name of the analysis file that shall be used for the report.
    :type analysis_file: str or `eos.AnalysisFile`
    :param template_file: The name of the Jinja2 template file that shall be used for the report.
    :type template_file: str
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    :param convert_to_pdf: The flag that enables the conversion of the intermediate file to PDF. Defaults to `True`.
    :type convert_to_pdf: bool, optional
    """
    import jinja2

    eos.inprogress(f'Generating report from template \'{template_file}\' for analysis file \'{analysis_file}\'')

    resource_directory = os.path.dirname(template_file)

    intermediate_file, template_ext = os.path.splitext(template_file)
    if template_ext not in ['.jinja', '.jinja2']:
        raise ValueError(f'Template file is expected to have an extension of \'.jinja\' or \'.jinja2\', but got \'{template_ext}\'')

    eos.info(f'Processing template \'{template_file}\' to intermediate \'{intermediate_file}\'')

    basename, intermediate_ext = os.path.splitext(intermediate_file)
    if intermediate_ext not in ['.md']:
        raise ValueError(f'Intermediate file is expected to have an extension of \'.md\', but got \'{intermediate_ext}\'')

    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(searchpath='./'),
        undefined=jinja2.StrictUndefined, # error for undefined variables in
                                          # template
        trim_blocks=True,
        lstrip_blocks=True)
    template = env.get_template(template_file)
    result = template.render(
        eos={ 'version': eos.__version__ },
        analysis_file=analysis_file,
        metadata=asdict(analysis_file._metadata),
        analyses={posterior: analysis_file.analysis(posterior) for posterior in analysis_file.posteriors},
        base_directory=base_directory,
        len=len,
        zip=zip,
        modes=lambda posterior: _get_modes(posterior, base_directory=base_directory),
        references=_get_references(analysis_file),
    )

    with open(intermediate_file, 'w') as f:
        f.write(result)

    if not generate_pdf:
        eos.completed(f'... report generation finished')
        return

    output_file = basename + '.pdf'

    eos.info(f'Processing intermediate file \'{intermediate_file}\' to final output \'{output_file}\'')

    import pypandoc

    pypandoc.convert_file(
        source_file=intermediate_file,
        outputfile=output_file,
        to='pdf',
        extra_args=[
            f'--resource-path=./:{resource_directory}/:{eos._pkg_data_dir}/',
            f'--variable=eoslogo:{eos._pkg_data_dir}/report-logo.pdf',
            f'--variable=eosversion:{eos.__version__}',
            '--template=report-template.tex',
            '--pdf-engine=xelatex'
        ]
    )

    eos.completed(f'... report generation finished')


# Draw figures
@task('draw-figure', 'figures', mode=lambda **kwargs: 'a')
def draw_figure(analysis_file:str, figure_name:str, base_directory:str='./', format:str|list[str]='pdf'):
    """
    Draws figures from the analysis file.

    The output files will be stored in EOS_BASE_DIRECTORY/figures/.

    :param analysis_file: The name of the analysis file that describes the named posterior, or an object of class `eos.AnalysisFile`.
    :type analysis_file: str or `eos.AnalysisFile`
    :param figure_name: The name of the figure to draw.
    :type figure_name: str
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    :param format: The file extension of the data files. Can also be a list of file extensions.
    :type format: str or list of str, optional
    """
    if figure_name not in analysis_file._figures.keys():
        raise ValueError(f'Figure with name \'{figure_name}\' not found in analysis file')

    if isinstance(format, str):
        format = [format]

    SUPPORTED_FORMATS = ['pdf', 'png', 'svg']
    for fmt in format:
        if fmt not in SUPPORTED_FORMATS:
            raise ValueError(f'Figure format \'{fmt}\' not supported. Supported formats are: {",".join(SUPPORTED_FORMATS)}')

    eos.inprogress(f'Drawing figure {figure_name}')
    context = eos.analysis_file_description.AnalysisFileContext(base_directory=base_directory)
    figure = analysis_file._figures[figure_name]
    for fmt in format:
        eos.info(f'Drawing figure {figure_name} in \'{fmt}\' format')
        figure.draw(context, output=os.path.join(base_directory, 'figures', f'{figure_name}.{fmt}'))
    eos.completed(f'... drawing finished')


# Corner plot
@task('corner-plot', 'data/{posterior}/plots')
def corner_plot(analysis_file:str, posterior:str, base_directory:str='./', format:str='pdf', distribution:str='posterior', begin:int=0, end:int=None, mask_name=None):
    """
    Generates a corner plot of the 1-D and 2-D marginalized posteriors.

    The input files are expected in EOS_BASE_DIRECTORY/data/POSTERIOR/samples.
    The output files will be stored in EOS_BASE_DIRECTORY/data/POSTERIOR/plots.

    :param analysis_file: The name of the analysis file that describes the named posterior, or an object of class `eos.AnalysisFile`.
    :type analysis_file: str or `eos.AnalysisFile`
    :param posterior: The name of the posterior.
    :type posterior: str
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    :param format: The file extension of the data files. Can also be a list of file extensions.
    :type format: str or list of str, optional
    :param distribution: The distribution to plot. Can be either 'posterior' or the name of a prediction.
    :type distribution: str, optional
    :param begin: The index of the first parameter to plot. Defaults to 0.
    :type begin: int, optional
    :param end: The index beyond the last parameter to plot. Defaults to ``None``.
    :type end: int or NoneType, optional
    :param mask_name: The label of the mask to apply to the observables. Defaults to None.
    :type mask_name: str, optional
    """
    import matplotlib.pyplot as _plt

    analysis = analysis_file.analysis(posterior)

    if mask_name is not None and (begin != 0 or end is not None):
        raise ValueError('The arguments mask-name and begin or end are mutually exclusive')
    if mask_name is not None:
        mask = eos.data.SampleMask(os.path.join(base_directory, 'data', posterior, f'mask-{mask_name}')).mask

    # Provide file with distribution samples and LaTeX labels for either parameters or observables
    if distribution == 'posterior':
        f = eos.data.ImportanceSamples(os.path.join(base_directory, 'data', posterior, 'samples'))
        _check_varied_parameters_match(analysis, f)
        labels = [p.latex() for p in analysis.varied_parameters]

    elif isinstance(distribution, str):
        try:
            prediction = analysis_file.predictions[distribution]
        except KeyError:
            raise RuntimeError(f"No prediction with name '{distribution}' exists in the analysis file")

        f = eos.data.Prediction(os.path.join(base_directory, 'data', posterior, f'pred-{distribution}'))

        observables = eos.Observables()
        qualified_names = [obs.name for obs in prediction.observables]
        labels = [f'${observables[qn].latex()}$' for qn in qualified_names]

    else:
        raise RuntimeError(f"Argument 'distribution' must be one of {['posterior',] + list(analysis_file.predictions.keys())}")

    # Apply slicing according to begin and end for samples and labels
    samples = f.samples[:, begin:end]
    weights = f.weights[:]
    labels = labels[begin:end]
    # Apply mask if specified
    if mask_name is not None:
        samples = samples[mask]
        weights = weights[mask]
    size = samples.shape[-1]

    fig, axes = _plt.subplots(size, size, figsize=(3.0 * size, 3.0 * size), dpi=100)

    for i in range(size):
        # diagonal
        ax = axes[i, i]
        variable_samples = samples[:, i] # For a single variable

        ax.hist(variable_samples, weights=weights, alpha=0.5, bins=100, density=True, stacked=True, color='C1')
        xmin = _np.min(variable_samples)
        xmax = _np.max(variable_samples)
        ax.set_xlim((xmin, xmax))
        ax.set_xlabel(labels[i])
        ax.set_ylabel(labels[i])
        ax.set_aspect(_np.diff((xmin, xmax))[0] / _np.diff(ax.get_ylim())[0])

        for j in range(0, size):
            # off-diagonal
            if j < i:
                axes[i, j].set_axis_off()

            if j <= i:
                continue

            ax = axes[i, j]

            # Samples for two single variables as x and y data
            xsamples = samples[:, j]
            ysamples = samples[:, i]

            xmin = _np.min(xsamples)
            xmax = _np.max(xsamples)
            ymin = _np.min(ysamples)
            ymax = _np.max(ysamples)
            ax.hist2d(xsamples, ysamples, weights=weights, alpha=1.0, bins=100, cmap='Greys', rasterized=True)
            ax.set_xlim((xmin, xmax))
            ax.set_ylim((ymin, ymax))
            ax.set_aspect(_np.diff((xmin, xmax))[0] / _np.diff((ymin, ymax))[0])

    fig.tight_layout()

    _format = _copy.copy(format)
    if isinstance(_format, str):
        _format = [ _format ]
    filename = 'corner-plot'
    if mask_name is not None:
        filename += f'-mask-{mask_name}'
    for f in _format:
        fig.savefig(os.path.join(base_directory, 'data', posterior, 'plots', f'{filename}.{f}'))


@task('validate', '', logfile=False)
def validate(analysis_file:str):
    """
    Validates the analysis file by checking that all posteriors and all prediction sets can be created.

    :param analysis_file: The name of the analysis file that describes the named posterior, or an object of class `eos.AnalysisFile`.
    :type analysis_file: str or :class:`eos.AnalysisFile`
    """
    analysis_file.validate()


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
    _analysis = analysis_file.analysis(posterior)
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
    eos.data.SampleMask.create(os.path.join(base_directory, 'data', posterior, f'mask-{mask_name}'), mask, observables)
    return mask


@task('list-figures', '', logfile=False)
def list_figures(analysis_file:str):
    """
    Lists all figures that can be created from the analysis file.

    :param analysis_file: The name of the analysis file, or an object of class `eos.AnalysisFile`.
    :type analysis_file: str or :class:`eos.AnalysisFile`
    """
    return [name for name in analysis_file._figures.keys()]


@task('list-steps', '', logfile=False)
def list_steps(analysis_file:str):
    """
    Lists the steps in the analysis file.

    :param analysis_file: The name of the analysis file that shall be inspected`.
    :type analysis_file: str or :class:`eos.AnalysisFile`
    """
    return analysis_file._steps.keys()


@task('list-step-dependencies', '', logfile=False)
def list_step_dependencies(analysis_file:str, id:str):
    """
    Lists all steps required to be completed before the given step can be executed.

    :param analysis_file: The name of the analysis file that shall be inspected`.
    :type analysis_file: str or :class:`eos.AnalysisFile`
    :param id: The id of the step to inspect.
    :type id: str
    """
    steps = analysis_file._steps

    if id not in steps:
        raise ValueError(f'Step with id \'{id}\' not found in analysis file')

    return steps[id].depends_on
