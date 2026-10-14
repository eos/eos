# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2021-2026 Danny van Dyk
# Copyright (c) 2021-2022 Méril Reboud
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

from ._framework import task, _staged


@task('sample-mcmc', 'data/{posterior}/mcmc-{chain:04}')
def sample_mcmc(analysis_file:str, posterior:str, chain:int, base_directory:str='./', pre_N:int=150, preruns:int=3, N:int=1000, stride:int=5, cov_scale:float=0.1, start_point:list=None):
    """
    Samples from a named posterior PDF using Markov Chain Monte Carlo (MCMC) methods.

    The output file will be stored in EOS_BASE_DIRECTORY/data/POSTERIOR/mcmc-CHAIN.

    :param analysis_file: The name of the analysis file that describes the named posterior, or an object of class `eos.AnalysisFile`.
    :type analysis_file: str or `eos.AnalysisFile`
    :param posterior: The name of the posterior PDF from which to draw the samples.
    :type posterior: str
    :param chain: The index assigned to the Markov chain, which must be non-negative. This value is used to seed the RNG for a reproducible analysis.
    :type chain: int
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
    :type start_point: collections.abc.Sequence, optional
    """

    eos.inprogress('Beginning sampling...')

    analysis = analysis_file.analysis(posterior, base_directory=base_directory)
    rng = _np.random.mtrand.RandomState(int(chain) + 1701)
    try:
        samples, usamples, weights = analysis.sample(N=N, stride=stride, pre_N=pre_N, preruns=preruns, rng=rng, cov_scale=cov_scale, start_point=start_point, return_uspace=True)
        eos.data.MarkovChain.create(_staged(os.path.join(base_directory, 'data', posterior, f'mcmc-{chain:04}')), analysis.varied_parameters, samples, usamples, weights)
    except RuntimeError as e:
        eos.error(f'encountered run time error ({e}) in parameter point:')
        for p in analysis.varied_parameters:
            eos.error(f' - {p.name()}: {p.evaluate()}')
        raise
    eos.completed('...finished!')
    eos.info(f'Generated {N} samples from posterior {posterior}.')


@task('find-clusters', 'data/{posterior}/clusters', modules=['pypmc'])
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
    :param threshold: The R value threshold, which must be larger than 1.0. If two sample subsets have an R value larger than this threshold, they will be treated as two distinct clusters. Defaults to 2.0.
    :type threshold: float, optional
    :param K_g: The number of mixture components per cluster, which must be at least 1. Default to 1.
    :type K_g: int, optional
    """

    import pathlib
    input_paths = [str(p) for p in pathlib.Path(os.path.join(base_directory, 'data', posterior)).glob('mcmc-*')]
    inputs    = [eos.data.MarkovChain(path) for path in input_paths]
    chains    = [mcmc.usamples for mcmc in inputs]
    parameters = inputs[0].varied_parameters
    n = len(chains[0])
    for mcmc, chain in zip(inputs, chains):
        assert len(chain) == n, 'Every chains must contain the same number of samples'
        assert mcmc.varied_parameters == parameters, 'Every chain must vary the same parameters, in the same order'

    groups = pypmc.mix_adapt.r_value.r_group([_np.mean(chain.T, axis=1) for chain in chains],
                           [_np.var (chain.T, axis=1, ddof=1) for chain in chains],
                           n, threshold)
    eos.info(f'Found {len(groups)} groups using an R value threshold of {threshold}')
    density   = pypmc.mix_adapt.r_value.make_r_gaussmix(chains, K_g=K_g, critical_r=threshold)
    eos.success(f'Created mixture density with {len(density.components)} components')
    eos.data.MixtureDensity.create(_staged(os.path.join(base_directory, 'data', posterior, 'clusters')), density, parameters)


@task('mixture-product', 'data/{posterior}/product', modules=['pypmc'])
def mixture_product(posterior:str, posteriors:list, base_directory:str='./', analysis_file:str=None):
    """
    Compute the cartesian product of the densities listed in posteriors. Note that this product is not commutative.

    The input densities are read from EOS_BASE_DIRECTORY/data/POSTERIOR_i/pmc, where POSTERIOR_i is listed in posteriors.
    The output density will be stored in EOS_BASE_DIRECTORY/data/POSTERIOR/product.

    :param posterior: The name of the posterior.
    :type posterior: str
    :param posteriors: The list of names of the posteriors whose mixture densities will be concatenated.
    :type posteriors: collections.abc.Iterable of str
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    """

    eos.inprogress('Beginning computation.')
    samplers = [eos.data.PMCSampler(os.path.join(base_directory, 'data', p, 'pmc')) for p in posteriors]
    densities = [sampler.density() for sampler in samplers]
    # the cartesian product concatenates the components in the order of the input densities
    parameters = [p for sampler in samplers for p in sampler.varied_parameters]
    output_path = _staged(os.path.join(base_directory, 'data', posterior, 'product'))
    eos.data.MixtureDensity.create(output_path, eos.data.MixtureDensity.cartesian_product(densities), parameters)
    eos.completed('...finished!')

# Sample PMC
@task('sample-pmc', 'data/{posterior}/pmc', mode=lambda initial_proposal, **kwargs: 'a' if initial_proposal != 'clusters' else 'a', modules=['pypmc'])
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
    :param step_N: The number of samples to be used in each adaptation step, which must be positive. These samples will be discarded. Defaults to 500.
    :type step_N: int, optional
    :param steps: The number of adaptation steps, which must be positive, and which are used to adapt the PMC proposal to the posterior. Defaults to 10.
    :type steps: int, optional
    :param final_N: The number of samples to be stored in the output file, which must be positive. Defaults to 5000,
    :type final_N: int, optional
    :param perplexity_threshold: The threshold for the perplexity in the last step after which further adaptation steps are to be skipped. Must be larger than 0.0 and at most 1.0. Defaults to 1.0.
    :type perplexity_threshold: float, optional
    :param weight_threshold: Mixture components with a weight smaller than this threshold are pruned. Must be larger than 0.0 and at most 1.0.
    :type weight_threshold: float, optional
    :param sigma_test_stat: If provided, the inverse CDF of -2*log(PDF) will be evaluated, using the provided values as the respective significance.
    :type sigma_test_stat: list or collections.abc.Iterable
    :param initial_proposal: Specify where the initial proposal should be taken from:

     * ``clusters``: use the proposal obtained using `find-clusters` (default)

     * ``product``: use the proposal obtained from `mixture_product`

     * ``pmc``: continue sampling from the previous `sample-pmc` results.

    :type initial_proposal: str, optional
    :param pmc_iterations: Maximum number of update of the PMC, which must be positive. This is an advanced option; changing this value may make the update unstable.
    :type pmc_iterations: int, optional
    :param pmc_rel_tol: Relative tolerance of the PMC, which must be positive. This is an advanced option. If two consecutive values of the current density log-likelihood are relatively smaller than this value, the convergence is declared.
    :type pmc_rel_tol: float, optional
    :param pmc_abs_tol: Absolute tolerance of the PMC, which must be positive. This is an advanced option. If two consecutive values of the current density log-likelihood are smaller than this value, the convergence is declared.
    :type pmc_abs_tol: float, optional
    :param pmc_lookback: Use reweighted samples from the previous update steps when adjusting the mixture density. The parameter determines the number of update steps to "look back", and must be non-negative. The default value of 1 disables this feature, a value of 0 means that all previous steps are used.
    :type pmc_lookback: int, optional
    """

    analysis = analysis_file.analysis(posterior, base_directory=base_directory)
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

    eos.data.PMCSampler.create(_staged(os.path.join(base_directory, 'data', posterior, 'pmc')), analysis.varied_parameters, proposal,
                               sigma_test_stat=sigma_test_stat, samples=samples, weights=weights)
    eos.data.ImportanceSamples.create(_staged(os.path.join(base_directory, 'data', posterior, 'samples')), analysis.varied_parameters,
                                      samples, weights, posterior_values=posterior_values)
    eos.completed('...finished!')
    eos.info(f'Finished sampling with {len(samples)} samples.')
