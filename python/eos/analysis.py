#!/usr/bin/python
# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2018, 2019, 2020 Danny van Dyk
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
import copy as _cp
from logging import info, warn, debug
import numpy as np
import scipy

class BestFitPoint:
    """
    Represents the best-fit point of a Bayesian analysis undertaken with the :class:`Analysis <eos.Analysis>` class.
    """
    def __init__(self, analysis, point):
        self.analysis = analysis
        self.point = point


    def _repr_html_(self):
        result = '<table>\n'
        result += '<tr><th>parameter</th><th>value</th></tr>\n'
        for p, v in zip(self.analysis.varied_parameters, self.point):
            name = p.name()
            latex = p.latex()
            name = latex if latex else name
            result += '<tr><td>{n}</td><td>{v:6.4f}</td></tr>'.format(n=name, v=v)
        result += '</table>'

        return(result)



class Analysis:
    """Represents a statistical analysis.

    Describes a Bayesian analysis in terms of a set of parameters, a log(likelihood),
    and a set containing one or more log(prior)s.

    :param global_options: The options as (key, value) pairs that shall be forwarded to all theory predictions.
    :type global_options: dict, optional
    :param priors: The priors for this analysis as a list of prior descriptions. See :ref:`below <eos-Analysis-prior-descriptions>` for what consitutes a valid prior description.
    :type priors: iterable
    :param likelihood: The likelihood as a list of individual constraints from the internal data base of experimental and theoretical constraints; cf. `the complete list of constraints <../constraints.html>`_.
    :type likelihood: iterable
    """

    def __init__(self, priors, likelihood, global_options=None):
        """Constructor."""
        self.init_args = { 'priors': priors, 'likelihood': likelihood, 'global_options': global_options }
        self.parameters = eos.Parameters.Defaults()
        self.global_options = eos.Options()
        self.log_likelihood = eos.LogLikelihood(self.parameters)
        self.log_posterior = eos.LogPosterior(self.log_likelihood)
        self.varied_parameters = []
        self.bounds = []

        # collect the global options
        if global_options:
            for key, value in global_options.items():
                self.global_options.set(key, value)

        # create the priors
        for prior in priors:
            parameter = prior['parameter']
            minv = prior['min']
            maxv = prior['max']
            prior_type = prior['type'] if 'type' in prior else 'uniform'
            if 'uniform' == prior_type or 'flat' == prior_type:
                self.log_posterior.add(eos.LogPrior.Flat(self.parameters, parameter, eos.ParameterRange(minv, maxv)), False)
            elif 'gauss' == prior_type or 'gaussian' == prior_type:
                central = prior['central']
                sigma = prior['sigma']
                if type(sigma) is list or type(sigma) is tuple:
                    sigma_lo = sigma[0]
                    sigma_hi = sigma[1]
                else:
                    sigma_lo = sigma
                    sigma_hi = sigma
                self.log_posterior.add(
                    eos.LogPrior.Gauss(
                        self.parameters, parameter, eos.ParameterRange(minv, maxv),
                        central - sigma_lo, central, central + sigma_hi
                    ),
                    False)
            else:
                raise ValueError('Unknown prior type \'{}\''.format(prior_type))

            self.bounds.append((minv, maxv))
            p = self.parameters[parameter]
            p.set_min(minv)
            p.set_max(maxv)
            self.varied_parameters.append(p)

        # create the likelihood
        for constraint_name in likelihood:
            constraint = eos.Constraint.make(constraint_name, self.global_options)
            self.log_likelihood.add(constraint)

        # perform some sanity checks
        varied_parameter_names = set([p.name() for p in self.varied_parameters])
        used_parameter_names = set()
        for observable in self.log_likelihood.observable_cache():
            for i in observable.used_parameter_ids():
                used_parameter_names.add(self.parameters.by_id(i).name())

        used_but_unvaried = used_parameter_names - varied_parameter_names
        if (len(used_but_unvaried) > 0):
            info('likelihood probably depends on {} parameter(s) that do not appear in the prior; check prior?'.format(len(used_but_unvaried)))
        for n in used_but_unvaried:
            debug('used, but not included in any prior: \'{}\''.format(n))
        for n in varied_parameter_names - used_parameter_names:
            warn('likelihood does not depend on parameter \'{}\'; remove from prior or check options!'.format(n))

    def clone(self):
        """Returns an independent instance of eos.Analysis."""
        return eos.Analysis(**self.init_args)

    def goodness_of_fit(self):
        """Returns a goodness-of-fit summary for the current parameter point."""
        return eos.GoodnessOfFit(self.log_posterior)


    def optimize(self, start_point=None, **kwargs):
        """
        Optimize the log(posterior) and returns a best-fit-point summary.

        :param start_point: Parameter point from which to start the optimization, with the elements in the same order as in eos.Analysis.varied_parameters. If not specified, optimization starts at the current parameter point.
        :param start_point: iterable, optional
        """
        from logging import info, warn

        if start_point == None:
            start_point = [float(p) for p in self.varied_parameters]

        default_kwargs = { 'method': 'SLSQP', 'options': { 'ftol': 1.0e-13 } }
        if kwargs is None:
            kwargs = default_kwargs

        res = scipy.optimize.minimize(
            self.negative_log_pdf,
            start_point,
            args=None,
            bounds=self.bounds,
            **kwargs)

        if not res.success:
            warn('Optimization did not succeed')
            warn('  optimizer'' message reas: {}'.format(res.message))
        else:
            info('Optimization goal achieved after {nfev} function evaluations'.format(nfev=res.nfev))

        for p, v in zip(self.varied_parameters, res.x):
            p.set(v)

        return eos.BestFitPoint(self, res.x)


    def log_pdf(self, x, *args):
        """
        Adapter for use with external optimization software (e.g. pypmc) to aid when optimizing the log(posterior).

        :param x: Parameter point, with the elements in the same order as in eos.Analysis.varied_parameters.
        :type x: iterable
        :param args: Dummy parameter (ignored)
        :type args: optional
        """
        for p, v in zip(self.varied_parameters, x):
            p.set(v)

        return(self.log_posterior.evaluate())


    def negative_log_pdf(self, x, *args):
        """
        Adapter for use with external optimization software (e.g. scipy.optimize.minimize) to aid when optimizing the log(posterior).

        :param x: Parameter point, with the elements in the same order as in eos.Analysis.varied_parameters.
        :type x: iterable
        :param args: Dummy parameter (ignored)
        :type args: optional
        """
        for p, v in zip(self.varied_parameters, x):
            p.set(v)

        return(-self.log_posterior.evaluate())


    def sample(self, N=1000, stride=5, pre_N=150, preruns=3, cov_scale=0.1, observables=None, start_point=None, rng=np.random.mtrand):
        """
        Return samples of the parameters, log(weights), and optionally posterior-predictive samples for a sequence of observables.

        Obtains random samples of the log(posterior) using an adaptive Markov Chain Monte Carlo with PyPMC.
        A prerun with adaptations is carried out first and its samples are discarded.

        :param N: Number of samples that shall be returned
        :param stride: Stride, i.e., the number by which the actual amount of samples shall be thinned to return N samples.
        :param pre_N: Number of samples in each prerun.
        :param preruns: Number of preruns.
        :param cov_scale: Scale factor for the initial guess of the covariance matrix.
        :param observables: Observables for which posterior-predictive samples shall be obtained.
        :type observables: list-like, optional
        :param start_point: Optional starting point for the chain
        :type start_point: list-like, optional
        :param rng: Optional random number generator (must be compatible with the requirements of pypmc.sampler.markov_chain.MarkovChain)

        :return: A tuple of the parameters as array of size N, the logarithmic weights as array of size N, and optionally the posterior-predictive samples of the observables as array of size N x len(observables).

        .. note::
           This method requiries the PyPMC python module, which can be installed from PyPI.
        """
        import logging
        import pypmc
        try:
            from tqdm import tqdm
            progressbar = tqdm
        except ImportError:
            progressbar = lambda x: x

        ind_lower = np.array([bound[0] for bound in self.bounds])
        ind_upper = np.array([bound[1] for bound in self.bounds])
        ind = pypmc.tools.indicator.hyperrectangle(ind_lower, ind_upper)

        log_target = pypmc.tools.indicator.merge_function_with_indicator(self.log_pdf, ind, -np.inf)

        # create initial covariance
        sigma = np.diag([np.square(bound[1] - bound[0]) / 12 * cov_scale for bound in self.bounds])
        log_proposal = pypmc.density.gauss.LocalGauss(sigma)

        # create start point, if not provided
        if start_point is None:
            u = np.array([rng.uniform(0.0, 1.0) for j in range(0, len(ind_lower))])
            ubar = 1.0 - u
            start_point = ubar * ind_upper + u * ind_lower

        # create MC sampler
        sampler = pypmc.sampler.markov_chain.AdaptiveMarkovChain(log_target, log_proposal, start_point, save_target_values=True, rng=rng)

        # pre run to adapt markov chains
        for i in range(0, preruns):
            logging.info('Prerun {} out of {}'.format(i, preruns))
            accept_count = sampler.run(pre_N)
            accept_rate  = accept_count / pre_N * 100
            logging.info('Prerun {}: acceptance rate is {:3.0f}%'.format(i, accept_rate))
            sampler.adapt()
        sampler.clear()

        # obtain final samples
        logging.info('Main run: started ...')
        sample_total  = N * stride
        sample_chunk  = sample_total // 100
        sample_chunks = [sample_chunk for i in range(0, 99)]
        sample_chunks.append(sample_total - 99 * sample_chunk)
        for current_chunk in progressbar(sample_chunks):
            accept_count = accept_count + sampler.run(current_chunk)
        accept_rate  = accept_count / (N * stride) * 100
        logging.info('Main run: acceptance rate is {:3.0f}%'.format(accept_rate))

        parameter_samples = sampler.samples[:][::stride]
        weights = sampler.target_values[:][::stride, 0]

        if not observables:
            return(parameter_samples, weights)
        else:
            observable_samples = []
            for parameters in parameter_samples:
                for p, v in zip(self.varied_parameters, parameters):
                    p.set(v)

                observable_samples.append([o.evaluate() for o in observables])

            return(parameter_samples, weights, np.array(observable_samples))


    def sample_pmc(self, log_proposal, step_N=1000, steps=10, final_N=5000, rng=np.random.mtrand):
        """
        Return samples of the parameters and log(weights)

        Obtains random samples of the log(posterior) using adaptive importance sampling following
        the Popoulation Monte Carlo approach with PyPMC.

        :param step_N: Number of samples that shall be drawn in each adaptation step.
        :param steps: Number of adaptation steps.
        :param final_N: Number of samples that shall be drawn after all adaptation steps.
        :param rng: Optional random number generator (must be compatible with the requirements of pypmc.sampler.importance_sampler.ImportancSampler)

        :return: A tuple of the parameters as array of length N = pre_N * steps + final_N, the (linear) weights as array of length N, and the
            final proposal function as pypmc.density.mixture.MixtureDensity.

        .. note::
           This method requires the PyPMC python module, which can be installed from PyPI.
        """
        import logging
        import pypmc
        try:
            from tqdm import tqdm
            progressbar = tqdm
        except ImportError:
            progressbar = lambda x: x

        ind_lower = np.array([bound[0] for bound in self.bounds])
        ind_upper = np.array([bound[1] for bound in self.bounds])
        ind = pypmc.tools.indicator.hyperrectangle(ind_lower, ind_upper)

        log_target = pypmc.tools.indicator.merge_function_with_indicator(self.log_pdf, ind, -np.inf)

        # create PMC sampler
        sampler = pypmc.sampler.importance_sampling.ImportanceSampler(log_target, log_proposal, save_target_values=True)#, rng=rng)
        generating_components = []

        # carry out adaptions
        for step in progressbar(range(steps)):
            origins = sampler.run(step_N, trace_sort=True)
            generating_components.append(origins)
            samples = sampler.samples[:]
            weights = sampler.weights[:][:, 0]
            normalized_weights = np.ma.masked_where(weights <= 0, weights) / np.sum(weights)
            entropy = -1.0 * np.dot(np.log(normalized_weights), normalized_weights)
            perplexity = np.exp(entropy) / len(normalized_weights)
            debug('Perplexity after sampling in step {}: {}'.format(step, perplexity))
            pypmc.mix_adapt.pmc.gaussian_pmc(samples, sampler.proposal, weights, mincount=0, rb=True, copy=False)
            sampler.proposal.normalize()

        # draw final samples
        origins = sampler.run(final_N, trace_sort=True)
        generating_components.append(origins)
        samples = sampler.samples[:]
        weights = sampler.weights[:][:, 0]
        normalized_weights = np.ma.masked_where(weights <= 0, weights) / np.sum(weights)
        entropy = -1.0 * np.dot(np.log(normalized_weights), normalized_weights)
        perplexity = np.exp(entropy) / len(normalized_weights)
        info('Perplexity after final samples: {}'.format(perplexity))

        return samples, weights, sampler.proposal
