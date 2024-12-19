#!/usr/bin/python
# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2018-2024 Danny van Dyk
# Copyright (c) 2024 Carolina Bolognani
# Copyright (c) 2023 Lorenz Gärtner
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
import numpy as np
import scipy
import pypmc

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
            result += f'<tr><td>{name}</td><td>{v:6.4f}</td></tr>'
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
    :param likelihood: The likelihood as a list of individual constraints from the internal data base of experimental and theoretical constraints; cf. `the complete list of constraints <../reference/constraints.html>`_.
    :type likelihood: iterable
    :param external_likelihood: The external likelihood blocks as a list or iterable of objects returned by :py:meth:`eos.LogLikelihoodBlock.External`.
    :type external_likelihood: list or iterable of :py:class:`eos.LogLikelihoodBlock`.
    :param manual_constraints: Additional manually-specified constraints that shall be added to the log(likelihood).
    :type manual_constraints: dict, optional
    :param fixed_parameters: Values of parameters that are set when the analysis is defined.
    :type fixed_parameters: dict, optional
    :param parameters: The optional set of parameters that shall be used for this analysis. Defaults to `None` which means that a new instance of :class:`eos.Parameters` is created.
    :type parameters: :class:`eos.Parameters` or None, optional
    """

    def __init__(self, priors, likelihood, external_likelihood=[], global_options={}, manual_constraints={}, fixed_parameters={}, parameters=None):
        """Constructor."""
        self.init_args = { 'priors': priors, 'likelihood': likelihood, 'external_likelihood': external_likelihood, 'global_options': global_options, 'manual_constraints': manual_constraints, 'fixed_parameters':fixed_parameters }
        self.parameters = parameters if parameters else eos.Parameters.Defaults()
        """The set of parameters used for this analysis."""
        self.global_options = eos.Options()
        self._constraint_names = []
        self._log_likelihood = eos.LogLikelihood(self.parameters)
        self._log_posterior = eos.LogPosterior(self._log_likelihood)
        self.varied_parameters = []
        self.varied_parameter_names = []

        eos.info('Creating analysis with {nprior} priors, {nconst} EOS-wide constraints, {nopts} global options, {nmanual} manually-entered constraints and {nparams} fixed parameters.'.format(
            nprior=len(priors), nconst=len(likelihood), nopts=len(global_options), nmanual=len(manual_constraints), nparams=len(fixed_parameters)))
        eos.debug('priors:')
        for p in priors:
            if 'parameter' in p:
                eos.debug(' - {name} ({type})'.format(name=p['parameter'], type=p['type']))
            elif 'constraint' in p:
                eos.debug(' - {name} (constraint)'.format(name=p['constraint']))
            elif 'parameters' in p:
                eos.debug(' - {names} ({type})'.format(names=p['parameters'], type=p['type']))
        eos.debug('constraints:')
        for cn in likelihood:
            eos.debug(f' - {cn}')
        eos.debug('manual_constraints:')
        for cn, ce in manual_constraints.items():
            eos.debug(f' - {cn}')
        eos.debug('fixed_parameters:')
        for pn, pe in fixed_parameters.items():
            eos.debug(f' - {pn}')

        # collect the global options
        for key, value in global_options.items():
            self.global_options.declare(key, value)

        # Fix specified parameters
        for param, value in fixed_parameters.items():
            self.parameters.set(param, value)

        # create the priors
        for prior in priors:
            if 'parameter' in prior and 'constraint' in prior:
                raise ValueError('Prior specification must not contain both a parameter and a constraint')

            if 'parameter' in prior and 'parameters' in prior:
                raise ValueError('Prior specification must not contain both a parameter and parameters')

            if 'parameters' in prior and 'constraint' in prior:
                raise ValueError('Prior specification must not contain both parameters and a constraint')

            if 'parameter' in prior:
                prior_type = prior['type'] if 'type' in prior else 'uniform'
                parameter = prior['parameter']

                if prior_type in ['uniform', 'flat', 'scale']: # min / max is mandatory
                    minv = float(prior['min'])
                    maxv = float(prior['max'])

                    if 'uniform' == prior_type or 'flat' == prior_type:
                        self._log_posterior.add(eos.LogPrior.Flat(self.parameters, parameter, minv, maxv), False)
                    elif 'scale' == prior_type:
                        mu_0 = prior['mu_0']
                        lambda_scale = prior['lambda']
                        self._log_posterior.add(eos.LogPrior.Scale(self.parameters,
                            parameter, minv, maxv, mu_0, lambda_scale), False)
                elif prior_type in ['gauss', 'gaussian']: # min / max is optional
                    if ('min' in prior) != ('max' in prior):
                        raise ValueError('Prior specification must contain both min and max or neither')

                    curtailed = 'min' in prior

                    central = prior['central']
                    sigma = prior['sigma']

                    if curtailed:
                        minv = float(prior['min'])
                        maxv = float(prior['max'])
                        if type(sigma) is list or type(sigma) is tuple:
                            sigma_lo = sigma[0]
                            sigma_hi = sigma[1]
                        else:
                            sigma_lo = sigma
                            sigma_hi = sigma
                        self._log_posterior.add(
                            eos.LogPrior.CurtailedGauss(
                                self.parameters, parameter, minv, maxv,
                                central - sigma_lo, central, central + sigma_hi
                            ),
                            False)
                    else:
                        if type(sigma) is list or type(sigma) is tuple:
                            raise ValueError('Prior specification must contain a scalar sigma value for a gaussian prior with infinite support')
                        self._log_posterior.add(
                            eos.LogPrior.Gaussian(
                                self.parameters, parameter, central, sigma
                            ),
                            False)
                elif prior_type == 'poisson':
                    k = float(prior['k'])
                    self._log_posterior.add(
                        eos.LogPrior.Poisson(self.parameters, parameter, k),
                        False)
                else:
                    raise ValueError(f'Unknown prior type \'{prior_type}\'')

                p = self.parameters[parameter]
                if p.name() in self.varied_parameter_names:
                    raise ValueError(f'Parameter {p} is repeated in a univariate prior.')
                self.varied_parameters.append(p)
                self.varied_parameter_names.append(p.name())
            elif 'parameters' in prior:
                prior_type = prior['type']
                parameters = prior['parameters']
                for p in parameters:
                    if p in self.varied_parameter_names:
                        raise ValueError(f'Parameter {p} is repeated in a multivariate prior.')
                    self.varied_parameters.append(self.parameters[p])
                    self.varied_parameter_names.append(p)

                if prior_type=='transform':
                    minv = list[float](prior['min'])
                    maxv = list[float](prior['max'])
                    shift = list[float](prior['shift'])
                    transform = list[list[float]](prior['transform'])
                    self._log_posterior.add(eos.LogPrior.Transform(self.parameters, parameters, shift, transform, minv, maxv), False)
                else:
                    raise ValueError(f'Unknown prior type \'{prior_type}\'')
            elif 'constraint' in prior:
                constraint_name = eos.QualifiedName(prior['constraint'])
                constraint_entry = eos.Constraints()[constraint_name]
                log_prior = constraint_entry.make_prior(self.parameters, constraint_name.options_part())
                self._log_posterior.add(log_prior, False)
                for p in log_prior.varied_parameters():
                    if p.name() in self.varied_parameter_names:
                        raise ValueError(f'Parameter {p} is repeated in a prior using constraint {constraint_name}')
                    self.varied_parameters.append(p)
                    self.varied_parameter_names.append(p.name())
            else:
                raise ValueError('Prior specification must contains either \'parameter\', \'parameters\', or \'constraint\'')

        # check for duplicate entries in the likelihood
        set_likelihood = set(likelihood)
        if len(set_likelihood) != len(likelihood):
            raise ValueError(f'The likelihood contains duplicate entries')
        set_manual_constraints = set(manual_constraints.keys())
        if len(set_manual_constraints) != len(manual_constraints):
            raise ValueError(f'The manual constraints contain duplicate entries')
        set_intersection = set_likelihood.intersection(set_manual_constraints)
        if len(set_intersection) > 0:
            raise ValueError(f'The likelihood contains constraint names also present in the manual_constraints: {set_intersection}')

        # record all constraints that comprise the likelihood
        self._constraint_names = list(set_likelihood.union(set_manual_constraints))

        # create the likelihood
        for constraint_name in self._constraint_names:
            if constraint_name in manual_constraints.keys():
                import yaml
                yaml_string = yaml.dump(self._sanitize_manual_input(manual_constraints[constraint_name]))
                constraint_entry = eos.ConstraintEntry.deserialize(constraint_name, yaml_string)
                constraint = constraint_entry.make(constraint_name, self.global_options)
            else:
                constraint = eos.Constraint.make(constraint_name, self.global_options)
            self._log_likelihood.add(constraint)

        # add external likelihood
        for likelihood in external_likelihood:
            if not isinstance(likelihood, eos.LogLikelihoodBlock):
                raise ValueError('`external_likelihood` must be a list of eos.LogLikelihoodBlock')
            self._log_likelihood.add(likelihood)

        # perform some sanity checks
        varied_parameter_names = {p.name() for p in self.varied_parameters}
        used_parameter_names = set()
        fixed_parameter_names = set(fixed_parameters.keys())
        for observable in self._log_likelihood.observable_cache():
            for i in observable.used_parameter_ids():
                used_parameter_names.add(self.parameters.by_id(i).name())

        used_but_unvaried = used_parameter_names - varied_parameter_names - fixed_parameter_names
        if (len(used_but_unvaried) > 0):
            eos.info(f'likelihood probably depends on {len(used_but_unvaried)} parameter(s) that do not appear in the prior; check prior?')
        for n in used_but_unvaried:
            eos.debug(f'used, but not included in any prior: \'{n}\'')
        for n in varied_parameter_names - used_parameter_names:
            eos.warn(f'likelihood does not depend on parameter \'{n}\'; remove from prior or check options!')


    def _u_to_par(self, u):
        """Internal function that uses the inverse prior transform to translate from u ∈ [0, 1)^D to the parameter space"""
        for p, uv in zip(self.varied_parameters, u):
            p.set_generator(uv)
        for prior in self._log_posterior.log_priors():
            prior.sample()
        return np.array([p.evaluate() for p in self.varied_parameters])


    def _par_to_u(self, par):
        """Internal function that used the CDF to translate from parameter space to u ∈ [0, 1)^D."""
        for p, pv in zip(self.varied_parameters, par):
            p.set(pv)
        for prior in self._log_posterior.log_priors():
            prior.compute_cdf()
        return np.array([p.evaluate_generator() for p in self.varied_parameters])


    @staticmethod
    def _sanitize_manual_input(data):
        """Helper function that converts all entries of a manual_constraint from numpy types to basic python types"""
        if np.issubdtype(type(data), int):
            return int(data)
        if np.issubdtype(type(data), float):
            return float(data)
        if np.issubdtype(type(data), str):
            return str(data)

        if type(data) is dict:
            return { k: Analysis._sanitize_manual_input(v) for k, v in data.items() }
        if np.issubdtype(type(data), list):  # data of type dict also matches but they are covered before
            return list(map(Analysis._sanitize_manual_input, data))

        # all valid cases are covered above
        raise ValueError(f"Unexpected entry type {type(data)} in manual_constraint")


    @staticmethod
    def _perplexity(weights):
        """Helper function that computes the perplexity of an array of weights.
           Non positive and NaN weights are neglected in the calculation."""
        # sum positive finite weights only
        weights_sum = np.sum(weights,
            where = np.logical_or(weights > 0, np.isfinite(weights)))
        if weights_sum <= 0:
            return 0.
        else:
            normalized_weights = weights / weights_sum
            # mask negative and nan weights
            normalized_weights = np.ma.MaskedArray(normalized_weights, copy=False,
                mask=(np.logical_or(normalized_weights <= 0, np.isnan(normalized_weights))))
            entropy = - np.sum(normalized_weights * np.log(normalized_weights.filled(1.0)))
            perplexity = np.exp(entropy) / len(normalized_weights)
            return perplexity

    @staticmethod
    def _ess(weights):
        """Helper function that computes the effective sample size of an array of weights"""
        return pypmc.tools.convergence.ess(weights)


    def clone(self):
        """Returns an independent instance of eos.Analysis."""
        return eos.Analysis(**self.init_args)


    def reset_parameters(self):
        """Resets the analysis parameters to their default values."""
        for p in eos.Parameters():
            self.parameters.set(p.name(), p.evaluate())


    def goodness_of_fit(self):
        """Returns a :class:`GoodnessOfFit` object that summarizes the quality of the fit for the current parameter point."""
        return eos.GoodnessOfFit(self._log_posterior)


    def optimize(self, start_point=None, rng=np.random.mtrand, **kwargs):
        r"""
        Optimize the log(posterior) and returns a best-fit-point summary.

        :param start_point: Parameter point from which to start the optimization, with the elements in the same order as in eos.Analysis.varied_parameters.
                            If set to "random", optimization starts at the random point in the space of the priors.
                            If not specified, optimization starts at the current parameter point.
        :type start_point: iterable, optional
        :param rng: Optional random number generator
        :param \**kwargs: Are passed to `scipy.optimize.minimize`

        """
        if str(start_point) == "random":
            # generate random uniform probabilities and store them in the generator values
            for param in self.varied_parameters:
                param.set_generator(rng.uniform())
            # sample the priors by using an inverse transform sampling based on the previously provided generator values
            for prior in self._log_posterior.log_priors():
                prior.sample()
            # use the sampled parameter point
            _start_point = [float(p) for p in self.varied_parameters]
        elif start_point is None:
            # use the current parameter point
            _start_point = [float(p) for p in self.varied_parameters]
        else:
            _start_point = np.array(start_point)

        scipy_opt_kwargs = { 'method': 'SLSQP', 'options': { 'ftol': 1.0e-13 } }
        # Update default values. If no keyword arguments are passed, kwargs is an empty dict
        scipy_opt_kwargs.update(kwargs)

        res = scipy.optimize.minimize(
            self.negative_log_pdf,
            self._par_to_u(_start_point),
            args=None,
            bounds=[(0.0, 1.0) for _ in self.varied_parameters],
            **scipy_opt_kwargs)

        if not res.success:
            eos.error('Optimization did not succeed')
            eos.error('  optimizer'' message reads: {}'.format(res.message))
            if str(start_point) == "random":
                eos.error('  random start point was used, perhaps try another method?')
        else:
            eos.success(f'Optimization goal achieved after {res.nfev} function evaluations')

        bfp = self._u_to_par(res.x)

        for p, v in zip(self.varied_parameters, bfp):
            p.set(v)

        return eos.BestFitPoint(self, bfp)


    def log_pdf(self, u, *args):
        """
        Adapter for use with external optimization software (e.g. pypmc) to aid when optimizing the log(posterior).

        :param u: Parameter point in u space, with the elements in the same order as in eos.Analysis.varied_parameters.
        :type u: iterable
        :param args: Dummy parameter (ignored)
        :type args: optional
        """
        self._u_to_par(u)

        try:
            return(self._log_posterior.evaluate())
        except RuntimeError as e:
            eos.error(f'encountered run time error ({e}) when evaluating log(posterior) in parameter point:')
            for p in self.varied_parameters:
                eos.error(f' - {p.name()}: {p.evaluate()}')
            return(-np.inf)


    def negative_log_pdf(self, u, *args):
        """
        Adapter for use with external optimization software (e.g. scipy.optimize.minimize) to aid when optimizing the log(posterior).

        :param u: Parameter point in u space, with the elements in the same order as in eos.Analysis.varied_parameters.
        :type u: iterable
        :param args: Dummy parameter (ignored)
        :type args: optional
        """
        return -self.log_pdf(u, *args)


    def sample_prior(self, N=1000, rng=np.random.mtrand):
        """
        Return prior samples of the parameters.

        Obtains random samples of the parameters based on their prior distributions. The code uses inverse transform sampling.

        :param N: Number of samples that shall be returned
        :type N: int
        :param rng: Optional random number generator

        :return: An iterable of the parameter samples of size N.
        """
        samples = []
        for _ in range(N):
            u_samples = rng.uniform(0.0, 1.0, len(self.varied_parameters))
            samples.append(self._u_to_par(u_samples))
        return np.array(samples)


    def sample(self, N=1000, stride=5, pre_N=150, preruns=3, cov_scale=0.1, observables=None, start_point=None, rng=np.random.mtrand,
               return_uspace=False):
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
        try:
            from tqdm.auto import tqdm
            progressbar = tqdm
        except ImportError:
            progressbar = lambda x, **kw: x

        ind_lower = np.array([ 0.0 for _ in self.varied_parameters])
        ind_upper = np.array([+1.0 for _ in self.varied_parameters])
        ind = pypmc.tools.indicator.hyperrectangle(ind_lower, ind_upper)

        log_target = pypmc.tools.indicator.merge_function_with_indicator(self.log_pdf, ind, -np.inf)

        # create initial covariance, assuming that each parameter's u-space value is uniformly distributed on [0, 1)
        sigma = np.diag([1.0 / 12.0 * cov_scale for _ in self.varied_parameters])   # 1 / 12 is the vairance U(0, 1)
        log_proposal = pypmc.density.gauss.LocalGauss(sigma)

        # create start point, if not provided or transform a provided start point to u space
        if start_point is None:
            start_point = np.array([rng.uniform(0.0, 1.0) for _ in self.varied_parameters])
        else:
            start_point = self._par_to_u(start_point)

        # create MC sampler
        sampler = pypmc.sampler.markov_chain.AdaptiveMarkovChain(log_target, log_proposal, start_point, save_target_values=True, rng=rng)

        # pre run to adapt markov chains
        eos.inprogress('Beginning preruns ...')
        for i in progressbar(range(0, preruns), desc="Preruns", leave=False):
            eos.info(f'Prerun {i} out of {preruns}')
            accept_count = sampler.run(pre_N)
            accept_rate  = accept_count / pre_N * 100
            eos.info(f'Prerun {i}: acceptance rate is {accept_rate:3.0f}%')
            sampler.adapt()
        sampler.clear()
        eos.completed(f'... completed {preruns} preruns')

        # obtain final samples
        eos.inprogress('Beginning main run ...')
        sample_total  = N * stride
        sample_chunk  = sample_total // 100
        sample_chunks = [sample_chunk for i in range(0, 99)]
        sample_chunks.append(sample_total - 99 * sample_chunk)
        for current_chunk in progressbar(sample_chunks, desc="Main run", leave=False):
            accept_count = accept_count + sampler.run(current_chunk)
        accept_rate  = accept_count / (N * stride) * 100
        eos.completed(f'... completed main run with acceptance rate {accept_rate:3.0f}%')

        # Transform from generator values in u space to the parameter values
        u_samples = sampler.samples[:][::stride]
        parameter_samples = np.apply_along_axis(self._u_to_par, 1, u_samples)
        weights = sampler.target_values[:][::stride, 0]

        if not observables:
            if return_uspace:
                return(parameter_samples, u_samples, weights)
            else:
                return(parameter_samples, weights)
        else:
            observable_samples = []
            for parameters in parameter_samples:
                for p, v in zip(self.varied_parameters, parameters):
                    p.set(v)

                observable_samples.append([o.evaluate() for o in observables])

            return(parameter_samples, weights, np.array(observable_samples))


    def sample_pmc(self, log_proposal, step_N=1000, steps=10, final_N=5000, rng=np.random.mtrand,
                    return_final_only=True, final_perplexity_threshold=1.0, weight_threshold=1e-10,
                    pmc_iterations=1, pmc_rel_tol=1e-10, pmc_abs_tol=1e-05, pmc_lookback=1):
        """
        Return samples of the parameters and log(weights), and a mixture density adapted to the posterior.

        Obtains random samples of the log(posterior) using adaptive importance sampling following
        the Population Monte Carlo approach with PyPMC.

        :param log_proposal: Initial gaussian mixture density that shall be adapted to the posterior density.
        :type log_proposal: pypmc.density.mixture.MixtureDensity
        :param step_N: Number of samples that shall be drawn in each adaptation step.
        :type step_N: int
        :param steps: Number of adaptation steps.
        :type steps: int
        :param final_N: Number of samples that shall be drawn after all adaptation steps.
        :type final_N: int
        :param rng: Optional random number generator (must be compatible with the requirements of pypmc.sampler.importance_sampler.ImportanceSampler)
        :param return_final_only: If set to True, only returns the samples and weights of the final sampling step, after all adaptations have finished.
        :param final_perplexity_threshold: Adaptations are stopped if the perplexity of the last adaptation step is above this threshold value.
        :param weight_threshold: Mixture components with a weight smaller than this threshold are pruned.
        :param pmc_iterations: (advanced) Maximum number of update of the PMC, changing this value may make the update unstable.
        :param pmc_rel_tol: (advanced) Relative tolerance of the PMC. If two consecutive values of the current density log-likelihood are relatively smaller than this value, the convergence is declared.
        :param pmc_abs_tol: (advanced) Absolute tolerance of the PMC. If two consecutive values of the current density log-likelihood are smaller than this value, the convergence is declared.
        :param pmc_lookback: (advanced) Use reweighted samples from the previous update steps when adjusting the mixture density.
            The parameter determines the number of update steps to "look back".
            The default value of 1 disables this feature, a value of 0 means that all previous steps are used.

        :return: A tuple of the parameters as array of length N = step_N * steps + final_N, the (linear) weights as array of length N, the posterior values as array of length N, and the
            final proposal function as pypmc.density.mixture.MixtureDensity.

        This method should be called after obtaining approximate samples of the
        log(posterior) by other means, e.g., by using :meth:`eos.Analysis.sample`.
        A possible (incomplete) example could look as follows:

        .. code-block:: python3

           from pypmc.mix_adapt.r_value import make_r_gaussmix
           chains = []
           for i in range(10):
               # run Markov Chains for your problem
               chain, _ = analysis.sample(...)
               chains.append(chain)

           # please consult the pypmc documentation for details on the call below
           proposal_density = make_r_gaussmix(chains, K_g=3, critical_r=1.1)

           # adapt the proposal to the posterior and obtain high-quality samples
           analysis.sample_pmc(proposal_density, ...)


        .. note::
           This method requires the PyPMC python module, which can be installed from PyPI.
        """
        try:
            from tqdm.auto import tqdm
            progressbar = tqdm
        except ImportError:
            progressbar = lambda x, **kw: x

        # create log_target
        ind_lower = np.array([ 0.0 for _ in self.varied_parameters])
        ind_upper = np.array([+1.0 for _ in self.varied_parameters])
        ind = pypmc.tools.indicator.hyperrectangle(ind_lower, ind_upper)

        log_target = pypmc.tools.indicator.merge_function_with_indicator(self.log_pdf, ind, -np.inf)

        # create PMC sampler
        sampler = pypmc.sampler.importance_sampling.ImportanceSampler(log_target, log_proposal, save_target_values=True, rng=rng)
        generating_components = []

        # list of proposals used to generate the samples. These proposals are not modified by `combine_weights`
        proposals = [sampler.proposal]

        # carry out adaptions
        eos.inprogress('Beggning PMC adaptations ...')
        for step in progressbar(range(steps), desc="Adaptations", leave=False):
            origins = sampler.run(step_N, trace_sort=True)
            generating_components.append(origins)

            # Compute the indicators for the current step
            last_weights = np.copy(sampler.weights[-1][:, 0])
            last_perplexity = self._perplexity(last_weights)
            last_ess = self._ess(last_weights)
            eos.info(f'Convergence diagnostics of the last samples after sampling in step {step}: '
                     f'perplexity = {last_perplexity}, ESS = {last_ess}')
            if last_perplexity < 0.05:
                eos.warn("Last step's perplexity is very low. This could possibly be improved by running "
                         "the markov chains that are used to form the initial PDF for a bit longer")

            # Use the samples of the last pmc_lookback steps to update the mixture
            samples = sampler.samples[-pmc_lookback:]
            weights = sampler.weights[-pmc_lookback:][:, 0]
            eos.info(f'Convergence diagnostics of all previous samples after sampling in step {step}: '
                     f'perplexity = {self._perplexity(weights)}, ESS = {self._ess(weights)}')

            # Reevaluate the weights of the previous pmc_lookback steps
            reevaluated_weights = pypmc.sampler.importance_sampling.combine_weights(
                 samples.reshape(-1, step_N, len(self.varied_parameters)),
                 weights.reshape(-1, step_N),
                 proposals[-pmc_lookback:]
                )[:][:,0]

            pmc = pypmc.mix_adapt.pmc.PMC(samples, sampler.proposal, reevaluated_weights, mincount=0, rb=True)
            # Update the proposal. Components with small weights are only pruned after the updates, this may slower the procedure but ensures that small weights are not removed to early.
            pmc.run(iterations=pmc_iterations, prune=0.0, rel_tol=pmc_rel_tol, abs_tol=pmc_abs_tol)
            sampler.proposal = pmc.density
            proposals.append(sampler.proposal)

            # Normalize the weights and remove components with a weight smaller than weight_threshold
            sampler.proposal.normalize()
            sampler.proposal.prune(threshold = weight_threshold)

            # stop adaptation if the perplexity of the last step is larger than the threshold
            if last_perplexity > final_perplexity_threshold:
                break
        eos.completed(f'... completed adaptations after {step} steps(s) with perplexity = {last_perplexity}')

        # draw final samples
        eos.inprogress(f'Beginning the final sampling ...')
        origins = sampler.run(final_N, trace_sort=True)
        generating_components.append(origins)

        # transform the samples back from u space to parameter space
        if return_final_only:
            # only returns the final_N final samples
            samples = np.apply_along_axis(self._u_to_par, 1, sampler.samples[:][-final_N:])
            weights = sampler.weights[:][-final_N:, 0]
            posterior_values = sampler.target_values[:][-final_N:]
        else:
            # returns all samples
            samples = np.apply_along_axis(self._u_to_par, 1, sampler.samples[:])
            weights = sampler.weights[:][:, 0]
            posterior_values = sampler.target_values[:]
        perplexity = self._perplexity(np.copy(weights))
        ess = self._ess(np.copy(weights))
        eos.completed(f'... completed final sampling with perplexity = {perplexity} and ESS = {ess}')

        return samples, weights, posterior_values, sampler.proposal


    def log_likelihood(self, p, *args):
        """
        Adapter for use with external sampling software (e.g. dynesty) to aid when sampling from the log(likelihood).

        :param p: Parameter point, with the elements in the same order as in eos.Analysis.varied_parameters.
        :type p: iterable
        :param args: Dummy parameter (ignored)
        :type args: optional
        """
        for p, pv in zip(self.varied_parameters, p):
            p.set(pv)

        try:
            return(self._log_likelihood.evaluate())
        except RuntimeError as e:
            eos.error(f'encountered run time error ({e}) when evaluating log(posterior) in parameter point:')
            for p in self.varied_parameters:
                eos.error(f' - {p.name()}: {p.evaluate()}')
            return(-np.inf)


    def _prior_transform(self, u):
        """
        Adapter for use with external sampling software to aid when sampling from the log(prior).

        :param u: The input probability point on the hypercube [0, 1)^D
        :type u: iterable
        """
        return self._u_to_par(u)


    def sample_nested(self, bound='multi', nlive=250, dlogz=1.0, maxiter=None, miniter=0, print_progress=True, print_function=None, seed=10, sample='auto'):
        """
        Return samples of the parameters.

        Obtains random samples of log(likelihood) using dynamic nested sampling with dynesty.

        :param bound: The option for bounding the target distribution. For valid values, see the dynesty documentation. Defaults to 'multi'.
        :type bound: str, optional
        :param nlive: The number of live points.
        :type nlive: int, optional
        :param dlogz: The relative tolerance for the remaining evidence. Defaults to 1.0.
        :type dlogz: float, optional
        :param maxiter: The maximum number of iterations. Iterations may stop earlier if the termination condition is reached.
        :type maxiter: int, optional
        :param miniter: The minimum number of iterations. Defaults to 0. Samples will be added until miniter is reached, even if the termination condition is reached earlier.
        :type miniter: int, optional
        :param print_function: The function used to print progress messages. Defaults to using a dynesty-based function.
        :type print_function: callable, optional
        :param seed: The seed used to initialize the Mersenne Twister pseudo-random number generator.
        :type seed: {None, int, array_like[ints], SeedSequence}, optional
        :param sample: The method used for sampling within the likelihood constraints. For valid values, see dynesty documentation. Defaults to 'auto'.
        :type sample: str, optional

        .. note::
           This method requires the dynesty python module, which can be installed from PyPI.
        """
        import dynesty, tqdm
        from functools import partial

        if print_function is None:
            print_function = partial(dynesty.results.print_fn, pbar=tqdm.tqdm())

        sampler = dynesty.DynamicNestedSampler(self.log_likelihood, self._prior_transform, len(self.varied_parameters), bound=bound, nlive=nlive, rstate = np.random.Generator(np.random.MT19937(seed)), sample=sample)
        sampler.run_nested(dlogz_init=dlogz, maxiter=maxiter, print_progress=print_progress, print_func=print_function)
        while sampler.results['niter'] < miniter:
            # using mode='full' ensures sampling from the entire posterior
            sampler.add_batch(mode='full', dlogz_init=dlogz, maxiter=maxiter, print_progress=print_progress, print_func=print_function)
        return sampler.results


    def _repr_html_(self):
        result = r'''
        <table>
            <colgroup>
                <col width="50%" id="qn"    style="min-width: 200px">
                <col width="20%" id="min"   style="min-width: 100px">
                <col width="20%" id="max"   style="min-width: 100px">
                <col width="30%" id="type"  style="min-width: 100px">
            </colgroup>
            <thead>
                <tr>
                    <th colspan="4">priors</th>
                </tr>
                <tr>
                    <th>qualified name</th>
                    <th>min</th>
                    <th>max</th>
                    <th>type</th>
                </tr>
            </thead>
            <tbody>'''

        for p in self.init_args['priors']:
            result += fr'''
                <tr>
                    <td><tt>{p['parameter']}</tt></td>
                    <td>{p['min'] if 'min' in p else None}</td>
                    <td>{p['max'] if 'min' in p else None}</td>
                    <td>{p['type']}</td>
                </tr>
            '''
        result += r'''
            </tbody>
            </table>
            <table>
            <thead>
                <tr>
                    <th colspan="2">likelihood</th>
                </tr>
                <tr>
                    <th>qualified name</th>
                    <th>&num; obs.</th>
                </tr>
            </thead>
            <tbody>
        '''
        for c in self._log_likelihood:
            result += fr'''
                <tr>
                    <td><tt>{c.name()}</tt></td>
                    <td>{len(list(c.observables()))}</td>
                </tr>
            '''
        result += r'''
            </tbody>
            </table>
        '''

        return(result)
