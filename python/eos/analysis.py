#!/usr/bin/python
# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2018 Danny van Dyk
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
import scipy

class BestFitPoint:
    """
    Represents the best-fit point of a Bayesian analysis undertaken with the eos.Analysis class.
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
    :param priors: The priors for this analysis as a list of prior descriptions. See :ref priors: for what consitutes a valid prior description.
    :type priors: iterable
    :param likelihood: The likelihood as a list of individual constraints from the internal data base of experimental and theoretical constraints. See :ref constraints: for a complete list of constraints.
    :type likelihood: iterable
    """

    def __init__(self, priors, likelihood, global_options=None):
        """Constructor."""
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
            prior_type = prior['type']
            if 'uniform' == prior_type or 'flat' == prior_type:
                self.log_posterior.add(eos.LogPrior.Flat(self.parameters, parameter, eos.ParameterRange(minv, maxv)), False)
            elif 'gauss' == prior_type or 'gaussian' == prior_type:
                central = prior['central']
                sigma = prior['sigma']
                self.log_posterior.add(
                    eos.LogPrior.Gauss(
                        self.parameters, parameter, eos.ParameterRange(minv, maxv),
                        central - sigma, central, central + sigma
                    ),
                    False)
            else:
                raise ValueError('Unknown prior type \'{}\''.format(prior_type))

            self.bounds.append((minv, maxv))
            self.varied_parameters.append(self.parameters[parameter])

        # create the likelihood
        for constraint_name in likelihood:
            constraint = eos.Constraint.make(constraint_name, self.global_options)
            self.log_likelihood.add(constraint)


    def goodness_of_fit(self):
        """Returns a goodness-of-fit summary for the current parameter point."""
        return eos.GoodnessOfFit(self.log_posterior)


    def optimize(self, start_point=None):
        """
        Optimize the log(posterior) and returns a best-fit-point summary.

        :param start_point: Parameter point from which to start the optimization, with the elements in the same order as in eos.Analysis.varied_parameters. If not specified, optimization starts at the current parameter point.
        :param start_point: iterable, optional
        """

        if start_point == None:
            start_point = [float(p) for p in self.varied_parameters]

        res = scipy.optimize.minimize(
            self.negative_log_pdf,
            start_point,
            args=None,
            bounds=self.bounds)

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


    def sample(self, N=1000, stride=5, pre_N=150, pre_runs=3, observables=None):
        import pypmc
        import numpy as np
        import random

        ind_lower = np.array([bound[0] for bound in self.bounds])
        ind_upper = np.array([bound[1] for bound in self.bounds])
        ind = pypmc.tools.indicator.hyperrectangle(ind_lower, ind_upper)

        log_target = pypmc.tools.indicator.merge_function_with_indicator(self.log_pdf, ind, -np.inf)

        # create initial covariance
        sigma = np.diag([np.square(bound[1] - bound[0]) / 12 for bound in self.bounds])
        log_proposal = pypmc.density.gauss.LocalGauss(sigma)

        # create start point
        u = np.array([random.uniform(0.0, 1.0) for j in range(0, len(ind_lower))])
        ubar = 1.0 - u
        start_point = ubar * ind_upper + u * ind_lower

        # create MC sampler
        sampler = pypmc.sampler.markov_chain.AdaptiveMarkovChain(log_target, log_proposal, start_point, save_target_values=True)

        # pre run to adapt markov chains
        for i in range(0, pre_runs):
            sampler.run(pre_N)
            sampler.adapt()
        sampler.clear()

        # obtain final samples
        sampler.run(N * stride)

        parameter_samples = sampler.samples[-1][::stride]
        weights = sampler.target_values[-1][::stride, 0]

        if not observables:
            return(parameter_samples, weights)
        else:
            observable_samples = []
            for parameters in parameter_samples:
                for p, v in zip(self.varied_parameters, parameters):
                    p.set(v)

                observable_samples.append([o.evaluate() for o in observables])

            return(parameter_samples, weights, np.array(observable_samples))


