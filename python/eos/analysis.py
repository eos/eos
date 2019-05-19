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
            result += '<tr><td>{n}</td><td>{v:4.2g}</td></tr>'.format(n=name, v=v)
        result += '</table>'

        return(result)


class Analysis:
    """
    Describes a Bayesian analysis in terms of a set of parameters, a log(likelihood),
    and potentially one or more log(prior)s.
    """
    def __init__(self, **kwargs):
        self.parameters = eos.Parameters.Defaults()
        self.global_options = eos.Options()
        self.log_likelihood = eos.LogLikelihood(self.parameters)
        self.log_posterior = eos.LogPosterior(self.log_likelihood)
        self.varied_parameters = []

        if 'global-options' in kwargs:
            global_options = kwargs['global-options']
            for key, value in global_options.items():
                self.global_options.set(key, value)

        if 'priors' in kwargs:
            priors = kwargs['priors']
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

                self.varied_parameters.append(self.parameters[parameter])

        if 'likelihood' in kwargs:
            constraints = kwargs['likelihood']
            for constraint_name in constraints:
                constraint = eos.Constraint.make(constraint_name, self.global_options)
                self.log_likelihood.add(constraint)

    def goodness_of_fit(self):
        return eos.GoodnessOfFit(self.log_posterior)

    def optimize(self, start_point=None):
        if start_point == None:
            start_point = [float(p) for p in self.varied_parameters]

        res = scipy.optimize.minimize(
            eos.Analysis.target_function, 
            start_point,
            args=(self,))

        for p, v in zip(self.varied_parameters, res.x):
            p.set(v)

        return eos.BestFitPoint(self, res.x)

    """
    Adapter for use with SciPy to aid when optimizing the log(posterior).
    """
    @staticmethod
    def target_function(x, *args):
        if type(args[0]) != eos.Analysis:
            raise ValueError('Expecting first elements of \'args\' to be of type \'eos.Analysis\', got \'{}\''.format(str(type(args[0]))))

        analysis = args[0]

        if len(x) != len(analysis.varied_parameters):
            raise ValueError('Analysis expects {na} parameters, but only {nx} parameters provided'.format(na=len(analysis.varied_parameters),nx=len(x)))

        for p, v in zip(analysis.varied_parameters, x):
            p.set(v)

        return -analysis.log_posterior.evaluate()


