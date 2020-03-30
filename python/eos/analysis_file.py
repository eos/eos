#!/usr/bin/python
# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2020 Danny van Dyk
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
from logging import info, warn
import os
import yaml

class AnalysisFile:
    """Represents a collection of statistical analyses and their building blocks.

    :param analysis_file: The path to the file to be parsed.
    :type analysis_file: str
    """

    def __init__(self, analysis_file):
        """Constructor."""

        self.analysis_file = analysis_file

        if not os.path.exists(analysis_file):
            raise RuntimeError('Cannot load analysis file: \'{}\' does not exist'.format(analysis_file))

        if not os.path.isfile(analysis_file):
            raise RuntimeError('Cannot load analysis file: \'{}\' is not a file'.format(analysis_file))

        instructions = None
        with open(analysis_file) as input_file:
            input_data = yaml.load(input_file, Loader=yaml.SafeLoader)

        if 'priors' not in input_data:
            raise RuntimeError('Cannot load analysis file: need at least one prior')

        self._priors = { prior['name'] : dict(prior) for prior in input_data['priors'] }

        if 'likelihoods' not in input_data:
            self._likelihoods = []
        else:
            self._likelihoods = { lh['name'] : dict(lh) for lh in input_data['likelihoods'] }

        if 'posteriors' not in input_data:
            raise RuntimeError('Cannot load analysis file: need at least one posterior')

        self._posteriors = { posterior['name'] : dict(posterior) for posterior in input_data['posteriors'] }

        if 'predictions' not in input_data:
            self._predictions = []
        else:
            self._predictions = { pred['name'] : dict(pred) for pred in input_data['predictions'] }


    def analysis(self, _posterior):
        """Create an eos.Analysis object for the named posterior."""
        if _posterior not in self._posteriors:
            raise RuntimeError('Cannot create analysis for unknown posterior: \'{}\''.format(_posterior))

        posterior = self._posteriors[_posterior]

        prior = []
        for p in posterior['prior']:
            prior.extend(self._priors[p]['parameters'])

        likelihood = []
        for lh in posterior['likelihood']:
            likelihood.extend(self._likelihoods[lh]['constraints'])

        global_options = posterior['global_options'] if 'global_options' in posterior else None

        return eos.Analysis(prior, likelihood, global_options)


    def observables(self, _prediction, parameters):
        """Creates a list of eos.Observable objects for the named set of predictions."""
        if _prediction not in self.predictions:
            raise RuntimeError('Cannot create observables for unknown set of predictions: \'{}\''.format(_prediction))

        prediction = self.predictions[_prediction]
        options = eos.Options(**prediction['global_options'])
        return [eos.Observable.make(
            o['name'],
            parameters,
            eos.Kinematics(**(o['kinematics'] if 'kinematics' in o else {})),
            options
        ) for o in prediction['observables']]


    @property
    def priors(self):
        """Returns a list of all priors recorded in the analysis file."""

        return self._priors


    @property
    def likelihoods(self):
        """Returns a list of all likelihoods recorded in the analysis file."""

        return self._likelihoods


    @property
    def posteriors(self):
        """Returns a list of all posteriors recorded in the analysis file."""

        return self._posteriors

    @property
    def predictions(self):
        """Returns a list of all predictions recorded in the analysis file."""

        return self._predictions
