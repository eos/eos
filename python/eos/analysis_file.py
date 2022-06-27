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

        if 'steps' not in input_data:
            self._steps = []
        else:
            self._steps = input_data['steps']


    def analysis(self, _posterior):
        """Create an eos.Analysis object for the named posterior."""
        if _posterior not in self._posteriors:
            raise RuntimeError('Cannot create analysis for unknown posterior: \'{}\''.format(_posterior))

        posterior = self._posteriors[_posterior]

        prior = []
        for p in posterior['prior']:
            prior.extend(self._priors[p]['parameters'])

        likelihood = []
        manual_constraints = {}
        for lh in posterior['likelihood']:
            LH_ALLOWED_KEYS = { 'name', 'constraints', 'manual_constraints' }
            for key in self._likelihoods[lh]:
                if key not in LH_ALLOWED_KEYS:
                    raise KeyError(f"Unsupported key in 'likelihoods['{lh}']': {key}")

            if 'constraints' not in self._likelihoods[lh] and 'manual_constraints' not in self._likelihoods[lh]:
                raise KeyError(f'Missing entry in \'likelihoods[\'{lh}\']\': neither \'constraints\' nor \'manual_constraints\' is provided')

            likelihood.extend(self._likelihoods[lh]['constraints'] if 'constraints' in self._likelihoods[lh] else [])
            manual_constraints.update(self._likelihoods[lh]['manual_constraints'] if 'manual_constraints' in self._likelihoods[lh] else {})

        global_options = posterior['global_options'] if 'global_options' in posterior else {}
        fixed_parameters = posterior['fixed_parameters'] if 'fixed_parameters' in posterior else {}

        return eos.Analysis(prior, likelihood, global_options,
                            manual_constraints=manual_constraints,
                            fixed_parameters=fixed_parameters)


    def observables(self, _prediction, parameters):
        """Creates a list of eos.Observable objects for the named set of predictions."""
        if _prediction not in self.predictions:
            raise RuntimeError('Cannot create observables for unknown set of predictions: \'{}\''.format(_prediction))

        prediction = self.predictions[_prediction]
        options = eos.Options(**prediction['global_options'])
        observables = [eos.Observable.make(
            o['name'],
            parameters,
            eos.Kinematics(**(o['kinematics'] if 'kinematics' in o else {})),
            options
        ) for o in prediction['observables']]

        if None in observables:
            unknown_observables = set()
            for p, o in zip(prediction['observables'], observables):
                if o is None:
                    unknown_observables.add(p['name'])
            raise RuntimeError('Prediction \'{}\' contains unknown observable names: {}'.format(_prediction, unknown_observables))

        return observables


    @staticmethod
    def _sanitize_params(command, params):
        """Helper functions to sanitize parameters for individual steps loaded from a raw YAML file."""
        general_params_maps = {
            'b': 'base_directory'
        }
        specific_params_map = {
            # find-clusters
            ('find-clusters', 't'): 'threshold',
            ('find-clusters', 'c'): 'K_g', ('find-clusters', 'clusters-per-group'): 'K_g',
            # predict-observables
            ('predict-observables', 'B'): 'begin',
            ('predict-observables', 'E'): 'end',
            # sample-mcmc
            ('sample-mcmc', 'n'): 'pre_N', ('sample-mcmc', 'prerun-samples'): 'pre_N',
            ('sample-mcmc', 'p'): 'preruns',
            ('sample-mcmc', 'S'): 'stride',
            # sample-pmc
            ('sample-pmc', 'n'): 'step_N', ('sample-pmc', 'step-samples'): 'step_N',
            ('sample-pmc', 's'): 'steps',
            ('sample-pmc', 'N'): 'final_N', ('sample-pmc', 'final-samples'): 'final_N'

        }
        return {
            specific_params_map[(command, k)] if (command, k) in specific_params_map
            else general_params_map[k] if k in general_params_map
            else k
            for k, v in params.items()
        }

    def steps(self, base_directory=None):
        """Runs predefined analysis steps recorded in the analysis file."""
        from .tasks import _tasks
        from collections import ChainMap
        from inspect import signature
        from itertools import chain, product
        import networkx as nx

        def _expand_arguments(step_name, task, arguments):
            task_sig = signature(_tasks[task])

            task_required_args = set()
            for n, p in task_sig.parameters.items():
                if p.default != p.empty:
                    continue
                task_required_args.add(n)
            for n in task_required_args:
                if n in arguments.keys():
                    continue
                raise ValueError(f'Task "{task}" in step "{step_name}" requires mandatory argument \'{n}\', which is not provided')

            outer_list = []
            for key, value in arguments.items():
                if key not in task_sig.parameters:
                    eos.warn(f'Task "{task}" in step "{step_name}" does not expect argument "{key}", possibly misspelled?')
                    continue

                param_type = task_sig.parameters[key].annotation

                # handle special cases
                if param_type is int:
                    # integer range?
                    if type(value) is str:
                        values = list(chain.from_iterable(range(int(v.split("-")[0]),int(v.split("-")[-1])+1) for v in value.split(",")))
                        outer_list.append([{key: v} for v in values])
                    else:
                        # append
                        outer_list.append([{key: value}])
                # otherwise append w/o modification
                else:
                    # append
                    outer_list.append([{key: value}])

            result = []
            for r in list(product(*outer_list)):
                _arguments = {}
                for d in r:
                    _arguments.update(d)

                # replace format strings
                arguments = {}
                for key, value in _arguments.items():
                    if type(value) is str:
                        arguments.update({key: value.format(**_arguments)})
                    else:
                        arguments.update({key: value})

                result.append(arguments)

            return result

        def _handle_single_step(step_desc, step_arguments):
            step_name = step_desc['name']
            tasks     = step_desc['tasks']

            for task_desc in tasks:
                task = task_desc['task']
                if task not in _tasks.keys():
                    raise ValueError(f'Unknown task \"{task}\" encountered in step {step_name}')

                arguments = { 'analysis_file': self }
                for key, value in step_arguments.items():
                    if key in _tasks.keys():
                        continue

                    arguments[key] = value

                if 'arguments' in task_desc:
                    arguments.update(task_desc['arguments'])

                if task in step_arguments:
                    arguments.update(step_arguments[task])

                if base_directory:
                    arguments.update({ 'base_directory': base_directory })

                return [(step_name, step_desc, task, a) for a in _expand_arguments(step_name, task, arguments)]

        # main part starts here

        # check inputs
        for step_idx, step_desc in enumerate(self._steps):
            if type(step_desc) is not dict:
                raise ValueError(f'Description of step #{step_idx} is not a key/value map')

            if 'name' not in step_desc:
                raise ValueError(f'Description of step #{step_idx} requires a name')

            if 'tasks' not in step_desc:
                raise ValueError('Description of step #{step_idx} ({step["name"]}) requires a list of tasks')

            if 'iterations' in step_desc:
                if type(step_desc['iterations']) is not list:
                    raise ValueError(f'Description of step #{step_idx} ({step["name"]} expects a list of iterations')

        # resolve dependencies
        graph = nx.MultiDiGraph()
        graph.add_nodes_from([step['name'] for step in self._steps])
        for step in self._steps:
            if not 'depends-on' in step:
                continue

            for dep in step['depends-on']:
                graph.add_edge(dep, step['name'])
        steps = list(nx.topological_sort(graph))

        # return steps and theirs tasks in the order imposed by the dependency graph
        result = []
        for step in steps:
            step_desc = { sd['name']: sd for sd in self._steps }[step]
            iterations = [{}]
            if 'iterations' in step_desc:
                iterations = step_desc['iterations']

            for iteration_arguments in iterations:
                result.extend(_handle_single_step(step_desc, iteration_arguments))

        return result


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


    def _repr_html_(self):
        result = r'''
        <table>
            <colgroup>
                <col width="100%" id="posterior"    style="min-width: 200px">
            </colgroup>
            <thead>
                <tr>
                    <th>posteriors</th>
                </tr>
            </thead>
            <tbody>
        '''
        for p in self._posteriors:
            result += fr'''
                <tr>
                    <td>{p}</td>
                </tr>
            '''

        result += r'''
            </tbody>
        </table>
        '''

        return result
