#!/usr/bin/python
# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2020-2023 Danny van Dyk
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
            raise RuntimeError(f'Cannot load analysis file: \'{analysis_file}\' does not exist')

        if not os.path.isfile(analysis_file):
            raise RuntimeError(f'Cannot load analysis file: \'{analysis_file}\' is not a file')

        instructions = None
        with open(analysis_file) as input_file:
            input_data = yaml.load(input_file, Loader=yaml.SafeLoader)

        if 'priors' not in input_data:
            raise RuntimeError('Cannot load analysis file: need at least one prior component')

        self._priors = { prior['name'] : dict(prior) for prior in input_data['priors'] }

        if 'likelihoods' not in input_data:
            raise RuntimeError('Cannot load analysis file: need at least one likelihood component')

        self._likelihoods = { lh['name'] : dict(lh) for lh in input_data['likelihoods'] }

        if 'posteriors' not in input_data:
            raise RuntimeError('Cannot load analysis file: need at least one posterior')

        self._posteriors = { posterior['name'] : dict(posterior) for posterior in input_data['posteriors'] }

        # Optional: provide a list of observables for posterior prediction
        if 'predictions' not in input_data:
            self._predictions = []
        else:
            self._predictions = { pred['name'] : dict(pred) for pred in input_data['predictions'] }

        # Optional: insert custom observables using the expression parser
        if 'observables' in input_data:
            for name in input_data['observables']:

                required_keys = {'latex', 'unit', 'options', 'expression'}
                optional_keys = {'alias_of'}
                provided_keys = input_data['observables'][name].keys()

                missing_keys  = required_keys - provided_keys
                if missing_keys:
                    raise KeyError(f'Missing keys for observable { name }: { missing_keys }')

                ignored_keys = provided_keys - required_keys - optional_keys
                if ignored_keys:
                    eos.warn(f'Ignoring unknown keys for observable { name }: { ignored_keys }')

                obs = input_data['observables'][name]
                options = eos.Options(**obs['options'])
                eos.Observables().insert(name, obs['latex'], eos.Unit(obs['unit']), options, obs['expression'])
                eos.info(f'Successfully inserted observable: { name }')

        # Optional: declare custom parameters before creating an analysis
        if 'parameters' in input_data:
            for name in input_data['parameters']:

                required_keys = {'latex', 'unit', 'central', 'min', 'max'}
                provided_keys = input_data['parameters'][name].keys()

                missing_keys  = required_keys - provided_keys
                if missing_keys:
                    raise KeyError(f'Missing keys for parameter { name }: { missing_keys }')

                ignored_keys = provided_keys - required_keys
                if ignored_keys:
                    eos.warn(f'Ignoring unknown keys for parameter { name }: { ignored_keys }')

                param = input_data['parameters'][name]
                try:
                    qn = eos.QualifiedName(name)
                    unit = eos.Unit(param['unit'])
                    central, min, max = float(param['central']), float(param['min']), float(param['max'])
                    id = eos.Parameters.declare(qn, param['latex'], unit, central, min, max)
                    eos.info(f'Successfully declared parameter: {qn}')
                    if 'alias_of' in param:
                        for aliased_qn in param['alias_of']:
                            eos.Parameters.redirect(eos.QualifiedName(aliased_qn), id)
                            eos.info(f'Successfully declared alias: {aliased_qn} -> {qn}')
                except eos.Exception as e:
                    raise ValueError(f'Unexpected value encountered in description of parameter \'{name}\': {e}')

        if 'steps' not in input_data:
            self._steps = []
        else:
            self._steps = input_data['steps']


    def analysis(self, _posterior):
        """Create an eos.Analysis object for the named posterior."""
        if _posterior not in self._posteriors:
            raise RuntimeError(f'Cannot create analysis for unknown posterior: \'{_posterior}\'')

        posterior = self._posteriors[_posterior]

        prior = []
        for p in posterior['prior']:
            descriptions = []

            if 'descriptions' in self._priors[p]:
                if 'parameters' in self._priors[p]:
                    eos.error(f'Both \'descriptions\' and \'parameters\' are provided for prior component \'{p}\', ignoring legacy support for \'parameters\'')

                descriptions = self._priors[p]['descriptions']
            elif 'parameters' in self._priors[p]:
                eos.warn(f'\'parameters\' is in the description of prior component \'{p}\', use \'descriptions\' instead')
                descriptions = self._priors[p]['parameters']
            else:
                raise KeyError(f'Missing entry in prior component \'{p}\': \'descriptions\' is not provided')

            prior.extend(descriptions)

        parameters = eos.Parameters()

        likelihood = []
        manual_constraints = {}
        external_likelihood = []
        for lh in posterior['likelihood']:
            LH_ALLOWED_KEYS = { 'name', 'constraints', 'manual_constraints', 'pyhf'}
            for key in self._likelihoods[lh]:
                if key not in LH_ALLOWED_KEYS:
                    raise KeyError(f"Unsupported key in 'likelihoods['{lh}']': {key}")

            if 'constraints' not in self._likelihoods[lh] and 'manual_constraints' not in self._likelihoods[lh] and 'pyhf' not in self._likelihoods[lh]:
                raise KeyError(f'Missing entry in \'likelihoods[\'{lh}\']\': need exactly one of: \'constraints\', \'manual_constraints\' or \'pyhf\'')

            likelihood.extend(self._likelihoods[lh]['constraints'] if 'constraints' in self._likelihoods[lh] else [])
            manual_constraints.update(self._likelihoods[lh]['manual_constraints'] if 'manual_constraints' in self._likelihoods[lh] else {})

            # create a pyhf likelihood
            if 'pyhf' in self._likelihoods[lh]:
                workspace = self._likelihoods[lh]['pyhf']['file']
                parameter_map = self._likelihoods[lh]['pyhf']['parameter_map'] if 'parameter_map' in self._likelihoods[lh]['pyhf'] else []
                cache = eos.ObservableCache(parameters)

                # extend prior
                pyhf_priors = eos.PyhfLogLikelihood.factory(cache, workspace, parameter_map).priors()
                existing_priors = [prior['parameter'] for prior in prior]
                # only add missing pyhf priors
                for pyhf_prior in pyhf_priors:
                    if pyhf_prior['parameter'] not in existing_priors:
                        eos.info(f'pyhf workspace parameter {pyhf_prior["parameter"]} added to prior; manually specify this prior to overwrite settings')
                        prior.append(pyhf_prior)

                # create likelihood block
                llh_block = eos.LogLikelihoodBlock.External(
                                cache,
                                lambda cache: eos.PyhfLogLikelihood.factory(cache, workspace, parameter_map)
                                )

                external_likelihood.extend([llh_block])

        global_options = posterior['global_options'] if 'global_options' in posterior else {}
        fixed_parameters = posterior['fixed_parameters'] if 'fixed_parameters' in posterior else {}

        return eos.Analysis(prior, likelihood,
                            external_likelihood=external_likelihood,
                            global_options=global_options,
                            manual_constraints=manual_constraints,
                            fixed_parameters=fixed_parameters,
                            parameters=parameters)


    def observables(self, _posterior, _prediction, parameters):
        """Creates a list of eos.Observable objects for the named set of posterior and predictions."""
        if _posterior not in self._posteriors:
            raise RuntimeError(f'Cannot create observables for unknown posterior: \'{_prediction}\'')
        if _prediction not in self.predictions:
            raise RuntimeError(f'Cannot create observables for unknown set of predictions: \'{_prediction}\'')

        posterior = self._posteriors[_posterior]
        prediction = self.predictions[_prediction]

        global_options = posterior['global_options'] if 'global_options' in posterior else {}
        fixed_parameters = posterior['fixed_parameters'] if 'fixed_parameters' in posterior else {}

        if 'global_options' in prediction:
            global_options.update(prediction['global_options'])

        if 'fixed_parameters' in prediction:
            fixed_parameters.update(prediction['fixed_parameters'])

        options = eos.Options(**global_options)
        for (p, v) in fixed_parameters.items():
            parameters.set(p, v)

        for o in prediction['observables']:
            options_part = eos.QualifiedName(o['name']).options_part()
            for key, value in options_part:
                if key in global_options and global_options[key] != value:
                    eos.error(f'Global option {key}={global_options[key]} overrides option part specification {key}={value} for observable {o["name"]} in prediction {_prediction}.')

        observables = []
        for o in prediction['observables']:
            if 'kinematics' not in o:
                kinematics = [{}]
            elif type(o['kinematics']) == dict:
                kinematics = [o['kinematics']]
            else:
                kinematics = [k for k in o['kinematics']]

            for k in kinematics:
                observables.append(eos.Observable.make(
                    o['name'],
                    parameters,
                    eos.Kinematics(k),
                    options
                ))

        if None in observables:
            unknown_observables = set()
            for p, o in zip(prediction['observables'], observables):
                if o is None:
                    unknown_observables.add(p['name'])
            raise RuntimeError(f'Prediction \'{_prediction}\' contains unknown observable names: {unknown_observables}')

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
                    <th style="text-align: center">priors</th>
                </tr>
            </thead>
            <tbody>
        '''
        for p in self._priors:
            result += fr'''
                <tr>
                    <td style="text-align: left">{p}</td>
                </tr>
            '''

        result += r'''
            </tbody>
            <thead>
                <tr>
                    <th style="text-align: center">likelihoods</th>
                </tr>
            </thead>
            <tbody>
        '''
        for l in self._likelihoods:
            result += fr'''
                <tr>
                    <td style="text-align: left">{l}</td>
                </tr>
            '''

        result += r'''
            </tbody>
            <thead>
                <tr>
                    <th style="text-align: center">posteriors</th>
                </tr>
            </thead>
            <tbody>
        '''

        for p in self._posteriors:
            result += fr'''
                <tr>
                    <td style="text-align: left">{p}</td>
                </tr>
            '''

        result += r'''
            </tbody>
        </table>
        '''

        return result
