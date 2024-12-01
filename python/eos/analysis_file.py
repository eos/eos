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
import sys
import yaml
from dataclasses import asdict
from eos.analysis_file_description import PriorComponent, LikelihoodComponent, PosteriorDescription, \
                                       PredictionDescription, ObservableComponent, ParameterComponent, \
                                       StepComponent, PriorDescription

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

        with open(analysis_file) as input_file:
            self.input_data = yaml.safe_load(input_file)

        if 'priors' not in self.input_data:
            raise RuntimeError('Cannot load analysis file: need at least one prior component')

        self._priors = { p["name"]: PriorComponent.from_dict(**p) for p in self.input_data['priors'] }

        if 'likelihoods' not in self.input_data:
            eos.warn('No likelihood components found in analysis file')

        self._likelihoods = { ll["name"]: LikelihoodComponent.from_dict(**ll) for ll in self.input_data['likelihoods'] }

        if 'posteriors' not in self.input_data:
            raise RuntimeError('Cannot load analysis file: need at least one posterior')

        self._posteriors = { p["name"]: PosteriorDescription.from_dict(**p) for p in self.input_data['posteriors'] }
        # Check that the priors and likelihoods referenced by the posteriors are defined
        for pc in self._posteriors.values():
            for p in pc.prior:
                if p not in self._priors:
                    raise RuntimeError(f'Posterior \'{pc.name}\' references prior \'{p}\' which is not defined')
            for l in pc.likelihood:
                if l not in self._likelihoods:
                    raise RuntimeError(f'Posterior \'{pc.name}\' references likelihood \'{l}\' which is not defined')

        # Optional: provide a list of observables for posterior prediction
        if 'predictions' not in self.input_data:
            self._predictions = []
        else:
            self._predictions = { p["name"]: PredictionDescription.from_dict(**p) for p in self.input_data['predictions'] }

        # Optional: insert custom observables using the expression parser
        if 'observables' not in self.input_data:
            self._obs = []
        else:
            eos.inprogress('Inserting custom observables ...')
            self._obs = [ObservableComponent.from_dict(name=n, **d) for n,d in self.input_data['observables'].items()]
            for o in self._obs:
                eos.Observables().insert(o.name, o.latex, eos.Unit(o.unit), eos.Options(**o.options), o.expression)
                eos.info(f'Inserted observable: { o.name }')
            eos.completed(f'... finished inserting {len(self._obs)} custom observables')

        if 'parameters' not in self.input_data:
            self._params = []
        else:
            eos.inprogress('Declaring custom parameters ...')
            self._params = [ParameterComponent.from_dict(name=n, **d) for n, d in self.input_data["parameters"].items()]
            for p in self._params:
                try:
                    qn = eos.QualifiedName(p.name)
                    id = eos.Parameters.declare(qn, p.latex, eos.Unit(p.unit), p.central, p.min, p.max)
                    eos.info(f'Declared parameter: {qn}')
                    if p.alias_of:
                        for aliased_qn in p.alias_of:
                            eos.Parameters.redirect(eos.QualifiedName(aliased_qn), id)
                            eos.info(f'Created parameter alias: {aliased_qn} -> {qn}')
                except eos.Exception as e:
                    raise ValueError(f'Unexpected value encountered in description of parameter \'{p.name}\': {e}')
            eos.completed(f'... finished declaring {len(self._params)} custom parameters')

        if 'steps' not in self.input_data:
            self._steps = []
        else:
            self._steps = [StepComponent.from_dict(**s) for s in self.input_data['steps']]

    def analysis(self, _posterior):
        """Create an eos.Analysis object for the named posterior."""
        if _posterior not in self._posteriors:
            raise RuntimeError(f'Cannot create analysis for unknown posterior: \'{_posterior}\'')

        posterior = self._posteriors[_posterior]

        prior = []
        for p in posterior.prior:
            prior.extend(self._priors[p].descriptions)

        parameters = eos.Parameters()

        likelihood = []
        manual_constraints = []
        external_likelihood = []
        for lh in posterior.likelihood:
            likelihood.extend(self._likelihoods[lh].constraints)
            manual_constraints.extend(self._likelihoods[lh].manual_constraints)

            # create a pyhf likelihood
            if self._likelihoods[lh].pyhf:
                workspace = self._likelihoods[lh].pyhf.file # TODO: Convert this path relative to analysis_file?
                parameter_map = self._likelihoods[lh].pyhf.parameter_map
                cache = eos.ObservableCache(parameters)

                # extend prior
                pyhf_priors = eos.PyhfLogLikelihood.factory(cache, workspace, parameter_map).priors()
                existing_priors = [p.parameter for p in prior]
                # only add missing pyhf priors
                for pyhf_prior in pyhf_priors:
                    if pyhf_prior['parameter'] not in existing_priors:
                        eos.info(f'pyhf workspace parameter {pyhf_prior["parameter"]} added to prior; manually specify this prior to overwrite settings')
                        prior.append(PriorDescription.from_dict(**pyhf_prior))

                # create likelihood block
                llh_block = eos.LogLikelihoodBlock.External(
                                cache,
                                lambda cache: eos.PyhfLogLikelihood.factory(cache, workspace, parameter_map)
                                )

                external_likelihood.extend([llh_block])

        global_options = posterior.global_options
        fixed_parameters = posterior.fixed_parameters

        #  Convert back to dictionaries (for now) as that is what Analysis expects
        prior = [ asdict(pc) for pc in prior]
        likelihood = [ d["constraint"] for d in (asdict(lc) for lc in likelihood) ]
        manual_constraints = { d["name"]: d["info"] for d in (asdict(mc) for mc in manual_constraints) }

        return eos.Analysis(prior, likelihood,
                            external_likelihood=external_likelihood,
                            global_options=global_options,
                            manual_constraints=manual_constraints,
                            fixed_parameters=fixed_parameters,
                            parameters=parameters)


    def observables(self, _posterior, _prediction, parameters):
        """Creates a list of eos.Observable objects for the named set of posterior and predictions."""
        if _posterior not in self._posteriors:
            raise RuntimeError(f'Cannot create observables for unknown posterior: \'{_posterior}\'')
        if _prediction not in self.predictions:
            raise RuntimeError(f'Cannot create observables for unknown set of predictions: \'{_prediction}\'')

        posterior = self._posteriors[_posterior]
        prediction = self.predictions[_prediction]

        global_options = posterior.global_options
        fixed_parameters = posterior.fixed_parameters

        global_options.update(prediction.global_options)
        fixed_parameters.update(prediction.fixed_parameters)

        for (p, v) in fixed_parameters.items():
            parameters.set(p, v)

        for o in prediction.observables:
            options_part = eos.QualifiedName(o.name).options_part()
            for key, value in options_part:
                if key in global_options and global_options[key] != value:
                    eos.error(f'Global option {key}={global_options[key]} overrides option part specification {key}={value} for observable {o.name} in prediction {_prediction} when using posterior {_posterior}.')
                if key in o.options and o.options[key] != value:
                    eos.error(f'Local option {key}={o.options[key]} overrides option part specification {key}={value} for observable {o.name} in prediction {_prediction}.')
            for key, global_value in global_options.items():
                if key in o.options and o.options[key] != global_value:
                    eos.warn(f'Local option {key}={o.options[key]} overrides global option {key}={global_value} for observable {o.name} in prediction {_prediction} when using posterior {_posterior}.')


        observables = []
        for o in prediction.observables:
            # Update global options with any options specified for this particular observable in the AnalysisFile
            local_options = eos.Options(**(global_options | o.options))

            if type(o.kinematics) == dict:
                kinematics = [o.kinematics]
            else:
                kinematics = [k for k in o.kinematics]

            for k in kinematics:
                observables.append(eos.Observable.make(
                    o.name,
                    parameters,
                    eos.Kinematics(k),
                    local_options
                ))

        if None in observables:
            unknown_observables = set()
            for p, o in zip(prediction.observables, observables):
                if o is None:
                    unknown_observables.add(p.name)
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
            step_name = step_desc.name
            tasks     = step_desc.tasks

            for task_desc in tasks:
                task = task_desc.task
                if task not in _tasks.keys():
                    raise ValueError(f'Unknown task \"{task}\" encountered in step {step_name}')

                arguments = { 'analysis_file': self }
                for key, value in step_arguments.items():
                    if key in _tasks.keys():
                        continue

                    arguments[key] = value

                if task_desc.arguments:
                    arguments.update(task_desc.arguments)

                if task in step_arguments:
                    arguments.update(step_arguments[task])

                if base_directory:
                    arguments.update({ 'base_directory': base_directory })

                return [(step_name, step_desc, task, a) for a in _expand_arguments(step_name, task, arguments)]

        # main part starts here

        # check inputs
        for step_idx, step_desc in enumerate(self._steps):
            if step_desc.iterations:
                if type(step_desc.iterations) is not list:
                    raise ValueError(f'Description of step #{step_idx} ({step_desc.name} expects a list of iterations')

        # resolve dependencies
        graph = nx.MultiDiGraph()
        graph.add_nodes_from([step.name for step in self._steps])
        for step in self._steps:
            if not step.depends_on:
                continue

            for dep in step.depends_on:
                graph.add_edge(dep, step.name)
        steps = list(nx.topological_sort(graph))

        # return steps and their tasks in the order imposed by the dependency graph
        result = []
        for step in steps:
            step_desc = { sd.name: sd for sd in self._steps }[step]
            iterations = [{}]
            if step_desc.iterations:
                iterations = step_desc.iterations

            for iteration_arguments in iterations:
                result.extend(_handle_single_step(step_desc, iteration_arguments))

        return result


    def validate(self):
        """Validates the analysis file."""
        messages = []
        # Check that all priors are known to EOS
        known_params = eos.Parameters()
        for p_name, pc in self._priors.items():
            for description in pc.descriptions:
                try:
                    known_params[description.parameter]
                except RuntimeError:
                    messages.append(f"Error in prior {p_name}: Parameter '{description.parameter}' not known to EOS")
        # Check that all likelihoods contain valid constraints and manual constraints
        known_constraints = eos.Constraints()
        for l_name, lc in self._likelihoods.items():
            for description in lc.constraints:
                try:
                    known_constraints[description.constraint]
                except RuntimeError:
                    messages.append(f"Error in likelihood {l_name}: Constraint '{description.constraint}'  not known to EOS")
            for manual in lc.manual_constraints:
                try:
                    known_constraints[manual.name]
                    messages.append(f"Error in likelihood {l_name}: The manual constraint named '{manual.name}' matches an already defined constraint name")
                except RuntimeError:
                    pass
        # Check that all observables in a prediction are known to EOS
        known_observables = eos.Observables()
        for pred_name, pd in self._predictions.items():
            for o in pd.observables:
                try:
                    known_observables[o.name]
                except RuntimeError:
                    messages.append(f"Error in prediction {pred_name}: Observable '{o.name}' not known to EOS")
        # Check that any parameters that are fixed (in posteriors or predictions) are known to EOS
        for p_name, posterior in self._posteriors.items():
            for param in posterior.fixed_parameters:
                try:
                    known_params[param]
                except RuntimeError:
                    messages.append(f"Error in posterior {p_name}: Fixed parameter '{param}' not known to EOS")
        for p_name, prediction in self._predictions.items():
            for param in prediction.fixed_parameters:
                try:
                    known_params[param]
                except RuntimeError:
                    messages.append(f"Error in prediction {p_name}: Fixed parameter '{param}' not known to EOS")

        # Check all the posteriors can be initialised, and used for the predictions specified in the analysis file
        # This will (hopefully) act as a catch all for any errors not spotted above
        for posterior in self._posteriors:
            try:
                self.analysis(posterior)
                eos.info(f'Successfully created analysis for posterior \'{posterior}\'')

                for prediction in self._predictions:
                    try:
                        self.observables(posterior, prediction, eos.Parameters())
                        eos.info(f'Successfully created prediction \'{prediction}\' set for posterior \'{posterior}\'')
                    except Exception as e:
                        messages.append(f'Error encountered when creating observables for prediction \'{prediction}\' of posterior \'{posterior}\': {e}')

            except Exception as e:
                messages.append(f'Error encountered when creating posterior \'{posterior}\': {e}')
        for m in messages:
            print(m)
        return messages


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

    def dump(self):
        """Dumps the contents of the analysis file in YAML format."""

        yaml.safe_dump(
            self.input_data,
            sys.stdout,
            width=100,
            allow_unicode=True,
            default_flow_style=False
            )

    def _repr_html_(self):
        result = r'''
            <table style="width: 80%">
                <colgroup>
                    <col width="40%" style="min-width: 160px">
                    <col width="15%" style="min-width:  80px">
                    <col width="15%" style="min-width:  80px">
                    <col width="30%" style="min-width:  80px">
                </colgroup>
                <thead>
                    <tr>
                        <th colspan="4" style="text-align: center">PRIORS</th>
                    </tr>
                </thead>
            '''
        for name, priors in self._priors.items():
            result += fr'''
                <thead>
                    <tr>
                        <th colspan="4" style="text-align: left">{name}</th>
                    </tr>
                    <tr>
                        <th style="text-align: center">qualified name</th>
                        <th style="text-align: center">min</th>
                        <th style="text-align: center">max</th>
                        <th style="text-align: center">type</th>
                    </tr>
                </thead>
                <tbody>
                '''
            for p in priors.descriptions:
                p_min, p_max = None, None
                if p.type == 'uniform':
                    p_min, p_max = p.min, p.max
                result += fr'''
                    <tr>
                        <td style="text-align: center"><code>{p.parameter}</code></td>
                        <td style="text-align: center">{p_min}</td>
                        <td style="text-align: center">{p_max}</td>
                        <td style="text-align: center">{p.type}</td>
                    </tr>
                '''

        result += r'''
            </tbody>
            </table>
            <br>
            <table style="width: 80%">
            <thead>
                <tr>
                    <th style="text-align: center">LIKELIHOODS</th>
                </tr>
            </thead>
        '''
        for name, likelihoods in self._likelihoods.items():
            result += rf'''
                <thead>
                    <tr>
                        <th style="text-align: left">{name}</th>
                    </tr>
                </thead>
                <tbody>
            '''
            for constraint in likelihoods.constraints:
                result += fr'''
                    <tr>
                        <td style="text-align: left"><code>{constraint.constraint}</code></td>
                    </tr>
                '''

        result += r'''
            </tbody>
            </table>
            <br>
            <table style="width: 80%">
            <colgroup>
                <col width="25%" id="priors"            style="min-width: 100px">
                <col width="25%" id="likelihoods"       style="min-width: 100px">
                <col width="25%" id="global_options"    style="min-width: 100px">
                <col width="25%" id="fixed_parameters"  style="min-width: 100px">
            </colgroup>
            <thead>
                <tr>
                    <th colspan="4" style="text-align: center">POSTERIORS</th>
                </tr>
            </thead>
        '''

        for name, posterior in self._posteriors.items():
            pp = ', '.join(posterior.prior)
            pl = ', '.join(posterior.likelihood)
            pg = ', '.join(f'{k} = {v}' for k, v in posterior.global_options.items())
            pf = ', '.join(f'{k} = {v}' for k, v in posterior.fixed_parameters.items())
            result += fr'''
                <thead>
                    <tr>
                        <th colspan="4" style="text-align: left">{name}</th>
                    </tr>
                    <tr>
                        <th style="text-align: center">priors</th>
                        <th style="text-align: center">likelihoods</th>
                        <th style="text-align: center">global options</th>
                        <th style="text-align: center">fixed parameters</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td style="text-align: center">{pp}</td>
                        <td style="text-align: center">{pl}</td>
                        <td style="text-align: center">{pg}</td>
                        <td style="text-align: center">{pf}</td>
                    </tr>
                </tbody>
                '''
        if self._predictions:
            result += r'''
                </tbody>
                </table>
                <br>
                <table style="width: 80%">
                <colgroup>
                    <col width="33%" id="observables"       style="min-width: 133px">
                    <col width="33%" id="global_options"    style="min-width: 133px">
                    <col width="33%" id="fixed_parameters"  style="min-width: 133px">
                </colgroup>
                <thead>
                    <tr>
                        <th colspan="3" style="text-align: center">PREDICTIONS</th>
                    </tr>
                </thead>
            '''

            for name, prediction in self._predictions.items():
                po = ', '.join(list({obs.name for obs in prediction.observables}))
                pg = ', '.join(f'{k} = {v}' for k, v in prediction.global_options.items())
                pf = ', '.join(f'{k} = {v}' for k, v in prediction.fixed_parameters.items())
                result += fr'''
                    <thead>
                        <tr>
                            <th colspan="3" style="text-align: left">{name}</th>
                        </tr>
                        <tr>
                            <th style="text-align: center">observables</th>
                            <th style="text-align: center">global options</th>
                            <th style="text-align: center">fixed parameters</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align: center">{po}</td>
                            <td style="text-align: center">{pg}</td>
                            <td style="text-align: center">{pf}</td>
                        </tr>
                    </tbody>
                    '''

        if self._obs:
            result += r'''
                </tbody>
                </table>
                <br>
                <table style="width: 80%">
                <colgroup>
                    <col width="25%" style="min-width:  100px">
                    <col width="25%" style="min-width:  100px">
                    <col width="25%" style="min-width:  100px">
                    <col width="25%" style="min-width:  100px">
                </colgroup>
                <thead>
                    <tr>
                        <th colspan="4" style="text-align: center">OBSERVABLES</th>
                    </tr>
                </thead>
            '''

            for obs in self._obs:
                oo = ', '.join(f'{k} = {v}' for k, v in obs.options.items())
                result += fr'''
                    <thead>
                        <tr>
                            <th colspan="1" style="text-align: left">{obs.name}</th>
                            <th colspan="2" style="text-align: left">{obs.latex}</th>
                            <th colspan="1" style="text-align: left">unit: {obs.unit}</th>
                        </tr>
                        <tr>
                            <th colspan="2" style="text-align: center">expression</th>
                            <th colspan="2" style="text-align: center">options</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td colspan="2" style="text-align: center">{obs.expression}</td>
                            <td colspan="2" style="text-align: center">{oo}</td>
                        </tr>
                    </tbody>
                    '''

        if self._params:
            result += r'''
                </tbody>
                </table>
                <br>
                <table style="width: 80%">
                <colgroup>
                    <col width="20%" style="min-width:  80px">
                    <col width="20%" style="min-width:  80px">
                    <col width="20%" style="min-width:  80px">
                    <col width="40%" style="min-width: 160px">
                </colgroup>
                <thead>
                    <tr>
                        <th colspan="4" style="text-align: center">PARAMETERS</th>
                    </tr>
                </thead>
            '''

            for param in self._params:
                po = ', '.join(param.alias_of)
                result += fr'''
                    <thead>
                        <tr>
                            <th colspan="1" style="text-align: left">{param.name}</th>
                            <th colspan="2" style="text-align: left">{param.latex}</th>
                            <th colspan="1" style="text-align: left">unit: {param.unit}</th>
                        </tr>
                        <tr>
                            <th style="text-align: center">central</th>
                            <th style="text-align: center">min</th>
                            <th style="text-align: center">max</th>
                            <th style="text-align: center">alias</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td style="text-align: center">{param.central}</td>
                            <td style="text-align: center">{param.min}</td>
                            <td style="text-align: center">{param.max}</td>
                            <td style="text-align: center">{po}</td>
                        </tr>
                    </tbody>
                    '''

        result += r'''
        </table>
        '''

        return result
