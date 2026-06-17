#!/usr/bin/python
# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2023 Lorenz Gärtner
# Copyright (c) 2023 Danny van Dyk
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
import json

class PyhfLogLikelihood:
    r"""Wraps a HistFactory likelihood, specified as a pyhf workspace, for use within EOS.

    On construction, the pyhf model's modifiers (its parameters) are mapped onto EOS observables or
    parameters: each pyhf parameter is associated with the EOS observable or parameter named by
    ``parameter_map``, defaulting to a parameter ``pyhf::<name>`` that is declared if it does not yet
    exist. The corresponding observables are registered with the supplied observable cache, so that
    :meth:`evaluate` can compute the main term of the pyhf log-likelihood from the current EOS parameter
    values. The constraint terms of the pyhf likelihood are not evaluated here; they are instead exposed
    as EOS priors via :meth:`priors`.

    Using this class requires the optional ``pyhf`` module.

    :param cache: The EOS observable cache to which the observables backing the pyhf parameters are added.
    :type cache: eos.ObservableCache
    :param workspace: A pyhf workspace, or the path to a JSON file specifying one.
    :type workspace: pyhf.workspace.Workspace | str
    :param parameter_map: An optional mapping from pyhf parameter names to EOS observables or parameters.
        Each value is either the qualified name of an EOS observable/parameter, or a dictionary with a
        ``name`` key and optional ``kinematics`` and ``options`` keys. Parameters not listed default to an
        EOS parameter named ``pyhf::<name>``.
    :type parameter_map: dict | None
    :raises RuntimeError: If the ``pyhf`` module is not installed.
    :raises ValueError: If ``workspace`` is neither a pyhf workspace nor a path to a valid JSON file, or
        if a ``parameter_map`` value is neither a string nor a dictionary.
    """

    def __init__(self, cache, workspace, parameter_map=None):
        if parameter_map is None:
            parameter_map = {}
        try:
            import pyhf
        except ModuleNotFoundError:
            raise RuntimeError('Attempting to use a PyHF likelihood without the `pyhf` module installed')
        self.cache = cache
        self.parameters = cache.parameters()
        self.parameter_map = parameter_map

        if isinstance(workspace, pyhf.workspace.Workspace):
            self.workspace = workspace
        else:
            try:
                with open(workspace) as f:
                    f = json.load(f)
                self.workspace = pyhf.Workspace(f)
            except Exception:
                raise ValueError('`workspace` must be a pyhf workspace or json file specifying one. ')

        self.model = self.workspace.model()
        self.data = np.array(self.workspace.data(self.model, include_auxdata=False))

        self._pyhf_parameters = []
        for parameter_name, init, bounds in zip(
            self.model.config.par_names,
            self.model.config.suggested_init(),
            self.model.config.suggested_bounds()
            ):

            kinematics = eos.Kinematics()
            options = eos.Options()

            if parameter_name in parameter_map:
                if isinstance(parameter_map[parameter_name], str):
                    qualified_name = eos.QualifiedName(
                        parameter_map[parameter_name]
                        )
                elif isinstance(parameter_map[parameter_name], dict):
                    qualified_name = eos.QualifiedName(
                        parameter_map[parameter_name]['name']
                        )
                    kinematics = eos.Kinematics(parameter_map[parameter_name].get('kinematics', {}))
                    options = eos.Options(parameter_map[parameter_name].get('options', {}))
                else:
                    raise ValueError('parameter_map values must be either strings or dictionaries.')
            else:
                qualified_name = eos.QualifiedName('pyhf::' + parameter_name)

            # check if qualified name corresponds to an existing observable
            observables = eos.Observables()
            if qualified_name in observables:
                pyhf_obs = eos.Observable.make(qualified_name, self.parameters, kinematics, options)
                self._pyhf_parameters.append(pyhf_obs)
            else:
                # check if qualified name corresponds to an existing parameter
                if self.parameters.has(qualified_name):
                    pyhf_param = self.parameters[qualified_name]
                # otherwise, create a new parameter
                else:
                    pyhf_param = self.parameters.declare_and_insert(
                        qualified_name,
                        parameter_name,
                        eos.Unit.Undefined(),
                        init,
                        bounds[0],
                        bounds[1]
                        )
                pyhf_param.set(init)
                self._pyhf_parameters.append(pyhf_param)
                pyhf_obs = eos.Observable.make(qualified_name, self.parameters, kinematics, options)

            self.cache.add(pyhf_obs)

        self.number_of_observations = int(1)

    def evaluate(self):
        """Evaluate the main term of the pyhf log-likelihood at the current parameter values.

        The pyhf parameters are read from their associated EOS observables/parameters, and the main
        (non-constraint) term of the HistFactory log-likelihood is returned. Constraint terms are handled
        separately within EOS via the priors returned by :meth:`priors`.

        :returns: The value of the main log-likelihood term.
        :rtype: float
        """
        parameter_values = np.array([p.evaluate() for p in self._pyhf_parameters])
        # only main term in pyhf - constraints are handled in EOS
        return self.model.mainlogpdf(self.data, parameter_values).item()

    @staticmethod
    def factory(cache, workspace, parameter_map=None):
        """Construct a :class:`PyhfLogLikelihood`.

        Convenience factory with the same signature as the constructor; see the class documentation for a
        description of the arguments.

        :returns: The constructed likelihood.
        :rtype: PyhfLogLikelihood
        """
        return PyhfLogLikelihood(cache, workspace, parameter_map)

    def priors(self):
        """Return the EOS prior descriptions corresponding to the pyhf constraint terms.

        Each modifier (parameter) of the pyhf model is translated into an EOS prior description:
        unconstrained parameters yield a uniform prior over their suggested bounds, while parameters with
        a normal or Poisson constraint yield Gaussian or Poisson priors, respectively. Parameters mapped
        to an EOS observable (a dictionary entry in ``parameter_map``) are skipped, since they are not
        free parameters.

        :returns: A list of prior descriptions, each a dictionary with a ``parameter`` name, a ``type``,
            and the type-specific keys (e.g. ``min``/``max``, ``central``/``sigma``, or ``k``).
        :rtype: list[dict]
        """
        priors = []
        for name, _ in self.model.config.modifiers:
            param_name = 'pyhf::' + name
            if name in self.parameter_map:
                if isinstance(self.parameter_map[name], str):
                    param_name = self.parameter_map[name]
                elif isinstance(self.parameter_map[name], dict):
                    continue
            param_set = self.model.config.param_set(name)

            if not param_set.constrained:
                prior = {
                    'parameter': param_name,
                    'min': param_set.suggested_bounds[0][0],
                    'max': param_set.suggested_bounds[0][1],
                    'type': 'uniform',
                }
                priors.append(prior)
                continue

            if param_set.pdf_type == 'normal':
                for i in range(param_set.n_parameters):
                    prior = {
                        'parameter': param_name + ('' if param_set.is_scalar else f'[{i}]'),
                        'type': 'gauss',
                        'central': param_set.auxdata[i],
                        'sigma': param_set.width()[i]
                    }
                    priors.append(prior)

            elif param_set.pdf_type == 'poisson':
                for i in range(param_set.n_parameters):
                    prior = {
                        'parameter': param_name + ('' if param_set.is_scalar else f'[{i}]'),
                        'type': 'poisson',
                        'k': param_set.auxdata[i],
                    }
                    priors.append(prior)
        return priors
