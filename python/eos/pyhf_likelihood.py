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
    def __init__(self, cache, workspace, parameter_map={}):
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
            except:
                raise ValueError('`workspace` must be a pyhf workspace or json file specifying one. ')

        self.model = self.workspace.model()
        self.data = np.array(self.workspace.data(self.model, include_auxdata=False))

        self.varied_parameters = []
        for pn, v, bounds in zip(self.model.config.par_names,
                                 self.model.config.suggested_init(),
                                 self.model.config.suggested_bounds()):

            if pn in parameter_map:
                qn = eos.QualifiedName(parameter_map[pn])
            else:
                qn = eos.QualifiedName('pyhf::' + pn)

            if self.parameters.has(qn):
                p = self.parameters[qn]
            else:
                p = self.parameters.declare_and_insert(qn, pn, eos.Unit.Undefined(), v, bounds[0], bounds[1])

            p.set(v)
            self.varied_parameters.append(p)

            o = eos.Observable.make(qn, self.parameters, eos.Kinematics(), eos.Options())
            self.cache.add(o)

        self.number_of_observations = int(1)

    def evaluate(self):
        parameter_values = np.array([p.evaluate() for p in self.varied_parameters])
        # only main term in pyhf - constraints are handled in EOS
        return self.model.mainlogpdf(self.data, parameter_values).item()

    @staticmethod
    def factory(cache, workspace, parameter_map={}):
        return PyhfLogLikelihood(cache, workspace, parameter_map)

    def priors(self):
        priors = []
        for name, _ in self.model.config.modifiers:
            param_name = self.parameter_map[name] if name in self.parameter_map else 'pyhf::' + name
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
