# Copyright (c) 2023 Danny van Dyk
# Copyright (c) 2023 Philip Lueghausen
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

class Item:
    "Can display itself using a Matplotlib ax object"

    def __init__(self):
        "Setup of the item attributes and validation of all supplied arguments"
        pass

    def prepare(self):
        "Setup of the item such as computing data points of an observable"
        pass

    def draw(self, ax):
        """
        Draw the item using a set of axis

        :param ax: Matplotlib axis object
        :type ax: matplotlib.axes
        """

        self.prepare()
        # Call to matplotlib method on ax here
        pass


class Observable(Item):
    r"""Represents an observable as a function of a variable.
    The variable can be either a kinematic variable or a parameter.


    :param observable: EOS qualified name of an observable
    :type observable: str
    :param variable: Name of a kinematic variable or of a parameter
    :type variable: str
    :param x_values: Values to be used as sampling points
    :type x_values: list of float numbers
    :param fixed_kineamtics: Names and values of fixed kinematic variables
    :type fixed_kineamtics: dict
    :param fixed_parameters: Names and values of fixed parameters
    :type fixed_parameters: dict
    :param fixed_parameters_from_file: Path to file that contains names and values of fixed parameters in the YAML format
    :type fixed_parameters_from_file: str
    :param options: Names and values of options passed to the EOS observable
    :type options: dict
    :param label: Label used in the legend
    :type label: str
    :param \**kwargs: Additional keyword arguments to pass to ``matplotlib.axes.Axes.plot``
    """

    def __init__(self,
            observable,
            variable:str,
            x_values:list,
            fixed_kinematics:dict=None,
            fixed_parameters:dict=None,
            fixed_parameters_from_file:str=None,
            options:dict=None,
            label:str=None,
            **kwargs
        ):
        super().__init__()

        obs_entry = eos.Observables._get_obs_entry(observable)
        valid_kin_vars = [kv for kv in obs_entry.kinematic_variables()]
        eos.info(f'Plotting EOS observable "{observable}"')

        self.label = label


        # Create kinematics
        eos_kinematics = eos.Kinematics()
        if fixed_kinematics:
            for k, v in fixed_kinematics.items():
                if k not in valid_kin_vars:
                    raise ValueError("Kinematic variable '" + k + "' does not " +
                    "match any of the declared kinematic variables '" + observable + "': " + valid_kin_vars.__repr__())
                eos_kinematics.declare(k, v)

        # Create parameters
        eos_parameters = eos.Parameters.Defaults()
        if fixed_parameters and fixed_parameters_from_file:
            eos.warn('Overriding values read from \'parameters-from-file\' with explicit values in \'parameters\'')
            raise RuntimeError("Providing 'fixed_parameters_from_file' is not yet implemented")

        elif fixed_parameters_from_file:
            eos.warn('Overriding parameters from file')
            parameters.override_from_file(fixed_parameters_from_file)

        elif fixed_parameters:
            for key, value in fixed_parameters.items():
                eos_parameters.set(key, value)


        # [To do] Check: variable is also specified as kinematic or parameter?
        # if 'parameters' in item and item['variable'] in item['parameters']:
        #     val = item['parameters'].get(item['variable'])
        #     raise ValueError("Variable '" + item['variable'] + "' of observable '" + oname + "' is " +
        #             "also specified as a fix parameter with value " + str(val))
        # if 'kinematics' in item and item['variable'] in item['kinematics']:
        #     val = item['kinematics'].get(item['variable'])
        #     raise ValueError("Variable '" + item['variable'] + "' of observable '" + oname + "' is " +
        #             "also specified as a fix kinematic with value " + str(val))


        # Declare variable that is either parameter or kinematic
        self.var = None

        # Does the variable correspond to one of the kinematic variables?
        if variable in valid_kin_vars:
            self.var = eos_kinematics.declare(variable, np.nan)
        else:
            # Is the variable name a QualifiedName?
            try:
                qn = eos.QualifiedName(variable)
                # Continues only if no failure occures
                self.var = eos_parameters.declare(qn, np.nan)
            except RuntimeError:
                raise ValueError("Value of 'variable' for observable '" + observable +
                    "' is neither a valid kinematic variable nor parameter: '" + variable + "'")


        # Create options
        eos_options = eos.Options()
        if options:
            for key, value in options.items():
                options.declare(key, value)

        self.eos_observable = eos.Observable.make(observable, eos_parameters, eos_kinematics, eos_options)

        self.x_values = x_values
        self.kwargs = kwargs


    def prepare(self):
        "Evaluate the observable at the sample points"
        self.y_values = np.empty((len(self.x_values),))
        for i, x in enumerate(self.x_values):
            self.var.set(x)
            self.y_values[i] = self.eos_observable.evaluate()


    def draw(self, ax):
        "Draw a line plot of the observable"
        super().draw(ax)
        ax.plot(self.x_values, self.y_values, label=self.label, **self.kwargs)


class ItemFactory:
    registry = {
        'observable': Observable,
    }


    @staticmethod
    def make(ax, item_type:str, **kwargs):
        if item_type not in ItemFactory.registry:
            raise ValueError(f'Unknown content item type: {item_type}')

        return ItemFactory.registry[item_type](**kwargs)
